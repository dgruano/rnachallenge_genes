"""Download or stage a single Phytozome GFF3 from a manifest entry.

Driven by the ``download_phytozome_gtf`` rule (Snakemake ``script:`` directive).

The numeric portal ``file_id`` stored in manifests is NOT a JGI download id;
the real download URL needs the Mongo ``_id`` from the file-list API. We resolve
it from the genome_id (config source of truth) via jgi_phytozome_lookup, pinning
the manifest's ``portal_file_name`` when present.
"""

import json
import os
import shutil
import sys
import urllib.request
from pathlib import Path

from snakemake.exceptions import WorkflowError

sys.path.insert(0, str(Path(__file__).resolve().parent))
from jgi_phytozome_lookup import request_restore, resolve_annotation


def load_token():
    """JGI bearer token (no "Bearer " prefix) from the environment or .env."""
    token = (
        os.environ.get("JGI_SESSION_TOKEN") or os.environ.get("PHYTOZOME_BEARER") or ""
    ).strip()
    if not token and Path(".env").exists():
        # ponytail: hand-rolled .env read, swap for python-dotenv if more keys appear.
        for line in Path(".env").read_text().splitlines():
            line = line.strip()
            if line.startswith(("JGI_SESSION_TOKEN=", "PHYTOZOME_BEARER=")):
                candidate = line.split("=", 1)[1].strip().strip("\"'")
                if candidate:
                    token = candidate
                    break
    # The value may already carry a "Bearer " prefix; strip it so the
    # f"Bearer {token}" call sites don't double it (which 401s).
    if token[:7].lower() == "bearer ":
        token = token[7:].strip()
    return token


def main(snakemake):
    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        _run(snakemake, log_path)
    except Exception as exc:
        # Never fail silently: the empty-log failures (config-key/basename
        # mismatch) were a diagnosis trap. Record the error before re-raising.
        if not log_path.exists() or not log_path.read_text().strip():
            log_path.write_text(f"ERROR resolving Phytozome GFF3: {exc}\n")
        raise


def _run(snakemake, log_path):
    species = snakemake.wildcards.species
    manifest_path = Path(snakemake.input.manifest)
    output_path = Path(snakemake.output[0])

    output_path.parent.mkdir(parents=True, exist_ok=True)

    manifest = json.loads(manifest_path.read_text())
    entry = manifest.get(species)
    if entry is None and isinstance(manifest.get("species"), dict):
        entry = manifest["species"].get(species)
    # The manifest only pins portal_file_name; genome_id in config is the source
    # of truth. A missing entry is fine — fall through to the config lookup below.
    if entry is None:
        entry = {}

    # The manifest's "status" is a stale snapshot. Don't gate on it here: the
    # live JGI resolution below checks the real current status and fires a
    # restore for PURGED files, so trusting the frozen field would block reruns
    # forever. We only fall back to it when there's no way to check live.
    stale_status = str(entry.get("status", "RESTORED")).upper()

    # Local/staged source short-circuit.
    source_path = entry.get("local_path") or entry.get("path") or entry.get("gtf")
    if source_path:
        source = Path(source_path)
        if not source.exists():
            raise WorkflowError(
                f"Configured/local Phytozome source missing for {species}: {source}"
            )
        if source.resolve() != output_path.resolve():
            shutil.copyfile(source, output_path)
        else:
            output_path.touch()
        log_path.write_text(f"staged {species} from {source}\n")
        return

    token = load_token()
    if not token:
        raise WorkflowError(
            f"Downloading Phytozome GFF3 for {species} needs a JGI bearer token. "
            f"Add JGI_SESSION_TOKEN (or PHYTOZOME_BEARER) to .env or the environment."
        )

    sources = snakemake.config.get("phytozome_gtf_sources", snakemake.config)
    config_source = sources.get(species) or {}

    genome_id = entry.get("genome_id")
    if genome_id is None:
        genome_id = config_source.get("genome_id")

    # Pin the exact JGI file. Config is the tracked source of truth (the manifest
    # is gitignored); fall back to a manifest pin when config doesn't set one.
    # Needed because JGI's selection heuristic prefers gene_exons over .gene.
    prefer_name = config_source.get("portal_file_name") or entry.get("portal_file_name")

    download_url = None
    if genome_id is not None:
        info = resolve_annotation(genome_id, f"Bearer {token}", prefer_name=prefer_name)
        if info is None:
            raise WorkflowError(
                f"JGI returned no gene annotation for {species} (genome_id={genome_id})"
            )
        if str(info.get("file_status", "")).upper() == "PURGED":
            # On tape: fire a restore request (safe to repeat) and fail with a
            # clear message. Restoration can take up to 24h, so we don't block
            # the job — rerun once JGI has staged the file to disk.
            resp = request_restore(info["_id"], f"Bearer {token}")
            log_path.write_text(
                f"requested JGI restore for {species} ({info['_id']}): {resp}\n"
            )
            raise WorkflowError(
                f"Phytozome GFF3 for {species} is on cold storage (PURGED). "
                f"Requested restoration from JGI (response={resp}); this can take "
                f"up to 24h. Rerun once it completes."
            )
        download_url = info["download_url"]

    # No genome_id to check live status with: honor the stale PURGED snapshot
    # so we emit a clear "restore it" message instead of a raw download error.
    if genome_id is None and stale_status in {
        "PURGED",
        "COLD",
        "COLD_STORAGE",
        "ARCHIVED",
    }:
        raise WorkflowError(
            f"Phytozome file for {species} is not downloadable (manifest "
            f"status={stale_status}) and has no genome_id to resolve a live "
            f"restore. Add a genome_id in config or restore it via JGI and rerun."
        )

    # Legacy fallback: explicit URL in the manifest (no genome_id available).
    if not download_url:
        download_url = entry.get("download_url") or entry.get("url")
    if not download_url:
        raise WorkflowError(
            f"Manifest entry for {species} needs a genome_id (via config) "
            f"or an explicit download_url/url"
        )

    request = urllib.request.Request(
        download_url, headers={"Authorization": f"Bearer {token}"}
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response, output_path.open(
            "wb"
        ) as out_fh:
            shutil.copyfileobj(response, out_fh)
    except Exception as exc:
        raise WorkflowError(
            f"Failed to download Phytozome GFF3 for {species}: {exc}"
        ) from exc

    log_path.write_text(f"downloaded {species} from {download_url}\n")


if __name__ == "__main__":
    main(snakemake)  # noqa: F821  (injected by Snakemake)
