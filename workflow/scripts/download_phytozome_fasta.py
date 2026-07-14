"""Download + cache a single Phytozome genome FASTA for extraction.

Driven by the ``download_phytozome_fasta`` rule (Snakemake ``script:``). Mirrors
``download_phytozome_gtf`` for auth/PURGED handling, but resolves the *assembly*
FASTA of the same Phytozome version as the resolved gene annotation (so the
coordinates and the sliced sequence come from the same assembly), decompresses
it, and indexes it where ``extract_sequences`` looks:
``resources/cache/phytozome_<species>/genome.fasta{,.fai}``.

Fault tolerance (like ``download_assembly``, unlike the blocking GFF3 rule):
this always exits 0 and writes a ``.download_done`` sentinel (``ok`` /
``failed: <reason>``). A cold-storage (PURGED) or missing FASTA then surfaces
downstream as extraction ``assembly_not_cached`` instead of blocking the run —
the FASTA is only needed for extraction, not for resolution.
"""

import gzip
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from jgi_phytozome_lookup import request_restore, resolve_sequence


def load_token():
    """JGI bearer token (no "Bearer " prefix) from the environment or .env."""
    token = (
        os.environ.get("JGI_SESSION_TOKEN") or os.environ.get("PHYTOZOME_BEARER") or ""
    ).strip()
    if not token and Path(".env").exists():
        # ponytail: hand-rolled .env read, mirrors download_phytozome_gtf.py.
        for line in Path(".env").read_text().splitlines():
            line = line.strip()
            if line.startswith(("JGI_SESSION_TOKEN=", "PHYTOZOME_BEARER=")):
                candidate = line.split("=", 1)[1].strip().strip("\"'")
                if candidate:
                    token = candidate
                    break
    if token[:7].lower() == "bearer ":
        token = token[7:].strip()
    return token


def _gunzip_and_index(gz_path, fasta_out):
    """Decompress .fa.gz -> fasta_out and build a .fai index.

    ponytail: ~same as download_assembly.py's decompress+faidx; two call sites,
    not worth a shared module yet.
    """
    with gzip.open(gz_path, "rb") as src, fasta_out.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    subprocess.run(["samtools", "faidx", str(fasta_out)], check=True)


def _run(snakemake, log_path):
    species = snakemake.wildcards.species
    status_path = Path(snakemake.output.status)
    cache_dir = status_path.parent
    cache_dir.mkdir(parents=True, exist_ok=True)
    fasta_out = cache_dir / "genome.fasta"
    fai_out = cache_dir / "genome.fasta.fai"

    def done(status):
        status_path.write_text(status + "\n")
        log_path.write_text(f"[{species}] {status}\n")

    # Idempotent cache guard: already downloaded + indexed.
    if fasta_out.exists() and fai_out.exists():
        done("ok")
        return

    token = load_token()
    if not token:
        done("failed: no JGI token (set JGI_SESSION_TOKEN/PHYTOZOME_BEARER in .env)")
        return

    sources = snakemake.config.get("phytozome_gtf_sources", snakemake.config)
    cfg = sources.get(species) or {}
    genome_id = cfg.get("genome_id")
    if genome_id is None:
        done(f"failed: no genome_id in config for {species}")
        return
    # Pin the annotation the FASTA must mate (so the assembly version matches the
    # GFF3 used for resolution); the config portal_file_name is the tracked pin.
    prefer_name = cfg.get("portal_file_name")

    info = resolve_sequence(genome_id, f"Bearer {token}", prefer_name=prefer_name)
    if info is None:
        done(
            f"failed: JGI returned no assembly FASTA for {species} (genome_id={genome_id})"
        )
        return

    if str(info.get("file_status", "")).upper() == "PURGED":
        # On tape: fire a restore (safe to repeat) and report failed. Rerun once
        # JGI stages it (≤24h). Non-blocking so the rest of the run proceeds.
        resp = request_restore(info["_id"], f"Bearer {token}")
        done(
            f"failed: restore requested ({info['_id']}); rerun after JGI stages it: {resp}"
        )
        return

    gz_path = cache_dir / "assembly.fa.gz"
    try:
        request = urllib.request.Request(
            info["download_url"], headers={"Authorization": f"Bearer {token}"}
        )
        with urllib.request.urlopen(request, timeout=600) as response, gz_path.open(
            "wb"
        ) as out_fh:
            shutil.copyfileobj(response, out_fh)
        _gunzip_and_index(gz_path, fasta_out)
    except Exception as exc:
        # Clean partial outputs so a rerun retries from scratch.
        for p in (gz_path, fasta_out, fai_out):
            p.unlink(missing_ok=True)
        done(f"failed: {exc}")
        return
    finally:
        gz_path.unlink(missing_ok=True)

    done(f"ok: {info['file_name']}")


def main(snakemake):
    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        _run(snakemake, log_path)
    except Exception as exc:  # never fail the job — keep the sentinel non-blocking
        Path(snakemake.output.status).write_text(f"failed: {exc}\n")
        if not log_path.exists() or not log_path.read_text().strip():
            log_path.write_text(f"ERROR downloading Phytozome FASTA: {exc}\n")


if __name__ == "__main__":
    main(snakemake)  # noqa: F821  (injected by Snakemake)
