# ============================================================
# Rules: download_phytozome_gtf  +  resolve_phytozome_gtf
# ============================================================
# GFF3-based fallback for plant transcript IDs backed by
# Phytozome annotations.
#
# Species metadata are loaded from config/phytozome_gtf_sources.yaml.
# The config may be either:
#   phytozome_gtf_sources:
#     <species>: {...}
# or a flat top-level mapping:
#   <species>: {...}
#
# The download step consumes a manifest produced elsewhere at
# resources/phytozome/manifest.json and supports either:
#   - a local source path (gtf/path/local_path), or
#   - a JGI file id / direct download URL.
# ============================================================

from pathlib import Path


def _phytozome_sources():
    nested = config.get("phytozome_gtf_sources")
    if isinstance(nested, dict):
        return nested

    return {
        key: value
        for key, value in config.items()
        if isinstance(value, dict)
        and any(field in value for field in ("species_query", "genome_id", "gtf", "phytozome_version"))
    }


_PHYTOZOME_GTF_SOURCES = _phytozome_sources()


# ── Download ─────────────────────────────────────────────────

rule download_phytozome_gtf:
    """Download or stage a single Phytozome GFF3 from a manifest entry."""
    input:
        manifest = "resources/phytozome/manifest.json",
    output:
        protected("resources/phytozome/{species}.gff3.gz"),
    log:
        f"{LOGS}/download_phytozome_gtf/{{species}}.log",
    benchmark:
        f"{BENCHMARKS}/download_phytozome_gtf/{{species}}.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 60,
        mem_mb          = 1024,
        cpus_per_task   = 1,
    run:
        import json
        import os
        import shutil
        import urllib.request
        from snakemake.exceptions import WorkflowError

        species = wildcards.species
        manifest_path = Path(input.manifest)
        output_path = Path(output[0])
        log_path = Path(log[0])

        output_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        with manifest_path.open() as fh:
            manifest = json.load(fh)

        entry = manifest.get(species)
        if entry is None and isinstance(manifest.get("species"), dict):
            entry = manifest["species"].get(species)
        if entry is None:
            raise WorkflowError(f"Phytozome manifest entry not found for species={species!r}")

        status = str(entry.get("status", "RESTORED")).upper()
        if status in {"PURGED", "COLD", "COLD_STORAGE", "ARCHIVED"}:
            raise WorkflowError(
                f"Phytozome file for {species} is not downloadable yet "
                f"(status={status}). Restore it via the JGI workflow and rerun."
            )

        source_path = entry.get("local_path") or entry.get("path") or entry.get("gtf")
        if source_path:
            source = Path(source_path)
            if not source.exists():
                raise WorkflowError(f"Configured/local Phytozome source missing for {species}: {source}")
            if source.resolve() != output_path.resolve():
                shutil.copyfile(source, output_path)
            else:
                output_path.touch()
            with log_path.open("w") as fh:
                fh.write(f"staged {species} from {source}\n")
            return

        download_url = entry.get("download_url") or entry.get("url")
        if not download_url:
            file_id = entry.get("file_id") or entry.get("_id")
            if file_id:
                download_url = f"https://files-download.jgi.doe.gov/download_files/{file_id}/"

        if not download_url:
            raise WorkflowError(
                f"Manifest entry for {species} must provide one of "
                f"local_path/path/gtf/url/download_url/file_id"
            )

        headers = {}
        token = os.environ.get("JGI_SESSION_TOKEN", "").strip()
        if token:
            headers["Authorization"] = f"Bearer {token}"

        request = urllib.request.Request(download_url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=120) as response, output_path.open("wb") as out_fh:
                shutil.copyfileobj(response, out_fh)
        except Exception as exc:
            raise WorkflowError(f"Failed to download Phytozome GFF3 for {species}: {exc}") from exc

        with log_path.open("w") as fh:
            fh.write(f"downloaded {species} from {download_url}\n")


# ── Resolve ──────────────────────────────────────────────────

def _phytozome_gtf_inputs_bak(wildcards):
    species_list = list(_PHYTOZOME_GTF_SOURCES.keys())
    return {
        "classified": f"{RESULTS}/classified_ids.tsv",
        "gff_files": expand(
            "resources/phytozome/{species}.gff3.gz",
            species=species_list,
        ),
    }


def _phytozome_gtf_inputs(wildcards):
    species_list = list(_PHYTOZOME_GTF_SOURCES.keys())
    gff_files = []
    phytozome_sources = config.get("phytozome_gtf_sources", {})
    if not phytozome_sources:
        raise ValueError("No phytozome_gtf_sources found in config")
    else:
        missing_species = [s for s in species_list if s not in phytozome_sources]
        if missing_species:
            raise ValueError(f"Species {missing_species} missing from phytozome_gtf_sources config")
        else:
            for s in species_list:
                source = phytozome_sources[s]
                #if not all(field in source for field in ("species_query", "genome_id", "gtf")):
                if "gtf" not in source:
                    raise ValueError(f"Species {s} missing 'gtf' field in phytozome_gtf_sources config")
                else:
                    gtf_path = Path(source["gtf"])
                    if not gtf_path.exists():
                        raise ValueError(f"GTF file for species {s} not found at {gtf_path}")
                    else:
                        gff_files.append(str(gtf_path))
    return {
        "classified": f"{RESULTS}/classified_ids.tsv",
        "gff_files": gff_files,
    }


rule resolve_phytozome_gtf:
    """Resolve configured Phytozome-backed plant IDs against GFF3 annotations."""
    input:
        unpack(_phytozome_gtf_inputs),
    output:
        resolved   = f"{RESULTS}/phytozome_gtf_resolved.tsv",
        unresolved = f"{RESULTS}/phytozome_gtf_unresolved.tsv",
    log:
        f"{LOGS}/resolve_phytozome_gtf.log",
    benchmark:
        f"{BENCHMARKS}/resolve_phytozome_gtf.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 60,
        mem_mb          = 8192,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_phytozome_gtf.py"
