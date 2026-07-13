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
    """Download or stage a single Phytozome GFF3 from a manifest entry.

    Layout: resources/phytozome/<species>/<source_file_name>.gff3.gz. The
    {species} wildcard is the folder — it equals the config key, so the script
    can look the genome_id up by it. The inner file keeps the JGI source name
    for traceability. The constraint stops {species} swallowing the '/'.
    """
    wildcard_constraints:
        species = r"[^/]+",
    input:
        manifest = "resources/phytozome/manifest.json",
    output:
        protected("resources/phytozome/{species}/{gff}"),
    log:
        f"{LOGS}/download_phytozome_gtf/{{species}}/{{gff}}.log",
    benchmark:
        f"{BENCHMARKS}/download_phytozome_gtf/{{species}}/{{gff}}.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 60,
        mem_mb          = 1024,
        cpus_per_task   = 1,
    script:
        "../scripts/download_phytozome_gtf.py"


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
                if "gtf" not in source:
                    raise ValueError(f"Species {s} missing 'gtf' field in phytozome_gtf_sources config")
    return {
        "classified": f"{RESULTS}/classified_ids.tsv",
        "gff_files": [phytozome_sources[s]["gtf"] for s in species_list],
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
