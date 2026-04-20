# ============================================================
# Rules: download_plant_gtf  +  resolve_plant_gtf
# ============================================================
# GTF-based fallback for plant transcript IDs that BioMart
# could not resolve.
#
# Species and GTF URLs are configured in
# config/plant_gtf_sources.yaml (merged into config by the
# top-level "configfile:" directive in the Snakefile).
#
# Chain:
#   parse_ids           → classified_ids.tsv
#           ↓
#   resolve_plant_gtf   → plant_gtf_resolved.tsv
#                         plant_gtf_unresolved.tsv
#           ↓
#   gramene_resolver    → gramene_resolved.tsv
# ============================================================

_PLANT_GTF_SOURCES = config.get("plant_gtf_sources", {})


# ── Download ─────────────────────────────────────────────────

rule download_plant_gtf:
    """Download a single species GTF from Ensembl Plants FTP."""
    output:
        "resources/plant_gtf/{species}.gtf.gz",
    params:
        url = lambda wc: _PLANT_GTF_SOURCES[wc.species]["url"],
    log:
        f"{LOGS}/download_plant_gtf/{{species}}.log",
    benchmark:
        f"{BENCHMARKS}/download_plant_gtf/{{species}}.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 60,
        mem_mb          = 1024,
        cpus_per_task   = 1,
    shell:
        "wget --quiet --tries=5 --timeout=120 -O {output} {params.url} 2> {log}"


# ── Resolve ──────────────────────────────────────────────────

def _plant_gtf_inputs(wildcards):
    """Collect all configured GTF files as named inputs."""
    species_list = list(_PLANT_GTF_SOURCES.keys())
    return {
        "classified": f"{RESULTS}/classified_ids.tsv",
        "gtf_files": expand(
            "resources/plant_gtf/{species}.gtf.gz",
            species=species_list,
        ),
    }


rule resolve_plant_gtf:
    """Check BioMart-unresolved plant IDs against downloaded GTF files."""
    input:
        unpack(_plant_gtf_inputs),
    output:
        resolved   = f"{RESULTS}/plant_gtf_resolved.tsv",
        unresolved = f"{RESULTS}/plant_gtf_unresolved.tsv",
    log:
        f"{LOGS}/resolve_plant_gtf.log",
    benchmark:
        f"{BENCHMARKS}/resolve_plant_gtf.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 60,
        mem_mb          = 8192,   # GTF indices can be large for maize / soybean
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_plant_gtf.py"
