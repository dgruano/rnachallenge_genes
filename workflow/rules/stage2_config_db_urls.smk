# ============================================================
# Rules: Stage 2 URL discovery for config-backed DBs
# ============================================================
# Reads annotation source YAMLs (plant_gtf_sources,
# metazoa_gtf_sources, yeast_gtf_sources) and emits a
# per-assembly URL table used by merge_resolved (Stage 3).
#
# Each rule writes a TSV with columns:
#   db_source, assembly_name, assembly_accession,
#   fasta_url, gtf_url, gtf_format, organism
# ============================================================


# ── Plant (EnsemblPlants) ─────────────────────────────────

rule discover_plant_urls:
    """Build per-assembly URL table for EnsemblPlants-backed species."""
    output:
        urls = f"{RESULTS}/plant_assembly_urls.tsv",
    params:
        config_key        = "plant_gtf_sources",
        db_source_override = "ensembl_plants",
    log:
        f"{LOGS}/discover_plant_urls.log",
    benchmark:
        f"{BENCHMARKS}/discover_plant_urls.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 5,
        mem_mb          = 256,
        cpus_per_task   = 1,
    script:
        "../scripts/discover_config_db_urls.py"


# ── Metazoa (wormbase + flybase) ──────────────────────────

rule discover_metazoa_urls:
    """Build per-assembly URL table for wormbase and flybase species."""
    output:
        urls = f"{RESULTS}/metazoa_assembly_urls.tsv",
    params:
        config_key        = "metazoa_gtf_sources",
        db_source_override = "",   # use entry key ("wormbase", "flybase", …)
    log:
        f"{LOGS}/discover_metazoa_urls.log",
    benchmark:
        f"{BENCHMARKS}/discover_metazoa_urls.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 5,
        mem_mb          = 256,
        cpus_per_task   = 1,
    script:
        "../scripts/discover_config_db_urls.py"


# ── Yeast (SGD) ───────────────────────────────────────────

rule discover_yeast_urls:
    """Build per-assembly URL table for SGD S. cerevisiae."""
    output:
        urls = f"{RESULTS}/yeast_assembly_urls.tsv",
    params:
        config_key        = "yeast_gtf_sources",
        db_source_override = "sgd",
    log:
        f"{LOGS}/discover_yeast_urls.log",
    benchmark:
        f"{BENCHMARKS}/discover_yeast_urls.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 5,
        mem_mb          = 256,
        cpus_per_task   = 1,
    script:
        "../scripts/discover_config_db_urls.py"
