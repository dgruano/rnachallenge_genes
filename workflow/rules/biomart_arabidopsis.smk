configfile: "config.yaml"
RESULTS = "results"
LOGS    = "logs"

# ============================================================
# Rule: list_biomart_plants_dbs
# ============================================================
# Diagnostic rule: lists available BioMart databases on the
# Ensembl Plants host and on the default vertebrate host.
# Run with: snakemake list_biomart_plants_dbs --profile profiles/default
# ============================================================

rule list_biomart_plants_dbs:
    output:
        f"{RESULTS}/biomart_plants_available_dbs.txt",
    log:
        f"{LOGS}/list_biomart_plants_dbs.log",
    resources:
        slurm_partition = "compute",
        runtime         = 10,
        mem_mb          = 1024,
        cpus_per_task   = 1,
    script:
        "../scripts/list_biomart_plants_dbs.R"


# ============================================================
# Rule: biomart_arabidopsis_table
# ============================================================
# Downloads Arabidopsis thaliana transcript→gene annotation
# from Ensembl Plants BioMart using the official Snakemake
# wrapper (v6.0.2).
#
# Output is cached by Snakemake (cache: "omit-software") so
# re-runs with the same release don't re-download.
# ============================================================

rule biomart_arabidopsis_table:
    output:
        table = f"{RESULTS}/biomart/arabidopsis_thaliana.tsv.gz",
    log:
        f"{LOGS}/biomart_arabidopsis_table.log",
    benchmark:
        f"{BENCHMARKS}/biomart_arabidopsis_table.tsv",
    params:
        release = config["ensembl_release"],
    cache: "omit-software"
    resources:
        slurm_partition = "compute",
        runtime         = 30,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/biomart_arabidopsis_download.py"
