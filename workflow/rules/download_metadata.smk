# ============================================================
# Rule: download_metadata_table
# ============================================================
# Downloads transcript metadata tables from Ensembl Plants BioMart
# for fast local ID resolution (avoids REST API rate limits).
# ============================================================

rule download_metadata_table:
    output:
        "resources/metadata/{species}.tsv"
    log:
        "logs/download_metadata_{species}.log"
    benchmark:
        "benchmarks/download_metadata_{species}.tsv"
    resources:
        slurm_partition = "compute",
        runtime         = 30,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/download_metadata_table.py"


# ============================================================
# Rule: download_all_metadata_tables
# ============================================================
# Downloads all configured metadata tables
# ============================================================

rule download_all_metadata_tables:
    input:
        expand(
            "resources/metadata/{species}.tsv",
            species=config.get("external_metadata_tables", {}).keys()
        )
    output:
        touch("resources/metadata/.download_complete")
    log:
        "logs/download_all_metadata_tables.log"
