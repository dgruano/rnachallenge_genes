# ============================================================
# Rule: download_id_mappings
# ============================================================
# Download legacy ID mapping files (Rice MSU, Maize GRMZM)
# ============================================================

rule download_id_mappings:
    output:
        touch("resources/id_mappings/.download_complete"),
    log:
        "logs/download_id_mappings.log",
    benchmark:
        "benchmarks/download_id_mappings.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 15,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/download_id_mappings.py"
