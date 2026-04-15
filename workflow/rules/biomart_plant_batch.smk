# ============================================================
# Rule: biomart_plant_batch
# ============================================================
# Batch query Ensembl Plants BioMart for unresolved plant IDs.
# Much more efficient than REST API for bulk queries.
# ============================================================

rule biomart_plant_batch:
    input:
        unresolved = f"{RESULTS}/external_unresolved.tsv",
    output:
        resolved   = f"{RESULTS}/biomart_resolved.tsv",
        unresolved = f"{RESULTS}/biomart_unresolved.tsv",
    log:
        f"{LOGS}/biomart_plant_batch.log",
    benchmark:
        f"{BENCHMARKS}/biomart_plant_batch.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 30,
        mem_mb          = 4096,
        cpus_per_task   = 1,
    script:
        "../scripts/biomart_plant_batch.py"
