# ============================================================
# Rule: gramene_resolver
# ============================================================
# Resolve legacy plant IDs using Gramene search API
# ============================================================

rule gramene_resolver:
    input:
        f"{RESULTS}/plant_gtf_unresolved.tsv",
    output:
        resolved   = "results/gramene_resolved.tsv",
        unresolved = "results/gramene_unresolved.tsv",
    log:
        "logs/gramene_resolver.log",
    benchmark:
        "benchmarks/gramene_resolver.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 120,  # 2 hours for rate-limited API calls
        mem_mb          = 4096,
        cpus_per_task   = 1,
    script:
        "../scripts/gramene_resolver.py"
