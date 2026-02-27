# ============================================================
# Rule: resolve_external_ids
# ============================================================
# Resolves non-NCBI/Ensembl/UCSC transcript IDs currently
# classified as "plant" or "wormbase". Uses Ensembl REST for
# plant IDs and parses WormBase-encoded headers when available.
# ============================================================

rule resolve_external_ids:
    input:
        classified = f"{RESULTS}/classified_ids.tsv",
    output:
        resolved   = f"{RESULTS}/external_resolved.tsv",
        ambiguous  = f"{RESULTS}/external_ambiguous.tsv",
        unresolved = f"{RESULTS}/external_unresolved.tsv",
    log:
        f"{LOGS}/resolve_external_ids.log",
    benchmark:
        f"{BENCHMARKS}/resolve_external_ids.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 60,
        mem_mb          = 4096,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_external_ids.py"
