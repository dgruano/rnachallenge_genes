# ============================================================
# Rule: resolve_ids  (NCBI + UCSC only)
# ============================================================
# Resolves NCBI RefSeq and UCSC transcript IDs to gene
# coordinates via their respective APIs.
#
# Ensembl IDs are NO LONGER handled here — they are resolved
# via the BioMart wrapper sub-DAG:
#   detect_ensembl_species → biomart_lookup → join_ensembl_results
#
# Outputs named ncbi_ucsc_* to distinguish from Ensembl outputs
# before they are merged in merge_resolved.
# ============================================================

rule resolve_ids:
    input:
        classified = f"{RESULTS}/classified_ids.tsv",
    output:
        resolved  = f"{RESULTS}/ncbi_ucsc_resolved.tsv",
        ambiguous = f"{RESULTS}/ncbi_ucsc_ambiguous.tsv",
    log:
        f"{LOGS}/resolve_ids.log",
    benchmark:
        f"{BENCHMARKS}/resolve_ids.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 120,
        mem_mb          = 4096,
        cpus_per_task   = 2,
    script:
        "../scripts/resolve_ids.py"
