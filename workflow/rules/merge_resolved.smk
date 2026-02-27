# ============================================================
# Rule: merge_resolved
# ============================================================
# Concatenates the three resolution streams into the single
# resolved_ids.tsv that downstream rules consume, and merges
# the ambiguous records from all sources into one file.
#
# Inputs:
#   ncbi_ucsc_resolved  — from resolve_ids (NCBI + UCSC path)
#   ensembl_resolved    — from join_ensembl_results
#   ncbi_ucsc_ambiguous — from resolve_ids
#   ensembl_ambiguous   — from join_ensembl_results
# ============================================================

rule merge_resolved:
    input:
        ncbi_ucsc_resolved  = f"{RESULTS}/ncbi_ucsc_resolved.tsv",
        ensembl_resolved    = f"{RESULTS}/ensembl_resolved.tsv",
        external_resolved   = f"{RESULTS}/external_resolved.tsv",
        ncbi_ucsc_ambiguous = f"{RESULTS}/ncbi_ucsc_ambiguous.tsv",
        ensembl_ambiguous   = f"{RESULTS}/ensembl_ambiguous.tsv",
        external_ambiguous  = f"{RESULTS}/external_ambiguous.tsv",
        unknown_ids         = f"{RESULTS}/unknown_ids.tsv",
        external_unresolved = f"{RESULTS}/external_unresolved.tsv",
    output:
        resolved  = f"{RESULTS}/resolved_ids.tsv",
        ambiguous = f"{RESULTS}/ambiguous.tsv",
        unresolved = f"{RESULTS}/unresolved.tsv",
    log:
        f"{LOGS}/merge_resolved.log",
    benchmark:
        f"{BENCHMARKS}/merge_resolved.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 15,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/merge_resolved.py"
