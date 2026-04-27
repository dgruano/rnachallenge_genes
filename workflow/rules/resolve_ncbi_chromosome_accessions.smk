# ============================================================
# Rule: resolve_ncbi_chromosome_accessions
# ============================================================
# Post-merge resolver — maps NC_/NT_/NW_ chromosomal accessions
# in resolved_ids.tsv to their parent GCF_/GCA_ assemblies via
# NCBI elink (batch-efficient: no full record downloads).
#
# Sits between merge_resolved and download_assemblies so that
# download_assemblies always receives GCF_/GCA_-level accessions.
#
# Input:  results/resolved_ids.tsv
# Output: results/ncbi_chromosome_resolved.tsv   (patched full TSV)
#         results/ncbi_chromosome_unresolved.tsv  (unmappable rows)
# ============================================================

rule resolve_ncbi_chromosome_accessions:
    input:
        resolved = f"{RESULTS}/resolved_ids.tsv",
    output:
        resolved   = f"{RESULTS}/ncbi_chromosome_resolved.tsv",
        unresolved = f"{RESULTS}/ncbi_chromosome_unresolved.tsv",
    log:
        f"{LOGS}/resolve_ncbi_chromosome_accessions.log",
    benchmark:
        f"{BENCHMARKS}/resolve_ncbi_chromosome_accessions.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 120,   # 2 hours for NCBI API calls
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_ncbi_chromosome_accessions.py"
