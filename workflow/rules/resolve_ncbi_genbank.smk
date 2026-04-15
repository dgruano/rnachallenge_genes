# ============================================================
# Rule: resolve_ncbi_genbank
# ============================================================
# Second-pass resolver for NCBI IDs that were NOT resolved by
# the primary esearch → elink → esummary pipeline.
#
# Uses the EPost → EFetch → GenBank parsing strategy
# (ncbi_genbank_fetcher.NCBIGenBankFetcher) which retrieves full
# GenBank records including suppressed / deleted accessions.
# A downstream esummary on the gene DB fills in genomic coords
# for any gene ID discovered from the GenBank record.
#
# Input:  results/ncbi_ucsc_unresolved.tsv  (from resolve_ids)
# Output: results/ncbi_genbank_resolved.tsv   (same schema as
#         ncbi_ucsc_resolved.tsv — feeds into merge_resolved)
#         results/ncbi_genbank_unresolved.tsv  (still unresolved)
# ============================================================

rule resolve_ncbi_genbank:
    input:
        unresolved = f"{RESULTS}/ncbi_ucsc_unresolved.tsv",
    output:
        resolved   = f"{RESULTS}/ncbi_genbank_resolved.tsv",
        unresolved = f"{RESULTS}/ncbi_genbank_unresolved.tsv",
    log:
        f"{LOGS}/resolve_ncbi_genbank.log",
    benchmark:
        f"{BENCHMARKS}/resolve_ncbi_genbank.tsv",
    threads: 1
    resources:
        slurm_partition = "compute",
        runtime         = 240,
        mem_mb          = 4096,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_ncbi_genbank.py"
