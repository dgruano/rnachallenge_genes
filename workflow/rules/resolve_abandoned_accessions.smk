# ============================================================
# Rule: resolve_abandoned_accessions
# ============================================================
# Third-pass resolver for NCBI transcript accessions that are
# withdrawn or suppressed (i.e. remain unresolved after the
# resolve_ncbi_genbank step).
#
# Strategy
# --------
# 1. esummary(nuccore) → recover AssemblyAcc from metadata.
#    NCBI retains this field even after the sequence record has
#    been suppressed or withdrawn.
# 2. elink(nuccore→assembly) fallback for accessions that still
#    have no AssemblyAcc.
# 3. Download assembly GTF from NCBI FTP; cache under
#    resources/cache/<assembly_accession>/genomic.gtf.gz.
# 4. Extract the gene and transcript annotation from the GTF.
# 5. Write resolved / unresolved TSVs.
#
# Input:   results/ncbi_genbank_unresolved.tsv   (resolve_ncbi_genbank)
# Output:  results/abandoned_resolved.tsv        (same schema as
#              ncbi_genbank_resolved.tsv → feeds merge_resolved)
#          results/abandoned_unresolved.tsv      (remaining failures)
# ============================================================

rule resolve_abandoned_accessions:
    input:
        unresolved = f"{RESULTS}/ncbi_genbank_unresolved.tsv",
    output:
        resolved   = f"{RESULTS}/abandoned_resolved.tsv",
        unresolved = f"{RESULTS}/abandoned_unresolved.tsv",
    log:
        f"{LOGS}/resolve_abandoned_accessions.log",
    benchmark:
        f"{BENCHMARKS}/resolve_abandoned_accessions.tsv",
    threads: 1
    resources:
        slurm_partition = "compute",
        runtime         = 480,   # 8 h — GTF downloads can be slow
        mem_mb          = 8192,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_abandoned_accessions.py"
