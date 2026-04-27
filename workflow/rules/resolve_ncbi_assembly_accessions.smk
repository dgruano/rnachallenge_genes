# ============================================================
# Rule: resolve_ncbi_assembly_accessions
# ============================================================
# Stage 2c — NCBI Assembly Accession Resolution
#
# For NCBI transcripts with NC_/NW_ sequence accessions but
# no parent assembly accession (GCF_/GCA_), resolve them to
# their parent assemblies and download the GTF.
#
# Input: resolved NCBI transcripts from Stage 2a
#        (ncbi_ucsc_resolved.tsv)
#
# Processing:
#   1. For each NC_/NW_ accession, fetch its parent GCF_/GCA_
#   2. Batch-download assembly GTF files
#   3. Extract coordinates from GTF (enriching empty chrom fields)
#   4. Update assembly_accession from NC_ to GCF_
#
# Output: enriched NCBI transcript table with complete assembly info
#
# Cache layout:
#   resources/cache/
#     <assembly_accession>/
#       genomic.gtf.gz
# ============================================================

checkpoint resolve_ncbi_assembly_accessions:
    input:
        resolved = f"{RESULTS}/ncbi_ucsc_resolved.tsv",
    output:
        resolved   = f"{RESULTS}/ncbi_assembly_resolved.tsv",
        unresolved = f"{RESULTS}/ncbi_assembly_unresolved.tsv",
        ambiguous  = f"{RESULTS}/ncbi_assembly_ambiguous.tsv",
    log:
        f"{LOGS}/resolve_ncbi_assembly_accessions.log",
    benchmark:
        f"{BENCHMARKS}/resolve_ncbi_assembly_accessions.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 240,  # 4 hours for NCBI API calls
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_ncbi_assembly_accessions.py"
