# ============================================================
# Rule: resolve_ensembl_assembly_accessions
# ============================================================
# Stage 2d — Ensembl Assembly Accession Resolution
#
# Enriches Ensembl-resolved rows by mapping assembly build names
# (e.g. GRCh38, GRCz11) to downloadable NCBI assembly accessions
# (GCF_/GCA_), then optionally fills missing coordinates from
# cached assembly GTF files.
#
# Input:  ensembl_resolved.tsv + ensembl_unresolved.tsv
# Output: ensembl_assembly_resolved.tsv + ensembl_assembly_unresolved.tsv
# ============================================================

checkpoint resolve_ensembl_assembly_accessions:
    input:
        resolved = f"{RESULTS}/ensembl_resolved.tsv",
    output:
        resolved = f"{RESULTS}/ensembl_assembly_resolved.tsv",
        unresolved = f"{RESULTS}/ensembl_assembly_unresolved.tsv",
    log:
        f"{LOGS}/resolve_ensembl_assembly_accessions.log",
    benchmark:
        f"{BENCHMARKS}/resolve_ensembl_assembly_accessions.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 180,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_ensembl_assembly_accessions.py"