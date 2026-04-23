# ============================================================
# Rule: resolve_noncode_assembly_accessions
# ============================================================
# Stage 3 — NONCODE Assembly Accession Resolution
#
# Enriches NONCODE-resolved rows by mapping UCSC genome build
# names (e.g. tair10, ce10, dm6) to downloadable NCBI assembly
# accessions (GCF_/GCA_).
#
# Input:  noncode_resolved.tsv
# Output: noncode_assembly_resolved.tsv + noncode_assembly_unresolved.tsv
# ============================================================

checkpoint resolve_noncode_assembly_accessions:
    input:
        resolved = f"{RESULTS}/noncode_resolved.tsv",
    output:
        resolved = f"{RESULTS}/noncode_assembly_resolved.tsv",
        unresolved = f"{RESULTS}/noncode_assembly_unresolved.tsv",
    log:
        f"{LOGS}/resolve_noncode_assembly_accessions.log",
    benchmark:
        f"{BENCHMARKS}/resolve_noncode_assembly_accessions.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 180,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_noncode_assembly_accessions.py"
