# ============================================================
# Rule: resolve_noncode_2016
# ============================================================
# Third-tier fallback for NONCODE IDs not found in NONCODEv5 BED/FASTA
# or the NONCODEv4 BED archive.  Uses the combined NONCODE2016.fa as an
# existence check only — resolved rows receive NA genomic coordinates.
#
# NONCODE2016 is an older combined release that retains many lncRNA IDs
# (especially Drosophila melanogaster and C. elegans) that were dropped
# from the species-specific NONCODEv5 per-species FASTA files.
#
# Resolution order per ID:
#   1. Full versioned ID in NONCODE2016.fa headers (e.g. NONCELT024012.1)
#   2. Base ID (version stripped)               (e.g. NONCELT024012)
#   Gene-ID: NONCODEv5_Transcript2Gene with full then base-ID probe.
#
# Inputs:
#   noncode_v4_unresolved — noncode_v4_unresolved.tsv from resolve_noncode_v4
#   transcript2gene       — NONCODEv5_Transcript2Gene
#
# Params:
#   nc2016_fa — path to extracted NONCODE2016.fa
#
# Outputs:
#   resolved   — noncode_2016_resolved.tsv   (RESOLVED_COLS; coords NA)
#   unresolved — noncode_2016_unresolved.tsv (reason: not_found_in_any_noncode)
# ============================================================

rule resolve_noncode_2016:
    input:
        noncode_v4_unresolved = f"{RESULTS}/noncode_v4_unresolved.tsv",
        transcript2gene       = "resources/NONCODEv5_Transcript2Gene",
    params:
        nc2016_fa = "resources/NONCODE2016.fa",
    output:
        resolved   = f"{RESULTS}/noncode_2016_resolved.tsv",
        unresolved = f"{RESULTS}/noncode_2016_unresolved.tsv",
    log:
        f"{LOGS}/resolve_noncode_2016.log",
    benchmark:
        f"{BENCHMARKS}/resolve_noncode_2016.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 10,
        mem_mb          = 4096,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_noncode_2016.py"
