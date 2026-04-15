# ============================================================
# Rule: resolve_noncode_v4
# ============================================================
# Fallback resolver for NONCODE IDs that were not found in NONCODEv5.
# Uses NONCODEv4_wholeSpecies_lncAndGene_bed.zip as the annotation source.
#
# V4 BED coverage per species (assemblies differ from v5):
#   C. elegans   → ce10   (v5 also ce10)
#   Fruit fly    → dm3    (v5: dm6)
#   Zebrafish    → danRer7 (v5: danRer10)
#   Chicken      → galGal3 (v5: galGal4)
#   Cow          → bosTau6 (v5: bosTau6)
#   Human        → hg19   (v5: hg38)
#   Mouse        → mm9    (v5: mm10)
#
# Resolution order per ID:
#   1. Full versioned ID in v4 BED  (always misses; v4 has no version suffix)
#   2. Base ID (version stripped) in v4 BED
#   Gene-ID: NONCODEv5_Transcript2Gene tried with full ID then base ID.
#
# Inputs:
#   noncode_unresolved — noncode_unresolved.tsv from resolve_noncode (v5)
#   transcript2gene    — NONCODEv5_Transcript2Gene
#
# Params:
#   v4_bed_zip — path to NONCODEv4_wholeSpecies_lncAndGene_bed.zip
#
# Outputs:
#   resolved   — noncode_v4_resolved.tsv  (RESOLVED_COLS schema)
#   unresolved — noncode_v4_unresolved.tsv
# ============================================================

rule resolve_noncode_v4:
    input:
        noncode_unresolved = f"{RESULTS}/noncode_unresolved.tsv",
        transcript2gene    = "resources/NONCODEv5_Transcript2Gene",
    params:
        v4_bed_zip = "resources/NONCODEv4_wholeSpecies_lncAndGene_bed.zip",
    output:
        resolved   = f"{RESULTS}/noncode_v4_resolved.tsv",
        unresolved = f"{RESULTS}/noncode_v4_unresolved.tsv",
    log:
        f"{LOGS}/resolve_noncode_v4.log",
    benchmark:
        f"{BENCHMARKS}/resolve_noncode_v4.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 15,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_noncode_v4.py"
