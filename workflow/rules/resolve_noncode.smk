# ============================================================
# Rule: resolve_noncode
# ============================================================
# Resolves NONCODE transcript and gene IDs (db_source == "noncode")
# classified by parse_ids into genomic coordinates and gene-level
# annotations using a cascade fallback:
#
#   1. NONCODEv5 BED + Transcript2Gene  (primary)
#   2. NONCODEv4 whole-species BED zip  (v4 coordinates; IDs without version)
#   3. NONCODE2016 FASTA                (existence only; NA coordinates)
#
# Inputs:
#   classified       — classified_ids.tsv from parse_ids
#   transcript2gene  — NONCODEv5_Transcript2Gene (2-col, space-delimited)
#
# Params:
#   bed_dir    — directory holding NONCODEv5_*.lncAndGene.bed.gz files
#   fa_dir     — directory holding NONCODEv5_*.fa.gz files
#   v4_bed_zip — path to NONCODEv4_wholeSpecies_lncAndGene_bed.zip
#   nc2016_fa  — path to extracted NONCODE2016.fa
#
# Outputs:
#   resolved   — noncode_resolved.tsv (RESOLVED_COLS schema)
#   unresolved — noncode_unresolved.tsv
# ============================================================

rule resolve_noncode:
    input:
        classified      = f"{RESULTS}/classified_ids.tsv",
        transcript2gene = "resources/NONCODEv5_Transcript2Gene",
    params:
        bed_dir = "resources/noncode_beds",
        fa_dir  = "resources/noncode_fa",
    output:
        resolved   = f"{RESULTS}/noncode_resolved.tsv",
        unresolved = f"{RESULTS}/noncode_unresolved.tsv",
    log:
        f"{LOGS}/resolve_noncode.log",
    benchmark:
        f"{BENCHMARKS}/resolve_noncode.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 30,
        mem_mb          = 4096,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_noncode.py"
