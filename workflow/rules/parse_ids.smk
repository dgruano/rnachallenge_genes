# ============================================================
# Rule: parse_ids
# ============================================================
# Stage 1 — Parse & Classify Transcript IDs
#
# Reads input FASTA file(s), extracts transcript IDs from headers,
# and classifies each ID into one of: ncbi | ensembl | ucsc | unknown.
#
# Outputs:
#   - classified_ids.tsv: transcript_id, db_source, raw_header, source_file
#   - unknown_ids.tsv:    transcript_id, raw_header, source_file, reason
# ============================================================

rule parse_ids:
    input:
        fastas = config["input_fastas"],
    output:
        classified = f"{RESULTS}/classified_ids.tsv",
        unknown    = f"{RESULTS}/unknown_ids.tsv",
    log:
        f"{LOGS}/parse_ids.log",
    benchmark:
        f"{BENCHMARKS}/parse_ids.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 30,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/parse_ids.py"
