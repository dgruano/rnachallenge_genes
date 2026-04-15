# ============================================================
# Rule: parse_ids
# ============================================================
# Stage 1 — Species-First Transcript Routing
#
# Reads input FASTA file(s), extracts transcript IDs from headers,
# and classifies each ID into:
#   - species_hint
#   - source_hint
#   - database route (db_source)
# using regex pattern matching on transcript names.
#
# Outputs:
#   - classified_ids.tsv: transcript_id, db_source, species_hint, source_hint,
#                         assembly_hint, raw_header, source_file
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
