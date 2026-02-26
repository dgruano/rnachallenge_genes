# ============================================================
# Rule: detect_ensembl_species  (checkpoint)
# ============================================================
# Reads classified_ids.tsv, finds all Ensembl transcript IDs,
# and maps each to a species using the prefix→species table in
# config. Writes:
#   results/ensembl_species_map.tsv  — transcript_id | prefix | species | build
#   results/ensembl_unmatched.tsv   — transcript_ids whose prefix is not in config
#
# Declared as a checkpoint so that downstream rules (biomart_lookup,
# join_ensembl_results) can fan out dynamically over the detected species.
# ============================================================

checkpoint detect_ensembl_species:
    input:
        classified = f"{RESULTS}/classified_ids.tsv",
    output:
        species_map = f"{RESULTS}/ensembl_species_map.tsv",
        unmatched   = f"{RESULTS}/ensembl_unmatched_prefix.tsv",
    log:
        f"{LOGS}/detect_ensembl_species.log",
    benchmark:
        f"{BENCHMARKS}/detect_ensembl_species.tsv",
    resources:
        slurm_partition = "short",
        runtime         = 15,
        mem_mb          = 1024,
        cpus_per_task   = 1,
    script:
        "../scripts/detect_ensembl_species.py"
