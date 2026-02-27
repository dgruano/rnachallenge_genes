# ============================================================
# Rule: detect_ensembl_species  (checkpoint)
# ============================================================
# Extracts unique Ensembl transcript ID prefixes from the input
# (equivalent to: grep -oE 'ENS[A-Z]*T' | sort | uniq), matches
# each against a built-in species reference table, and fans out
# one biomart_lookup job per detected species.
#
# If any prefix is unrecognised, the pipeline stops with a clear
# actionable error message and instructions for config.yaml.
# ============================================================

checkpoint detect_ensembl_species:
    input:
        classified = f"{RESULTS}/classified_ids.tsv",
    output:
        species_map = f"{RESULTS}/ensembl_species_map.tsv",
        unmatched   = f"{RESULTS}/ensembl_unknown_prefixes.tsv",
    log:
        f"{LOGS}/detect_ensembl_species.log",
    benchmark:
        f"{BENCHMARKS}/detect_ensembl_species.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 15,
        mem_mb          = 1024,
        cpus_per_task   = 1,
    script:
        "../scripts/detect_ensembl_species.py"
