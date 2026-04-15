# ============================================================
# Rule: biomart_lookup
# ============================================================
# Downloads full annotation table for each species via Ensembl
# BioMart REST API using a custom Python script.
#
# The wildcard {species} is populated dynamically from the
# detect_ensembl_species checkpoint output.
#
# Output is cached by Snakemake (cache: "omit-software") so
# re-runs with the same release don't re-download.
#
# Note: The previous R-wrapper implementation had variable name
# bugs in the R code. This custom Python implementation:
#   - Directly queries the Ensembl BioMart REST API
#   - Has proper error handling and retry logic
#   - Provides detailed logging
#   - Uses less memory and CPU resources
# ============================================================


def _biomart_input(wildcards):
    """Declare checkpoint output as input to wait for checkpoint completion."""
    return checkpoints.detect_ensembl_species.get().output.species_map


rule biomart_lookup:
    input:
        _biomart_input,
    output:
        table = f"{RESULTS}/biomart/{{species}}.tsv.gz",
    log:
        f"{LOGS}/biomart_lookup/{{species}}.log",
    params:
        species = lambda wc: wc.species,
        release = config["ensembl_release"],
    cache:
        "omit-software"    # Reuse across runs with same release; invalidated if release changes
    resources:
        slurm_partition = "compute",
        runtime         = 120,    # Slightly increased for network variability
        mem_mb          = 4096,
        cpus_per_task   = 1,
    script:
        "../scripts/biomart_download.py"
