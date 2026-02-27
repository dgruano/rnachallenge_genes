# ============================================================
# Rule: biomart_lookup
# ============================================================
# Invokes the official snakemake-wrappers BioMart wrapper once
# per species detected in the input. Each invocation downloads
# the full annotation table for that species (all transcripts),
# which is then filtered to our IDs in join_ensembl_results.
#
# The wildcard {species} is populated dynamically from the
# detect_ensembl_species checkpoint output.
#
# Output is cached by Snakemake (cache: "omit-software") so
# re-runs with the same release don't re-download.
# ============================================================

def _biomart_build(wildcards):
    """Look up the genome build for a given species from config."""
    for prefix, info in config["ensembl_species"].items():
        if info["species"] == wildcards.species:
            return info["build"]
    raise ValueError(f"No build configured for species: {wildcards.species}")


rule biomart_lookup:
    output:
        table = f"{RESULTS}/biomart/{{species}}.tsv.gz",
    log:
        f"{LOGS}/biomart_lookup/{{species}}.log",
    params:
        biomart    = "genes",
        species    = lambda wc: wc.species,
        build      = _biomart_build,
        release    = config["ensembl_release"],
        attributes = [
            "ensembl_transcript_id",
            "ensembl_transcript_id_version",
            "ensembl_gene_id",
            "external_gene_name",
            "chromosome_name",
            "start_position",
            "end_position",
            "strand",
        ],
        # No filters — we download the full species table and filter locally.
        # This maximises Snakemake caching value across different input sets.
    cache:
        "omit-software"    # Reuse across runs with same release; invalidated if release changes
    resources:
        slurm_partition = "compute",
        runtime         = 60,
        mem_mb          = 8192,
        cpus_per_task   = 2,
    wrapper:
        "v5.1.0/bio/reference/ensembl-biomart-table"
