# ============================================================
# Rule: join_ensembl_results
# ============================================================
# After all per-species BioMart tables are downloaded, this rule:
#   1. Reads all BioMart TSVs
#   2. Concatenates them into one lookup table
#   3. Left-joins our Ensembl transcript IDs against it
#   4. Normalises columns to the pipeline's unified schema
#   5. Flags and records any ambiguous transcript→gene mappings
#
# Output feeds directly into merge_resolved (alongside NCBI/UCSC).
# ============================================================

def _all_biomart_tables(wildcards):
    """
    Collect all per-species BioMart output files by reading the
    detect_ensembl_species checkpoint output.
    """
    checkpoints.detect_ensembl_species.get()
    import pandas as pd
    species_map = pd.read_csv(
        checkpoints.detect_ensembl_species.get().output.species_map,
        sep="\t",
    )
    species_list = species_map["species"].unique().tolist()
    return expand(f"{RESULTS}/biomart/{{species}}.tsv.gz", species=species_list)


rule join_ensembl_results:
    input:
        classified   = f"{RESULTS}/classified_ids.tsv",
        species_map  = f"{RESULTS}/ensembl_species_map.tsv",
        biomart_tables = _all_biomart_tables,
    output:
        resolved  = f"{RESULTS}/ensembl_resolved.tsv",
        ambiguous = f"{RESULTS}/ensembl_ambiguous.tsv",
        unresolved = f"{RESULTS}/ensembl_unresolved.tsv",
    log:
        f"{LOGS}/join_ensembl_results.log",
    benchmark:
        f"{BENCHMARKS}/join_ensembl_results.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 30,
        mem_mb          = 8192,
        cpus_per_task   = 2,
    script:
        "../scripts/join_ensembl_results.py"
