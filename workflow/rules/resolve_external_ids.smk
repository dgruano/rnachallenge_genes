# ============================================================
# Rule: resolve_external_ids
# ============================================================
# Resolves non-NCBI/Ensembl/UCSC transcript IDs currently
# classified as "plant" or "wormbase". Uses Ensembl REST for
# plant IDs and parses WormBase-encoded headers when available.
# Can use local metadata tables for faster resolution.
# ============================================================

def get_external_unresolved_input(wildcards):
    """Get unresolved input file if specified and avoid cyclic dependency."""
    unresolved_input = config.get("external_unresolved_input", "")
    if unresolved_input and unresolved_input != f"{RESULTS}/external_unresolved.tsv":
        return unresolved_input
    return []

def get_metadata_tables_input(wildcards):
    """Get metadata table files if configured."""
    tables = config.get("external_metadata_tables", {})
    if tables:
        return [f"resources/metadata/{species}.tsv" for species in tables.keys()]
    return []

rule resolve_external_ids:
    input:
        classified = f"{RESULTS}/classified_ids.tsv",
        unresolved = get_external_unresolved_input,
        metadata   = get_metadata_tables_input,
    output:
        resolved   = f"{RESULTS}/external_resolved.tsv",
        ambiguous  = f"{RESULTS}/external_ambiguous.tsv",
        unresolved = f"{RESULTS}/external_unresolved.tsv",
    log:
        f"{LOGS}/resolve_external_ids.log",
    benchmark:
        f"{BENCHMARKS}/resolve_external_ids.tsv",
    threads: 2
    resources:
        slurm_partition = "compute",
        runtime         = 60,
        mem_mb          = 4096,
        cpus_per_task   = 2,
    script:
        "../scripts/resolve_external_ids.py"
