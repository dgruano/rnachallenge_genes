# ============================================================
# Rule: merge_resolved
# ============================================================
# Concatenates the three resolution streams into the single
# resolved_ids.tsv that downstream rules consume, and merges
# the ambiguous records from all sources into one file.
#
# Inputs:
#   ncbi_ucsc_resolved   — from resolve_ids
#   ensembl_resolved     — from join_ensembl_results
#   external_resolved    — from resolve_external_ids (plant)
#   biomart_resolved     — from biomart_plant_batch
#   worm_gtf_resolved    — from resolve_worm_gtf
#   fly_gtf_resolved     — from resolve_fly_gtf
#   yeast_gtf_resolved   — from resolve_yeast_gtf (SGD transcript fallback)
#   gramene_resolved     — from gramene_resolver
#
# Outputs include explicit unresolved classes:
#   pattern_unmatched.tsv   — IDs with no pattern match in parse_ids
#   matched_not_found.tsv   — IDs routed to a DB but unresolved by that DB path
# ============================================================

rule merge_resolved:
    input:
        ncbi_ucsc_resolved    = f"{RESULTS}/ncbi_ucsc_resolved.tsv",
        ensembl_resolved      = f"{RESULTS}/ensembl_resolved.tsv",
        external_resolved     = f"{RESULTS}/external_resolved.tsv",
        biomart_resolved      = f"{RESULTS}/biomart_resolved.tsv",
        plant_gtf_resolved    = f"{RESULTS}/plant_gtf_resolved.tsv",
        worm_gtf_resolved     = f"{RESULTS}/worm_gtf_resolved.tsv",
        worm_gtf_unresolved   = f"{RESULTS}/worm_gtf_unresolved.tsv",
        fly_gtf_resolved      = f"{RESULTS}/fly_gtf_resolved.tsv",
        fly_gtf_unresolved    = f"{RESULTS}/fly_gtf_unresolved.tsv",
        yeast_gtf_resolved    = f"{RESULTS}/yeast_gtf_resolved.tsv",
        yeast_gtf_unresolved  = f"{RESULTS}/yeast_gtf_unresolved.tsv",
        gramene_resolved      = f"{RESULTS}/gramene_resolved.tsv",
        noncode_resolved      = f"{RESULTS}/noncode_resolved.tsv",
        noncode_v4_resolved   = f"{RESULTS}/noncode_v4_resolved.tsv",
        noncode_2016_resolved = f"{RESULTS}/noncode_2016_resolved.tsv",
        abandoned_resolved    = f"{RESULTS}/abandoned_resolved.tsv",
        ncbi_ucsc_ambiguous   = f"{RESULTS}/ncbi_ucsc_ambiguous.tsv",
        ensembl_ambiguous   = f"{RESULTS}/ensembl_ambiguous.tsv",
        external_ambiguous  = f"{RESULTS}/external_ambiguous.tsv",
        unknown_ids         = f"{RESULTS}/unknown_ids.tsv",
        ncbi_ucsc_unresolved = f"{RESULTS}/ncbi_ucsc_unresolved.tsv",
        ensembl_unresolved   = f"{RESULTS}/ensembl_unresolved.tsv",
        gramene_unresolved   = f"{RESULTS}/gramene_unresolved.tsv",
        noncode_v4_unresolved   = f"{RESULTS}/noncode_v4_unresolved.tsv",
        noncode_2016_unresolved = f"{RESULTS}/noncode_2016_unresolved.tsv",
    output:
        resolved  = f"{RESULTS}/resolved_ids.tsv",
        ambiguous = f"{RESULTS}/ambiguous.tsv",
        unresolved = f"{RESULTS}/unresolved.tsv",
        unmatched = f"{RESULTS}/pattern_unmatched.tsv",
        not_found = f"{RESULTS}/matched_not_found.tsv",
    log:
        f"{LOGS}/merge_resolved.log",
    benchmark:
        f"{BENCHMARKS}/merge_resolved.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 15,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/merge_resolved.py"
