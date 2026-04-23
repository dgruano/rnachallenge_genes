# ============================================================
# Rule: merge_resolved
# ============================================================
# Concatenates the three resolution streams into the single
# resolved_ids.tsv that downstream rules consume, and merges
# the ambiguous records from all sources into one file.
#
# Inputs:
#   ncbi_assembly_resolved   — from resolve_ncbi_assembly_accessions (enriched ncbi_ucsc)
#   ensembl_assembly_resolved — from resolve_ensembl_assembly_accessions
#   noncode_assembly_resolved — from resolve_noncode_assembly_accessions
#   external_resolved    — from resolve_external_ids (plant)
#   biomart_resolved     — from biomart_plant_batch
#   worm_gtf_resolved    — from resolve_worm_gtf
#   fly_gtf_resolved     — from resolve_fly_gtf
#   yeast_gtf_resolved   — from resolve_yeast_gtf (SGD transcript fallback)
#   phytozome_gtf_resolved — from resolve_phytozome_gtf
#   gramene_resolved     — from gramene_resolver
#
# Outputs include explicit unresolved classes:
#   pattern_unmatched.tsv   — IDs with no pattern match in parse_ids
#   matched_not_found.tsv   — IDs routed to a DB but unresolved by that DB path
# ============================================================

rule merge_resolved:
    input:
        ncbi_assembly_resolved = f"{RESULTS}/ncbi_assembly_resolved.tsv",
        ncbi_assembly_unresolved = f"{RESULTS}/ncbi_assembly_unresolved.tsv",
        ensembl_assembly_resolved = f"{RESULTS}/ensembl_assembly_resolved.tsv",
        external_resolved     = f"{RESULTS}/external_resolved.tsv",
        biomart_resolved      = f"{RESULTS}/biomart_resolved.tsv",
        plant_gtf_resolved    = f"{RESULTS}/plant_gtf_resolved.tsv",
        phytozome_gtf_resolved = f"{RESULTS}/phytozome_gtf_resolved.tsv",
        phytozome_gtf_unresolved = f"{RESULTS}/phytozome_gtf_unresolved.tsv",
        worm_gtf_resolved     = f"{RESULTS}/worm_gtf_resolved.tsv",
        worm_gtf_unresolved   = f"{RESULTS}/worm_gtf_unresolved.tsv",
        fly_gtf_resolved      = f"{RESULTS}/fly_gtf_resolved.tsv",
        fly_gtf_unresolved    = f"{RESULTS}/fly_gtf_unresolved.tsv",
        yeast_gtf_resolved    = f"{RESULTS}/yeast_gtf_resolved.tsv",
        yeast_gtf_unresolved  = f"{RESULTS}/yeast_gtf_unresolved.tsv",
        gramene_resolved      = f"{RESULTS}/gramene_resolved.tsv",
        noncode_assembly_resolved = f"{RESULTS}/noncode_assembly_resolved.tsv",
        noncode_assembly_unresolved = f"{RESULTS}/noncode_assembly_unresolved.tsv",
        noncode_v4_resolved   = f"{RESULTS}/noncode_v4_resolved.tsv",
        noncode_2016_resolved = f"{RESULTS}/noncode_2016_resolved.tsv",
        abandoned_resolved    = f"{RESULTS}/abandoned_resolved.tsv",
        ncbi_assembly_ambiguous = f"{RESULTS}/ncbi_assembly_ambiguous.tsv",
        ensembl_ambiguous   = f"{RESULTS}/ensembl_ambiguous.tsv",
        external_ambiguous  = f"{RESULTS}/external_ambiguous.tsv",
        unknown_ids         = f"{RESULTS}/unknown_ids.tsv",
        ensembl_assembly_unresolved = f"{RESULTS}/ensembl_assembly_unresolved.tsv",
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
