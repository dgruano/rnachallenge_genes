# ============================================================
# Rule: resolution_report
# Intermediate report — resolution stage only.
# Does not require downloaded assemblies or extracted sequences.
# Run with:  snakemake results/resolution_report.html
# ============================================================
rule resolution_report:
    input:
        classified       = f"{RESULTS}/classified_ids.tsv",
        resolved         = f"{RESULTS}/resolved_ids.tsv",
        unresolved       = f"{RESULTS}/unresolved.tsv",
        unmatched        = f"{RESULTS}/pattern_unmatched.tsv",
        not_found        = f"{RESULTS}/matched_not_found.tsv",
        ambiguous        = f"{RESULTS}/ambiguous.tsv",
        species_map      = f"{RESULTS}/ensembl_species_map.tsv",
        unknown_prefixes              = f"{RESULTS}/ensembl_unknown_prefixes.tsv",
        ncbi_assembly_resolved       = f"{RESULTS}/ncbi_assembly_resolved.tsv",
        ncbi_assembly_unresolved     = f"{RESULTS}/ncbi_assembly_unresolved.tsv",
        ensembl_assembly_resolved    = f"{RESULTS}/ensembl_assembly_resolved.tsv",
        ensembl_assembly_unresolved  = f"{RESULTS}/ensembl_assembly_unresolved.tsv",
        noncode_assembly_resolved    = f"{RESULTS}/noncode_assembly_resolved.tsv",
        noncode_assembly_unresolved  = f"{RESULTS}/noncode_assembly_unresolved.tsv",
        gramene_resolved           = f"{RESULTS}/gramene_resolved.tsv",
        gramene_unresolved         = f"{RESULTS}/gramene_unresolved.tsv",
        phytozome_resolved         = f"{RESULTS}/phytozome_gtf_resolved.tsv",
        phytozome_unresolved       = f"{RESULTS}/phytozome_gtf_unresolved.tsv",
        noncode_resolved           = f"{RESULTS}/noncode_resolved.tsv",
        noncode_unresolved         = f"{RESULTS}/noncode_unresolved.tsv",
        noncode_v4_resolved        = f"{RESULTS}/noncode_v4_resolved.tsv",
        noncode_v4_unresolved      = f"{RESULTS}/noncode_v4_unresolved.tsv",
        noncode_2016_resolved      = f"{RESULTS}/noncode_2016_resolved.tsv",
        noncode_2016_unresolved    = f"{RESULTS}/noncode_2016_unresolved.tsv",
        plant_gtf_resolved         = f"{RESULTS}/plant_gtf_resolved.tsv",
        plant_gtf_unresolved       = f"{RESULTS}/plant_gtf_unresolved.tsv",
        worm_gtf_resolved          = f"{RESULTS}/worm_gtf_resolved.tsv",
        worm_gtf_unresolved        = f"{RESULTS}/worm_gtf_unresolved.tsv",
        fly_gtf_resolved           = f"{RESULTS}/fly_gtf_resolved.tsv",
        fly_gtf_unresolved         = f"{RESULTS}/fly_gtf_unresolved.tsv",
        yeast_gtf_resolved         = f"{RESULTS}/yeast_gtf_resolved.tsv",
        yeast_gtf_unresolved       = f"{RESULTS}/yeast_gtf_unresolved.tsv",
        benchmarks       = expand(
            f"{BENCHMARKS}/{{rule}}.tsv",
            rule=[
                "parse_ids",
                "detect_ensembl_species",
                "resolve_ids",
                "resolve_external_ids",
                "biomart_plant_batch",
                "gramene_resolver",
                "join_ensembl_results",
                "merge_resolved",
                "resolve_ncbi_assembly_accessions",
                "resolve_ensembl_assembly_accessions",
                "resolve_noncode_assembly_accessions",
            ]
        ),
    output:
        html = f"{RESULTS}/resolution_report.html",
    log:
        f"{LOGS}/resolution_report.log",
    benchmark:
        f"{BENCHMARKS}/resolution_report.tsv",
    params:
        upstream   = config["upstream_bp"],
        downstream = config["downstream_bp"],
        release    = config["ensembl_release"],
    resources:
        slurm_partition = "compute",
        runtime         = 30,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/generate_resolution_report.py"


# ============================================================
# Rule: report
# ============================================================
rule report:
    input:
        classified       = f"{RESULTS}/classified_ids.tsv",
        resolved         = f"{RESULTS}/resolved_ids.tsv",
        unresolved       = f"{RESULTS}/unresolved.tsv",
        unmatched        = f"{RESULTS}/pattern_unmatched.tsv",
        not_found        = f"{RESULTS}/matched_not_found.tsv",
        ambiguous        = f"{RESULTS}/ambiguous.tsv",
        species_map      = f"{RESULTS}/ensembl_species_map.tsv",
        unknown_prefixes              = f"{RESULTS}/ensembl_unknown_prefixes.tsv",
        ncbi_assembly_resolved       = f"{RESULTS}/ncbi_assembly_resolved.tsv",
        ncbi_assembly_unresolved     = f"{RESULTS}/ncbi_assembly_unresolved.tsv",
        ensembl_assembly_resolved    = f"{RESULTS}/ensembl_assembly_resolved.tsv",
        ensembl_assembly_unresolved  = f"{RESULTS}/ensembl_assembly_unresolved.tsv",
        noncode_assembly_resolved    = f"{RESULTS}/noncode_assembly_resolved.tsv",
        noncode_assembly_unresolved  = f"{RESULTS}/noncode_assembly_unresolved.tsv",
        gramene_resolved           = f"{RESULTS}/gramene_resolved.tsv",
        gramene_unresolved         = f"{RESULTS}/gramene_unresolved.tsv",
        phytozome_resolved         = f"{RESULTS}/phytozome_gtf_resolved.tsv",
        phytozome_unresolved       = f"{RESULTS}/phytozome_gtf_unresolved.tsv",
        noncode_resolved           = f"{RESULTS}/noncode_resolved.tsv",
        noncode_unresolved         = f"{RESULTS}/noncode_unresolved.tsv",
        noncode_v4_resolved        = f"{RESULTS}/noncode_v4_resolved.tsv",
        noncode_v4_unresolved      = f"{RESULTS}/noncode_v4_unresolved.tsv",
        noncode_2016_resolved      = f"{RESULTS}/noncode_2016_resolved.tsv",
        noncode_2016_unresolved    = f"{RESULTS}/noncode_2016_unresolved.tsv",
        plant_gtf_resolved         = f"{RESULTS}/plant_gtf_resolved.tsv",
        plant_gtf_unresolved       = f"{RESULTS}/plant_gtf_unresolved.tsv",
        worm_gtf_resolved          = f"{RESULTS}/worm_gtf_resolved.tsv",
        worm_gtf_unresolved        = f"{RESULTS}/worm_gtf_unresolved.tsv",
        fly_gtf_resolved           = f"{RESULTS}/fly_gtf_resolved.tsv",
        fly_gtf_unresolved         = f"{RESULTS}/fly_gtf_unresolved.tsv",
        yeast_gtf_resolved         = f"{RESULTS}/yeast_gtf_resolved.tsv",
        yeast_gtf_unresolved       = f"{RESULTS}/yeast_gtf_unresolved.tsv",
        fasta            = f"{RESULTS}/output.fasta",
        bed              = f"{RESULTS}/output.bed",
        extraction_failed = f"{RESULTS}/extraction_failed.tsv",
        benchmarks       = expand(
            f"{BENCHMARKS}/{{rule}}.tsv",
            rule=[
                "parse_ids",
                "detect_ensembl_species",
                "resolve_ids",
                "resolve_external_ids",
                "biomart_plant_batch",
                "gramene_resolver",
                "join_ensembl_results",
                "merge_resolved",
                "resolve_ncbi_assembly_accessions",
                "resolve_ensembl_assembly_accessions",
                "resolve_noncode_assembly_accessions",
                "prepare_accession_list",
                "download_assemblies_done",
                "split_batches",
                "extract_sequences",
            ]
        ),
    output:
        html = f"{RESULTS}/report.html",
    log:
        f"{LOGS}/report.log",
    benchmark:
        f"{BENCHMARKS}/report.tsv",
    params:
        upstream   = config["upstream_bp"],
        downstream = config["downstream_bp"],
        release    = config["ensembl_release"],
    resources:
        slurm_partition = "compute",
        runtime         = 30,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/generate_report.py"
