# ============================================================
# Rule: report
# ============================================================
rule report:
    input:
        classified       = f"{RESULTS}/classified_ids.tsv",
        resolved         = f"{RESULTS}/resolved_ids.tsv",
        unresolved       = f"{RESULTS}/unresolved.tsv",
        ambiguous        = f"{RESULTS}/ambiguous.tsv",
        species_map      = f"{RESULTS}/ensembl_species_map.tsv",
        unknown_prefixes = f"{RESULTS}/ensembl_unknown_prefixes.tsv",
        fasta            = f"{RESULTS}/output.fasta",
        bed              = f"{RESULTS}/output.bed",
        benchmarks       = expand(
            f"{BENCHMARKS}/{{rule}}.tsv",
            rule=[
                "parse_ids",
                "detect_ensembl_species",
                "resolve_ids",
                "resolve_external_ids",
                "join_ensembl_results",
                "merge_resolved",
                "download_assemblies",
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
