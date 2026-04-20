# ============================================================
# Rules: download_metazoa_gtf + direct worm/fly resolvers
# ============================================================
# Direct annotation-based resolution for all WormBase and FlyBase
# transcript IDs, analogous to the direct SGD yeast resolver.
# ============================================================

_METAZOA_GTF_SOURCES = config.get("metazoa_gtf_sources", {})


rule download_metazoa_gtf:
    """Download a metazoa annotation file used for direct resolution."""
    output:
        "resources/metazoa_gtf/{source}.gtf.gz",
    params:
        url = lambda wc: _METAZOA_GTF_SOURCES[wc.source]["url"],
    log:
        f"{LOGS}/download_metazoa_gtf/{{source}}.log",
    benchmark:
        f"{BENCHMARKS}/download_metazoa_gtf/{{source}}.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 60,
        mem_mb          = 1024,
        cpus_per_task   = 1,
    shell:
        "wget --quiet --tries=5 --timeout=120 -O {output} {params.url} 2> {log}"


rule resolve_worm_gtf:
    """Resolve all WormBase IDs directly from the C. elegans annotation."""
    input:
        classified = f"{RESULTS}/classified_ids.tsv",
        gtf        = "resources/metazoa_gtf/wormbase.gtf.gz",
    output:
        resolved   = f"{RESULTS}/worm_gtf_resolved.tsv",
        unresolved = f"{RESULTS}/worm_gtf_unresolved.tsv",
    log:
        f"{LOGS}/resolve_worm_gtf.log",
    benchmark:
        f"{BENCHMARKS}/resolve_worm_gtf.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 30,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_worm_gtf.py"


rule resolve_fly_gtf:
    """Resolve all FlyBase IDs directly from the D. melanogaster annotation."""
    input:
        classified = f"{RESULTS}/classified_ids.tsv",
        gtf        = "resources/metazoa_gtf/flybase.gtf.gz",
    output:
        resolved   = f"{RESULTS}/fly_gtf_resolved.tsv",
        unresolved = f"{RESULTS}/fly_gtf_unresolved.tsv",
    log:
        f"{LOGS}/resolve_fly_gtf.log",
    benchmark:
        f"{BENCHMARKS}/resolve_fly_gtf.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 30,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_fly_gtf.py"

rule resolve_metazoa_gtf:
    input:
        worm_resolved = f"{RESULTS}/worm_gtf_resolved.tsv",
        fly_resolved  = f"{RESULTS}/fly_gtf_resolved.tsv",
