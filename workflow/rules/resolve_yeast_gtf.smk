# ============================================================
# Rules: download_yeast_gtf  +  resolve_yeast_gtf
# ============================================================
# GFF3-based resolver for all SGD (yeast) transcript IDs.
#
# The GFF3 annotation is downloaded once from yeastgenome.org.
# Version and URL are configured in config/yeast_gtf_sources.yaml
# (loaded as a top-level configfile in the Snakefile).
#
# Chain:
#   parse_ids             →  classified_ids.tsv
#                                    ↓
#   resolve_yeast_gtf     →  yeast_gtf_resolved.tsv
#                             yeast_gtf_unresolved.tsv
# ============================================================

_YEAST_GTF_CFG = config.get("yeast_gtf_sources", {}).get(
    "saccharomyces_cerevisiae", {}
)
_YEAST_GFF = "resources/yeast_gtf/saccharomyces_cerevisiae.gff.gz"


# ── Download ─────────────────────────────────────────────────

rule download_yeast_gtf:
    """Download the SGD S. cerevisiae GFF3 annotation."""
    output:
        protected(_YEAST_GFF),
    params:
        url = _YEAST_GTF_CFG.get("url", ""),
    log:
        f"{LOGS}/download_yeast_gtf.log",
    benchmark:
        f"{BENCHMARKS}/download_yeast_gtf.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 60,
        mem_mb          = 512,
        cpus_per_task   = 1,
    shell:
        "wget --quiet --tries=5 --timeout=120 -O {output} {params.url} 2> {log}"


# ── Resolve ──────────────────────────────────────────────────

rule resolve_yeast_gtf:
    """Resolve all SGD IDs directly from the SGD GFF3."""
    input:
        classified = f"{RESULTS}/classified_ids.tsv",
        gff        = _YEAST_GFF,
    output:
        resolved   = f"{RESULTS}/yeast_gtf_resolved.tsv",
        unresolved = f"{RESULTS}/yeast_gtf_unresolved.tsv",
    log:
        f"{LOGS}/resolve_yeast_gtf.log",
    benchmark:
        f"{BENCHMARKS}/resolve_yeast_gtf.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 30,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_yeast_gtf.py"
