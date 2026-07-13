# ============================================================
# Rule: annotate_with_tools
# ============================================================
# Add tool-source annotations to unresolved reports for triage.
# Supported bases:
#   - matched_not_found
#   - pattern_unmatched
# ============================================================

rule annotate_with_tools:
    input:
        records = f"{RESULTS}/{{base}}.tsv",
        tool_map = f"{RESULTS}/tool_source_map.tsv",
    output:
        annotated = f"{RESULTS}/{{base}}_with_tools.tsv",
    wildcard_constraints:
        base = "matched_not_found|pattern_unmatched",
    log:
        f"{LOGS}/annotate_with_tools_{{base}}.log",
    benchmark:
        f"{BENCHMARKS}/annotate_with_tools_{{base}}.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 10,
        mem_mb          = 1024,
        cpus_per_task   = 1,
    script:
        "../scripts/annotate_with_tools.py"
