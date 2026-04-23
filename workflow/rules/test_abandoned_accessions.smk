# ============================================================
# Rule: test_abandoned_accessions
# ============================================================
# Test harness for resolve_abandoned_accessions optimization.
# Creates a 100-transcript subset from ncbi_genbank_unresolved.tsv,
# runs the resolver, and collects metrics.

rule subset_unresolved_for_test:
    """Extract first 100 unresolved NCBI accessions for testing."""
    input:
        f"{RESULTS}/ncbi_genbank_unresolved.tsv",
    output:
        temp(f"{RESULTS}/.test_abandoned_100.tsv"),
    shell:
        "head -101 {input} > {output}"


rule test_abandoned_accessions:
    """Test resolve_abandoned_accessions with 100-transcript subset."""
    input:
        unresolved = rules.subset_unresolved_for_test.output,
    output:
        resolved   = f"{RESULTS}/test_abandoned_resolved.tsv",
        unresolved = f"{RESULTS}/test_abandoned_unresolved.tsv",
    log:
        f"{LOGS}/test_abandoned_accessions.log",
    benchmark:
        f"{BENCHMARKS}/test_abandoned_accessions.tsv",
    threads: 1
    resources:
        slurm_partition = "compute",
        runtime         = 120,  # 2 hours max for test
        mem_mb          = 8192,
        cpus_per_task   = 1,
    script:
        "../scripts/resolve_abandoned_accessions.py"


rule test_abandoned_report:
    """Summarize test results: runtime, resolution rate, API efficiency."""
    input:
        resolved   = rules.test_abandoned_accessions.output.resolved,
        unresolved = rules.test_abandoned_accessions.output.unresolved,
        benchmark  = rules.test_abandoned_accessions.benchmark,
        log        = rules.test_abandoned_accessions.log,
    output:
        report = f"{RESULTS}/test_abandoned_report.txt",
    run:
        import pandas as pd
        import json

        with open(input.resolved) as f_res, \
             open(input.unresolved) as f_unres, \
             open(input.benchmark) as f_bench, \
             open(input.log[0]) as f_log:

            df_res = pd.read_csv(f_res, sep="\t")
            df_unres = pd.read_csv(f_unres, sep="\t")

            bench_data = pd.read_csv(f_bench, sep="\t").iloc[0].to_dict()
            log_content = f_log.read()

        total_input = len(df_res) + len(df_unres)
        resolution_rate = len(df_res) / total_input * 100 if total_input > 0 else 0

        report_lines = [
            "=" * 70,
            "TEST RESULTS: resolve_abandoned_accessions (100-transcript subset)",
            "=" * 70,
            "",
            f"Input transcripts:        {total_input}",
            f"Resolved:                 {len(df_res)} ({resolution_rate:.1f}%)",
            f"Unresolved:               {len(df_unres)}",
            "",
            f"Runtime (s):              {bench_data.get('s', 'N/A')}",
            f"CPU time (s):             {bench_data.get('cpu_seconds', 'N/A')}",
            f"Max memory (MB):          {bench_data.get('max_rss', 'N/A')}",
            "",
            "Unresolved breakdown (if any):",
        ]

        if not df_unres.empty:
            for reason, grp in df_unres.groupby("reason"):
                report_lines.append(f"  {reason:<30}: {len(grp)}")
        else:
            report_lines.append("  (none)")

        report_lines.extend([
            "",
            "Key metrics (from log):",
            "  - Batched efetch calls vs sequential (see Step 1/2 counts)",
            "  - Rate-limit delays reduced from 0.11s to 0.02s",
            "  - GTF download + parsing time (Step 4/5)",
            "",
            "=" * 70,
        ])

        report_text = "\n".join(report_lines)

        with open(output.report, "w") as f:
            f.write(report_text)
