"""
scripts/generate_report.py
Stage 5 — Generate Summary Report
===================================
Reads all pipeline outputs and produces:
  - results/report.html   : human-readable HTML summary
Covers:
  - Total transcripts, breakdown by DB
  - Resolution success/failure/ambiguity
  - Per-species breakdown
  - Benchmark timing table (s/min per rule)
  - Links to output files
"""

import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("generate_report", snakemake.log[0])
out_html = snakemake.output.html
UPSTREAM = snakemake.params.upstream
DOWNSTREAM = snakemake.params.downstream

in_classified = snakemake.input.classified
in_resolved = snakemake.input.resolved
in_unresolved = snakemake.input.unresolved
in_ambiguous = snakemake.input.ambiguous
in_fasta = snakemake.input.fasta
in_bed = snakemake.input.bed
in_benchmarks = snakemake.input.benchmarks

log.info("Stage 5: Generating HTML summary report")

# ── Load data ─────────────────────────────────────────────────
df_cls = pd.read_csv(in_classified, sep="\t")
df_res = pd.read_csv(in_resolved, sep="\t")
df_unr = pd.read_csv(in_unresolved, sep="\t")
df_amb = pd.read_csv(in_ambiguous, sep="\t")

# Count FASTA output records
fasta_count = sum(1 for line in open(in_fasta) if line.startswith(">"))
bed_count = max(0, sum(1 for _ in open(in_bed)) - 1)  # minus header

# ── Benchmark table ───────────────────────────────────────────
bench_rows = []
for bpath in in_benchmarks:
    rule_name = Path(bpath).stem
    try:
        bdf = pd.read_csv(bpath, sep="\t")
        row = bdf.iloc[0]
        bench_rows.append(
            {
                "Rule": rule_name,
                "Wall time": f"{float(row.get('s', 0)):.1f} s",
                "CPU time": f"{float(row.get('cpu_time', row.get('s', 0))):.1f} s",
                "Max RSS (MB)": f"{float(row.get('max_rss', 0)):.0f}",
            }
        )
    except Exception as exc:
        log.warning(f"Could not parse benchmark {bpath}: {exc}")
        bench_rows.append(
            {
                "Rule": rule_name,
                "Wall time": "N/A",
                "CPU time": "N/A",
                "Max RSS (MB)": "N/A",
            }
        )

df_bench = pd.DataFrame(bench_rows)

# ── Per-DB stats ──────────────────────────────────────────────
db_counts = df_cls["db_source"].value_counts().to_dict()
db_resolved = df_res["db_source"].value_counts().to_dict()

# ── Per-species stats ─────────────────────────────────────────
if "organism" in df_res.columns:
    species_counts = df_res["organism"].value_counts().head(20)
else:
    species_counts = pd.Series(dtype=int)


# ── HTML helpers ──────────────────────────────────────────────
def df_to_html_table(df: pd.DataFrame, table_id: str = "", max_rows: int = 200) -> str:
    if df.empty:
        return "<p><em>No data.</em></p>"
    if len(df) > max_rows:
        df = df.head(max_rows)
        footer = f"<p><em>Showing first {max_rows} rows.</em></p>"
    else:
        footer = ""
    html = df.to_html(index=False, border=0, classes="data-table", table_id=table_id)
    return html + footer


def stat_card(label: str, value, color: str = "#3b82f6") -> str:
    return f"""
    <div class="stat-card">
      <div class="stat-value" style="color:{color}">{value}</div>
      <div class="stat-label">{label}</div>
    </div>"""


# ── HTML template ─────────────────────────────────────────────
total_input = len(df_cls) + len(df_unr)
total_cls = len(df_cls)
total_res = len(df_res)
total_unres = len(df_unr)
total_ambig = len(df_amb)
total_out = fasta_count

now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

db_rows_html = "".join(
    f"<tr><td>{db}</td><td>{db_counts.get(db,0)}</td><td>{db_resolved.get(db,0)}</td></tr>"
    for db in ["ncbi", "ensembl", "ucsc"]
)

species_rows_html = "".join(
    f"<tr><td>{sp}</td><td>{cnt}</td></tr>" for sp, cnt in species_counts.items()
)

HTML = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>RNA Flanking Pipeline Report</title>
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ font-family: 'Segoe UI', sans-serif; background: #f8fafc; color: #1e293b; }}
    header {{ background: #1e40af; color: white; padding: 24px 40px; }}
    header h1 {{ font-size: 1.8rem; }}
    header p  {{ opacity: 0.8; font-size: 0.9rem; margin-top: 4px; }}
    main {{ max-width: 1100px; margin: 0 auto; padding: 32px 20px; }}
    h2 {{ font-size: 1.2rem; margin: 32px 0 12px; color: #1e40af; border-bottom: 2px solid #bfdbfe; padding-bottom: 6px; }}
    h3 {{ font-size: 1rem; margin: 16px 0 8px; color: #475569; }}
    .stat-grid {{ display: flex; flex-wrap: wrap; gap: 16px; margin: 16px 0; }}
    .stat-card {{ background: white; border-radius: 10px; padding: 20px 24px; min-width: 160px;
                  box-shadow: 0 1px 3px rgba(0,0,0,.1); flex: 1; text-align: center; }}
    .stat-value {{ font-size: 2rem; font-weight: 700; }}
    .stat-label {{ font-size: 0.8rem; color: #64748b; margin-top: 4px; }}
    table.data-table {{ width: 100%; border-collapse: collapse; background: white;
                        border-radius: 8px; overflow: hidden;
                        box-shadow: 0 1px 3px rgba(0,0,0,.08); font-size: 0.875rem; }}
    table.data-table th {{ background: #e0e7ff; color: #1e40af; padding: 10px 14px;
                           text-align: left; font-weight: 600; }}
    table.data-table td {{ padding: 8px 14px; border-bottom: 1px solid #f1f5f9; }}
    table.data-table tr:hover td {{ background: #f8fafc; }}
    .badge {{ display: inline-block; padding: 2px 8px; border-radius: 99px;
              font-size: 0.75rem; font-weight: 600; }}
    .badge-green  {{ background: #dcfce7; color: #166534; }}
    .badge-yellow {{ background: #fef9c3; color: #854d0e; }}
    .badge-red    {{ background: #fee2e2; color: #991b1b; }}
    .info-box {{ background: #eff6ff; border-left: 4px solid #3b82f6;
                 padding: 12px 16px; border-radius: 4px; font-size: 0.875rem; margin: 12px 0; }}
    footer {{ text-align: center; color: #94a3b8; font-size: 0.8rem; padding: 24px; margin-top: 32px; }}
  </style>
</head>
<body>
<header>
  <h1>🧬 RNA Flanking Sequence Pipeline — Run Report</h1>
  <p>Generated: {now} &nbsp;|&nbsp; Upstream: {UPSTREAM} bp &nbsp;|&nbsp; Downstream: {DOWNSTREAM} bp</p>
</header>
<main>

  <h2>Overview</h2>
  <div class="stat-grid">
    {stat_card("Total Input IDs",       total_input,  "#1e40af")}
    {stat_card("Classified",            total_cls,    "#0369a1")}
    {stat_card("Resolved",              total_res,    "#059669")}
    {stat_card("Unresolved",            total_unres,  "#dc2626")}
    {stat_card("Ambiguous (alts)",      total_ambig,  "#d97706")}
    {stat_card("Output Sequences",      total_out,    "#7c3aed")}
  </div>

  <div class="info-box">
    Pipeline configuration: <strong>{UPSTREAM} bp upstream</strong> and
    <strong>{DOWNSTREAM} bp downstream</strong> of each gene.
    Output: <code>output.fasta</code> ({fasta_count} sequences) and
    <code>output.bed</code> ({bed_count} records).
  </div>

  <h2>By Database Source</h2>
  <table class="data-table">
    <thead><tr><th>Database</th><th>Input IDs</th><th>Resolved</th></tr></thead>
    <tbody>{db_rows_html}</tbody>
  </table>

  <h2>By Species (top 20)</h2>
  <table class="data-table">
    <thead><tr><th>Organism</th><th>Resolved Transcripts</th></tr></thead>
    <tbody>{species_rows_html if species_rows_html else "<tr><td colspan='2'>No species data available</td></tr>"}</tbody>
  </table>

  <h2>Unresolved IDs</h2>
  {df_to_html_table(df_unr, table_id="tbl-unresolved")}

  <h2>Ambiguous IDs (alternatives not chosen)</h2>
  {df_to_html_table(df_amb, table_id="tbl-ambiguous")}

  <h2>Benchmark — Rule Timing</h2>
  {df_to_html_table(df_bench, table_id="tbl-benchmarks")}

</main>
<footer>RNA Flanking Sequence Pipeline &nbsp;|&nbsp; Powered by Snakemake + BioPython</footer>
</body>
</html>"""

Path(out_html).write_text(HTML, encoding="utf-8")

log.info("=" * 60)
log.info(f"Total input IDs              : {total_input}")
log.info(f"Resolved                     : {total_res}")
log.info(f"Unresolved                   : {total_unres}")
log.info(f"Ambiguous alternatives       : {total_ambig}")
log.info(f"Output FASTA sequences       : {fasta_count}")
log.info(f"Written report → {out_html}")
log.info("Stage 5 complete. Pipeline finished.")
