"""
scripts/generate_report.py
Stage 5 — Generate Summary Report
===================================
Produces results/report.html covering:
  - Total transcripts, breakdown by DB
  - Ensembl: species detected, BioMart release used
  - Resolution success/failure/ambiguity per source
  - Per-species breakdown
  - Benchmark timing table (all rules)
"""

import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

log          = get_logger("generate_report", snakemake.log[0])
out_html     = snakemake.output.html
UPSTREAM     = snakemake.params.upstream
DOWNSTREAM   = snakemake.params.downstream
RELEASE      = snakemake.params.release

in_classified      = snakemake.input.classified
in_resolved        = snakemake.input.resolved
in_unresolved      = snakemake.input.unresolved
in_ambiguous       = snakemake.input.ambiguous
in_species_map     = snakemake.input.species_map
in_unknown_pfx     = snakemake.input.unknown_prefixes
in_fasta           = snakemake.input.fasta
in_bed             = snakemake.input.bed
in_benchmarks      = snakemake.input.benchmarks

log.info("Generating HTML summary report")

df_cls    = pd.read_csv(in_classified,  sep="\t")
df_res    = pd.read_csv(in_resolved,    sep="\t")
df_unr    = pd.read_csv(in_unresolved,  sep="\t")
df_amb    = pd.read_csv(in_ambiguous,   sep="\t")
df_spmap  = pd.read_csv(in_species_map, sep="\t")
df_unkpfx = pd.read_csv(in_unknown_pfx, sep="\t") if Path(in_unknown_pfx).exists() else pd.DataFrame()

fasta_count = sum(1 for line in open(in_fasta) if line.startswith(">"))
bed_count   = max(0, sum(1 for _ in open(in_bed)) - 1)

# ── Benchmark table ───────────────────────────────────────────
bench_rows = []
for bpath in in_benchmarks:
    rule_name = Path(bpath).stem
    try:
        bdf = pd.read_csv(bpath, sep="\t")
        row = bdf.iloc[0]
        bench_rows.append({
            "Rule":          rule_name,
            "Wall time (s)": f"{float(row.get('s', 0)):.1f}",
            "CPU time (s)":  f"{float(row.get('cpu_time', row.get('s', 0))):.1f}",
            "Max RSS (MB)":  f"{float(row.get('max_rss', 0)):.0f}",
        })
    except Exception as exc:
        log.warning(f"Could not parse benchmark {bpath}: {exc}")
        bench_rows.append({"Rule": rule_name, "Wall time (s)": "N/A", "CPU time (s)": "N/A", "Max RSS (MB)": "N/A"})

df_bench = pd.DataFrame(bench_rows)

# ── Per-DB stats ──────────────────────────────────────────────
db_input    = df_cls["db_source"].value_counts().to_dict()
db_resolved = df_res["db_source"].value_counts().to_dict() if "db_source" in df_res.columns else {}

# ── Ensembl species breakdown ─────────────────────────────────
ensembl_species_rows = ""
if not df_spmap.empty and "species" in df_spmap.columns:
    sp_counts = df_spmap.groupby(["species", "build"]).size().reset_index(name="n_transcripts")
    for _, r in sp_counts.iterrows():
        ensembl_species_rows += (
            f"<tr><td>{r['species']}</td><td>{r['build']}</td>"
            f"<td>{r['n_transcripts']}</td></tr>"
        )

# ── Per-species resolved breakdown ────────────────────────────
species_rows_html = ""
if "organism" in df_res.columns:
    for sp, cnt in df_res["organism"].value_counts().head(20).items():
        species_rows_html += f"<tr><td>{sp}</td><td>{cnt}</td></tr>"

# ── Helpers ───────────────────────────────────────────────────
def df_to_html_table(df: pd.DataFrame, max_rows: int = 200) -> str:
    if df.empty:
        return "<p><em>No data.</em></p>"
    if len(df) > max_rows:
        df = df.head(max_rows)
        footer = f"<p><em>Showing first {max_rows} rows.</em></p>"
    else:
        footer = ""
    return df.to_html(index=False, border=0, classes="data-table") + footer

def stat_card(label, value, color="#3b82f6"):
    return f"""<div class="stat-card">
      <div class="stat-value" style="color:{color}">{value}</div>
      <div class="stat-label">{label}</div></div>"""

# ── Stats ─────────────────────────────────────────────────────
total_input  = len(df_cls) + len(df_unr)
now          = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

db_table_rows = "".join(
    f"<tr><td>{db}</td><td>{db_input.get(db,0)}</td>"
    f"<td>{db_resolved.get(db,0)}</td>"
    f"<td>{'BioMart wrapper' if db == 'ensembl' else 'REST API'}</td></tr>"
    for db in ["ncbi", "ensembl", "ucsc"]
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
    header p  {{ opacity: 0.8; font-size: 0.9rem; margin-top: 6px; }}
    main {{ max-width: 1100px; margin: 0 auto; padding: 32px 20px; }}
    h2 {{ font-size: 1.2rem; margin: 32px 0 12px; color: #1e40af;
          border-bottom: 2px solid #bfdbfe; padding-bottom: 6px; }}
    .stat-grid {{ display: flex; flex-wrap: wrap; gap: 16px; margin: 16px 0; }}
    .stat-card {{ background: white; border-radius: 10px; padding: 20px 24px;
                  min-width: 150px; box-shadow: 0 1px 3px rgba(0,0,0,.1);
                  flex: 1; text-align: center; }}
    .stat-value {{ font-size: 2rem; font-weight: 700; }}
    .stat-label {{ font-size: 0.8rem; color: #64748b; margin-top: 4px; }}
    table.data-table {{ width: 100%; border-collapse: collapse; background: white;
                        border-radius: 8px; overflow: hidden;
                        box-shadow: 0 1px 3px rgba(0,0,0,.08); font-size: 0.875rem; }}
    table.data-table th {{ background: #e0e7ff; color: #1e40af; padding: 10px 14px;
                           text-align: left; font-weight: 600; }}
    table.data-table td {{ padding: 8px 14px; border-bottom: 1px solid #f1f5f9; }}
    table.data-table tr:hover td {{ background: #f8fafc; }}
    .info-box {{ background: #eff6ff; border-left: 4px solid #3b82f6;
                 padding: 12px 16px; border-radius: 4px; font-size: 0.875rem; margin: 12px 0; }}
    .biomart-box {{ background: #f0fdf4; border-left: 4px solid #22c55e;
                    padding: 12px 16px; border-radius: 4px; font-size: 0.875rem; margin: 12px 0; }}
    footer {{ text-align: center; color: #94a3b8; font-size: 0.8rem;
              padding: 24px; margin-top: 32px; }}
  </style>
</head>
<body>
<header>
  <h1>🧬 RNA Flanking Sequence Pipeline — Run Report</h1>
  <p>Generated: {now} &nbsp;|&nbsp;
     Upstream: {UPSTREAM} bp &nbsp;|&nbsp; Downstream: {DOWNSTREAM} bp &nbsp;|&nbsp;
     Ensembl release: {RELEASE}</p>
</header>
<main>

  <h2>Overview</h2>
  <div class="stat-grid">
    {stat_card("Total Input IDs",   total_input,    "#1e40af")}
    {stat_card("Classified",        len(df_cls),    "#0369a1")}
    {stat_card("Resolved",          len(df_res),    "#059669")}
    {stat_card("Unresolved",        len(df_unr),    "#dc2626")}
    {stat_card("Ambiguous (alts)",  len(df_amb),    "#d97706")}
    {stat_card("Output Sequences",  fasta_count,    "#7c3aed")}
  </div>

  <div class="info-box">
    Output: <code>output.fasta</code> ({fasta_count} sequences) and
    <code>output.bed</code> ({bed_count} records) with
    <strong>{UPSTREAM} bp upstream</strong> and <strong>{DOWNSTREAM} bp downstream</strong>
    flanking regions per gene.
  </div>

  <h2>Resolution by Database</h2>
  <div class="biomart-box">
    <strong>Ensembl IDs</strong> were resolved using the
    <strong>snakemake-wrappers BioMart wrapper</strong>
    (release <strong>{RELEASE}</strong>) — one bulk download per species,
    fully reproducible and cached across pipeline runs.
    NCBI and UCSC IDs were resolved via their respective REST APIs.
  </div>
  <table class="data-table">
    <thead><tr><th>Database</th><th>Input IDs</th><th>Resolved</th><th>Resolution method</th></tr></thead>
    <tbody>{db_table_rows}</tbody>
  </table>

  <h2>Ensembl Species Detected (BioMart runs)</h2>
  <table class="data-table">
    <thead><tr><th>Species</th><th>Genome Build</th><th>Transcripts</th></tr></thead>
    <tbody>{ensembl_species_rows if ensembl_species_rows else
            "<tr><td colspan='3'>No Ensembl species detected</td></tr>"}</tbody>
  </table>

  <h2>Unknown Ensembl Prefixes</h2>
  {"<p style='color:#059669'>✓ All Ensembl transcript prefixes were recognised automatically.</p>" if df_unkpfx.empty else
   "<div style='background:#fef2f2;border-left:4px solid #dc2626;padding:14px 16px;border-radius:4px;margin:12px 0'>"
   "<strong style='color:#dc2626'>⚠ Unknown prefixes detected during this run.</strong>"
   " The pipeline stopped and required manual intervention. "
   "These prefixes were not in the built-in species reference table — they were resolved "
   "by adding entries to <code>ensembl_species_overrides</code> in config.yaml.</div>"
   + df_to_html_table(df_unkpfx)}

  <h2>Resolved Transcripts by Organism</h2>
  <table class="data-table">
    <thead><tr><th>Organism</th><th>Resolved Transcripts</th></tr></thead>
    <tbody>{species_rows_html if species_rows_html else
            "<tr><td colspan='2'>No organism data available</td></tr>"}</tbody>
  </table>

  <h2>Unresolved IDs</h2>
  {df_to_html_table(df_unr)}

  <h2>Ambiguous IDs (alternatives not chosen)</h2>
  {df_to_html_table(df_amb)}

  <h2>Benchmark — Rule Timing</h2>
  {df_to_html_table(df_bench)}

</main>
<footer>RNA Flanking Sequence Pipeline &nbsp;|&nbsp;
Snakemake + BioPython + snakemake-wrappers BioMart</footer>
</body>
</html>"""

Path(out_html).write_text(HTML, encoding="utf-8")

log.info("=" * 60)
log.info(f"Total input          : {total_input}")
log.info(f"Resolved             : {len(df_res)}")
log.info(f"Unresolved           : {len(df_unr)}")
log.info(f"Ambiguous alts       : {len(df_amb)}")
log.info(f"Output sequences     : {fasta_count}")
log.info(f"Written report → {out_html}")
log.info("Stage 5 complete.")
