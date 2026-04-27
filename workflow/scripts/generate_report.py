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
in_extraction_failed = snakemake.input.extraction_failed

in_ncbi_asm_res   = snakemake.input.ncbi_assembly_resolved
in_ncbi_asm_unr   = snakemake.input.ncbi_assembly_unresolved
in_ensl_asm_res   = snakemake.input.ensembl_assembly_resolved
in_ensl_asm_unr   = snakemake.input.ensembl_assembly_unresolved
in_nonc_asm_res   = snakemake.input.noncode_assembly_resolved
in_nonc_asm_unr   = snakemake.input.noncode_assembly_unresolved

in_gramene_res     = snakemake.input.gramene_resolved
in_gramene_unr     = snakemake.input.gramene_unresolved
in_phyto_res       = snakemake.input.phytozome_resolved
in_phyto_unr       = snakemake.input.phytozome_unresolved
in_noncode_res        = snakemake.input.noncode_resolved
in_noncode_unr        = snakemake.input.noncode_unresolved
in_noncode_v4_res     = snakemake.input.noncode_v4_resolved
in_noncode_v4_unr     = snakemake.input.noncode_v4_unresolved
in_noncode_2016_res   = snakemake.input.noncode_2016_resolved
in_noncode_2016_unr   = snakemake.input.noncode_2016_unresolved
in_plant_gtf_res      = snakemake.input.plant_gtf_resolved
in_plant_gtf_unr      = snakemake.input.plant_gtf_unresolved
in_worm_gtf_res       = snakemake.input.worm_gtf_resolved
in_worm_gtf_unr       = snakemake.input.worm_gtf_unresolved
in_fly_gtf_res        = snakemake.input.fly_gtf_resolved
in_fly_gtf_unr        = snakemake.input.fly_gtf_unresolved
in_yeast_gtf_res      = snakemake.input.yeast_gtf_resolved
in_yeast_gtf_unr      = snakemake.input.yeast_gtf_unresolved

log.info("Generating HTML summary report")

df_cls    = pd.read_csv(in_classified,  sep="\t")
df_res    = pd.read_csv(in_resolved,    sep="\t")
df_unr    = pd.read_csv(in_unresolved,  sep="\t")
df_amb    = pd.read_csv(in_ambiguous,   sep="\t")
df_spmap  = pd.read_csv(in_species_map, sep="\t")
df_unkpfx = pd.read_csv(in_unknown_pfx, sep="\t") if Path(in_unknown_pfx).exists() else pd.DataFrame()

df_failed = pd.read_csv(in_extraction_failed, sep="\t") if Path(in_extraction_failed).exists() else pd.DataFrame(columns=["transcript_id", "assembly_accession", "chrom", "fail_reason"])

df_ncbi_asm_res  = pd.read_csv(in_ncbi_asm_res,  sep="\t") if Path(in_ncbi_asm_res).exists()  else pd.DataFrame()
df_ncbi_asm_unr  = pd.read_csv(in_ncbi_asm_unr,  sep="\t") if Path(in_ncbi_asm_unr).exists()  else pd.DataFrame()
df_ensl_asm_res  = pd.read_csv(in_ensl_asm_res,  sep="\t") if Path(in_ensl_asm_res).exists()  else pd.DataFrame()
df_ensl_asm_unr  = pd.read_csv(in_ensl_asm_unr,  sep="\t") if Path(in_ensl_asm_unr).exists()  else pd.DataFrame()
df_nonc_asm_res  = pd.read_csv(in_nonc_asm_res,  sep="\t") if Path(in_nonc_asm_res).exists()  else pd.DataFrame()
df_nonc_asm_unr  = pd.read_csv(in_nonc_asm_unr,  sep="\t") if Path(in_nonc_asm_unr).exists()  else pd.DataFrame()

def _read(p):
    return pd.read_csv(p, sep="\t") if Path(p).exists() else pd.DataFrame()

df_gramene_res   = _read(in_gramene_res)
df_gramene_unr   = _read(in_gramene_unr)
df_phyto_res     = _read(in_phyto_res)
df_phyto_unr     = _read(in_phyto_unr)
df_noncode_res      = _read(in_noncode_res)
df_noncode_unr      = _read(in_noncode_unr)
df_noncode_v4_res   = _read(in_noncode_v4_res)
df_noncode_v4_unr   = _read(in_noncode_v4_unr)
df_noncode_2016_res = _read(in_noncode_2016_res)
df_noncode_2016_unr = _read(in_noncode_2016_unr)
df_plant_gtf_res    = _read(in_plant_gtf_res)
df_plant_gtf_unr    = _read(in_plant_gtf_unr)
df_worm_gtf_res     = _read(in_worm_gtf_res)
df_worm_gtf_unr     = _read(in_worm_gtf_unr)
df_fly_gtf_res      = _read(in_fly_gtf_res)
df_fly_gtf_unr      = _read(in_fly_gtf_unr)
df_yeast_gtf_res    = _read(in_yeast_gtf_res)
df_yeast_gtf_unr    = _read(in_yeast_gtf_unr)

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

# Pipeline funnel
_asm = df_res["assembly_accession"] if "assembly_accession" in df_res.columns else pd.Series(dtype=str)
n_has_assembly = int((_asm.notna() & (_asm.str.strip() != "") & (_asm.str.lower() != "nan") & (_asm.str.lower() != "none")).sum())
n_has_coords   = int((df_res["start"].notna() & df_res["end"].notna()).sum()) if {"start","end"} <= set(df_res.columns) else 0
n_total_failed = len(df_failed)

def _pct(num, denom):
    return f"{100*num/denom:.0f}%" if denom else "—"

funnel_rows = [
    ("Input IDs",             total_input,      _pct(total_input,    total_input)),
    ("Classified",            len(df_cls),       _pct(len(df_cls),    total_input)),
    ("Resolved",              len(df_res),       _pct(len(df_res),    total_input)),
    ("Has assembly accession",n_has_assembly,    _pct(n_has_assembly, len(df_res))),
    ("Has coordinates",       n_has_coords,      _pct(n_has_coords,   len(df_res))),
    ("Sequence extracted",    fasta_count,       _pct(fasta_count,    len(df_res))),
    ("Extraction failed",     n_total_failed,    _pct(n_total_failed, len(df_res))),
]
funnel_html = "".join(
    f"<tr><td>{step}</td><td>{cnt}</td><td>{pct}</td></tr>"
    for step, cnt, pct in funnel_rows
)

# ── Per-resolver resolution funnel ────────────────────────────
RESOLVERS = [
    ("BioMart Ensembl",    "ensembl",      df_ensl_asm_res,      df_ensl_asm_unr,      "BioMart API"),
    ("BioMart Plants GTF", "plant_gtf",    df_plant_gtf_res,     df_plant_gtf_unr,     "Ensembl Plants GTF"),
    ("Gramene",            "gramene",       df_gramene_res,       df_gramene_unr,        "Gramene REST API"),
    ("Phytozome",          "phytozome",     df_phyto_res,         df_phyto_unr,          "JGI Phytozome GFF3"),
    ("NCBI GenBank",       "ncbi",          df_ncbi_asm_res,      df_ncbi_asm_unr,       "NCBI Entrez API"),
    ("NONCODE v5",         "noncode",       df_noncode_res,       df_noncode_unr,        "NONCODE BED/FASTA"),
    ("NONCODE v4",         "noncode_v4",    df_noncode_v4_res,    df_noncode_v4_unr,     "NONCODE v4 BED/FASTA"),
    ("NONCODE 2016",       "noncode_2016",  df_noncode_2016_res,  df_noncode_2016_unr,   "NONCODE 2016 BED/FASTA"),
    ("WormBase GTF",       "wormbase",      df_worm_gtf_res,      df_worm_gtf_unr,       "Ensembl Metazoa GTF"),
    ("FlyBase GTF",        "flybase",       df_fly_gtf_res,       df_fly_gtf_unr,        "Ensembl Metazoa GTF"),
    ("SGD/Yeast GTF",      "sgd",           df_yeast_gtf_res,     df_yeast_gtf_unr,      "SGD GTF"),
    ("Phytozome (plants)", "plant",         df_phyto_res,         df_phyto_unr,          "JGI Phytozome GFF3"),
]

resolver_rows_html = ""
for label, db, df_r, df_u, method in RESOLVERS:
    n_in  = db_input.get(db, 0)
    n_res = len(df_r)
    n_unr = len(df_u)
    rate  = f"{100*n_res/(n_res+n_unr):.0f}%" if (n_res + n_unr) > 0 else "—"
    resolver_rows_html += (
        f"<tr><td>{label}</td><td><code>{db}</code></td>"
        f"<td>{n_in}</td><td>{n_res}</td><td>{n_unr}</td>"
        f"<td>{rate}</td><td>{method}</td></tr>"
    )

# ── Gramene stats ─────────────────────────────────────────────
gramene_species_rows = ""
if not df_gramene_res.empty and "species" in df_gramene_res.columns:
    for sp, cnt in df_gramene_res["species"].value_counts().head(15).items():
        gramene_species_rows += f"<tr><td>{sp}</td><td>{cnt}</td></tr>"

gramene_biotype_rows = ""
if not df_gramene_res.empty and "biotype" in df_gramene_res.columns:
    for bt, cnt in df_gramene_res["biotype"].value_counts().head(10).items():
        gramene_biotype_rows += f"<tr><td>{bt}</td><td>{cnt}</td></tr>"

# ── Extraction failure breakdown by db_source ─────────────────
if not df_failed.empty and "db_source" in df_failed.columns:
    fail_by_db = (
        df_failed.groupby(["db_source", "fail_reason"])
        .size()
        .reset_index(name="count")
        .sort_values(["db_source", "count"], ascending=[True, False])
    )
    fail_detail_html = df_to_html_table(fail_by_db)
else:
    fail_detail_html = "<p><em>No failure data.</em></p>"

# ── NONCODE version breakdown ─────────────────────────────────
noncode_rows_html = ""
for label, df_r, df_u in [
    ("NONCODE v5",   df_noncode_res,     df_noncode_unr),
    ("NONCODE v4",   df_noncode_v4_res,  df_noncode_v4_unr),
    ("NONCODE 2016", df_noncode_2016_res, df_noncode_2016_unr),
]:
    noncode_rows_html += (
        f"<tr><td>{label}</td><td>{len(df_r)}</td><td>{len(df_u)}</td>"
        f"<td>{len(df_r)+len(df_u)}</td></tr>"
    )

asm_rows_html = ""
for resolver, df_r, df_u in [
    ("NCBI",    df_ncbi_asm_res, df_ncbi_asm_unr),
    ("Ensembl", df_ensl_asm_res, df_ensl_asm_unr),
    ("NONCODE", df_nonc_asm_res, df_nonc_asm_unr),
]:
    n_res = len(df_r)
    n_unr = len(df_u)
    asm_rows_html += (
        f"<tr><td>{resolver}</td><td>{n_res}</td><td>{n_unr}</td>"
        f"<td>{n_res + n_unr}</td></tr>"
    )

fail_counts = df_failed["fail_reason"].value_counts().to_dict() if not df_failed.empty else {}
n_asm_missing   = fail_counts.get("assembly_not_cached", 0)
n_coord_missing = fail_counts.get("missing_coordinates", 0)
n_chrom_missing = fail_counts.get("chrom_not_found", 0)
n_seq_error     = fail_counts.get("sequence_error", 0)
n_total_failed  = len(df_failed)

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
    {stat_card("Total Input IDs",        total_input,       "#1e40af")}
    {stat_card("Classified",             len(df_cls),       "#0369a1")}
    {stat_card("Resolved",               len(df_res),       "#059669")}
    {stat_card("Has Assembly Accession", n_has_assembly,    "#0891b2")}
    {stat_card("Has Coordinates",        n_has_coords,      "#0891b2")}
    {stat_card("Unresolved",             len(df_unr),       "#dc2626")}
    {stat_card("Ambiguous (alts)",       len(df_amb),       "#d97706")}
    {stat_card("Output Sequences",       fasta_count,       "#7c3aed")}
    {stat_card("Extraction Failures",    n_total_failed,    "#dc2626")}
  </div>

  <h2>Pipeline Funnel</h2>
  <table class="data-table">
    <thead><tr><th>Stage</th><th>Transcripts</th><th>% of resolved (or input)</th></tr></thead>
    <tbody>{funnel_html}</tbody>
  </table>

  <div class="info-box">
    Output: <code>output.fasta</code> ({fasta_count} sequences) and
    <code>output.bed</code> ({bed_count} records) with
    <strong>{UPSTREAM} bp upstream</strong> and <strong>{DOWNSTREAM} bp downstream</strong>
    flanking regions per gene.
  </div>

  <h2>Extraction Failures</h2>
  <table class="data-table">
    <thead><tr><th>Failure reason</th><th>Count</th></tr></thead>
    <tbody>
      <tr><td>Assembly not cached</td><td>{n_asm_missing}</td></tr>
      <tr><td>Missing coordinates (NaN)</td><td>{n_coord_missing}</td></tr>
      <tr><td>Chromosome not found in index</td><td>{n_chrom_missing}</td></tr>
      <tr><td>Sequence extraction error</td><td>{n_seq_error}</td></tr>
    </tbody>
  </table>

  <h2>Extraction Failures by Resolver and Reason</h2>
  {fail_detail_html}

  <h2>Resolution by Resolver</h2>
  <table class="data-table">
    <thead><tr>
      <th>Resolver</th><th>db_source</th><th>Input IDs</th>
      <th>Resolved</th><th>Unresolved</th><th>Rate</th><th>Method</th>
    </tr></thead>
    <tbody>{resolver_rows_html}</tbody>
  </table>

  <h2>Ensembl Species Detected (BioMart runs)</h2>
  <table class="data-table">
    <thead><tr><th>Species</th><th>Genome Build</th><th>Transcripts</th></tr></thead>
    <tbody>{ensembl_species_rows if ensembl_species_rows else
            "<tr><td colspan='3'>No Ensembl species detected</td></tr>"}</tbody>
  </table>

  <h2>Gramene Resolution Details</h2>
  <div style="display:flex;gap:24px;flex-wrap:wrap;">
    <div style="flex:1;min-width:280px;">
      <strong style="display:block;margin-bottom:8px;">Top species</strong>
      <table class="data-table">
        <thead><tr><th>Species</th><th>Transcripts</th></tr></thead>
        <tbody>{gramene_species_rows if gramene_species_rows else "<tr><td colspan='2'>No data</td></tr>"}</tbody>
      </table>
    </div>
    <div style="flex:1;min-width:280px;">
      <strong style="display:block;margin-bottom:8px;">Biotype breakdown</strong>
      <table class="data-table">
        <thead><tr><th>Biotype</th><th>Transcripts</th></tr></thead>
        <tbody>{gramene_biotype_rows if gramene_biotype_rows else "<tr><td colspan='2'>No data</td></tr>"}</tbody>
      </table>
    </div>
  </div>

  <h2>NONCODE Version Breakdown</h2>
  <table class="data-table">
    <thead><tr><th>Version</th><th>Resolved</th><th>Unresolved</th><th>Total Input</th></tr></thead>
    <tbody>{noncode_rows_html}</tbody>
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

  <h2>Assembly Accession Resolution</h2>
  <table class="data-table">
    <thead><tr><th>Resolver</th><th>Assembly Resolved</th><th>Assembly Unresolved</th><th>Total Input</th></tr></thead>
    <tbody>{asm_rows_html}</tbody>
  </table>

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
