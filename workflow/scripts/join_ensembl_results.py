"""
scripts/join_ensembl_results.py
Join BioMart Tables Against Ensembl Transcript IDs
===================================================
After all per-species BioMart TSVs are downloaded by the wrapper,
this script:
  1. Loads and concatenates all BioMart tables
  2. Normalises column names to the pipeline's unified schema
  3. Handles versioned IDs (e.g. ENST00000456328.2 → match on base ID)
  4. Left-joins our Ensembl transcript IDs against the combined table
  5. Detects and records ambiguous mappings (transcript → multiple genes)
  6. Converts BioMart strand (1/-1) → (+/-)
  7. Writes ensembl_resolved.tsv and ensembl_ambiguous.tsv

BioMart attributes used (set in biomart_lookup.smk):
  ensembl_transcript_id        — base ID without version
  ensembl_transcript_id_version — full versioned ID
  ensembl_gene_id
  external_gene_name
  chromosome_name
  start_position
  end_position
  strand                       — 1 or -1
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("join_ensembl_results", snakemake.log[0])
input_classified = snakemake.input.classified
input_species_map = snakemake.input.species_map
biomart_tables = snakemake.input.biomart_tables  # list of per-species .tsv.gz
out_resolved = snakemake.output.resolved
out_ambiguous = snakemake.output.ambiguous
out_unresolved = snakemake.output.unresolved
cfg = snakemake.config

# ── Unified output columns (same schema as NCBI/UCSC resolver) ─
RESOLVED_COLS = [
    "transcript_id",
    "db_source",
    "gene_id",
    "gene_symbol",
    "organism",
    "assembly_accession",
    "chrom",
    "start",
    "end",
    "strand",
    "is_ambiguous",
]
AMBIG_COLS = [
    "transcript_id",
    "db_source",
    "chosen_gene_id",
    "alternative_gene_id",
    "alternative_gene_symbol",
    "organism",
    "assembly_accession",
    "chrom",
    "start",
    "end",
    "strand",
]
UNRESOLVED_COLS = ["transcript_id", "db_source", "reason"]

log.info("join_ensembl_results: joining BioMart tables against Ensembl transcript IDs")

# ── Load our Ensembl transcript IDs ──────────────────────────
df_cls = pd.read_csv(input_classified, sep="\t")
df_ensembl_ids = df_cls[df_cls["db_source"] == "ensembl"][["transcript_id"]].copy()
log.info(f"Our Ensembl IDs to resolve: {len(df_ensembl_ids)}")

# ── Load species map (to get organism + build per transcript) ─
df_species_map = pd.read_csv(input_species_map, sep="\t")
# transcript_id | prefix | species | build

# ── Load and concatenate all BioMart tables ───────────────────
biomart_dfs: list[pd.DataFrame] = []

for tbl_path in biomart_tables:
    species = Path(tbl_path).name.replace(".tsv.gz", "")
    log.info(f"  Loading BioMart table: {tbl_path}")
    try:
        df_tbl = pd.read_csv(tbl_path, sep="\t", compression="gzip", low_memory=False)
    except Exception as exc:
        log.error(f"  Failed to load {tbl_path}: {exc}")
        continue

    log.info(f"  {species}: {len(df_tbl):,} rows, columns: {list(df_tbl.columns)}")

    # Attach species + build for later use
    df_tbl["_species"] = species
    # Get the build from species map
    build_match = df_species_map[df_species_map["species"] == species]["build"]
    df_tbl["_build"] = build_match.iloc[0] if len(build_match) > 0 else ""

    biomart_dfs.append(df_tbl)

if not biomart_dfs:
    log.error("No BioMart tables loaded — writing empty outputs")
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=AMBIG_COLS).to_csv(out_ambiguous, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    sys.exit(0)

df_biomart = pd.concat(biomart_dfs, ignore_index=True)
log.info(
    f"Combined BioMart table: {len(df_biomart):,} rows across {len(biomart_dfs)} species"
)

# ── Normalise column names ────────────────────────────────────
# BioMart REST TSV returns human-readable headers:
col_map = {
    "Transcript stable ID": "base_transcript_id",
    "Transcript stable ID version": "versioned_transcript_id",
    "Gene stable ID": "gene_id",
    "Gene name": "gene_symbol",
    "Chromosome/scaffold name": "chrom",
    "Gene start (bp)": "start",
    "Gene end (bp)": "end",
    "Strand": "strand_raw",
}
df_biomart = df_biomart.rename(
    columns={k: v for k, v in col_map.items() if k in df_biomart.columns}
)

# Convert strand: 1 → "+", -1 → "-"
if "strand_raw" in df_biomart.columns:
    df_biomart["strand"] = df_biomart["strand_raw"].apply(
        lambda s: "+" if int(s) == 1 else "-"
    )

# ── Match our IDs against the BioMart table ───────────────────
# IDs in our input may be versioned (ENST00000456328.2) or bare (ENST00000456328).
# BioMart provides both; we try versioned first, fall back to base.


def strip_version(tid: str) -> str:
    return tid.split(".")[0]


df_ensembl_ids["base_id"] = df_ensembl_ids["transcript_id"].apply(strip_version)

# Build lookup index on both versioned and base IDs for O(1) lookup
if "versioned_transcript_id" in df_biomart.columns:
    bm_versioned = df_biomart.set_index("versioned_transcript_id")
else:
    bm_versioned = pd.DataFrame()

bm_base = (
    df_biomart.set_index("base_transcript_id")
    if "base_transcript_id" in df_biomart.columns
    else pd.DataFrame()
)

# ── Join and build output rows ────────────────────────────────
resolved_rows: list[dict] = []
ambig_rows: list[dict] = []
missing: list[str] = []

for _, row in df_ensembl_ids.iterrows():
    tid = str(row["transcript_id"])
    base_id = str(row["base_id"])

    # Try versioned match first, then base ID match
    hits = pd.DataFrame()
    if not bm_versioned.empty and tid in bm_versioned.index:
        hits = bm_versioned.loc[[tid]].reset_index()
    elif not bm_base.empty and base_id in bm_base.index:
        hits = bm_base.loc[[base_id]].reset_index()

    if hits.empty:
        log.warning(f"  {tid!r} not found in any BioMart table")
        missing.append(tid)
        continue

    # Check ambiguity: same transcript → multiple genes
    unique_genes = hits["gene_id"].unique() if "gene_id" in hits.columns else []
    is_ambiguous = len(unique_genes) > 1

    if is_ambiguous:
        log.info(
            f"  Ambiguous: {tid} maps to {len(unique_genes)} genes — picking primary"
        )

    # Primary = first row (BioMart returns canonical first for most species)
    primary = hits.iloc[0]

    # Get organism from species map
    sp_match = df_species_map[df_species_map["transcript_id"] == tid]
    organism = (
        sp_match["species"].iloc[0].replace("_", " ") if len(sp_match) > 0 else ""
    )
    assembly = primary.get("_build", "")

    resolved_rows.append(
        {
            "transcript_id": tid,
            "db_source": "ensembl",
            "gene_id": primary.get("gene_id", ""),
            "gene_symbol": primary.get("gene_symbol", ""),
            "organism": organism,
            "assembly_accession": assembly,
            "chrom": str(primary.get("chrom", "")),
            "start": int(primary.get("start", 0)),
            "end": int(primary.get("end", 0)),
            "strand": primary.get("strand", "+"),
            "is_ambiguous": is_ambiguous,
        }
    )

    # Record all alternative gene mappings for traceability
    for _, alt in hits.iloc[1:].iterrows():
        if alt.get("gene_id", "") == primary.get("gene_id", ""):
            continue  # same gene, different transcript version — not truly ambiguous
        ambig_rows.append(
            {
                "transcript_id": tid,
                "db_source": "ensembl",
                "chosen_gene_id": primary.get("gene_id", ""),
                "alternative_gene_id": alt.get("gene_id", ""),
                "alternative_gene_symbol": alt.get("gene_symbol", ""),
                "organism": organism,
                "assembly_accession": assembly,
                "chrom": str(alt.get("chrom", "")),
                "start": int(alt.get("start", 0)),
                "end": int(alt.get("end", 0)),
                "strand": alt.get("strand", "+"),
            }
        )

df_resolved = pd.DataFrame(resolved_rows, columns=RESOLVED_COLS)
df_ambig = pd.DataFrame(ambig_rows, columns=AMBIG_COLS)
df_unresolved = pd.DataFrame(
    [
        {
            "transcript_id": tid,
            "db_source": "ensembl",
            "reason": "not_found_in_ensembl_biomart",
        }
        for tid in missing
    ],
    columns=UNRESOLVED_COLS,
)

df_resolved.to_csv(out_resolved, sep="\t", index=False)
df_ambig.to_csv(out_ambiguous, sep="\t", index=False)
df_unresolved.to_csv(out_unresolved, sep="\t", index=False)

# ── Summary ──────────────────────────────────────────────────
log.info("=" * 60)
log.info(f"Ensembl IDs input             : {len(df_ensembl_ids)}")
log.info(f"Resolved via BioMart join     : {len(df_resolved)}")
log.info(f"Not found in BioMart          : {len(missing)}")
log.info(f"Ambiguous alternatives logged : {len(df_ambig)}")
if missing:
    log.warning(f"Missing IDs (first 10): {missing[:10]}")
log.info(f"Written ensembl_resolved  → {out_resolved}")
log.info(f"Written ensembl_ambiguous → {out_ambiguous}")
log.info("join_ensembl_results complete.")
