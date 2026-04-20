"""
scripts/resolve_yeast_gtf.py
GFF3-based resolver for SGD (yeast) gene IDs
============================================
Consumes all classified SGD IDs directly rather than waiting for
resolve_external_ids to emit them as unresolved.

Reads the S. cerevisiae annotation GFF3 from the Saccharomyces Genome
Database (yeastgenome.org), builds a systematic-name → coordinates index,
and looks up each unresolved ID with normalised candidate variants.

GFF3 index keys
---------------
Each gene feature is indexed by:
  • ID attribute   (systematic name, e.g. YDL184C, Q0275)
  • Name attribute (gene symbol, e.g. SMC4) — if different from ID
  • Alias values   (comma-separated list in the Alias= field)

Candidate generation for each transcript_id
--------------------------------------------
  1. raw ID as-is
  2. _A / _B locus suffix stripped  (YLR264C_A  → YLR264C)
  3. underscore → dash              (YLR264C_A  → YLR264C-A)
  4. suffix-stripped + dash variant (YLR264C_A  → YLR264C-A stripped)
  5. trailing _\\d+ version removed  (YAL001C_1  → YAL001C)
  6. _cdna / _mRNA suffix stripped

Output columns match RESOLVED_COLS used throughout the pipeline so
merge_resolved can concatenate yeast_gtf_resolved.tsv with all other
resolution streams.

Output files
------------
  yeast_gtf_resolved.tsv   — SGD IDs successfully resolved from the GFF3
  yeast_gtf_unresolved.tsv — SGD IDs still unresolved (reason: sgd_gtf_not_resolved)
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from logging_utils import get_logger
from utils.annotation_resolver import (
    RESOLVED_COLS,
    UNRESOLVED_COLS,
    build_annotation_index,
    resolve_classified_ids,
    yeast_candidates,
)

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("resolve_yeast_gtf", snakemake.log[0])

input_tsv      = snakemake.input.classified   # classified_ids.tsv
gff_file       = snakemake.input.gff
out_resolved   = snakemake.output.resolved
out_unresolved = snakemake.output.unresolved

cfg = snakemake.config
_YEAST_SRC    = cfg.get("yeast_gtf_sources", {}).get("saccharomyces_cerevisiae", {})
ASSEMBLY_NAME = _YEAST_SRC.get("assembly_name", "R64")
ORGANISM      = "saccharomyces_cerevisiae"

# Feature types in the SGD GFF3 that represent gene-level entries
_GENE_TYPES = {
    "gene",
    "pseudogene",
    "transposable_element_gene",
    "blocked_reading_frame",
    "long_terminal_repeat",
}

# ── Main ─────────────────────────────────────────────────────

log.info("Loading classified SGD IDs from classified_ids.tsv")
classified_df = pd.read_csv(input_tsv, sep="\t")

sgd_df = classified_df[classified_df["db_source"].astype(str) == "sgd"].copy()
log.info(f"SGD IDs to resolve via GFF3: {len(sgd_df)}")

if sgd_df.empty:
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    log.info("No SGD IDs to resolve; wrote empty outputs.")
    sys.exit(0)

log.info(f"Building gene index from {gff_file}")
gene_index = build_annotation_index(
    gff_file,
    feature_types=_GENE_TYPES,
    transcript_fields=(),
    gene_id_fields=("ID",),
    gene_symbol_fields=("Name",),
    alias_fields=("Alias",),
    log=log,
)

res_df, unres_out = resolve_classified_ids(
    classified_df,
    db_source="sgd",
    organism=ORGANISM,
    assembly_name=ASSEMBLY_NAME,
    index=gene_index,
    candidates_fn=yeast_candidates,
    unresolved_reason="sgd_gtf_not_resolved",
)

res_df.to_csv(out_resolved,   sep="\t", index=False)
unres_out.to_csv(out_unresolved, sep="\t", index=False)

log.info("=" * 60)
log.info(f"SGD IDs resolved via GFF3  : {len(res_df)}")
log.info(f"SGD IDs still unresolved   : {len(unres_out)}")
log.info(f"Written yeast_gtf_resolved.tsv   → {out_resolved}")
log.info(f"Written yeast_gtf_unresolved.tsv → {out_unresolved}")
log.info("resolve_yeast_gtf complete.")
