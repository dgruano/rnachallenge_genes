"""
Direct GTF-based resolver for FlyBase transcript IDs.
"""

import re
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
    normalize_strand,
    generic_candidates,
    resolve_classified_ids,
)

log = get_logger("resolve_fly_gtf", snakemake.log[0])

input_tsv = snakemake.input.classified
gtf_file = snakemake.input.gtf
out_resolved = snakemake.output.resolved
out_unresolved = snakemake.output.unresolved

cfg = snakemake.config.get("metazoa_gtf_sources", {}).get("flybase", {})
ASSEMBLY_NAME = cfg.get("assembly_name", "BDGP6")
ORGANISM = cfg.get("organism", "drosophila_melanogaster")

FLYBASE_COORD_RE = re.compile(
    r"chromosome:(?P<assembly>[^:]+):(?P<chrom>[^:]+):(?P<start>\d+):(?P<end>\d+):(?P<strand>[-+_]?1)"
    r"\s+gene:(?P<gene_id>FBgn\d+).*?gene_symbol:(?P<gene_symbol>[^\s]+)",
    re.IGNORECASE,
)


def pre_resolve(row):
    transcript_id = str(row["transcript_id"])
    raw_header = str(row.get("raw_header", ""))
    match = FLYBASE_COORD_RE.search(raw_header)
    if not match:
        return None

    return {
        "transcript_id": transcript_id,
        "db_source": "flybase",
        "gene_id": match.group("gene_id"),
        "gene_symbol": match.group("gene_symbol"),
        "organism": ORGANISM,
        "assembly_accession": match.group("assembly"),
        "chrom": match.group("chrom"),
        "start": int(match.group("start")),
        "end": int(match.group("end")),
        "strand": normalize_strand(match.group("strand")),
        "is_ambiguous": False,
    }


log.info("Loading classified FlyBase IDs from classified_ids.tsv")
classified_df = pd.read_csv(input_tsv, sep="\t")
fly_df = classified_df[classified_df["db_source"].astype(str) == "flybase"].copy()
log.info(f"FlyBase IDs to resolve via direct annotation: {len(fly_df)}")

if fly_df.empty:
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    log.info("No FlyBase IDs to resolve; wrote empty outputs.")
    sys.exit(0)

log.info(f"Building FlyBase annotation index from {gtf_file}")
index = build_annotation_index(
    gtf_file,
    feature_types={"transcript"},
    transcript_fields=("transcript_id", "transcript_name"),
    gene_id_fields=("gene_id",),
    gene_symbol_fields=("gene_name",),
    log=log,
)

res_df, unres_df = resolve_classified_ids(
    classified_df,
    db_source="flybase",
    organism=ORGANISM,
    assembly_name=ASSEMBLY_NAME,
    index=index,
    candidates_fn=generic_candidates,
    unresolved_reason="fly_gtf_not_resolved",
    pre_resolve_fn=pre_resolve,
)

res_df.to_csv(out_resolved, sep="\t", index=False)
unres_df.to_csv(out_unresolved, sep="\t", index=False)

log.info("=" * 60)
log.info(f"FlyBase IDs resolved via annotation : {len(res_df)}")
log.info(f"FlyBase IDs still unresolved       : {len(unres_df)}")
log.info(f"Written fly_gtf_resolved.tsv   → {out_resolved}")
log.info(f"Written fly_gtf_unresolved.tsv → {out_unresolved}")
log.info("resolve_fly_gtf complete.")
