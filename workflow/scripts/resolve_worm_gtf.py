"""
Direct GTF-based resolver for WormBase transcript IDs.
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
    resolve_classified_ids,
    wormbase_candidates,
)

log = get_logger("resolve_worm_gtf", snakemake.log[0])

input_tsv = snakemake.input.classified
gtf_file = snakemake.input.gtf
out_resolved = snakemake.output.resolved
out_unresolved = snakemake.output.unresolved

cfg = snakemake.config.get("metazoa_gtf_sources", {}).get("wormbase", {})
ASSEMBLY_NAME = cfg.get("assembly_name", "WBcel235")
ASSEMBLY_ACCESSION = cfg.get("assembly_accession")
FASTA_URL = cfg.get("fasta_url")
GTF_URL = cfg.get("url")
GTF_FORMAT = cfg.get("gtf_format", "gtf")
ORGANISM = cfg.get("organism", "caenorhabditis_elegans")

WORMBASE_COORD_RE = re.compile(
    r"^(?P<transcript>[^_]+)_wormbase:known_chromosome:(?P<assembly>WBcel\d+):"
    r"(?P<chrom>[^:]+):(?P<start>\d+):(?P<end>\d+):.*gene:(?P<gene_id>WBGene\d+)",
    re.IGNORECASE,
)
TAIR_COORD_RE = re.compile(
    r"chromosome:(?P<assembly>TAIR10):(?P<chrom>[^:]+):(?P<start>\d+):(?P<end>\d+):(?P<strand>[-+_]?1)"
    r"\s+gene:(?P<gene_id>AT\dG\d+)",
    re.IGNORECASE,
)


def pre_resolve(row):
    transcript_id = str(row["transcript_id"])
    raw_header = str(row.get("raw_header", ""))

    if transcript_id.startswith("AT"):
        match = TAIR_COORD_RE.search(raw_header)
        if match:
            strand = normalize_strand(match.group("strand"))
            return {
                "transcript_id": transcript_id,
                "db_source": "plant",
                "gene_id": match.group("gene_id"),
                "gene_symbol": match.group("gene_id"),
                "organism": "arabidopsis_thaliana",
                "assembly_accession": match.group("assembly"),
                "chrom": match.group("chrom"),
                "start": int(match.group("start")),
                "end": int(match.group("end")),
                "strand": strand,
                "is_ambiguous": False,
            }

    match = WORMBASE_COORD_RE.match(transcript_id)
    if match:
        return {
            "transcript_id": transcript_id,
            "db_source": "wormbase",
            "gene_id": match.group("gene_id"),
            "gene_symbol": match.group("transcript"),
            "organism": ORGANISM,
            "assembly_accession": match.group("assembly"),
            "chrom": match.group("chrom"),
            "start": int(match.group("start")),
            "end": int(match.group("end")),
            "strand": ".",
            "is_ambiguous": False,
        }
    return None


log.info("Loading classified WormBase IDs from classified_ids.tsv")
classified_df = pd.read_csv(input_tsv, sep="\t")
worm_df = classified_df[classified_df["db_source"].astype(str) == "wormbase"].copy()
log.info(f"WormBase IDs to resolve via direct annotation: {len(worm_df)}")

if worm_df.empty:
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    log.info("No WormBase IDs to resolve; wrote empty outputs.")
    sys.exit(0)

log.info(f"Building WormBase annotation index from {gtf_file}")
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
    db_source="wormbase",
    organism=ORGANISM,
    assembly_name=ASSEMBLY_NAME,
    assembly_accession=ASSEMBLY_ACCESSION,
    fasta_url=FASTA_URL,
    gtf_url=GTF_URL,
    gtf_format=GTF_FORMAT,
    index=index,
    candidates_fn=wormbase_candidates,
    unresolved_reason="worm_gtf_not_resolved",
    pre_resolve_fn=pre_resolve,
)

res_df.to_csv(out_resolved, sep="\t", index=False)
unres_df.to_csv(out_unresolved, sep="\t", index=False)

log.info("=" * 60)
log.info(f"WormBase IDs resolved via annotation : {len(res_df)}")
log.info(f"WormBase IDs still unresolved       : {len(unres_df)}")
log.info(f"Written worm_gtf_resolved.tsv   → {out_resolved}")
log.info(f"Written worm_gtf_unresolved.tsv → {out_unresolved}")
log.info("resolve_worm_gtf complete.")
