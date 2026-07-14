"""
scripts/resolve_phytozome_gtf.py
GFF3-based resolver for Phytozome-backed plant transcript IDs
==============================================================
Consumes plant / phytozome-classified IDs from classified_ids.tsv,
indexes configured Phytozome GFF3 files by mRNA transcript ID, and
resolves transcript -> gene/coordinates per species.

Runs under Snakemake OR standalone via argparse:
  snakemake resolve_phytozome_gtf  (Snakemake mode)
  python resolve_phytozome_gtf.py --classified <tsv> --gff-files <paths...> ...  (CLI mode)
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from logging_utils import get_logger
from resolution_guard import check_match_rates
from utils.annotation_resolver import (
    UNRESOLVED_COLS,
    build_annotation_index,
    generic_candidates,
)

# ── Auto-detect mode: Snakemake vs CLI ────────────────────────
_is_snakemake = "snakemake" in dir()

if not _is_snakemake:
    parser = argparse.ArgumentParser(
        description="Resolve Phytozome GFF3-backed plant transcript IDs"
    )
    parser.add_argument(
        "--classified",
        required=True,
        help="Path to classified_ids.tsv input file",
    )
    parser.add_argument(
        "--gff-files",
        nargs="+",
        required=True,
        help="List of GFF3 file paths (order must match species in config)",
    )
    parser.add_argument(
        "--output-resolved",
        required=True,
        help="Output path for resolved IDs",
    )
    parser.add_argument(
        "--output-unresolved",
        required=True,
        help="Output path for unresolved IDs",
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to phytozome_gtf_sources config file (YAML or JSON)",
    )
    parser.add_argument(
        "--log",
        default=None,
        help="Log file path (optional)",
    )

    args = parser.parse_args()

    # Create a mock snakemake namespace from args
    class _Input:
        def __init__(self, classified, gff_files):
            self.classified = classified
            self.gff_files = gff_files

    class _Output:
        def __init__(self, resolved, unresolved):
            self.resolved = resolved
            self.unresolved = unresolved

        def __getitem__(self, idx):
            if idx == 0:
                return self.resolved
            raise IndexError(f"Output index {idx} out of range")

    class _Snakemake:
        def __init__(
            self,
            classified,
            gff_files,
            output_resolved,
            output_unresolved,
            config_path,
            log_path,
        ):
            self.input = _Input(classified, gff_files)
            self.output = _Output(output_resolved, output_unresolved)
            self.log = [log_path] if log_path else ["/dev/null"]

            # Load config from YAML or JSON
            with open(config_path) as f:
                if config_path.endswith(".json"):
                    data = json.load(f)
                else:
                    data = yaml.safe_load(f)
            self.config = data

    snakemake = _Snakemake(
        args.classified,
        args.gff_files,
        args.output_resolved,
        args.output_unresolved,
        args.config,
        args.log,
    )

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("resolve_phytozome_gtf", snakemake.log[0])

input_tsv = snakemake.input.classified
gff_files = snakemake.input.gff_files
out_resolved = snakemake.output.resolved
out_unresolved = snakemake.output.unresolved


def _phytozome_sources(cfg: dict) -> Dict[str, dict]:
    nested = cfg.get("phytozome_gtf_sources")
    if isinstance(nested, dict):
        return nested

    return {
        key: value
        for key, value in cfg.items()
        if isinstance(value, dict)
        and any(
            field in value
            for field in ("species_query", "genome_id", "gtf", "phytozome_version")
        )
    }


GTF_SOURCES: Dict[str, dict] = _phytozome_sources(snakemake.config)
SPECIES_ORDER = list(GTF_SOURCES.keys())
SPECIES_TO_GFF: Dict[str, str] = {
    sp: gff_files[i] for i, sp in enumerate(SPECIES_ORDER) if i < len(gff_files)
}

RESOLVED_COLS = [
    "transcript_id",
    "db_source",
    "gene_id",
    "gene_symbol",
    "organism",
    "assembly_accession",
    "assembly_name",
    "chrom",
    "start",
    "end",
    "strand",
    "phytozome_version",
    "is_ambiguous",
]


# ── Prefix → species ─────────────────────────────────────────
PREFIX_TO_SPECIES = {
    "orange": "citrus_sinensis",
    "Sb": "sorghum_bicolor",
    "Sobic": "sorghum_bicolor",
    "PGSC": "solanum_tuberosum",
    "Glyma": "glycine_max",
    "Bradi": "brachypodium_distachyon",
    "Cre": "chlamydomonas_reinhardtii",
    "Pp": "physcomitrella_patens",
    "Medtr": "medicago_truncatula",
    "Si": "setaria_italica",
    "Thecc1EG": "theobroma_cacao",
    "Potri": "populus_trichocarpa",
    "VIT_": "vitis_vinifera",
    "GSVIVT": "vitis_vinifera",
    "GTVIVG": "vitis_vinifera",
}


def infer_species(row: pd.Series) -> Optional[str]:
    species_hint = str(row.get("species_hint", "")).strip()
    if species_hint in GTF_SOURCES:
        return species_hint

    transcript_id = str(row.get("transcript_id", ""))
    for prefix, species in PREFIX_TO_SPECIES.items():
        if transcript_id.startswith(prefix):
            return species
    return None


def phytozome_candidates(tid: str) -> List[str]:
    candidates = list(generic_candidates(tid))

    if "." in tid:
        parts = tid.split(".")
        while len(parts) > 1:
            parts = parts[:-1]
            candidates.append(".".join(parts))

    if tid.endswith("_cdna"):
        candidates.append(tid[: -len("_cdna")])

    return [candidate for candidate in dict.fromkeys(candidates) if candidate]


def build_phytozome_index(gff_path: str) -> Dict[str, dict]:
    return build_annotation_index(
        gff_path,
        feature_types=("mRNA",),
        transcript_fields=("ID", "Name"),
        gene_id_fields=("Parent", "gene", "locusName"),
        gene_symbol_fields=("gene_name", "gene_symbol", "locusName"),
        alias_fields=(),
        log=log,
    )


# ── Main ─────────────────────────────────────────────────────
log.info("Loading classified IDs from classified_ids.tsv")
classified_df = pd.read_csv(input_tsv, sep="\t", low_memory=False)

configured_species = set(GTF_SOURCES)
df = classified_df[
    (classified_df["db_source"].astype(str).isin(["plant", "phytozome"]))
    & (classified_df["species_hint"].astype(str).isin(configured_species))
].copy()

if df.empty:
    fallback_df = classified_df[
        classified_df["db_source"].astype(str).isin(["plant", "phytozome"])
    ].copy()
    fallback_df["inferred_species"] = fallback_df.apply(infer_species, axis=1)
    df = fallback_df[
        fallback_df["inferred_species"].astype(str).isin(configured_species)
    ].copy()
else:
    df["inferred_species"] = df.apply(infer_species, axis=1)

log.info(f"Phytozome IDs to resolve via GFF3: {len(df)}")

if df.empty or not GTF_SOURCES:
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    log.info("Nothing to resolve; wrote empty outputs.")
    sys.exit(0)

needed_species = set(df["inferred_species"].dropna().unique()) & set(SPECIES_TO_GFF)
log.info(f"Building GFF3 indices for: {sorted(needed_species)}")

gff_indices: Dict[str, Dict[str, dict]] = {}
for species in sorted(needed_species):
    gff_path = SPECIES_TO_GFF[species]
    log.info(f"  Parsing Phytozome GFF3 for {species}: {gff_path}")
    gff_indices[species] = build_phytozome_index(gff_path)

resolved_rows: List[dict] = []
unresolved_rows: List[dict] = []

for _, row in df.iterrows():
    tid = str(row["transcript_id"])
    species = row.get("inferred_species")

    if not species or species not in GTF_SOURCES:
        unresolved_rows.append(
            {
                "transcript_id": tid,
                "raw_header": str(row.get("raw_header", "")),
                "source_file": str(row.get("source_file", "")),
                "db_source": "phytozome",
                "reason": "no_phytozome_species_match",
            }
        )
        continue

    if species not in gff_indices:
        unresolved_rows.append(
            {
                "transcript_id": tid,
                "raw_header": str(row.get("raw_header", "")),
                "source_file": str(row.get("source_file", "")),
                "db_source": "phytozome",
                "reason": f"no_phytozome_gff_for_{species}",
            }
        )
        continue

    index = gff_indices[species]
    src_meta = GTF_SOURCES[species]
    hit = None

    for candidate in phytozome_candidates(tid):
        if candidate in index:
            hit = index[candidate]
            break

    if hit is None:
        unresolved_rows.append(
            {
                "transcript_id": tid,
                "raw_header": str(row.get("raw_header", "")),
                "source_file": str(row.get("source_file", "")),
                "db_source": "phytozome",
                "reason": f"phytozome_gff_no_match_{species}",
            }
        )
        continue

    resolved_rows.append(
        {
            "transcript_id": tid,
            "db_source": "phytozome",
            "gene_id": hit.get("gene_id", ""),
            "gene_symbol": hit.get("gene_symbol", "") or hit.get("gene_id", ""),
            "organism": species,
            # Per-species cache key: extract maps assembly_accession -> cache dir
            # (resources/cache/<acc>/genome.fasta), so a shared "Phytozome" would
            # collide across species. download_phytozome_fasta caches the genome
            # FASTA under this same key.
            "assembly_accession": f"phytozome_{species}",
            "assembly_name": src_meta.get(
                "assembly_name", src_meta.get("assembly", "")
            ),
            "chrom": hit.get("chrom", ""),
            "start": hit.get("start", ""),
            "end": hit.get("end", ""),
            "strand": hit.get("strand", "."),
            "phytozome_version": str(
                src_meta.get("phytozome_version", src_meta.get("version", ""))
            ),
            "is_ambiguous": False,
        }
    )

res_df = pd.DataFrame(resolved_rows, columns=RESOLVED_COLS)
unres_df = pd.DataFrame(unresolved_rows, columns=UNRESOLVED_COLS)

res_df.to_csv(out_resolved, sep="\t", index=False)
unres_df.to_csv(out_unresolved, sep="\t", index=False)

log.info("Phytozome GFF3 resolution complete:")
log.info(f"  Resolved:         {len(res_df)}")
log.info(f"  Still unresolved: {len(unres_df)}")
if not res_df.empty:
    for species, group in res_df.groupby("organism"):
        src = GTF_SOURCES.get(species, {})
        log.info(
            f"  {species}: {len(group)} transcripts "
            f"[{src.get('assembly_name', '?')} / Phytozome {src.get('phytozome_version', '?')}]"
        )

# ── Presence verification (namespace→assembly correctness guard) ──
strict = bool(snakemake.config.get("plant_resolution_strict", True))
min_rate = float(snakemake.config.get("plant_resolution_min_match_rate", 0.02))
attempted_counts = (
    df[df["inferred_species"].isin(gff_indices)]
    .groupby("inferred_species")
    .size()
    .to_dict()
)
matched_counts = res_df.groupby("organism").size().to_dict() if not res_df.empty else {}
failures = check_match_rates(
    matched_counts, attempted_counts, min_rate=min_rate, log=log
)
if failures and strict:
    log.error(
        "Strict presence-check failed (phytozome); exiting non-zero. "
        "Set plant_resolution_strict: false in config to override."
    )
    sys.exit(1)
