"""
scripts/detect_ensembl_species.py
Detect Ensembl Species from Classified IDs
==========================================
Reads classified_ids.tsv, isolates all Ensembl transcript IDs,
and maps each to a (species, build) pair using the prefix→species
table defined in config["ensembl_species"].

This is the checkpoint that allows the BioMart lookup rules to
fan out dynamically — one rule invocation per detected species.

Outputs
-------
ensembl_species_map.tsv     : transcript_id | prefix | species | build
ensembl_unmatched_prefix.tsv: transcript_id | raw_header | detected_prefix | reason
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("detect_ensembl_species", snakemake.log[0])
input_tsv = snakemake.input.classified
out_map = snakemake.output.species_map
out_unmat = snakemake.output.unmatched
cfg = snakemake.config

# Config: {prefix: {species: ..., build: ...}}
PREFIX_MAP: dict[str, dict] = cfg.get("ensembl_species", {})

log.info("Detecting Ensembl species from classified transcript IDs")
log.info(f"Configured prefixes: {list(PREFIX_MAP.keys())}")

df = pd.read_csv(input_tsv, sep="\t")
df_ensembl = df[df["db_source"] == "ensembl"].copy()

log.info(f"Total classified IDs  : {len(df)}")
log.info(f"Ensembl IDs to map    : {len(df_ensembl)}")

if df_ensembl.empty:
    log.warning("No Ensembl IDs found — writing empty outputs")
    pd.DataFrame(columns=["transcript_id", "prefix", "species", "build"]).to_csv(
        out_map, sep="\t", index=False
    )
    pd.DataFrame(
        columns=["transcript_id", "raw_header", "detected_prefix", "reason"]
    ).to_csv(out_unmat, sep="\t", index=False)
    log.info("Stage: detect_ensembl_species complete (no Ensembl IDs)")
    sys.exit(0)


def detect_prefix(transcript_id: str) -> str | None:
    """
    Match a transcript ID against the configured prefix table.
    Tries longest prefix first to avoid ENST matching ENSMUST etc.
    Returns the matching prefix key, or None.
    """
    for prefix in sorted(PREFIX_MAP.keys(), key=len, reverse=True):
        if transcript_id.upper().startswith(prefix.upper()):
            return prefix
    return None


map_rows: list[dict] = []
unmat_rows: list[dict] = []

for _, row in df_ensembl.iterrows():
    tid = str(row["transcript_id"])
    header = str(row.get("raw_header", ""))

    prefix = detect_prefix(tid)

    if prefix is None:
        log.warning(f"  No configured prefix matches {tid!r}")
        unmat_rows.append(
            {
                "transcript_id": tid,
                "raw_header": header,
                "detected_prefix": tid[:8],  # show first 8 chars as hint
                "reason": (
                    "Transcript ID prefix not found in config ensembl_species. "
                    "Add the prefix to config/config.yaml to resolve this ID."
                ),
            }
        )
        continue

    info = PREFIX_MAP[prefix]
    log.debug(
        f"  {tid!r} → prefix={prefix!r} species={info['species']} build={info['build']}"
    )
    map_rows.append(
        {
            "transcript_id": tid,
            "prefix": prefix,
            "species": info["species"],
            "build": info["build"],
        }
    )

df_map = pd.DataFrame(map_rows, columns=["transcript_id", "prefix", "species", "build"])
df_unmat = pd.DataFrame(
    unmat_rows, columns=["transcript_id", "raw_header", "detected_prefix", "reason"]
)

df_map.to_csv(out_map, sep="\t", index=False)
df_unmat.to_csv(out_unmat, sep="\t", index=False)

# ── Summary ──────────────────────────────────────────────────
species_counts = df_map["species"].value_counts()
log.info("=" * 60)
log.info(f"Ensembl IDs mapped successfully  : {len(df_map)}")
log.info(f"Ensembl IDs with unknown prefix  : {len(df_unmat)}")
log.info("Detected species:")
for sp, cnt in species_counts.items():
    log.info(f"  {sp:<40} {cnt} transcripts")
log.info(f"BioMart lookup will run for {species_counts.shape[0]} species")
log.info(f"Written species_map  → {out_map}")
log.info(f"Written unmatched    → {out_unmat}")
log.info("detect_ensembl_species complete.")
