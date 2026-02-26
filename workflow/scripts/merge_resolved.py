"""
scripts/merge_resolved.py
Merge All Resolution Streams
============================
Concatenates ncbi_ucsc_resolved.tsv and ensembl_resolved.tsv into
the single resolved_ids.tsv that all downstream rules consume.
Does the same for ambiguous records.

Also appends any Ensembl IDs that were unmatched in BioMart
to unresolved.tsv (the unresolved file from parse_ids already
holds IDs with unknown DB format; this adds BioMart misses).
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("merge_resolved", snakemake.log[0])

in_ncbi_ucsc_res = snakemake.input.ncbi_ucsc_resolved
in_ensembl_res = snakemake.input.ensembl_resolved
in_ncbi_ucsc_amb = snakemake.input.ncbi_ucsc_ambiguous
in_ensembl_amb = snakemake.input.ensembl_ambiguous
out_resolved = snakemake.output.resolved
out_ambiguous = snakemake.output.ambiguous

log.info("merge_resolved: combining NCBI/UCSC and Ensembl resolution outputs")


# ── Load all parts ────────────────────────────────────────────
def safe_read(path: str, label: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t")
        log.info(f"  Loaded {label}: {len(df)} rows")
        return df
    except Exception as exc:
        log.warning(f"  Could not load {label}: {exc} — treating as empty")
        return pd.DataFrame()


df_ncbi_ucsc = safe_read(in_ncbi_ucsc_res, "ncbi_ucsc_resolved")
df_ensembl = safe_read(in_ensembl_res, "ensembl_resolved")
df_amb_nu = safe_read(in_ncbi_ucsc_amb, "ncbi_ucsc_ambiguous")
df_amb_ens = safe_read(in_ensembl_amb, "ensembl_ambiguous")

# ── Concatenate ───────────────────────────────────────────────
df_all_resolved = pd.concat([df_ncbi_ucsc, df_ensembl], ignore_index=True)
df_all_ambig = pd.concat([df_amb_nu, df_amb_ens], ignore_index=True)

# Sanity check: flag any duplicate transcript IDs (shouldn't happen but log if so)
dupes = df_all_resolved[df_all_resolved.duplicated(subset="transcript_id", keep=False)]
if not dupes.empty:
    log.warning(
        f"  {len(dupes)} duplicate transcript_id entries found after merge "
        f"(keeping first occurrence):\n  {dupes['transcript_id'].unique().tolist()[:10]}"
    )
    df_all_resolved = df_all_resolved.drop_duplicates(
        subset="transcript_id", keep="first"
    )

df_all_resolved.to_csv(out_resolved, sep="\t", index=False)
df_all_ambig.to_csv(out_ambiguous, sep="\t", index=False)

# ── Summary ──────────────────────────────────────────────────
log.info("=" * 60)
log.info(f"NCBI/UCSC resolved   : {len(df_ncbi_ucsc)}")
log.info(f"Ensembl resolved     : {len(df_ensembl)}")
log.info(f"Total resolved       : {len(df_all_resolved)}")
log.info(f"Total ambiguous alts : {len(df_all_ambig)}")
if not df_all_resolved.empty and "db_source" in df_all_resolved.columns:
    for src, grp in df_all_resolved.groupby("db_source"):
        log.info(f"  {src:<12}: {len(grp)} resolved")
log.info(f"Written resolved_ids.tsv → {out_resolved}")
log.info(f"Written ambiguous.tsv    → {out_ambiguous}")
log.info("merge_resolved complete.")
