"""
scripts/merge_resolved.py
Merge All Resolution Streams
============================
Concatenates ncbi_ucsc_resolved.tsv, ensembl_resolved.tsv and plant-specific
resolution streams into
the single resolved_ids.tsv that all downstream rules consume.
Does the same for ambiguous records.

Produces two explicit unresolved reports:
    - pattern_unmatched.tsv   (IDs that match no known pattern)
    - matched_not_found.tsv   (IDs matched to a route but not found in DB)
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
in_external_res = snakemake.input.external_resolved
in_biomart_res = snakemake.input.biomart_resolved
in_plant_gtf_res = snakemake.input.plant_gtf_resolved
in_gramene_res = snakemake.input.gramene_resolved
in_noncode_res = snakemake.input.noncode_resolved
in_noncode_v4_res = snakemake.input.noncode_v4_resolved
in_noncode_2016_res = snakemake.input.noncode_2016_resolved
in_abandoned_res = snakemake.input.abandoned_resolved
in_ncbi_ucsc_amb = snakemake.input.ncbi_ucsc_ambiguous
in_ensembl_amb = snakemake.input.ensembl_ambiguous
in_external_amb = snakemake.input.external_ambiguous
in_unknown = snakemake.input.unknown_ids
in_ncbi_ucsc_unres = snakemake.input.ncbi_ucsc_unresolved
in_ensembl_unres = snakemake.input.ensembl_unresolved
in_gramene_unres = snakemake.input.gramene_unresolved
in_noncode_v4_unres = snakemake.input.noncode_v4_unresolved
in_noncode_2016_unres = snakemake.input.noncode_2016_unresolved
in_noncode_unres = snakemake.input.noncode_2016_unresolved
out_resolved = snakemake.output.resolved
out_ambiguous = snakemake.output.ambiguous
out_unresolved = snakemake.output.unresolved
out_unmatched = snakemake.output.unmatched
out_not_found = snakemake.output.not_found

log.info("merge_resolved: combining RefSeq/Ensembl/plant resolution outputs")


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
df_external = safe_read(in_external_res, "external_resolved")
df_biomart    = safe_read(in_biomart_res,    "biomart_resolved")
df_plant_gtf  = safe_read(in_plant_gtf_res, "plant_gtf_resolved")
df_gramene    = safe_read(in_gramene_res,   "gramene_resolved")
df_noncode = safe_read(in_noncode_res, "noncode_resolved")
df_noncode_v4 = safe_read(in_noncode_v4_res, "noncode_v4_resolved")
df_noncode_2016 = safe_read(in_noncode_2016_res, "noncode_2016_resolved")
df_abandoned = safe_read(in_abandoned_res, "abandoned_resolved")
df_amb_nu = safe_read(in_ncbi_ucsc_amb, "ncbi_ucsc_ambiguous")
df_amb_ens = safe_read(in_ensembl_amb, "ensembl_ambiguous")
df_amb_ext = safe_read(in_external_amb, "external_ambiguous")
df_unknown = safe_read(in_unknown, "unknown_ids")
df_nu_unres = safe_read(in_ncbi_ucsc_unres, "ncbi_ucsc_unresolved")
df_ens_unres = safe_read(in_ensembl_unres, "ensembl_unresolved")
df_gram_unres = safe_read(in_gramene_unres, "gramene_unresolved")
df_noncode_v4_unres = safe_read(in_noncode_v4_unres, "noncode_v4_unresolved")
df_noncode_2016_unres = safe_read(in_noncode_2016_unres, "noncode_2016_unresolved")
df_noncode_unres = df_noncode_2016_unres  # final unresolved after all fallbacks

# ── Concatenate ───────────────────────────────────────────────
df_all_resolved = pd.concat(
    [df_ncbi_ucsc, df_ensembl, df_external, df_biomart, df_plant_gtf, df_gramene, df_noncode, df_noncode_v4, df_noncode_2016, df_abandoned],
    ignore_index=True,
)
df_all_ambig = pd.concat([df_amb_nu, df_amb_ens, df_amb_ext], ignore_index=True)
df_pattern_unmatched = df_unknown.copy()

# Normalize matched-not-found columns across resolvers
normalized_not_found = []
for frame in (df_nu_unres, df_ens_unres, df_gram_unres, df_noncode_unres):
    if frame.empty:
        continue
    cols = set(frame.columns)
    if {"transcript_id", "db_source", "reason"}.issubset(cols):
        normalized_not_found.append(frame[["transcript_id", "db_source", "reason"]].copy())
    elif {"transcript_id", "reason"}.issubset(cols):
        tmp = frame[["transcript_id", "reason"]].copy()
        tmp["db_source"] = "plant"
        normalized_not_found.append(tmp[["transcript_id", "db_source", "reason"]])
    elif {"transcript_id", "inferred_species", "reason"}.issubset(cols):
        tmp = frame[["transcript_id", "reason"]].copy()
        tmp["db_source"] = "plant"
        normalized_not_found.append(tmp[["transcript_id", "db_source", "reason"]])

df_matched_not_found = (
    pd.concat(normalized_not_found, ignore_index=True)
    if normalized_not_found
    else pd.DataFrame(columns=["transcript_id", "db_source", "reason"])
)

df_all_unresolved = pd.concat([df_pattern_unmatched, df_matched_not_found], ignore_index=True)

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
df_all_unresolved.to_csv(out_unresolved, sep="\t", index=False)
df_pattern_unmatched.to_csv(out_unmatched, sep="\t", index=False)
df_matched_not_found.to_csv(out_not_found, sep="\t", index=False)

# ── Summary ──────────────────────────────────────────────────
log.info("=" * 60)
log.info(f"NCBI/UCSC resolved   : {len(df_ncbi_ucsc)}")
log.info(f"Ensembl resolved     : {len(df_ensembl)}")
log.info(f"External resolved    : {len(df_external)}")
log.info(f"BioMart resolved     : {len(df_biomart)}")
log.info(f"Plant GTF resolved   : {len(df_plant_gtf)}")
log.info(f"Gramene resolved     : {len(df_gramene)}")
log.info(f"NONCODE resolved     : {len(df_noncode)}")
log.info(f"NONCODEv4 resolved   : {len(df_noncode_v4)}")
log.info(f"NONCODE2016 resolved : {len(df_noncode_2016)}")
log.info(f"Abandoned resolved   : {len(df_abandoned)}")
log.info(f"Total resolved       : {len(df_all_resolved)}")
log.info(f"Total ambiguous alts : {len(df_all_ambig)}")
log.info(f"Pattern unmatched    : {len(df_pattern_unmatched)}")
log.info(f"Matched not found    : {len(df_matched_not_found)}")
log.info(f"Total unresolved     : {len(df_all_unresolved)}")
if not df_all_resolved.empty and "db_source" in df_all_resolved.columns:
    for src, grp in df_all_resolved.groupby("db_source"):
        log.info(f"  {src:<12}: {len(grp)} resolved")
log.info(f"Written resolved_ids.tsv → {out_resolved}")
log.info(f"Written ambiguous.tsv    → {out_ambiguous}")
log.info(f"Written unresolved.tsv   → {out_unresolved}")
log.info(f"Written pattern_unmatched.tsv → {out_unmatched}")
log.info(f"Written matched_not_found.tsv → {out_not_found}")
log.info("merge_resolved complete.")
