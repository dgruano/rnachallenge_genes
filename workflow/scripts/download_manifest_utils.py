"""Build download manifest rows from resolved assembly metadata."""

from __future__ import annotations

import pandas as pd

from cache_key_utils import build_cache_key_from_url


def _clean_optional_string(series: pd.Series) -> pd.Series:
    cleaned = series.astype("string").str.strip()
    return cleaned.mask(cleaned.str.lower().isin(["", "nan", "none"]))


def _is_ncbi_assembly(series: pd.Series) -> pd.Series:
    return series.str.match(r"GC[FA]_\d+\.\d+", na=False)


def build_download_manifest(df: pd.DataFrame) -> pd.DataFrame:
    """Return ``cache_key,fasta_url`` rows for assemblies that can be downloaded."""
    if "assembly_accession" not in df.columns:
        raise ValueError("resolved dataframe missing 'assembly_accession' column")

    acc = _clean_optional_string(df["assembly_accession"])
    if "fasta_url" in df.columns:
        url = _clean_optional_string(df["fasta_url"])
    else:
        url = pd.Series(pd.NA, index=df.index, dtype="string")

    is_ncbi = _is_ncbi_assembly(acc)
    has_url = url.notna()
    eligible = df.loc[is_ncbi | has_url].copy()
    if eligible.empty:
        return pd.DataFrame(columns=["cache_key", "fasta_url"])

    eligible_acc = acc.loc[eligible.index]
    eligible_url = url.loc[eligible.index]
    eligible_is_ncbi = is_ncbi.loc[eligible.index]

    keys = []
    for idx in eligible.index:
        if bool(eligible_is_ncbi.loc[idx]):
            keys.append(str(eligible_acc.loc[idx]))
        else:
            keys.append(build_cache_key_from_url(str(eligible_url.loc[idx])))

    manifest = pd.DataFrame({
        "cache_key": keys,
        "fasta_url": eligible_url.values,
    })

    # Prefer entries that carry a direct URL when duplicate cache keys exist.
    manifest["_has_url"] = manifest["fasta_url"].notna()
    manifest = (
        manifest.sort_values("_has_url", ascending=False)
        .drop_duplicates(subset="cache_key", keep="first")
        .drop(columns=["_has_url"])
        .reset_index(drop=True)
    )
    return manifest
