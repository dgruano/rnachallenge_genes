"""Helpers for filling URL columns from per-assembly lookup tables."""

from __future__ import annotations

import pandas as pd


def _norm_text(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().str.lower()


def fill_urls_from_table(
    df: pd.DataFrame,
    url_table: pd.DataFrame,
    *,
    fill_cols: list[str],
    fallback_on_organism: bool = False,
) -> tuple[pd.DataFrame, int]:
    """Fill NA URL/accession fields from url_table.

    Primary join key is assembly_name. Optionally applies a second pass on
    organism for rows still missing values.
    """
    if url_table.empty or df.empty:
        return df, 0

    usable_fill_cols = [
        c for c in fill_cols if c in url_table.columns and c in df.columns
    ]
    if not usable_fill_cols:
        return df, 0

    merged = df.copy()
    before = merged["fasta_url"].notna().sum() if "fasta_url" in merged.columns else 0

    if "assembly_name" in merged.columns and "assembly_name" in url_table.columns:
        slim = (
            url_table[["assembly_name"] + usable_fill_cols]
            .dropna(subset=["assembly_name"])
            .assign(_asm_key=lambda x: _norm_text(x["assembly_name"]))
            .drop(columns=["assembly_name"])
            .drop_duplicates("_asm_key")
            .rename(columns={c: f"_fill_{c}" for c in usable_fill_cols})
        )
        merged = merged.assign(_asm_key=_norm_text(merged["assembly_name"]))
        merged = merged.merge(slim, on="_asm_key", how="left")
        for col in usable_fill_cols:
            tmp = f"_fill_{col}"
            if tmp in merged.columns:
                merged[col] = merged[col].combine_first(merged[tmp])
                merged = merged.drop(columns=[tmp])
        merged = merged.drop(columns=["_asm_key"])

    if (
        fallback_on_organism
        and "organism" in merged.columns
        and "organism" in url_table.columns
    ):
        slim = (
            url_table[["organism"] + usable_fill_cols]
            .dropna(subset=["organism"])
            .assign(_org_key=lambda x: _norm_text(x["organism"]))
            .drop(columns=["organism"])
            .drop_duplicates("_org_key")
            .rename(columns={c: f"_orgfill_{c}" for c in usable_fill_cols})
        )
        merged = merged.assign(_org_key=_norm_text(merged["organism"]))
        merged = merged.merge(slim, on="_org_key", how="left")
        for col in usable_fill_cols:
            tmp = f"_orgfill_{col}"
            if tmp in merged.columns:
                merged[col] = merged[col].combine_first(merged[tmp])
                merged = merged.drop(columns=[tmp])
        merged = merged.drop(columns=["_org_key"])

    after = (
        merged["fasta_url"].notna().sum() if "fasta_url" in merged.columns else before
    )
    return merged, int(after - before)
