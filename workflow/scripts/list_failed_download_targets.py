"""Write targeted download sentinel paths for failed assembly downloads.

This helper reads the unresolved assemblies table and emits one
resources/cache/<cache_key>/.download_done path per failed download row.
It is meant to support rerunning only the missing or failed sentinels,
not the entire download rule.

Usage:
    python workflow/scripts/list_failed_download_targets.py \
        results/unresolved_assemblies.tsv \
        /tmp/download_targets.txt

The output file is a plain newline-delimited target list suitable for
``xargs -a`` or similar tooling.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

TARGET_TEMPLATE = "resources/cache/{cache_key}/.download_done"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Write download sentinel targets for failed assemblies."
    )
    parser.add_argument(
        "input_tsv",
        help="Path to results/unresolved_assemblies.tsv or a compatible table.",
    )
    parser.add_argument(
        "output_file",
        help="Path to write newline-delimited .download_done targets.",
    )
    parser.add_argument(
        "--reason-prefix",
        default="failed:",
        help=(
            "Only rows whose reason starts with this prefix are targeted. "
            "Defaults to all failure rows."
        ),
    )
    parser.add_argument(
        "--reason-column",
        default="reason",
        help="Column containing the row classification / download reason.",
    )
    parser.add_argument(
        "--key-column",
        default="assembly_accession",
        help=(
            "Column containing the cache key / accession. Defaults to assembly_accession."
        ),
    )
    return parser


def _select_key_column(df: pd.DataFrame, preferred: str) -> str:
    if preferred in df.columns:
        return preferred
    for candidate in ("cache_key", "assembly_accession"):
        if candidate in df.columns:
            return candidate
    raise ValueError(
        "input table must contain a cache key column such as 'assembly_accession' or 'cache_key'"
    )


def build_targets(
    df: pd.DataFrame, *, reason_prefix: str, key_column: str, reason_column: str
) -> list[str]:
    if reason_column not in df.columns:
        raise ValueError(f"input table missing required '{reason_column}' column")

    selected_key_column = _select_key_column(df, key_column)

    reasons = df[reason_column].astype("string").fillna("").str.strip()
    mask = reasons.str.startswith(reason_prefix)
    keys = df.loc[mask, selected_key_column].dropna().astype(str).str.strip()
    targets = [TARGET_TEMPLATE.format(cache_key=key) for key in keys if key]
    return list(dict.fromkeys(targets))


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    input_path = Path(args.input_tsv)
    output_path = Path(args.output_file)

    df = pd.read_csv(input_path, sep="\t")
    targets = build_targets(
        df,
        reason_prefix=args.reason_prefix,
        key_column=args.key_column,
        reason_column=args.reason_column,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(targets) + ("\n" if targets else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
