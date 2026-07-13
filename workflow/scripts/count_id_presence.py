#!/usr/bin/env python3
"""Count species-specific IDs from a TSV that are present in a GFF3 or GTF file.

Usage:
    python count_id_presence.py <input.tsv> <species_name> <annotations.gff3|gtf[.gz]>

Expected TSV columns by default:
    - transcript_id
    - inferred_species
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def _clean_token(token: str) -> str:
    """Normalize a token extracted from GFF3 attributes."""
    return token.strip().strip('"').strip("'")


def _extract_attribute_values(item: str) -> list[str]:
    """Extract one or more values from a GFF3 or GTF attribute item."""
    item = item.strip()
    if not item:
        return []

    if "=" in item:
        _, value = item.split("=", 1)
    else:
        parts = item.split(None, 1)
        if len(parts) != 2:
            return []
        value = parts[1]

    return [_clean_token(token) for token in value.split(",") if _clean_token(token)]


def read_species_ids(
    tsv_path: Path,
    species_name: str,
    id_column: str,
    species_column: str,
) -> list[str]:
    """Read IDs from rows where species_column equals species_name."""
    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {tsv_path}")

        missing = [
            col
            for col in (id_column, species_column)
            if col not in set(reader.fieldnames)
        ]
        if missing:
            cols = ", ".join(missing)
            raise ValueError(f"Missing required TSV column(s): {cols}")

        ids = [
            row[id_column].strip()
            for row in reader
            if row[species_column].strip() == species_name and row[id_column].strip()
        ]

    return ids


def _open_text_maybe_gzip(path: Path):
    """Open text file, supporting .gz transparently."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def extract_gff_ids(gff_path: Path) -> set[str]:
    """Extract all raw attribute values from GFF3 or GTF and de-prefixed forms.

    For each attribute value token, both forms are collected:
    - token itself
    - token with any leading '<prefix>:' removed (e.g. 'transcript:ABC' -> 'ABC')
    """
    ids: set[str] = set()

    with _open_text_maybe_gzip(gff_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            attrs = parts[8]
            for item in attrs.split(";"):
                for token in _extract_attribute_values(item):
                    ids.add(token)

                    if ":" in token:
                        deprefixed = _clean_token(token.split(":", 1)[1])
                        if deprefixed:
                            ids.add(deprefixed)

    return ids


def count_presence(species_ids: list[str], gff_ids: set[str]) -> tuple[int, int, int]:
    """Return total, present and absent counts for species IDs."""
    total = len(species_ids)
    present = sum(1 for value in species_ids if value in gff_ids)
    absent = total - present
    return total, present, absent


def get_missing_ids(species_ids: list[str], gff_ids: set[str]) -> list[str]:
    """Return IDs from the TSV that are not present in the annotation file."""
    return [value for value in species_ids if value not in gff_ids]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Count how many IDs in a TSV (filtered by species) are present in a GFF3 or GTF file."
        )
    )
    parser.add_argument("tsv", type=Path, help="Input TSV path")
    parser.add_argument(
        "species", type=str, help="Species name to filter (exact match)"
    )
    parser.add_argument(
        "gff",
        type=Path,
        help="Input annotation path (.gff3, .gtf, .gff3.gz, or .gtf.gz)",
    )
    parser.add_argument(
        "--id-column",
        default="transcript_id",
        help="TSV column containing IDs (default: transcript_id)",
    )
    parser.add_argument(
        "--species-column",
        default="inferred_species",
        help="TSV column containing species labels (default: inferred_species)",
    )
    parser.add_argument(
        "--unique",
        action="store_true",
        help="Count unique IDs only (deduplicates TSV IDs before matching)",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run internal sanity checks and exit",
    )
    parser.add_argument(
        "--missing-out",
        type=Path,
        help="Optional output file for IDs absent from the annotation file",
    )
    return parser


def _self_test() -> None:
    sample_ids = ["A", "A", "B", "C"]
    sample_gff = {"A", "C", "X"}

    total, present, absent = count_presence(sample_ids, sample_gff)

    assert total == 4
    assert present == 3
    assert absent == 1

    assert _extract_attribute_values("ID=transcript:ABC123") == ["transcript:ABC123"]
    assert _extract_attribute_values("Parent=gene:XYZ,transcript:ABC123") == [
        "gene:XYZ",
        "transcript:ABC123",
    ]
    assert _extract_attribute_values('transcript_id "ABC123"') == ["ABC123"]
    assert _extract_attribute_values('gene_id "XYZ"') == ["XYZ"]
    assert get_missing_ids(sample_ids, sample_gff) == ["B"]


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.self_test:
        _self_test()
        print("self_test=ok")
        return

    if not args.tsv.exists():
        raise FileNotFoundError(f"TSV not found: {args.tsv}")
    if not args.gff.exists():
        raise FileNotFoundError(f"GFF not found: {args.gff}")

    species_ids = read_species_ids(
        tsv_path=args.tsv,
        species_name=args.species,
        id_column=args.id_column,
        species_column=args.species_column,
    )

    if args.unique:
        species_ids = list(dict.fromkeys(species_ids))

    gff_ids = extract_gff_ids(args.gff)
    missing_ids = get_missing_ids(species_ids, gff_ids)

    total, present, absent = count_presence(species_ids, gff_ids)

    if args.missing_out is not None:
        args.missing_out.parent.mkdir(parents=True, exist_ok=True)
        args.missing_out.write_text(
            "\n".join(missing_ids) + ("\n" if missing_ids else "")
        )

    print(f"species={args.species}")
    print(f"total_ids={total}")
    print(f"present_ids={present}")
    print(f"absent_ids={absent}")
    if args.missing_out is not None:
        print(f"missing_ids_file={args.missing_out}")


if __name__ == "__main__":
    main()
