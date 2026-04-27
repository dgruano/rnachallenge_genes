"""
tests/test_resolve_ncbi_chromosome_accessions.py

Unit tests for resolve_ncbi_chromosome_accessions script logic.
Tests the filtering, mapping, and output writing logic in isolation
(snakemake interface is excluded).
"""

import sys
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).parent.parent / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from resolve_ncbi_chromosome_accessions import (
    is_chromosomal_accession,
    resolve_chromosomal_rows,
)


SAMPLE_COLS = [
    "transcript_id", "db_source", "gene_id", "gene_symbol", "organism",
    "assembly_accession", "chrom", "start", "end", "strand", "is_ambiguous",
]


def _make_df(rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    for col in SAMPLE_COLS:
        if col not in df.columns:
            df[col] = pd.NA
    return df[SAMPLE_COLS]


class TestIsChromosomalAccession:

    def test_nc_prefix(self):
        assert is_chromosomal_accession("NC_000001.11") is True

    def test_nt_prefix(self):
        assert is_chromosomal_accession("NT_033779.5") is True

    def test_nw_prefix(self):
        assert is_chromosomal_accession("NW_006239xxx") is True

    def test_gcf_is_not_chromosomal(self):
        assert is_chromosomal_accession("GCF_000001405.40") is False

    def test_gca_is_not_chromosomal(self):
        assert is_chromosomal_accession("GCA_000001405.1") is False

    def test_none_is_not_chromosomal(self):
        assert is_chromosomal_accession(None) is False

    def test_nan_is_not_chromosomal(self):
        assert is_chromosomal_accession(float("nan")) is False

    def test_empty_string_is_not_chromosomal(self):
        assert is_chromosomal_accession("") is False


class TestResolveChromosomalRows:

    def test_nc_row_replaced_with_gcf(self):
        df = _make_df([
            {"transcript_id": "NM_001.1", "db_source": "ncbi",
             "assembly_accession": "NC_000001.11"},
        ])
        mapping = {"NC_000001.11": "GCF_000001405.40"}
        resolved, unresolved = resolve_chromosomal_rows(df, mapping)
        assert len(resolved) == 1
        assert resolved.iloc[0]["assembly_accession"] == "GCF_000001405.40"
        assert len(unresolved) == 0

    def test_unmapped_nc_goes_to_unresolved(self):
        df = _make_df([
            {"transcript_id": "NM_001.1", "db_source": "ncbi",
             "assembly_accession": "NT_033779.5"},
        ])
        mapping = {"NT_033779.5": None}
        resolved, unresolved = resolve_chromosomal_rows(df, mapping)
        assert len(resolved) == 0
        assert len(unresolved) == 1
        assert unresolved.iloc[0]["reason"] == "chromosomal_mapping_failed:NT_033779.5"

    def test_gcf_rows_pass_through_unchanged(self):
        df = _make_df([
            {"transcript_id": "TX1", "db_source": "ncbi",
             "assembly_accession": "GCF_000001405.40"},
        ])
        resolved, unresolved = resolve_chromosomal_rows(df, mapping={})
        assert len(resolved) == 1
        assert resolved.iloc[0]["assembly_accession"] == "GCF_000001405.40"
        assert len(unresolved) == 0

    def test_mixed_df_routes_correctly(self):
        df = _make_df([
            {"transcript_id": "TX1", "assembly_accession": "GCF_000001405.40"},
            {"transcript_id": "TX2", "assembly_accession": "NC_000001.11"},
            {"transcript_id": "TX3", "assembly_accession": "NT_033779.5"},
        ])
        mapping = {"NC_000001.11": "GCF_000001405.40", "NT_033779.5": None}
        resolved, unresolved = resolve_chromosomal_rows(df, mapping)
        assert len(resolved) == 2   # GCF_ passthrough + NC_ mapped
        assert len(unresolved) == 1  # NT_ unmapped
        gcf_vals = set(resolved["assembly_accession"])
        assert gcf_vals == {"GCF_000001405.40"}

    def test_empty_df(self):
        df = pd.DataFrame(columns=SAMPLE_COLS)
        resolved, unresolved = resolve_chromosomal_rows(df, mapping={})
        assert resolved.empty
        assert unresolved.empty
