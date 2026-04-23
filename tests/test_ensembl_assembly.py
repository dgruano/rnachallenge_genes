"""
tests/test_ensembl_assembly.py

Unit tests for Phase 2 Ensembl Assembly Accession Resolution

Tests the mapping of GRC* assembly names to GCF_/GCA_ accessions
and proper filtering of rows.
"""

import sys
import tempfile
from pathlib import Path
from typing import Optional

import pandas as pd
import pytest

# Add workflow scripts to path
SCRIPT_DIR = Path(__file__).parent.parent / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

# Import functions from the resolver script
from resolve_ensembl_assembly_accessions import (
    ASSEMBLY_NAME_MAPPING,
    normalize_build_name,
    is_grc_assembly_accession,
    map_grc_to_gcf,
    filter_and_resolve_ensembl,
)


# ============================================================================
# Test Cases
# ============================================================================


class TestGRCNameMapping:
    """Test GRC name to GCF_ mapping function."""

    def test_map_grch38_to_gcf(self):
        """Test mapping of GRCh38 to standard GCF accession."""
        result = map_grc_to_gcf("GRCh38")
        assert result == "GCF_000001405.40"

    def test_map_grch38_lowercase(self):
        """Test mapping is case-insensitive."""
        result = map_grc_to_gcf("grch38")
        assert result == "GCF_000001405.40"

    def test_map_grch38_with_spaces(self):
        """Test mapping handles whitespace."""
        result = map_grc_to_gcf("  GRCh38  ")
        assert result == "GCF_000001405.40"

    def test_map_grcz11(self):
        """Test mapping of zebrafish GRCz11."""
        result = map_grc_to_gcf("GRCz11")
        assert result == "GCF_000002035.6"

    def test_map_grcrh1(self):
        """Test mapping of rhesus macaque GRCRH1."""
        result = map_grc_to_gcf("GRCRH1")
        assert result == "GCF_000008735.2"

    def test_map_unknown_grc(self):
        """Test that unknown GRC versions return None."""
        result = map_grc_to_gcf("GRCX99")
        assert result is None

    def test_map_empty_string(self):
        """Test that empty string returns None."""
        result = map_grc_to_gcf("")
        assert result is None

    def test_map_none(self):
        """Test that None returns None."""
        result = map_grc_to_gcf(None)
        assert result is None


class TestGRCAccessionDetection:
    """Test detection of GRC* assembly accessions."""

    def test_detect_grch38(self):
        """Test detection of GRCh38 pattern."""
        assert is_grc_assembly_accession("GRCh38") is True

    def test_detect_grcz11(self):
        """Test detection of GRCz11 pattern."""
        assert is_grc_assembly_accession("GRCz11") is True

    def test_detect_grcrh1(self):
        """Test detection of GRCRH1 pattern."""
        assert is_grc_assembly_accession("GRCRH1") is True

    def test_detect_grc_lowercase(self):
        """Test detection is case-insensitive."""
        assert is_grc_assembly_accession("grch38") is True

    def test_detect_gcf_not_grc(self):
        """Test that GCF_ accessions are not detected as GRC."""
        assert is_grc_assembly_accession("GCF_000001405.40") is False

    def test_detect_ncbi_scaffold(self):
        """Test that NCBI scaffold accessions are not detected as GRC."""
        assert is_grc_assembly_accession("NW_017955610.1") is False

    def test_detect_empty(self):
        """Test that empty string returns False."""
        assert is_grc_assembly_accession("") is False

    def test_detect_none(self):
        """Test that None returns False."""
        assert is_grc_assembly_accession(None) is False

    def test_detect_nan(self):
        """Test that NaN returns False."""
        assert is_grc_assembly_accession(pd.NA) is False


class TestRowFiltering:
    """Test filtering of rows - only GRC* are mapped, others pass through."""

    def test_filter_only_grc_rows_processed(self):
        """Test that only rows with GRC* are flagged for mapping."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2", "TX3"],
            "db_source": ["ensembl", "ensembl", "ensembl"],
            "assembly_accession": ["GRCh38", "GCF_000001405.40", "unknown"],
            "organism": ["homo sapiens", "homo sapiens", "homo sapiens"],
        })

        resolved, unresolved = filter_and_resolve_ensembl(df)

        # All three should end up in resolved (none fail mapping)
        assert len(resolved) == 3
        assert len(unresolved) == 0

    def test_filter_grc_rows_are_mapped(self):
        """Test that GRC* rows are actually mapped."""
        df = pd.DataFrame({
            "transcript_id": ["TX1"],
            "db_source": ["ensembl"],
            "assembly_accession": ["GRCh38"],
            "organism": ["homo sapiens"],
        })

        resolved, _ = filter_and_resolve_ensembl(df)

        assert len(resolved) == 1
        assert resolved.iloc[0]["assembly_accession"] == "GCF_000001405.40"

    def test_filter_non_grc_rows_unchanged(self):
        """Test that non-GRC rows pass through unchanged."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2"],
            "db_source": ["ensembl", "ensembl"],
            "assembly_accession": ["GCF_000001405.40", "unknown_value"],
            "organism": ["homo sapiens", "canis lupus"],
        })

        resolved, _ = filter_and_resolve_ensembl(df)

        # Non-GRC rows should pass through exactly as-is
        assert len(resolved) == 2
        assert resolved.iloc[0]["assembly_accession"] == "GCF_000001405.40"
        assert resolved.iloc[1]["assembly_accession"] == "unknown_value"


class TestResolvedVsUnresolvedSplit:
    """Test splitting of rows into resolved and unresolved outputs."""

    def test_successful_mapping_goes_to_resolved(self):
        """Test that successfully mapped GRC* rows go to resolved output."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2"],
            "db_source": ["ensembl", "ensembl"],
            "assembly_accession": ["GRCh38", "GRCz11"],
            "organism": ["homo sapiens", "danio rerio"],
            "gene_id": ["G1", "G2"],
            "gene_symbol": ["GENE1", "GENE2"],
            "chrom": ["1", "1"],
            "start": [100, 200],
            "end": [200, 300],
            "strand": ["+", "-"],
        })

        resolved, unresolved = filter_and_resolve_ensembl(df)

        assert len(resolved) == 2
        assert len(unresolved) == 0
        assert resolved.iloc[0]["assembly_accession"] == "GCF_000001405.40"
        assert resolved.iloc[1]["assembly_accession"] == "GCF_000002035.6"

    def test_unmapped_grc_goes_to_unresolved(self):
        """Test that unmappable GRC* rows go to unresolved output."""
        df = pd.DataFrame({
            "transcript_id": ["TX1"],
            "db_source": ["ensembl"],
            "assembly_accession": ["GRCX99_unknown"],
            "organism": ["unknown_organism"],
        })

        resolved, unresolved = filter_and_resolve_ensembl(df)

        assert len(resolved) == 0
        assert len(unresolved) == 1
        assert unresolved.iloc[0]["transcript_id"] == "TX1"
        assert unresolved.iloc[0]["db_source"] == "ensembl"
        assert "grc_mapping_failed" in unresolved.iloc[0]["reason"]

    def test_mixed_rows_correct_split(self):
        """Test correct splitting when mix of mappable and non-GRC rows."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2", "TX3"],
            "db_source": ["ensembl", "ensembl", "ensembl"],
            "assembly_accession": ["GRCh38", "GCF_000001405.40", "GRCz11"],
            "organism": ["homo sapiens", "homo sapiens", "danio rerio"],
            "gene_id": ["G1", "G2", "G3"],
            "gene_symbol": ["GENE1", "GENE2", "GENE3"],
            "chrom": ["1", "1", "1"],
            "start": [100, 100, 100],
            "end": [200, 200, 200],
            "strand": ["+", "+", "+"],
        })

        resolved, unresolved = filter_and_resolve_ensembl(df)

        # All three should resolve (two are non-GRC, one GRC maps successfully)
        assert len(resolved) == 3
        assert len(unresolved) == 0


class TestNonGRCRowsPassThrough:
    """Test that non-GRC rows pass through completely unchanged."""

    def test_gcf_rows_unchanged(self):
        """Test that GCF_ accessions pass through without modification."""
        df = pd.DataFrame({
            "transcript_id": ["TX1"],
            "db_source": ["ensembl"],
            "assembly_accession": ["GCF_000001405.40"],
            "organism": ["homo sapiens"],
            "gene_id": ["GENE1"],
            "gene_symbol": ["SYM1"],
            "chrom": ["1"],
            "start": [100],
            "end": [200],
            "strand": ["+"],
            "is_ambiguous": [False],
        })

        resolved, _ = filter_and_resolve_ensembl(df)

        assert len(resolved) == 1
        assert resolved.iloc[0]["assembly_accession"] == "GCF_000001405.40"
        assert resolved.iloc[0]["transcript_id"] == "TX1"
        assert resolved.iloc[0]["gene_id"] == "GENE1"

    def test_unknown_accessions_unchanged(self):
        """Test that unknown accessions pass through unchanged."""
        df = pd.DataFrame({
            "transcript_id": ["TX1"],
            "db_source": ["ensembl"],
            "assembly_accession": ["some_unknown_value"],
            "organism": ["custom_organism"],
        })

        resolved, _ = filter_and_resolve_ensembl(df)

        assert len(resolved) == 1
        assert resolved.iloc[0]["assembly_accession"] == "some_unknown_value"

    def test_all_columns_preserved_for_passthrough(self):
        """Test that all columns are preserved for pass-through rows."""
        df = pd.DataFrame({
            "transcript_id": ["TX1"],
            "db_source": ["ensembl"],
            "assembly_accession": ["GCF_000001405.40"],
            "organism": ["homo sapiens"],
            "gene_id": ["G1"],
            "gene_symbol": ["SYM1"],
            "chrom": ["1"],
            "start": [100],
            "end": [200],
            "strand": ["+"],
            "is_ambiguous": [False],
            "custom_field": ["value1"],
        })

        resolved, _ = filter_and_resolve_ensembl(df)

        # All columns should be present in output
        assert "custom_field" in resolved.columns
        assert resolved.iloc[0]["custom_field"] == "value1"


class TestEmptyAndEdgeCases:
    """Test edge cases and empty dataframes."""

    def test_empty_dataframe(self):
        """Test handling of empty input dataframe."""
        df = pd.DataFrame()
        resolved, unresolved = filter_and_resolve_ensembl(df)

        assert len(resolved) == 0
        assert len(unresolved) == 0

    def test_dataframe_with_no_grc_rows(self):
        """Test dataframe containing no GRC* rows."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2"],
            "db_source": ["ensembl", "ensembl"],
            "assembly_accession": ["GCF_000001405.40", "GCA_000001405.1"],
            "organism": ["homo sapiens", "homo sapiens"],
        })

        resolved, unresolved = filter_and_resolve_ensembl(df)

        # All rows should pass through
        assert len(resolved) == 2
        assert len(unresolved) == 0

    def test_dataframe_all_grc_mappable(self):
        """Test dataframe where all rows are GRC and mappable."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2"],
            "db_source": ["ensembl", "ensembl"],
            "assembly_accession": ["GRCh38", "GRCz11"],
            "organism": ["homo sapiens", "danio rerio"],
        })

        resolved, unresolved = filter_and_resolve_ensembl(df)

        assert len(resolved) == 2
        assert len(unresolved) == 0

    def test_missing_assembly_accession_column(self):
        """Test handling of rows with missing assembly_accession column."""
        df = pd.DataFrame({
            "transcript_id": ["TX1"],
            "db_source": ["ensembl"],
            "organism": ["homo sapiens"],
        })

        # Should handle gracefully
        try:
            resolved, unresolved = filter_and_resolve_ensembl(df)
            # Rows without assembly_accession should pass through or be handled safely
            assert len(resolved) >= 0  # Just verify it doesn't crash
        except KeyError:
            pytest.fail("Should handle missing assembly_accession column gracefully")

    def test_null_assembly_accession_values(self):
        """Test handling of NaN/None in assembly_accession."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2"],
            "db_source": ["ensembl", "ensembl"],
            "assembly_accession": [None, pd.NA],
            "organism": ["homo sapiens", "homo sapiens"],
        })

        resolved, unresolved = filter_and_resolve_ensembl(df)

        # Null values should pass through (they're not GRC*)
        assert len(resolved) == 2
        assert len(unresolved) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
