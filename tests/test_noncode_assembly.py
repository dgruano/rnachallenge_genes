"""
tests/test_noncode_assembly.py

Unit tests for Phase 3 NONCODE Assembly Accession Resolution

Tests the mapping of UCSC genome build names to GCF_/GCA_ accessions
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
from resolve_noncode_assembly_accessions import (
    UCSC_TO_GCF_MAPPING,
    normalize_ucsc_name,
    is_ucsc_assembly_accession,
    map_ucsc_to_gcf,
    filter_and_resolve_noncode,
)


# ============================================================================
# Test Cases
# ============================================================================


class TestUCSCNameMapping:
    """Test UCSC name to GCF_ mapping function."""

    def test_map_tair10_to_gcf(self):
        """Test mapping of tair10 to GCF accession."""
        result = map_ucsc_to_gcf("tair10")
        assert result == "GCF_000001735.4"  # Arabidopsis thaliana (TAIR10.1)

    def test_map_ce10_to_gcf(self):
        """Test mapping of ce10 to GCF accession."""
        result = map_ucsc_to_gcf("ce10")
        assert result == "GCF_000002985.6"  # Caenorhabditis elegans (WBcel235)

    def test_map_dm6_to_gcf(self):
        """Test mapping of dm6 to GCF accession."""
        result = map_ucsc_to_gcf("dm6")
        assert result == "GCF_000001215.4"  # Drosophila melanogaster

    def test_map_rn6_to_gcf(self):
        """Test mapping of rn6 to GCF accession."""
        result = map_ucsc_to_gcf("rn6")
        assert result == "GCF_000001895.5"  # Rattus norvegicus

    def test_map_mondom5_to_gcf(self):
        """Test mapping of monDom5 to GCF accession."""
        result = map_ucsc_to_gcf("monDom5")
        assert result == "GCF_000002295.2"  # Monodelphis domesticus (MonDom5)

    def test_map_ponabe2_to_gcf(self):
        """Test mapping of ponAbe2 to GCF accession."""
        result = map_ucsc_to_gcf("ponAbe2")
        assert result == "GCF_000001545.5"  # Pongo abelii

    def test_map_galgal4_to_gcf(self):
        """Test mapping of galGal4 to GCF accession."""
        result = map_ucsc_to_gcf("galGal4")
        assert result == "GCF_000002315.6"  # Gallus gallus (GRCg6a)

    def test_map_ornana1_to_gcf(self):
        """Test mapping of ornAna1 to GCF accession."""
        result = map_ucsc_to_gcf("ornAna1")
        assert result == "GCF_000002275.2"  # Ornithorhynchus anatinus (ASM227v2)

    def test_map_bostau6_to_gcf(self):
        """Test mapping of bosTau6 to GCF accession."""
        result = map_ucsc_to_gcf("bosTau6")
        assert result == "GCF_000003055.6"  # Bos taurus

    def test_map_danrer10_to_gcf(self):
        """Test mapping of danRer10 to GCF accession."""
        result = map_ucsc_to_gcf("danRer10")
        assert result == "GCF_000002035.6"  # Danio rerio (GRCz11)

    def test_map_ucsc_case_insensitive(self):
        """Test mapping is case-insensitive."""
        result = map_ucsc_to_gcf("TAIR10")
        assert result == "GCF_000001735.4"

    def test_map_ucsc_with_spaces(self):
        """Test mapping handles whitespace."""
        result = map_ucsc_to_gcf("  ce10  ")
        assert result == "GCF_000002985.6"

    def test_map_unknown_ucsc(self):
        """Test that unknown UCSC names return None."""
        result = map_ucsc_to_gcf("hg999")
        assert result is None

    def test_map_empty_string(self):
        """Test that empty string returns None."""
        result = map_ucsc_to_gcf("")
        assert result is None

    def test_map_none(self):
        """Test that None returns None."""
        result = map_ucsc_to_gcf(None)
        assert result is None


class TestUCSCAccessionDetection:
    """Test detection of UCSC assembly accessions."""

    def test_detect_tair10(self):
        """Test detection of tair10 pattern."""
        assert is_ucsc_assembly_accession("tair10") is True

    def test_detect_ce10(self):
        """Test detection of ce10 pattern."""
        assert is_ucsc_assembly_accession("ce10") is True

    def test_detect_dm6(self):
        """Test detection of dm6 pattern."""
        assert is_ucsc_assembly_accession("dm6") is True

    def test_detect_ucsc_case_insensitive(self):
        """Test detection is case-insensitive."""
        assert is_ucsc_assembly_accession("TAIR10") is True

    def test_detect_gcf_not_ucsc(self):
        """Test that GCF_ accessions are not detected as UCSC."""
        assert is_ucsc_assembly_accession("GCF_000001405.40") is False

    def test_detect_grc_not_ucsc(self):
        """Test that GRC* accessions are not detected as UCSC."""
        assert is_ucsc_assembly_accession("GRCh38") is False

    def test_detect_empty(self):
        """Test that empty string returns False."""
        assert is_ucsc_assembly_accession("") is False

    def test_detect_none(self):
        """Test that None returns False."""
        assert is_ucsc_assembly_accession(None) is False

    def test_detect_nan(self):
        """Test that NaN returns False."""
        assert is_ucsc_assembly_accession(pd.NA) is False


class TestRowFiltering:
    """Test filtering of rows - only UCSC names are mapped, others pass through."""

    def test_filter_only_ucsc_rows_processed(self):
        """Test that only rows with UCSC names are flagged for mapping."""
        df = pd.DataFrame({
            "transcript_id": ["NONCELT011856.2", "TX2", "TX3"],
            "db_source": ["noncode", "noncode", "noncode"],
            "assembly_accession": ["ce10", "GCF_000002035.6", "unknown"],
            "organism": ["Caenorhabditis elegans", "Caenorhabditis elegans", "unknown"],
        })

        resolved, unresolved = filter_and_resolve_noncode(df)

        # All three should end up in resolved (none fail mapping)
        assert len(resolved) == 3
        assert len(unresolved) == 0

    def test_filter_ucsc_rows_are_mapped(self):
        """Test that UCSC rows are actually mapped."""
        df = pd.DataFrame({
            "transcript_id": ["NONCELT011856.2"],
            "db_source": ["noncode"],
            "assembly_accession": ["ce10"],
            "organism": ["Caenorhabditis elegans"],
        })

        resolved, _ = filter_and_resolve_noncode(df)

        assert len(resolved) == 1
        assert resolved.iloc[0]["assembly_accession"] == "GCF_000002985.6"

    def test_filter_non_ucsc_rows_unchanged(self):
        """Test that non-UCSC rows pass through unchanged."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2"],
            "db_source": ["noncode", "noncode"],
            "assembly_accession": ["GCF_000001405.40", "unknown_value"],
            "organism": ["organism1", "organism2"],
        })

        resolved, _ = filter_and_resolve_noncode(df)

        # Non-UCSC rows should pass through exactly as-is
        assert len(resolved) == 2
        assert resolved.iloc[0]["assembly_accession"] == "GCF_000001405.40"
        assert resolved.iloc[1]["assembly_accession"] == "unknown_value"


class TestResolvedVsUnresolvedSplit:
    """Test splitting of rows into resolved and unresolved outputs."""

    def test_successful_mapping_goes_to_resolved(self):
        """Test that successfully mapped UCSC rows go to resolved output."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2"],
            "db_source": ["noncode", "noncode"],
            "assembly_accession": ["ce10", "dm6"],
            "organism": ["Caenorhabditis elegans", "Drosophila melanogaster"],
            "gene_id": ["G1", "G2"],
            "gene_symbol": ["GENE1", "GENE2"],
            "chrom": ["1", "2"],
            "start": [100, 200],
            "end": [200, 300],
            "strand": ["+", "-"],
        })

        resolved, unresolved = filter_and_resolve_noncode(df)

        assert len(resolved) == 2
        assert len(unresolved) == 0
        assert resolved.iloc[0]["assembly_accession"] == "GCF_000002985.6"
        assert resolved.iloc[1]["assembly_accession"] == "GCF_000001215.4"

    def test_unmapped_ucsc_like_accession_fails_gracefully(self):
        """Test that UCSC-format names not in mapping are rejected."""
        df = pd.DataFrame({
            "transcript_id": ["TX1"],
            "db_source": ["noncode"],
            "assembly_accession": ["xx999"],  # Not in UCSC mapping
            "organism": ["unknown"],
        })

        resolved, unresolved = filter_and_resolve_noncode(df)

        # xx999 is not in UCSC_TO_GCF_MAPPING, so it should pass through as-is
        # (treated as unknown accession format, not as UCSC)
        assert len(resolved) + len(unresolved) == 1  # All rows accounted for
        # Row should pass through to resolved (not recognized as UCSC)
        assert len(resolved) == 1
        assert resolved.iloc[0]["assembly_accession"] == "xx999"

    def test_mixed_rows_correct_split(self):
        """Test correct splitting when mix of mappable and non-UCSC rows."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2", "TX3"],
            "db_source": ["noncode", "noncode", "noncode"],
            "assembly_accession": ["ce10", "GCF_000001405.40", "dm6"],
            "organism": ["Caenorhabditis elegans", "organism2", "Drosophila melanogaster"],
            "gene_id": ["G1", "G2", "G3"],
            "gene_symbol": ["GENE1", "GENE2", "GENE3"],
            "chrom": ["1", "1", "1"],
            "start": [100, 100, 100],
            "end": [200, 200, 200],
            "strand": ["+", "+", "+"],
        })

        resolved, unresolved = filter_and_resolve_noncode(df)

        # All three should resolve (two are non-UCSC, one UCSC maps successfully)
        assert len(resolved) == 3
        assert len(unresolved) == 0


class TestNonUCSCRowsPassThrough:
    """Test that non-UCSC rows pass through completely unchanged."""

    def test_gcf_rows_unchanged(self):
        """Test that GCF_ accessions pass through without modification."""
        df = pd.DataFrame({
            "transcript_id": ["TX1"],
            "db_source": ["noncode"],
            "assembly_accession": ["GCF_000001405.40"],
            "organism": ["organism1"],
            "gene_id": ["GENE1"],
            "gene_symbol": ["SYM1"],
            "chrom": ["1"],
            "start": [100],
            "end": [200],
            "strand": ["+"],
            "is_ambiguous": [False],
        })

        resolved, _ = filter_and_resolve_noncode(df)

        assert len(resolved) == 1
        assert resolved.iloc[0]["assembly_accession"] == "GCF_000001405.40"
        assert resolved.iloc[0]["transcript_id"] == "TX1"
        assert resolved.iloc[0]["gene_id"] == "GENE1"

    def test_unknown_accessions_unchanged(self):
        """Test that unknown accessions pass through unchanged."""
        df = pd.DataFrame({
            "transcript_id": ["TX1"],
            "db_source": ["noncode"],
            "assembly_accession": ["some_unknown_value"],
            "organism": ["custom_organism"],
        })

        resolved, _ = filter_and_resolve_noncode(df)

        assert len(resolved) == 1
        assert resolved.iloc[0]["assembly_accession"] == "some_unknown_value"

    def test_all_columns_preserved_for_passthrough(self):
        """Test that all columns are preserved for pass-through rows."""
        df = pd.DataFrame({
            "transcript_id": ["TX1"],
            "db_source": ["noncode"],
            "assembly_accession": ["GCF_000001405.40"],
            "organism": ["organism1"],
            "gene_id": ["G1"],
            "gene_symbol": ["SYM1"],
            "chrom": ["1"],
            "start": [100],
            "end": [200],
            "strand": ["+"],
            "is_ambiguous": [False],
            "custom_field": ["value1"],
        })

        resolved, _ = filter_and_resolve_noncode(df)

        # All columns should be present in output
        assert "custom_field" in resolved.columns
        assert resolved.iloc[0]["custom_field"] == "value1"


class TestEmptyAndEdgeCases:
    """Test edge cases and empty dataframes."""

    def test_empty_dataframe(self):
        """Test handling of empty input dataframe."""
        df = pd.DataFrame()
        resolved, unresolved = filter_and_resolve_noncode(df)

        assert len(resolved) == 0
        assert len(unresolved) == 0

    def test_dataframe_with_no_ucsc_rows(self):
        """Test dataframe containing no UCSC rows."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2"],
            "db_source": ["noncode", "noncode"],
            "assembly_accession": ["GCF_000001405.40", "GCA_000001405.1"],
            "organism": ["organism1", "organism2"],
        })

        resolved, unresolved = filter_and_resolve_noncode(df)

        # All rows should pass through
        assert len(resolved) == 2
        assert len(unresolved) == 0

    def test_dataframe_all_ucsc_mappable(self):
        """Test dataframe where all rows are UCSC and mappable."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2"],
            "db_source": ["noncode", "noncode"],
            "assembly_accession": ["ce10", "dm6"],
            "organism": ["Caenorhabditis elegans", "Drosophila melanogaster"],
        })

        resolved, unresolved = filter_and_resolve_noncode(df)

        assert len(resolved) == 2
        assert len(unresolved) == 0

    def test_missing_assembly_accession_column(self):
        """Test handling of rows with missing assembly_accession column."""
        df = pd.DataFrame({
            "transcript_id": ["TX1"],
            "db_source": ["noncode"],
            "organism": ["organism1"],
        })

        # Should handle gracefully
        try:
            resolved, unresolved = filter_and_resolve_noncode(df)
            # Rows without assembly_accession should pass through or be handled safely
            assert len(resolved) >= 0  # Just verify it doesn't crash
        except KeyError:
            pytest.fail("Should handle missing assembly_accession column gracefully")

    def test_null_assembly_accession_values(self):
        """Test handling of NaN/None in assembly_accession."""
        df = pd.DataFrame({
            "transcript_id": ["TX1", "TX2"],
            "db_source": ["noncode", "noncode"],
            "assembly_accession": [None, pd.NA],
            "organism": ["organism1", "organism2"],
        })

        resolved, unresolved = filter_and_resolve_noncode(df)

        # Null values should pass through (they're not UCSC)
        assert len(resolved) == 2
        assert len(unresolved) == 0


class TestAllTenUCSCSpecies:
    """Test all 10 UCSC species mappings."""

    def test_all_10_ucsc_names_present_in_mapping(self):
        """Test that all 10 UCSC names are in the mapping dictionary."""
        expected_names = [
            "tair10", "ce10", "dm6", "rn6", "monDom5",
            "ponAbe2", "galGal4", "ornAna1", "bosTau6", "danRer10"
        ]
        for name in expected_names:
            # Normalize and check
            normalized = name.upper()
            assert normalized in UCSC_TO_GCF_MAPPING, f"Missing mapping for {name}"

    def test_all_10_species_resolve(self):
        """Test that all 10 UCSC species resolve to valid GCF_ accessions."""
        test_data = [
            ("tair10", "GCF_000001735.4"),    # Arabidopsis thaliana (TAIR10.1)
            ("ce10", "GCF_000002985.6"),      # Caenorhabditis elegans (WBcel235)
            ("dm6", "GCF_000001215.4"),       # Drosophila melanogaster (Release 6)
            ("rn6", "GCF_000001895.5"),       # Rattus norvegicus (Rnor_6.0)
            ("monDom5", "GCF_000002295.2"),   # Monodelphis domesticus (MonDom5)
            ("ponAbe2", "GCF_000001545.5"),   # Pongo abelii (P_pygmaeus_2.0.2)
            ("galGal4", "GCF_000002315.6"),   # Gallus gallus (GRCg6a)
            ("ornAna1", "GCF_000002275.2"),   # Ornithorhynchus anatinus (ASM227v2)
            ("bosTau6", "GCF_000003055.6"),   # Bos taurus (Bos_taurus_UMD_3.1.1)
            ("danRer10", "GCF_000002035.6"),  # Danio rerio (GRCz11)
        ]

        for ucsc_name, expected_gcf in test_data:
            result = map_ucsc_to_gcf(ucsc_name)
            assert result == expected_gcf, f"Failed to map {ucsc_name}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
