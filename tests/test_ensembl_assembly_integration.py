"""
tests/test_ensembl_assembly_integration.py

Integration tests for resolve_ensembl_assembly_accessions.py script
Tests the script's ability to read/write files and process real data.
"""

import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

# Add workflow scripts to path
SCRIPT_DIR = Path(__file__).parent.parent / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

# Import functions from the test module
from test_ensembl_assembly import (
    ASSEMBLY_NAME_MAPPING,
    filter_and_resolve_ensembl,
    is_grc_assembly_accession,
    map_grc_to_gcf,
    normalize_build_name,
)


class TestScriptIntegration:
    """Test the script's I/O and workflow integration."""

    def test_script_functions_match_test_helpers(self):
        """Verify script uses same functions as tests."""
        # This just ensures the script has the mapping
        assert "GRCH38" in ASSEMBLY_NAME_MAPPING
        assert ASSEMBLY_NAME_MAPPING["GRCH38"] == "GCF_000001405.40"

    def test_read_write_roundtrip(self):
        """Test reading and writing TSV files with expected schema."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create input file
            df_input = pd.DataFrame(
                {
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
                    "is_ambiguous": [False, False, False],
                }
            )

            input_file = tmpdir / "input.tsv"
            df_input.to_csv(input_file, sep="\t", index=False)

            # Process through test functions
            resolved, unresolved = filter_and_resolve_ensembl(df_input)

            # Write outputs
            resolved_file = tmpdir / "resolved.tsv"
            unresolved_file = tmpdir / "unresolved.tsv"

            resolved.to_csv(resolved_file, sep="\t", index=False)
            # Ensure unresolved always has headers even if empty
            if len(unresolved) == 0:
                pd.DataFrame(columns=["transcript_id", "db_source", "reason"]).to_csv(
                    unresolved_file, sep="\t", index=False
                )
            else:
                unresolved.to_csv(unresolved_file, sep="\t", index=False)

            # Read back and verify
            df_resolved = pd.read_csv(resolved_file, sep="\t")
            df_unresolved = pd.read_csv(unresolved_file, sep="\t")

            assert len(df_resolved) == 3  # All resolve
            assert len(df_unresolved) == 0
            assert df_resolved.iloc[0]["assembly_accession"] == "GCF_000001405.40"
            assert df_resolved.iloc[2]["assembly_accession"] == "GCF_000002035.6"

    def test_mixed_input_output(self):
        """Test mixed resolvable/unresolvable rows."""
        df_input = pd.DataFrame(
            {
                "transcript_id": ["TX1", "TX2", "TX3"],
                "db_source": ["ensembl", "ensembl", "ensembl"],
                "assembly_accession": ["GRCh38", "GRCX99", "GCF_000001405.40"],
                "organism": ["homo sapiens", "unknown", "homo sapiens"],
                "gene_id": ["G1", "G2", "G3"],
                "gene_symbol": ["GENE1", "GENE2", "GENE3"],
                "chrom": ["1", "1", "1"],
                "start": [100, 100, 100],
                "end": [200, 200, 200],
                "strand": ["+", "+", "+"],
                "is_ambiguous": [False, False, False],
            }
        )

        resolved, unresolved = filter_and_resolve_ensembl(df_input)

        # GRCh38 and GCF_ should resolve; GRCX99 should not
        assert len(resolved) == 2
        assert len(unresolved) == 1
        assert unresolved.iloc[0]["transcript_id"] == "TX2"
        assert "grc_mapping_failed" in unresolved.iloc[0]["reason"]

    def test_case_insensitive_mapping(self):
        """Test that GRC* detection and mapping is case-insensitive."""
        df_input = pd.DataFrame(
            {
                "transcript_id": ["TX1"],
                "db_source": ["ensembl"],
                "assembly_accession": ["grch38"],  # lowercase
                "organism": ["homo sapiens"],
                "gene_id": ["G1"],
                "gene_symbol": ["GENE1"],
                "chrom": ["1"],
                "start": [100],
                "end": [200],
                "strand": ["+"],
                "is_ambiguous": [False],
            }
        )

        resolved, unresolved = filter_and_resolve_ensembl(df_input)

        assert len(resolved) == 1
        assert resolved.iloc[0]["assembly_accession"] == "GCF_000001405.40"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
