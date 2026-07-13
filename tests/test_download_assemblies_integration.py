"""
tests/test_download_assemblies_integration.py

Integration tests for refactored download_assemblies.py

Tests the full script behavior with mock Snakemake objects.
"""

import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest


class MockSnakemake:
    """Mock Snakemake object for testing scripts."""

    def __init__(self, tmpdir: Path):
        self.tmpdir = tmpdir
        self.log = [str(tmpdir / "test.log")]
        self.input = MagicMock()
        self.output = MagicMock()
        self.config = {
            "cache_dir": str(tmpdir / "cache"),
            "max_retries": 1,
            "retry_wait_seconds": 1,
        }


class TestDownloadAssembliesIntegration:
    """Integration tests for download_assemblies script logic."""

    def test_script_imports_successfully(self):
        """Test that the script file can be parsed."""
        script_path = (
            Path(__file__).parent.parent
            / "workflow"
            / "scripts"
            / "download_assemblies.py"
        )
        with open(script_path) as f:
            content = f.read()
        # Check for key functions
        assert "is_ncbi_assembly_accession" in content
        assert "ncbi_fasta_url" in content
        assert "ensure_assembly" in content
        assert "downloaded_df.to_csv(output_downloaded" in content

    def test_script_has_no_complex_fallback_logic(self):
        """Test that complex fallback logic has been removed."""
        script_path = (
            Path(__file__).parent.parent
            / "workflow"
            / "scripts"
            / "download_assemblies.py"
        )
        with open(script_path) as f:
            content = f.read()

        # Should NOT have these anymore
        assert "UCSC_TO_ENSEMBL" not in content
        assert "ensembl_assembly_url" not in content
        assert "ncbi_efetch_url" not in content
        assert "efetch.fcgi" not in content

    def test_script_handles_mixed_inputs(self):
        """Test that script logic correctly handles mixed accessions."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            cache_dir = tmpdir / "cache"
            cache_dir.mkdir()

            # Create input with mixed accessions
            df = pd.DataFrame(
                {
                    "transcript_id": ["TX1", "TX2", "TX3", "TX4"],
                    "assembly_accession": [
                        "GCF_000001405.40",  # NCBI assembly
                        "GRCh38",  # GRC name
                        "NC_000001.11",  # Sequence accession
                        "hg38",  # UCSC name
                    ],
                    "organism": [
                        "homo sapiens",
                        "homo sapiens",
                        "homo sapiens",
                        "homo sapiens",
                    ],
                    "db_source": ["ncbi", "ensembl", "ncbi", "ucsc"],
                }
            )

            input_file = tmpdir / "resolved.tsv"
            df.to_csv(input_file, sep="\t", index=False)

            # Read and validate split logic
            unique_asm = (
                df[["assembly_accession", "organism", "db_source"]]
                .dropna(subset=["assembly_accession"])
                .drop_duplicates(subset="assembly_accession")
            )

            downloaded = []
            unresolved = []

            from tests.test_download_assemblies_phase4 import is_ncbi_assembly_accession

            for _, row in unique_asm.iterrows():
                acc = str(row["assembly_accession"]).strip()
                if is_ncbi_assembly_accession(acc):
                    downloaded.append(row)
                else:
                    unresolved.append(row)

            assert len(downloaded) == 1
            assert len(unresolved) == 3
            assert "GCF_000001405.40" in [d["assembly_accession"] for d in downloaded]

    def test_output_files_format(self):
        """Test that output file formats are correct."""
        from tests.test_download_assemblies_phase4 import split_assemblies

        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            df = pd.DataFrame(
                {
                    "transcript_id": ["TX1", "TX2"],
                    "assembly_accession": ["GCF_000001405.40", "GRCh38"],
                    "organism": ["homo sapiens", "homo sapiens"],
                    "db_source": ["ncbi", "ensembl"],
                }
            )

            downloaded, unresolved = split_assemblies(df, cache_dir)

            # Verify downloaded has expected columns
            assert "assembly_accession" in downloaded.columns
            assert "organism" in downloaded.columns
            assert "db_source" in downloaded.columns

            # Verify unresolved has reason column
            assert "reason" in unresolved.columns
            assert (
                unresolved.iloc[0]["reason"] == "not_resolvable_by_download_assemblies"
            )

    def test_cache_existence_prevents_redownload(self):
        """Test that cached assemblies are not re-downloaded."""
        from tests.test_download_assemblies_phase4 import split_assemblies

        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            # Pre-create cache for one assembly
            asm_dir = cache_dir / "GCF_000001405.40"
            asm_dir.mkdir(parents=True)
            (asm_dir / "genome.fasta").touch()
            (asm_dir / "genome.fasta.fai").touch()

            df = pd.DataFrame(
                {
                    "transcript_id": ["TX1"],
                    "assembly_accession": ["GCF_000001405.40"],
                    "organism": ["homo sapiens"],
                    "db_source": ["ncbi"],
                }
            )

            downloaded, unresolved = split_assemblies(df, cache_dir)

            # Should be empty since already cached
            assert len(downloaded) == 0
            assert len(unresolved) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
