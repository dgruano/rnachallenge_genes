"""
tests/test_download_assemblies_phase4.py

Unit tests for Phase 4: Download Assemblies (Simplified Option B)

Tests the simplified download_assemblies that:
- Only handles GCF_/GCA_ accessions
- Validates cache existence
- Gracefully handles non-GCF_/GCA_ accessions by marking them unresolved
- Produces two outputs: downloaded_assemblies.tsv and unresolved_assemblies.tsv
"""

import sys
import tempfile
from pathlib import Path
from typing import Tuple

import pandas as pd
import pytest

# Add workflow scripts to path
SCRIPT_DIR = Path(__file__).parent.parent / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))


# ============================================================================
# Helper Functions (to be extracted into download_assemblies.py)
# ============================================================================

def is_ncbi_assembly_accession(accession: str) -> bool:
    """
    Check if accession is a downloadable NCBI assembly accession (GCF_/GCA_).
    """
    try:
        if pd.isna(accession):
            return False
    except (TypeError, ValueError):
        pass
    if not accession:
        return False
    acc_str = str(accession).strip()
    return acc_str.startswith(("GCF_", "GCA_"))


def should_download_assembly(row, cache_dir: Path) -> bool:
    """
    Determine if an assembly should be downloaded based on:
    - Must be GCF_/GCA_ accession
    - Must not be cached already
    """
    if not is_ncbi_assembly_accession(row["assembly_accession"]):
        return False

    cached = cache_dir / row["assembly_accession"] / "genomic.gtf.gz"
    return not cached.exists()


def split_assemblies(
    df: pd.DataFrame, cache_dir: Path
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split input dataframe into:
    - downloaded: GCF_/GCA_ accessions not in cache
    - unresolved: Non-GCF_/GCA_ accessions

    Returns (downloaded, unresolved)
    """
    # Handle empty or incomplete DataFrames
    if df.empty or "assembly_accession" not in df.columns:
        return pd.DataFrame(), pd.DataFrame()

    # Get unique assemblies
    unique_asm = (
        df[["assembly_accession", "organism", "db_source"]]
        .dropna(subset=["assembly_accession"])
        .drop_duplicates(subset="assembly_accession")
    )

    downloaded = []
    unresolved = []

    for _, row in unique_asm.iterrows():
        if is_ncbi_assembly_accession(row["assembly_accession"]):
            cached = cache_dir / row["assembly_accession"] / "genomic.gtf.gz"
            if not cached.exists():
                downloaded.append(row)
            # else: already cached, skip (not in either output)
        else:
            # Non-GCF_/GCA_ accession
            unresolved_row = row.copy()
            unresolved_row["reason"] = "not_resolvable_by_download_assemblies"
            unresolved.append(unresolved_row)

    downloaded_df = (
        pd.DataFrame(downloaded) if downloaded else pd.DataFrame()
    )
    unresolved_df = (
        pd.DataFrame(unresolved) if unresolved else pd.DataFrame()
    )

    return downloaded_df, unresolved_df


# ============================================================================
# Test Cases
# ============================================================================

class TestAccessionDetection:
    """Test detection of NCBI assembly accessions."""

    def test_detect_gcf_accession(self):
        """Test detection of GCF_ accession."""
        assert is_ncbi_assembly_accession("GCF_000001405.40") is True

    def test_detect_gca_accession(self):
        """Test detection of GCA_ accession."""
        assert is_ncbi_assembly_accession("GCA_000001405.1") is True

    def test_reject_nc_accession(self):
        """Test that NC_ sequence accessions are not detected."""
        assert is_ncbi_assembly_accession("NC_000001.11") is False

    def test_reject_nw_accession(self):
        """Test that NW_ scaffold accessions are not detected."""
        assert is_ncbi_assembly_accession("NW_017955610.1") is False

    def test_reject_ucsc_name(self):
        """Test that UCSC names are not detected."""
        assert is_ncbi_assembly_accession("hg38") is False

    def test_reject_grc_name(self):
        """Test that GRC assembly names are not detected."""
        assert is_ncbi_assembly_accession("GRCh38") is False

    def test_reject_empty_string(self):
        """Test that empty string returns False."""
        assert is_ncbi_assembly_accession("") is False

    def test_reject_none(self):
        """Test that None returns False."""
        assert is_ncbi_assembly_accession(None) is False

    def test_reject_nan(self):
        """Test that NaN/pd.NA returns False."""
        assert is_ncbi_assembly_accession(pd.NA) is False
        assert is_ncbi_assembly_accession(float("nan")) is False

    def test_reject_whitespace_only(self):
        """Test that whitespace-only string returns False."""
        assert is_ncbi_assembly_accession("   ") is False


class TestCacheChecking:
    """Test cache existence checking."""

    def test_should_download_when_not_cached(self):
        """Test that assembly is marked for download when not cached."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            row = {
                "assembly_accession": "GCF_000001405.40",
                "organism": "homo sapiens",
                "db_source": "ncbi",
            }
            assert should_download_assembly(row, cache_dir) is True

    def test_should_not_download_when_cached(self):
        """Test that assembly is skipped when already cached."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            asm_dir = cache_dir / "GCF_000001405.40"
            asm_dir.mkdir(parents=True, exist_ok=True)
            (asm_dir / "genomic.gtf.gz").touch()

            row = {
                "assembly_accession": "GCF_000001405.40",
                "organism": "homo sapiens",
                "db_source": "ncbi",
            }
            assert should_download_assembly(row, cache_dir) is False

    def test_should_not_download_non_ncbi_accession(self):
        """Test that non-NCBI accessions are never marked for download."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            row = {
                "assembly_accession": "GRCh38",
                "organism": "homo sapiens",
                "db_source": "ensembl",
            }
            assert should_download_assembly(row, cache_dir) is False

    def test_should_not_download_sequence_accession(self):
        """Test that sequence accessions are not marked for download."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            row = {
                "assembly_accession": "NC_000001.11",
                "organism": "homo sapiens",
                "db_source": "ncbi",
            }
            assert should_download_assembly(row, cache_dir) is False


class TestAssemblySplitting:
    """Test splitting assemblies into downloadable and unresolved."""

    def test_empty_dataframe(self):
        """Test handling of empty input."""
        with tempfile.TemporaryDirectory() as tmpdir:
            df = pd.DataFrame()
            downloaded, unresolved = split_assemblies(df, Path(tmpdir))
            assert len(downloaded) == 0
            assert len(unresolved) == 0

    def test_all_gcf_not_cached(self):
        """Test all GCF_ accessions go to downloaded when not cached."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            df = pd.DataFrame({
                "transcript_id": ["TX1", "TX2"],
                "assembly_accession": ["GCF_000001405.40", "GCF_000002035.6"],
                "organism": ["homo sapiens", "danio rerio"],
                "db_source": ["ncbi", "ncbi"],
            })
            downloaded, unresolved = split_assemblies(df, cache_dir)

            assert len(downloaded) == 2
            assert len(unresolved) == 0
            assert list(downloaded["assembly_accession"]) == [
                "GCF_000001405.40",
                "GCF_000002035.6",
            ]

    def test_all_gca_not_cached(self):
        """Test GCA_ accessions are handled same as GCF_."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            df = pd.DataFrame({
                "transcript_id": ["TX1"],
                "assembly_accession": ["GCA_000001405.1"],
                "organism": ["homo sapiens"],
                "db_source": ["ncbi"],
            })
            downloaded, unresolved = split_assemblies(df, cache_dir)

            assert len(downloaded) == 1
            assert len(unresolved) == 0

    def test_all_non_ncbi_go_to_unresolved(self):
        """Test non-NCBI accessions go to unresolved."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            df = pd.DataFrame({
                "transcript_id": ["TX1", "TX2", "TX3"],
                "assembly_accession": ["GRCh38", "NC_000001.11", "hg38"],
                "organism": ["homo sapiens", "homo sapiens", "homo sapiens"],
                "db_source": ["ensembl", "ncbi", "ucsc"],
            })
            downloaded, unresolved = split_assemblies(df, cache_dir)

            assert len(downloaded) == 0
            assert len(unresolved) == 3
            for row in unresolved.itertuples():
                assert row.reason == "not_resolvable_by_download_assemblies"

    def test_mixed_assemblies(self):
        """Test correct splitting of mixed accessions."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            df = pd.DataFrame({
                "transcript_id": ["TX1", "TX2", "TX3", "TX4"],
                "assembly_accession": [
                    "GCF_000001405.40",  # NCBI, not cached → download
                    "GRCh38",             # Non-NCBI → unresolved
                    "GCA_000001405.1",   # NCBI, not cached → download
                    "NC_000001.11",       # Sequence accession → unresolved
                ],
                "organism": ["homo sapiens", "homo sapiens", "homo sapiens", "homo sapiens"],
                "db_source": ["ncbi", "ensembl", "ncbi", "ncbi"],
            })
            downloaded, unresolved = split_assemblies(df, cache_dir)

            assert len(downloaded) == 2
            assert len(unresolved) == 2
            assert "GCF_000001405.40" in downloaded["assembly_accession"].values
            assert "GCA_000001405.1" in downloaded["assembly_accession"].values
            assert "GRCh38" in unresolved["assembly_accession"].values
            assert "NC_000001.11" in unresolved["assembly_accession"].values

    def test_cached_assemblies_excluded(self):
        """Test that cached NCBI assemblies are not in either output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            # Create cache for one assembly
            asm_dir = cache_dir / "GCF_000001405.40"
            asm_dir.mkdir(parents=True, exist_ok=True)
            (asm_dir / "genomic.gtf.gz").touch()

            df = pd.DataFrame({
                "transcript_id": ["TX1", "TX2"],
                "assembly_accession": ["GCF_000001405.40", "GCF_000002035.6"],
                "organism": ["homo sapiens", "danio rerio"],
                "db_source": ["ncbi", "ncbi"],
            })
            downloaded, unresolved = split_assemblies(df, cache_dir)

            # Only uncached assembly should be in downloaded
            assert len(downloaded) == 1
            assert len(unresolved) == 0
            assert downloaded.iloc[0]["assembly_accession"] == "GCF_000002035.6"

    def test_unresolved_rows_have_reason(self):
        """Test that unresolved rows include reason field."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            df = pd.DataFrame({
                "transcript_id": ["TX1"],
                "assembly_accession": ["GRCh38"],
                "organism": ["homo sapiens"],
                "db_source": ["ensembl"],
            })
            _, unresolved = split_assemblies(df, cache_dir)

            assert "reason" in unresolved.columns
            assert unresolved.iloc[0]["reason"] == "not_resolvable_by_download_assemblies"

    def test_duplicates_removed_before_split(self):
        """Test that duplicate assemblies are deduplicated."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            df = pd.DataFrame({
                "transcript_id": ["TX1", "TX2", "TX3"],
                "assembly_accession": [
                    "GCF_000001405.40",
                    "GCF_000001405.40",
                    "GCF_000001405.40",
                ],
                "organism": ["homo sapiens", "homo sapiens", "homo sapiens"],
                "db_source": ["ncbi", "ncbi", "ncbi"],
            })
            downloaded, unresolved = split_assemblies(df, cache_dir)

            # Should be only one unique assembly
            assert len(downloaded) == 1
            assert len(unresolved) == 0

    def test_missing_assembly_accession_column(self):
        """Test handling of missing assembly_accession."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            df = pd.DataFrame({
                "transcript_id": ["TX1"],
                "organism": ["homo sapiens"],
                "db_source": ["ncbi"],
            })
            # Should not crash, just have no assemblies
            downloaded, unresolved = split_assemblies(df, cache_dir)
            assert len(downloaded) == 0
            assert len(unresolved) == 0

    def test_null_assembly_accession_values(self):
        """Test handling of NaN/None in assembly_accession."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            df = pd.DataFrame({
                "transcript_id": ["TX1", "TX2"],
                "assembly_accession": [None, pd.NA],
                "organism": ["homo sapiens", "homo sapiens"],
                "db_source": ["ncbi", "ncbi"],
            })
            downloaded, unresolved = split_assemblies(df, cache_dir)

            # Null values should be dropped, so no output
            assert len(downloaded) == 0
            assert len(unresolved) == 0


class TestOutputFormat:
    """Test output TSV format for downloaded and unresolved assemblies."""

    def test_downloaded_tsv_has_required_columns(self):
        """Test that downloaded TSV has required columns."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            df = pd.DataFrame({
                "transcript_id": ["TX1"],
                "assembly_accession": ["GCF_000001405.40"],
                "organism": ["homo sapiens"],
                "db_source": ["ncbi"],
            })
            downloaded, _ = split_assemblies(df, cache_dir)

            required = ["assembly_accession", "organism", "db_source"]
            for col in required:
                assert col in downloaded.columns

    def test_unresolved_tsv_has_reason_column(self):
        """Test that unresolved TSV has reason column."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            df = pd.DataFrame({
                "transcript_id": ["TX1"],
                "assembly_accession": ["GRCh38"],
                "organism": ["homo sapiens"],
                "db_source": ["ensembl"],
            })
            _, unresolved = split_assemblies(df, cache_dir)

            assert "reason" in unresolved.columns


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
