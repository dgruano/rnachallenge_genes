"""Tests for chrom_translation.load_chrom_translation (audit breakage #1)."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.parent / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from chrom_translation import load_chrom_translation

FIXTURE = Path(__file__).parent / "data" / "GCF_000001405.40_assembly_report.txt"


class TestLoadChromTranslation:
    def test_sequence_name_alias(self):
        assert load_chrom_translation(FIXTURE)["1"] == "NC_000001.11"

    def test_ucsc_alias(self):
        assert load_chrom_translation(FIXTURE)["chr1"] == "NC_000001.11"

    def test_genbank_alias(self):
        assert load_chrom_translation(FIXTURE)["CM000663.2"] == "NC_000001.11"

    def test_refseq_identity(self):
        assert load_chrom_translation(FIXTURE)["NC_000001.11"] == "NC_000001.11"

    def test_na_refseq_rows_excluded(self):
        xlate = load_chrom_translation(FIXTURE)
        assert "HSCHR1_CTG1_UNLOCALIZED" not in xlate
        assert "KI270706.1" not in xlate

    def test_missing_file_returns_empty(self):
        assert load_chrom_translation(Path("/nonexistent/report.txt")) == {}
