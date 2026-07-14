"""Tests for chrom_translation.load_chrom_translation (audit breakage #1)."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.parent / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from chrom_translation import load_chrom_translation, resolve_chrom_key

FIXTURE = Path(__file__).parent / "data" / "GCF_000001405.40_assembly_report.txt"
MOLECULE_FIXTURE = Path(__file__).parent / "data" / "assigned_molecule_report.txt"


class TestLoadChromTranslation:
    def test_sequence_name_alias(self):
        assert "NC_000001.11" in load_chrom_translation(FIXTURE)["1"]

    def test_ucsc_alias(self):
        assert "NC_000001.11" in load_chrom_translation(FIXTURE)["chr1"]

    def test_genbank_alias(self):
        assert "NC_000001.11" in load_chrom_translation(FIXTURE)["CM000663.2"]

    def test_genbank_seqid_candidate(self):
        # Ensembl/GenBank toplevel FASTAs use the CM_ seqid, not RefSeq — the
        # friendly name must offer it as a candidate too (tomato SL3.0 case).
        assert "CM000663.2" in load_chrom_translation(FIXTURE)["1"]

    def test_refseq_identity(self):
        assert "NC_000001.11" in load_chrom_translation(FIXTURE)["NC_000001.11"]

    def test_na_refseq_rows_excluded(self):
        xlate = load_chrom_translation(FIXTURE)
        assert "HSCHR1_CTG1_UNLOCALIZED" not in xlate
        assert "KI270706.1" not in xlate

    def test_case_insensitive_lookup(self):
        # Map carries lowercased keys; callers lowercase their query to match.
        xlate = load_chrom_translation(FIXTURE)
        assert "NC_000001.11" in xlate["CHR1".lower()]

    def test_assigned_molecule_bare_integer(self):
        # Sequence-Name is "Chr1" but resolvers emit bare "1" (Assigned-Molecule).
        assert "NC_052532.1" in load_chrom_translation(MOLECULE_FIXTURE)["1"]

    def test_assigned_molecule_does_not_clobber_chromosome(self):
        # A scaffold shares Assigned-Molecule "1" but must not overwrite the
        # real chromosome's RefSeq mapping.
        xlate = load_chrom_translation(MOLECULE_FIXTURE)
        assert "NC_052532.1" in xlate["1"]  # not the scaffold's NW_020110099.1
        assert "NW_020110099.1" not in xlate["1"]

    def test_mito_case_insensitive(self):
        assert "NC_053523.1" in load_chrom_translation(MOLECULE_FIXTURE)["chrmt"]

    def test_missing_file_returns_empty(self):
        assert load_chrom_translation(Path("/nonexistent/report.txt")) == {}


class TestResolveChromKey:
    # Report maps bare '1'/'V' -> [RefSeq, GenBank]; .fai keys are the seqids.
    XLATE = {
        "1": ["NC_1", "CM_1"],
        "v": ["NC_5"],
        "V": ["NC_5"],
        "mt": ["NC_MT"],
        "MT": ["NC_MT"],
    }
    FAI = {"NC_1", "NC_5", "NC_MT", "5"}

    def test_genbank_seqid_fallback(self):
        # RefSeq 'NC_9' absent from .fai; GenBank 'CM_9' present → pick GenBank.
        assert resolve_chrom_key("1", {"1": ["NC_9", "CM_9"]}, {"CM_9"}) == "CM_9"

    def test_chr_prefixed_arabic_via_report(self):
        # 'chr1' -> strip 'chr' -> '1' -> report -> NC_1 (the .fai seqid).
        assert resolve_chrom_key("chr1", self.XLATE, self.FAI) == "NC_1"

    def test_chr_prefixed_roman_via_report(self):
        # 'chrV' -> 'V' -> report -> NC_5. This is the C. elegans case.
        assert resolve_chrom_key("chrV", self.XLATE, self.FAI) == "NC_5"

    def test_case_insensitive_mito(self):
        # Resolver emits bare 'mt'; report Assigned-Molecule is 'MT'.
        assert resolve_chrom_key("mt", self.XLATE, self.FAI) == "NC_MT"

    def test_raw_name_in_fai_no_report(self):
        # Ensembl-style FASTA: empty map, seqid '5' is literally in the .fai.
        assert resolve_chrom_key("chr5", {}, self.FAI) == "5"

    def test_unresolvable_returns_none(self):
        assert resolve_chrom_key("chrZ", self.XLATE, self.FAI) is None
