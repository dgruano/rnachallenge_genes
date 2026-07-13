"""Tests for SGD Source-tagged ID canonicalization (ROI strategy #6).

Input IDs of the form ``Source:SGD;Acc:S000028522`` never matched the SGD
GFF3 index because (a) candidate generation didn't extract the SGDID and
(b) the index didn't key on the ``dbxref`` attribute (``SGD:S000028522``).
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow"))
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from utils.annotation_resolver import build_annotation_index, yeast_candidates

# One real gene line from the SGD GFF3 (YJL127C-B / MCO6 / SGD:S000028522).
_GFF_LINE = (
    "chrX\tSGD\tgene\t181275\t183844\t.\t-\t.\t"
    "ID=YJL127C-B;Name=YJL127C-B;gene=MCO6;Alias=MCO6,YJL127C-A;"
    "dbxref=SGD:S000028522;curie=SGD:S000028522\n"
)


class TestYeastCandidates:
    def test_extracts_sgdid_from_source_tag(self):
        cands = yeast_candidates("Source:SGD;Acc:S000028522")
        assert "S000028522" in cands
        assert "SGD:S000028522" in cands

    def test_plain_systematic_name_untouched(self):
        cands = yeast_candidates("YDL184C")
        assert cands[0] == "YDL184C"
        assert not any(c.startswith("SGD:") for c in cands)


class TestSourceTagResolvesAgainstIndex:
    def test_source_tag_hits_dbxref_key(self, tmp_path):
        gff = tmp_path / "sgd.gff"
        gff.write_text(_GFF_LINE)
        index = build_annotation_index(
            str(gff),
            feature_types={"gene"},
            transcript_fields=(),
            gene_id_fields=("ID",),
            gene_symbol_fields=("Name",),
            alias_fields=("Alias", "dbxref"),
        )
        # The canonicalized candidate must land on a real record.
        hit = next(
            (
                index[c]
                for c in yeast_candidates("Source:SGD;Acc:S000028522")
                if c in index
            ),
            None,
        )
        assert hit is not None
        assert hit["chrom"] == "chrX"
        assert hit["start"] == 181275
        assert hit["gene_symbol"] == "YJL127C-B"
