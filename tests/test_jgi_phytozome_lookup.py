"""Unit tests for the pure selection logic in jgi_phytozome_lookup.

No network: assert on canned in-memory JGI file dicts modeled on real shapes.
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from jgi_phytozome_lookup import select_files


def _f(name, ftype, _id="deadbeef", status="RESTORED"):
    return {
        "_id": _id,
        "file_name": name,
        "file_status": status,
        "metadata": {"type": ftype},
    }


CANNED = [
    # Potato annotations: unversioned gene, unversioned gene_exons,
    # versioned gene_exons, plus a repeatmasked one to exclude.
    _f("Stuberosum_206_gene.gff3.gz", "annotation/gene"),
    _f("Stuberosum_206_gene_exons.gff3.gz", "annotation/gene"),
    _f("Stuberosum_206_v3.4.gene_exons.gff3.gz", "annotation/gene", _id="winner_ann"),
    _f("Stuberosum_206_repeatmasked.gff3.gz", "annotation/repeat"),
    # Assemblies: softmasked + hardmasked + unmasked.
    _f("Stuberosum_206_softmasked.fa.gz", "assembly"),
    _f("Stuberosum_206_hardmasked.fa.gz", "assembly"),
    _f("Stuberosum_206.fa.gz", "assembly", _id="winner_seq"),
]


def test_gene_exons_versioned_chosen():
    result = select_files(CANNED)
    assert result["annotation"]["_id"] == "winner_ann"
    assert result["annotation"]["file_name"] == "Stuberosum_206_v3.4.gene_exons.gff3.gz"


def test_gene_exons_preferred_over_plain_gene():
    files = [
        _f("Stuberosum_206_gene.gff3.gz", "annotation/gene"),
        _f("Stuberosum_206_gene_exons.gff3.gz", "annotation/gene", _id="exons"),
    ]
    assert select_files(files)["annotation"]["_id"] == "exons"


def test_versioned_preferred_over_unversioned():
    files = [
        _f("Stuberosum_206_gene_exons.gff3.gz", "annotation/gene"),
        _f("Stuberosum_206_v3.4.gene_exons.gff3.gz", "annotation/gene", _id="ver"),
    ]
    assert select_files(files)["annotation"]["_id"] == "ver"


def test_repeatmasked_excluded():
    files = [_f("Stuberosum_206_repeatmasked.gff3.gz", "annotation/repeat")]
    assert select_files(files)["annotation"] is None


def test_unmasked_assembly_chosen_over_masked():
    result = select_files(CANNED)
    assert result["sequence"]["_id"] == "winner_seq"
    assert "masked" not in result["sequence"]["file_name"]


def test_prefer_name_overrides_ranking():
    # Even though the versioned gene_exons ranks highest by default, an
    # explicit prefer_name pins the exact portal file the manifest names.
    result = select_files(CANNED, prefer_name="Stuberosum_206_gene.gff3.gz")
    assert result["annotation"]["file_name"] == "Stuberosum_206_gene.gff3.gz"


def test_prefer_name_case_insensitive_and_falls_back_when_absent():
    # Unknown prefer_name -> default ranking still wins (no crash).
    result = select_files(CANNED, prefer_name="does_not_exist.gff3.gz")
    assert result["annotation"]["_id"] == "winner_ann"


def test_none_when_category_absent():
    # Only an assembly -> annotation is None.
    files = [_f("Stuberosum_206.fa.gz", "assembly")]
    result = select_files(files)
    assert result["annotation"] is None
    assert result["sequence"] is not None

    # Only an annotation -> sequence is None.
    files = [_f("Stuberosum_206_v3.4.gene_exons.gff3.gz", "annotation/gene")]
    result = select_files(files)
    assert result["sequence"] is None
    assert result["annotation"] is not None
