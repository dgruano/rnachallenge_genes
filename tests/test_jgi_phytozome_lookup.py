"""Unit tests for the pure selection logic in jgi_phytozome_lookup.

No network: assert on canned in-memory JGI file dicts modeled on real shapes.
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from jgi_phytozome_lookup import annotation_stem, select_assembly_for, select_files


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


# ── version-matched assembly selection ───────────────────────


def test_annotation_stem_strips_variants():
    assert annotation_stem("Csinensis_154_v1.1.gene.gff3.gz") == "Csinensis_154_v1.1"
    assert (
        annotation_stem("Stuberosum_206_v3.4.gene_exons.gff3.gz")
        == "Stuberosum_206_v3.4"
    )
    assert annotation_stem("amborella.gff3.gz") == "amborella"
    assert annotation_stem(None) == ""


def test_select_assembly_matches_annotation_version_unmasked():
    # Two versions on disk; the assembly must mate the resolved annotation's
    # version (v1.1), and prefer the unmasked file over the softmasked one.
    files = [
        _f("Csinensis_154_v1.0.fa.gz", "assembly", _id="wrong_version"),
        _f("Csinensis_154_v1.1.softmasked.fa.gz", "assembly", _id="masked"),
        _f("Csinensis_154_v1.1.fa.gz", "assembly", _id="right"),
    ]
    picked = select_assembly_for(files, "Csinensis_154_v1.1.gene.gff3.gz")
    assert picked["_id"] == "right"


def test_select_assembly_version_boundary_not_prefix():
    # stem 'v1.1' must not swallow 'v1.10'.
    files = [
        _f("X_1_v1.10.fa.gz", "assembly", _id="v110"),
        _f("X_1_v1.1.fa.gz", "assembly", _id="v11"),
    ]
    assert select_assembly_for(files, "X_1_v1.1.gene.gff3.gz")["_id"] == "v11"


def test_select_assembly_falls_back_when_no_stem_match():
    # Annotation version has no assembly mate -> best-guess unmasked assembly.
    files = [
        _f("X_1_v9.9.softmasked.fa.gz", "assembly", _id="masked"),
        _f("X_1_v9.9.fa.gz", "assembly", _id="best"),
    ]
    picked = select_assembly_for(files, "X_1_v1.1.gene.gff3.gz")
    assert picked["_id"] == "best"


def test_select_assembly_none_when_no_assembly():
    files = [_f("X_1_v1.1.gene.gff3.gz", "annotation/gene")]
    assert select_assembly_for(files, "X_1_v1.1.gene.gff3.gz") is None
