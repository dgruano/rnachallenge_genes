"""
scripts/chrom_translation.py

Translate friendly chromosome names to the RefSeq accessions used as seqids in
NCBI genome FASTAs (e.g. "1"/"chr1"/"CM000663.2" -> "NC_000001.11").

NCBI ships a `*_assembly_report.txt` next to every genome; download_assembly.py
caches it as `resources/cache/<accession>/assembly_report.txt`. This parses it
into an alias -> RefSeq-Accn map. Gated to NCBI GCF_/GCA_ assemblies only —
Ensembl-style FASTAs (worm/yeast/fly/plant) already carry friendly seqids.
"""

from pathlib import Path

# assembly_report.txt tab-separated columns (0-based):
#   0 Sequence-Name  4 GenBank-Accn  6 RefSeq-Accn  9 UCSC-style-name
_ALIAS_COLS = (0, 4, 9, 6)  # 6 included for RefSeq -> RefSeq identity
_REFSEQ_COL = 6


def load_chrom_translation(report_path) -> dict[str, str]:
    """Parse an NCBI *_assembly_report.txt into {alias -> RefSeq-Accn}.

    Keys include Sequence-Name, GenBank-Accn, UCSC-style-name, and the
    RefSeq-Accn itself (identity). Skips '#' comments and rows whose RefSeq-Accn
    is 'na'. Returns {} if the file is missing so callers fall back to the
    existing chr-prefix toggle.
    """
    report_path = Path(report_path)
    if not report_path.exists():
        return {}

    xlate: dict[str, str] = {}
    with open(report_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= _REFSEQ_COL:
                continue
            refseq = cols[_REFSEQ_COL].strip()
            if not refseq or refseq == "na":
                continue
            for i in _ALIAS_COLS:
                alias = cols[i].strip()
                if alias and alias != "na":
                    xlate[alias] = refseq
    return xlate
