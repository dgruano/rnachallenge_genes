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
from typing import Optional

# assembly_report.txt tab-separated columns (0-based):
#   0 Sequence-Name  1 Sequence-Role  2 Assigned-Molecule
#   4 GenBank-Accn   6 RefSeq-Accn    9 UCSC-style-name
_ALIAS_COLS = (0, 4, 9, 6)  # 6 included for RefSeq -> RefSeq identity
_REFSEQ_COL = 6
_GENBANK_COL = 4
_ROLE_COL = 1
_MOLECULE_COL = 2  # bare "1"/"MT"; only trustworthy on assembled-molecule rows


def load_chrom_translation(report_path) -> dict[str, list[str]]:
    """Parse an NCBI *_assembly_report.txt into {alias -> [candidate seqids]}.

    Keys include Sequence-Name, GenBank-Accn, UCSC-style-name, and the
    accessions themselves (identity). Each alias maps to BOTH the RefSeq-Accn
    and the GenBank-Accn: NCBI RefSeq FASTAs use the NC_/NW_ seqids, but an
    Ensembl/GenBank toplevel FASTA of the same GCA_ assembly uses the CM_/GL_
    GenBank seqids — resolve_chrom_key picks whichever is actually in the .fai.
    Skips '#' comments and rows whose RefSeq-Accn is 'na'. Returns {} if the
    file is missing so callers fall back to the existing chr-prefix toggle.
    """
    report_path = Path(report_path)
    if not report_path.exists():
        return {}

    xlate: dict[str, list[str]] = {}
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
            genbank = cols[_GENBANK_COL].strip()
            targets = [t for t in (refseq, genbank) if t and t != "na"]
            aliases = [cols[i].strip() for i in _ALIAS_COLS]
            # Assigned-Molecule ("1"/"MT") only on assembled-molecule rows —
            # on scaffolds it repeats the chromosome and would clobber it.
            if cols[_ROLE_COL].strip() == "assembled-molecule":
                aliases.append(cols[_MOLECULE_COL].strip())
            for alias in aliases:
                if alias and alias != "na":
                    xlate[alias] = targets
                    xlate[alias.lower()] = targets  # case-insensitive lookup
    return xlate


def resolve_chrom_key(chrom: str, xlate: dict, chrom_lengths) -> Optional[str]:
    """Resolve a friendly chrom name to an actual .fai seqid, or None.

    Tries the name as-is and with the 'chr' prefix toggled; for each, the
    report map (friendly -> [RefSeq, GenBank]) is consulted first, returning the
    first candidate seqid present in `chrom_lengths`, then the raw name is
    checked against `chrom_lengths`. Catches e.g. 'chr5'/'chrV' whose bare form
    ('5'/'V') is in the report but whose 'chr' prefix hides it from both.
    """
    candidates = [chrom, chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"]
    for cand in candidates:
        mapped = xlate.get(cand) or xlate.get(cand.lower()) or ()
        for seqid in mapped:
            if seqid in chrom_lengths:
                return seqid
        if cand in chrom_lengths:
            return cand
    return None
