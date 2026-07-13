#!/usr/bin/env python3
"""Audit files in results/ against outputs declared by current Snakemake rules.

Modes:
- list:   print orphan files (not produced by current rules)
- move:   move orphan files to old_results/
- remove: delete orphan files
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS_DIR = ROOT / "results"
OLD_RESULTS_DIR = ROOT / "old results"


def gather_declared_outputs() -> tuple[set[str], bool]:
    """Return declared output patterns under results/ and whether species wildcard exists.

    This parser is intentionally lightweight and pattern-based, matching the
    conventions used in current .smk files.
    """
    files = sorted((ROOT / "workflow" / "rules").glob("*.smk")) + [ROOT / "Snakefile"]

    # f"{RESULTS}/foo.tsv"
    pat_results_f = re.compile(r'f"\{RESULTS\}/([^"]+)"')
    # "results/foo.tsv" literals
    pat_results_lit = re.compile(r'"results/([^"]+)"')

    declared: set[str] = set()
    has_species_wildcard = False

    for path in files:
        text = path.read_text(encoding="utf-8")
        for m in pat_results_f.finditer(text):
            rel = m.group(1)
            declared.add(rel)
            if "{species}" in rel or "{{species}}" in rel:
                has_species_wildcard = True
        for m in pat_results_lit.finditer(text):
            rel = m.group(1)
            declared.add(rel)
            if "{species}" in rel or "{{species}}" in rel:
                has_species_wildcard = True

    return declared, has_species_wildcard


def is_declared(
    rel_path: str, declared_patterns: set[str], has_species_wildcard: bool
) -> bool:
    """Check whether a results-relative file path is declared by current rules."""
    if rel_path in declared_patterns:
        return True

    # Handle biomart/{species}.tsv.gz wildcard without full Snakemake parsing.
    if has_species_wildcard:
        if rel_path.startswith("biomart/") and rel_path.endswith(".tsv.gz"):
            return True

    # Side-effect file produced by current script (not declared in rule output).
    if rel_path == "ensembl_unknown_prefixes_ACTION_REQUIRED.txt":
        return True

    return False


def orphan_files(results_dir: Path) -> list[Path]:
    """Return files in results/ that are not declared outputs."""
    declared, has_species_wildcard = gather_declared_outputs()
    files = sorted(p for p in results_dir.rglob("*") if p.is_file())

    orphans: list[Path] = []
    for p in files:
        rel = p.relative_to(results_dir).as_posix()
        if not is_declared(rel, declared, has_species_wildcard):
            orphans.append(p)
    return orphans


def cmd_list(files: list[Path], results_dir: Path) -> int:
    if not files:
        print("No orphan files found in results/.")
        return 0
    print(f"Found {len(files)} orphan file(s):")
    for p in files:
        print(f"- results/{p.relative_to(results_dir).as_posix()}")
    return 0


def cmd_move(files: list[Path], results_dir: Path, old_results_dir: Path) -> int:
    if not files:
        print("No orphan files to move.")
        return 0

    moved = 0
    for src in files:
        rel = src.relative_to(results_dir)
        dst = old_results_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dst))
        moved += 1
        print(f"moved: results/{rel.as_posix()} -> old_results/{rel.as_posix()}")

    print(f"Moved {moved} file(s) to old_results/.")
    return 0


def cmd_remove(files: list[Path], results_dir: Path) -> int:
    if not files:
        print("No orphan files to remove.")
        return 0

    removed = 0
    for p in files:
        rel = p.relative_to(results_dir)
        p.unlink()
        removed += 1
        print(f"removed: results/{rel.as_posix()}")

    print(f"Removed {removed} file(s).")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit results/ for files not produced by current rule outputs."
    )
    parser.add_argument(
        "action",
        choices=["list", "move", "remove"],
        help="What to do with orphan files",
    )
    parser.add_argument(
        "--results-dir",
        default=str(RESULTS_DIR),
        help="Path to results directory (default: ./results)",
    )
    parser.add_argument(
        "--old-results-dir",
        default=str(OLD_RESULTS_DIR),
        help="Path to old results directory for move action (default: ./old_results)",
    )

    args = parser.parse_args()

    results_dir = Path(args.results_dir).resolve()
    old_results_dir = Path(args.old_results_dir).resolve()

    if not results_dir.exists() or not results_dir.is_dir():
        raise SystemExit(f"results dir not found: {results_dir}")

    files = orphan_files(results_dir)

    if args.action == "list":
        return cmd_list(files, results_dir)
    if args.action == "move":
        return cmd_move(files, results_dir, old_results_dir)
    return cmd_remove(files, results_dir)


if __name__ == "__main__":
    raise SystemExit(main())
