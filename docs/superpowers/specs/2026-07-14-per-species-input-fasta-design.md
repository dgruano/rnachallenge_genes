# Per-species input FASTA split — design

**Date:** 2026-07-14
**Status:** Approved, pre-implementation

## Goal

Produce one FASTA file per species containing the **original input transcript
sequences** (verbatim from `config["input_fastas"]`), grouped by the species
that resolution assigned to each transcript.

This is deliberately built from the pipeline **input**, not `output.fasta`:
the input records are the sequences of interest, and sourcing species from the
resolution table captures every transcript that resolved to a species —
including those that later failed flanking-sequence extraction (which
`output.fasta` drops).

## Data flow

```
config["input_fastas"]  ─┐
                         ├─►  split_input_by_species  ─►  results/sequences_by_species/<species>.fasta
results/ncbi_chromosome_resolved.tsv (transcript_id → organism) ─┘
```

No dependency on the download/extract stages — this runs straight off the
resolution table, so per-species files are buildable without the multi-hour
assembly download.

## Species key

The `organism` column (col 5) is the only populated species field (the
`species` column is empty). It is inconsistent across sources:

- casing: `Homo sapiens` vs `homo sapiens`
- separators: `Arabidopsis thaliana` vs `arabidopsis_thaliana`
- extra tokens: `Oryza sativa Japonica Group`, `Drosophila yakuba (flies)`

**Normalization** (applied to build the file key):
1. `strip().lower()`
2. collapse any run of non-alphanumeric characters to a single `_`
3. strip leading/trailing `_`

Examples:
- `Arabidopsis thaliana`, `arabidopsis_thaliana` → `arabidopsis_thaliana`
- `Homo sapiens`, `homo sapiens` → `homo_sapiens`
- `Danio rerio`, `danio rerio` → `danio_rerio`
- `Drosophila yakuba (flies)` → `drosophila_yakuba_flies`

Keys are filename-safe by construction.

## Rule

`split_input_by_species`, appended to
`workflow/rules/extract_sequences.smk`.

- **input:**
  - `resolved = f"{RESULTS}/ncbi_chromosome_resolved.tsv"`
  - `inputs = config["input_fastas"]` (a list)
- **output:** `directory(f"{RESULTS}/sequences_by_species")`
  (dynamic species set → `directory()` output, no checkpoint since nothing
  downstream consumes the individual files as Snakemake targets)
- **wired into** `rule all` in `Snakefile` as an additional target.
- Modest resources (single core, minutes) — it is a streaming text pass.

`results/output.fasta` is untouched (additive change). `report.smk` is
unaffected.

## Script

`workflow/scripts/split_input_by_species.py` (~35 lines).

1. Read the resolved TSV; build `transcript_id → normalized_species`.
   First organism seen wins on duplicate transcript IDs (ambiguous rows).
2. Ensure output directory exists.
3. Stream each input FASTA (Biopython `SeqIO.parse`); for each record, look up
   the species by `record.id`. Append the record to that bucket's file.
4. IDs not present in the resolved table → `_unresolved.fasta` (nothing is
   silently dropped).
5. Log per-species counts and the unresolved count.

Open file handles are kept in a dict keyed by species and closed at the end
(species count is small — ~90 — so this is safe and avoids reopen/append
churn).

## Testing

`__main__` self-check (assert-based, no framework) that:
- the normalization function maps the known variant pairs above to a single
  key each;
- `Drosophila yakuba (flies)` → `drosophila_yakuba_flies`.

## Explicitly out of scope (YAGNI)

- Checkpoint / fan-out over per-species files — add only if a later stage needs
  to iterate them as Snakemake targets.
- A config toggle for output location or unresolved handling — no value varies.
- Re-deriving species from `output.fasta` headers — rejected in favor of the
  resolution table (more complete, extraction-independent).
