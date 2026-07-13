# Pipeline Audit — RNA Flanking Sequence Pipeline

_Date: 2026-07-12 · Branch: `main` · HEAD: `bab60c9` (WIP: merge resolved)_

> **📌 This document is the source of truth for pipeline health.**
> Every agent (or person) who investigates or fixes anything **must finish by updating this file** — add findings to the right section (Things That Break / Overcomplicated Patterns / Legitimate Data Issues), tick or renumber the Actionables, and refresh Next Steps. Leave the audit consistent with reality before ending your turn. If a prior conclusion is wrong, mark it superseded rather than deleting the history.

## TL;DR — UPDATED 2026-07-12 POST-FIX

**Status: Extraction working, four blockers/refinements fixed.**
- **Blocker #1 (chrom translation)** ✅ DONE — assembly reports fetched, chromosome names translated
- **Blocker #2 (extract file wiring)** ✅ DONE — extract now reads remapped `ncbi_chromosome_resolved.tsv`
- **Blocker #3 (backwards coordinates)** ✅ DONE — defensive swap added, sequence_error collapsed from 1,658 → 1
- **Cleanup #5 (coordinate-less resolvers)** ✅ DONE — gramene now returns coords, NONCODE2016 reclassified to matched_not_found, merge guards invalid resolved rows

**Current extraction:** **13,199 sequences** from 15,544 resolved IDs (~85% extraction rate) after landing bottleneck fix #10 (manifest-driven URL downloads) and rerunning. Up from 12,163.

Remaining failures are now mostly legitimate data gaps. The only sizable unresolved download bucket is **phytozome (203)**, previously deferred due missing genome FASTA URLs/auth constraints — **JGI auth is now implemented (2026-07-12), see [Phytozome/JGI unblock](#phytozomejgi-unblock-2026-07-12).**

---

## Phytozome/JGI unblock (2026-07-12)

The phytozome deferral had two distinct blockers; the auth one is now closed:

1. **Annotation (GFF3) download — auth. ✅ UNBLOCKED.** `download_phytozome_gtf` ([resolve_phytozome_gtf.smk](workflow/rules/resolve_phytozome_gtf.smk)) now loads a JGI bearer token from `.env` (`JGI_SESSION_TOKEN` or `PHYTOZOME_BEARER`) and fails loudly if absent. Five new species wired into [config/phytozome_gtf_sources.yaml](config/phytozome_gtf_sources.yaml) + [resources/phytozome/manifest.json](resources/phytozome/manifest.json): amborella, chlamydomonas, physcomitrella, **vitis_vinifera** (`VIT_`/`GSVIVT`/`GTVIVG`), **solanum_tuberosum** (`PGSC`). `parse_ids.py` + `PREFIX_TO_SPECIES` route the new prefixes. This expands **transcript→coordinate resolution**, not extraction.
2. **Genome FASTA for extraction — still open.** The 203 `assembly_not_cached` phytozome rows resolve to coordinates but can't be sliced: the manifest downloads GFF3 annotations, not genome FASTAs. To extract, add genome-FASTA entries (same JGI token mechanism) to the manifest/config, then re-run download+extract. Until then these remain resolved-but-not-extracted.

**Net:** the token fix grows the *resolved* set (new plant/algae/moss/grape/potato transcripts get coordinates); closing the extraction gap for phytozome is a separate, now-tractable follow-up (JGI FASTA URLs behind the same auth).

### Download id + restore fix (2026-07-12, later — makes downloads actually work)

The token wiring above was necessary but not sufficient: `download_phytozome_gtf` still `401`/`404`'d because it built the URL from the manifest's numeric portal `file_id`, which is **not** a JGI download id. Root-caused and fixed:

- **Correct id resolution.** JGI's `download_files/{_id}/` needs the Mongo `_id` from the file-list API, not `file_id`. The rule now resolves it at download time from the `genome_id` in [config/phytozome_gtf_sources.yaml](config/phytozome_gtf_sources.yaml) via `resolve_annotation()`, pinning the manifest's `portal_file_name`. New helpers in [jgi_phytozome_lookup.py](workflow/scripts/jgi_phytozome_lookup.py); `per_page` capped at 50 (JGI 400s above that). Verified: Amborella downloads a valid 3.5 MB GFF3.
- **Run block → script.** The rule's inline `run:` block moved to [workflow/scripts/download_phytozome_gtf.py](workflow/scripts/download_phytozome_gtf.py) (`script:` directive, repo convention).
- **PURGED → auto restore.** On-tape files now POST to `request_archived_files/` and fail with the request id + rerun instructions (safe to repeat). Verified live: physcomitrella → `request_id 652368`. Of the 4 manifest species, amborella + chlamydomonas are `RESTORED`; **physcomitrella + vitis_vinifera are `PURGED`** — restore requested, rerun after ≤24 h.
- Tests: 8/8 in [tests/test_jgi_phytozome_lookup.py](tests/test_jgi_phytozome_lookup.py) (2 new for `prefer_name` selection). Token/PURGED workflow documented in [CLAUDE.md](CLAUDE.md) → "Phytozome (JGI) access".

Run the two RESTORED species (Amborella output is `protected()` + already present, so clear it first to regenerate):
```bash
conda activate rnachallenge_genes
chmod u+w resources/phytozome/amborella_trichopoda.gff3.gz 2>/dev/null; rm -f resources/phytozome/amborella_trichopoda.gff3.gz
snakemake --profile profiles/default --rerun-triggers=mtime \
  resources/phytozome/amborella_trichopoda.gff3.gz \
  resources/phytozome/chlamydomonas_reinhardtii.gff3.gz
```

---

## Evidence (post-fix run, 2026-07-12 ~12:00 and reruns after bottleneck fixes)

| Signal | Before fixes | After fixes 1–3 | After cleanup #5 | Status |
|---|---|---|---|---|
| Classified IDs | 22,249 | 22,249 | 22,249 | — |
| Unresolved / pattern-unmatched | ~8,844 | ~8,844 | ~8,844 | — |
| **Resolved to coordinates** | 15,549 | 15,549 | 15,544 | ✅ stable |
| **Sequences actually extracted** | 795 | 10,402 | 10,177 → 12,163 → **13,199** (after #10) | ✅✅ FIXED |
| Extraction failures | sequence_error 6,860+ · assembly_not_cached 6,803 | assembly_not_cached 2,187 · chrom_not_found 1,874 · missing_coordinates 1,080 · sequence_error 1 | **assembly_not_cached 314 · chrom_not_found 82 · sequence_error 1** (post-#10 rerun) | ✅ data bugs gone |

**Post-#10 run (2026-07-12, after manifest + URL-slug fan-out):** extracted **13,199**. `assembly_not_cached` **1,322 → 314**. Plant + plant_gtf rows are now downloaded through URL-backed cache keys; residual is dominated by deferred phytozome:

| db_source | count | note |
|---|---|---|
| phytozome | 203 | No genome FASTA URL source in-config; deferred by decision |
| ensembl | 53 | non-plant Ensembl stragglers |
| flybase | 33 | metazoa GTF genomes not cached |
| noncode(+v4) | 25 | NONCODE genomes not cached |

**Resolved-table URL/key coverage now matches the expected findings:**
- `plant_gtf`: 485 rows, `fasta_url` 485/485, `assembly_accession` missing 0 (467 URL-slug cache keys assigned)
- `plant`: 541 rows, `fasta_url` 541/541, `assembly_accession` missing 0 (organism fallback join in merge)
- `phytozome`: 203 rows, `fasta_url` 163/203 from table fill, but only local annotation-backed rows are currently cacheable; 203 still fail extraction as `assembly_not_cached`

Download manifest now has 112 cache keys: 107 NCBI accessions + 5 URL-slug keys (non-NCBI).

**Discovery: 2,941 rows had backwards coordinates** (start > end) in resolved TSV. This alone caused ~1,657 sequence extraction failures. Defensive fix in extract_sequences.py swaps them before samtools faidx.

**Current bottlenecks (all legitimate):**
- **assembly_not_cached (314)** — mostly deferred phytozome plus smaller ensembl/flybase/noncode URL gaps
- **chrom_not_found (82)** — residual naming/assembly-edge gaps
- **sequence_error (1)** — DEBUGGED 2026-07-12: out-of-bounds coordinate, not a bug (see [Edge-case debug](#edge-case-debug-2026-07-12))

**Real failure total: 397** (314 + 82 + 1).

---

## Edge-case debug (2026-07-12)

The two smallest failure buckets (`fail_reason 3`, `sequence_error 1`) were investigated. **Neither is a pipeline bug.**

### "fail_reason (3)" — counting artifact, not failures
The verification command counted the last column across all four `batch_*.failed.tsv` files with `awk 'NR>1'`. `NR` is cumulative, so it only skips the header of the *first* file; the other three files' header rows (`transcript_id … fail_reason`) leak in as data and get tallied as three phantom "fail_reason" failures.
- **Fix:** use `FNR>1` (per-file line number) instead of `NR>1`. Zero code change to the pipeline.

### "sequence_error (1)" — out-of-bounds coordinate (assembly-version mismatch)
- **Record:** `NONDRET005521.2` · `noncode_v4` · `GCF_000002035.5` (zebrafish GRCz10) · `chr2:60,205,149-60,216,252` (− strand).
- **What happens:** `chr2` translates correctly to `NC_007113.6`, whose length is **59,543,403 bp**. The requested start (60,205,149) is **past the end of the chromosome**. `samtools faidx` returns exit 0 but prints `[faidx] Zero length sequence` and emits a header with no bases. `faidx_extract_seq()` ([extract_sequences.py:70-72](workflow/scripts/extract_sequences.py#L70-L72)) sees empty output → returns `None` → row is tagged `sequence_error` at [extract_sequences.py:176](workflow/scripts/extract_sequences.py#L176).
- **Root cause:** NONCODE v4 zebrafish coordinates sit on an **older assembly** (Zv9/GRCz9-era) than the GRCz10 genome the pipeline downloaded. Chromosome names match after translation, but the coordinates don't fit the newer assembly's sequence lengths. This is a data-provenance mismatch, not a wiring bug.
- **Scope:** exactly 1 record. Not worth remapping NONCODE v4 to its source assembly.

**Optional diagnostic improvement (low value, ~5 lines):** in `faidx_extract_seq`, distinguish the empty-output / "Zero length sequence" case from a genuine samtools error and tag it `coord_out_of_bounds` instead of `sequence_error`. Purely for clearer failure reports; does not recover the sequence.

---

## Legitimate bottleneck analysis (2026-07-12)

The two real failure buckets were traced to root causes. **Both are fixable; neither is "just non-NCBI genomes we can't get."** Evidence gathered from `unresolved_assemblies.tsv`, per-accession `.download_done` sentinels, `logs/download_assembly/*.log`, cached `assembly_report.txt` files, and the "Extraction Failures by Resolver and Reason" table in `report.html`.

### chrom_not_found (1,633) — mostly a one-line translation bug

| Source | Count | chrom values | Root cause |
|---|---|---|---|
| ncbi | 1,025 | `1`, `2`, `3`… `X` | `chrom_translation.py` `_ALIAS_COLS = (0,4,9,6)` = Sequence-Name/GenBank/UCSC/RefSeq. **Omits column 2 (Assigned-Molecule)** — the exact column holding bare `1`/`2`/`MT`. Resolvers emit `1`; report row is `Chr1  assembled-molecule  1  …  NC_053035.3`, so `1` never enters the alias map. |
| noncode / noncode_v4 | 592 | `chrV`, `chr5`, `chrIV`, `chrII` | UCSC-ish names. `chr5`→`5` recovered once col 2 is added (via existing chr-prefix toggle). Roman (`chrIV` vs arabic `4`) still needs roman/case handling. |
| sgd | 16 | `chrmt` | Mito alias gap — yeast genome seqid is `NC_001224`; assembly report molecule is `MT`. Needs a `chrmt`/`mt`→`MT` alias + case-insensitive match. |

**Fix (do first, near-zero risk):** add `2` to `_ALIAS_COLS` in [chrom_translation.py:17](workflow/scripts/chrom_translation.py#L17). Recovers ~1,025 ncbi immediately, plus the arabic noncode chroms. Make the alias map case-insensitive to also catch `chrmt` and mixed-case UCSC names. Extend `tests/test_chrom_translation.py` with a bare-integer and a `chrmt` case.

**✅ DONE (2026-07-12)** — Assigned-Molecule column (col 2) now feeds the alias map in [chrom_translation.py](workflow/scripts/chrom_translation.py), **but gated to `assembled-molecule` rows only.** Reason: unlocalized/unplaced scaffolds repeat the parent chromosome's Assigned-Molecule (e.g. `X` with its own `NW_...` RefSeq); since the map is built top-to-bottom and later keys overwrite earlier ones, an unguarded col 2 would clobber `X`→`NC_000023.11` with a junk scaffold accession. The `Sequence-Role == "assembled-molecule"` guard keeps col 2 trustworthy. Confirmed against real reports (`GCF_000001215.4` has `X`/`NW_007931105.1` scaffold rows).

Lookups are now case-insensitive: the map stores lowercased keys and the caller lowercases its query. This picks up mixed-case UCSC names and mito `chrMT`/`chrmt` variants for NCBI assemblies. Tests: `tests/test_chrom_translation.py` + new fixture `tests/data/assigned_molecule_report.txt` covering the `Chr1`/bare-`1` split, the scaffold no-clobber guard, and `chrmt`. **The sgd `chrmt` (16) is NOT recovered here** — yeast uses Ensembl-style FASTAs with no `assembly_report.txt`, so translation returns `{}` and never runs; that subset needs a separate seqid-naming fix, not this map.

**After first rerun: 1,633 → 674.** The remaining 674 split into three *different* causes (evidence: `results/extraction_failed.tsv` grouped by assembly + cached `assembly_report.txt` + `.fai`):

| Sub-cause | Count | Assembly / db | Root cause |
|---|---|---|---|
| `chr`-prefix hid the bare name | ~564 | Arabidopsis `GCF_000001735.4` (`chr1`–`chr5`), C. elegans `GCF_000002985.6` (`chrI`–`chrV`, roman) | The report has the **bare** name (`5`, `V`) but the resolver emits `chr5`/`chrV`. The old `chr`-toggle stripped `chr`→`5`/`V` and checked the **raw `.fai`** (whose seqids are `NC_...`) — it never fed the stripped name back through the translation map. |
| assembly report never cached | 76 | `GCF_036512215.1` (ncbi, bare `1`–`X`) | Genome FASTA was cached on a **prior run (Apr 28)**; the `.download_done` sentinel then short-circuits the rule, and report fetch is **coupled to FASTA-URL resolution** which hit a 429. → `xlate={}`, bare `2` never maps. **Download-side coupling bug, not translation.** See [download_assembly.py](workflow/scripts/download_assembly.py) — decouple report fetch from URL resolution and don't gate it on the `.download_done` sentinel. |
| genuine naming gap | ~34 | sgd `chrmt` (16), Arabidopsis `chrM`/`chrC` (2), noncode_v4 `Contig*`/`Zv9_NA*` scaffolds | sgd: no report (Ensembl FASTA). `chrM`↔`MT`, `chrC`↔`Pltd` are organelle-name mismatches the report can't bridge. Legacy noncode scaffold names absent from the modern assembly. Low ROI. |

**✅ Translation fix (2026-07-12, 2nd pass):** the `chr`-toggle now runs each candidate (as-is + prefix-toggled) **through the report map before** the raw `.fai`, in a new pure `resolve_chrom_key()` in [chrom_translation.py](workflow/scripts/chrom_translation.py) (extracted from `extract_sequences.py` so it's unit-testable — 5 added cases).

**After 2nd rerun: 674 → 90.** Recovered the ~564 Arabidopsis + C. elegans rows *and* the sgd `chrmt` bucket cleared itself. Confirmed residual (`results/extraction_failed.tsv`):

| Residual | Count | Cause |
|---|---|---|
| `GCF_036512215.1` (ncbi, bare `1`–`X`) | 76 | Report-cache coupling — **fixed below.** |
| Small assemblies + organelle/legacy names | ~14 | `GCF_041296265.1` (5), `GCF_000002275.2` (5), Arabidopsis `chrM`/`chrC` (2), noncode_v4 `Contig*`/`Zv9_NA*` scaffolds. Genuine data gaps, low ROI — left as-is. |

**✅ Report-cache coupling fix (2026-07-12, 3rd pass):** report fetch was resolving its URL through `ncbi_fasta_url()` → the **datasets API**, which 429'd on rerun; already-cached genomes (FASTA from a prior run) therefore never got a report → bare `1`–`X` never mapped → `chrom_not_found`. [download_assembly.py](workflow/scripts/download_assembly.py) now resolves the report straight from the **FTP directory listing** via new `ftp_assembly_folder()` (no API, no rate limit); the fiddly FTP path-math is extracted to pure `ncbi_ftp_species_dir()` in [ncbi_assembly_utils.py](workflow/scripts/ncbi_assembly_utils.py) (de-dupes two inline copies, unit-tested in `tests/test_ncbi_ftp_species_dir.py`). `ncbi_fasta_url()` now uses the same helper for its FTP fallback. Live-verified: `GCF_036512215.1`'s report URL returns HTTP 200 with no API call. Expected to recover the 76 on rerun. **Leaves ~14 genuine data gaps.**

### assembly_not_cached (1,857) — download failures, NOT ungettable genomes

The audit's earlier framing ("non-NCBI genomes not cached") is **only half right**. Breakdown by accession/fail_detail:

| Category | Transcripts | Recoverable? | How |
|---|---|---|---|
| **NCBI GCF, HTTP 429 rate-limit** | ~494 | ✅ Yes | `ncbi_fasta_url()` hit `429 Too Many Requests`; API works now, but `.download_done` cached "failed" so Snakemake won't retry. **The correct FTP URL is already in the `fasta_url` column.** |
| **Non-NCBI (plant/phytozome/ensembl/flybase)** | ~1,315 | ✅ Yes | `download_assembly.py` only builds NCBI FTP URLs. These carry (or can carry) direct `fasta_url` to ensemblgenomes/phytozome FTP. |
| **Genuinely dead (404 / suppressed old assemblies)** | ~48 | ⚠️ Partial | `GCF_034140825.1`, `GCF_000001215.3` → 404; `GCF_000001545.5` (ponAbe2) → no FTP folder. Old assemblies pulled from NCBI FTP. Recoverable only by bumping to a current assembly version. |

**Fix (the strategic one): make `download_assembly.py` consume `fasta_url` from the resolved table instead of re-deriving it via the NCBI datasets API.** This single change:
1. Eliminates the 429 problem — no API call, the URL is already resolved (~494 NCBI transcripts).
2. Enables non-NCBI genome downloads — plant/ensembl/phytozome/flybase FTP URLs download directly (~1,315 transcripts).
3. Makes the "dead plumbing" from section A **earn its keep** as the actual download mechanism.

> ⚠️ **This reverses actionable #6.** The `fasta_url` threading should be **finished, not reverted.** Evidence: `unresolved_assemblies.tsv` already shows correct `ftp://…_genomic.fna.gz` URLs in the `fasta_url` column for the failed NCBI rows, and `ensemblgenomes.ebi.ac.uk` URLs for plant rows. The plumbing has a consumer now — it's the fix.

**Open work to actually land it:**
- The `download_assembly` fan-out wildcard is `accession = GC[FA]_\d+\.\d+`. Non-NCBI rows have `accession` = `nan`/`Phytozome`/`ARS-UCD2.0`/`BDGP6`, which don't match. Re-key the fan-out on a URL-derived slug (or `assembly_name`) so non-NCBI genomes get a download job.
- Rows where `fasta_url` is empty (plant 541, phytozome 203, plain `nan`) need the resolver to populate `fasta_url` from the per-source config (`plant_gtf_sources.yaml`, `phytozome_gtf_sources.yaml`). Without a URL there's nothing to download.
- Keep a bounded NCBI-API fallback (with retry/backoff + `ncbi_api_key` + `ncbi_connections` throttle) only for rows that carry an accession but no `fasta_url`.

### Quick win available right now (no code): retry the 429 casualties
The ~494 NCBI-429 failures can be recovered immediately by clearing their failed sentinels and re-running, now that the rate-limit window has passed:
```bash
# Remove only the failed sentinels, then re-run download (API works now)
grep -rl '^failed' resources/cache/*/.download_done | xargs rm -f
snakemake --cores 4 --resources ncbi_connections=1 results/.assemblies_ready
```
`ponytail:` this is the band-aid; consuming `fasta_url` is the durable fix that also handles non-NCBI.

**✅ Retry executed 2026-07-12 (post FTP-URL refactor):** cleared the 14 `failed` sentinels in the current 107-accession run and re-ran downloads with `--rerun-triggers=mtime` (so the `download_assembly.py` code-change didn't force a redownload of the 93 good genomes). **12 of 14 recovered.** Procedure that works:
```bash
# 1. list + clear ONLY failed sentinels (and their logs)
for d in resources/cache/*/.download_done; do
  [ "$(head -c4 "$d")" = fail ] && a=$(basename "$(dirname "$d")") && rm -f "$d" "logs/download_assembly/$a.log" && echo "$a"
done > /tmp/failed_acc.txt
# 2. re-download just those (mtime avoids redownloading the good ones after a code edit)
snakemake --profile profiles/default --rerun-triggers=mtime \
  $(sed 's#^#resources/cache/#; s#$#/.download_done#' /tmp/failed_acc.txt)
# 3. re-aggregate + re-extract + report
snakemake --profile profiles/default --rerun-triggers=mtime
```
The 2 that stayed dead are genuine data gaps (see **L2**): `GCF_000001215.3` (HTTP 404, superseded by `.4`) and `GCF_000001545.5` (whole FTP folder removed).

---

## Things That Break (ranked by impact)

### 1. No chromosome-name reconciliation — costs ~6,860 extractions ✅ DONE
_Fixed 2026-07-12 — see [reports/2026-07-12_extract_blockers.md](reports/2026-07-12_extract_blockers.md)._
`extract_sequences.py` only tries toggling the `chr` prefix (lines 126–139). NCBI genome FASTAs name sequences by **RefSeq accession** (`NC_000001.11`); resolvers emit **friendly names** (`1`, `2L`, `chrIV`, `Chr10`, `II`). These never match.

**Fix (laziest robust option):** NCBI ships a tiny `*_assembly_report.txt` next to every genome with columns `Sequence-Name / GenBank-Accn / RefSeq-Accn / UCSC-style-name`. At download time fetch it (a few KB) into `resources/cache/<acc>/assembly_report.txt`; at extract time build `{seq-name, ucsc-name, genbank} → refseq-accn` and translate `chrom` before `faidx`. ~30 lines, one helper. This single change plausibly takes extraction from 795 → ~12–13k.
- `ponytail:` don't reinvent — the report file is authoritative; no API guessing.

### 2. Extract reads the wrong file — costs most of the 6,803 "assembly_not_cached" ✅ DONE
_Fixed 2026-07-12 — see [reports/2026-07-12_extract_blockers.md](reports/2026-07-12_extract_blockers.md)._
- `download` path reads `ncbi_chromosome_resolved.tsv` (NC_/NT_/NW_ already remapped → GCF_).
- `extract_sequences.smk:30` reads **`resolved_ids.tsv`** (pre-remap; still has sequence-level `NT_…` accessions and `nan`).

So genomes are cached under `GCF_…/` but extract looks up `CACHE/NT_033779.5/` → `assembly_not_cached`. Commit `b7d0149` fixed this for `download` but not for `extract`.

**Fix:** point `extract_sequences.smk` input at `ncbi_chromosome_resolved.tsv` (one line). Verify no other rule still consumes `resolved_ids.tsv` as the "final" table.

### 2.5. Resolver outputs contain backwards coordinates — costs ~1,657 extractions ✅ DONE
_Discovered during post-fix validation; fixed 2026-07-12._

**Problem:** 2,941 rows in `ncbi_chromosome_resolved.tsv` have `start > end`. Example: `chr 8210574-8204730`. This breaks `samtools faidx`, which expects `start < end`.

**Root cause:** One or more resolvers (likely early NCBI or Ensembl stages) are emitting reversed coordinates without normalization.

**Fix (defensive):** Added coordinate swap in `extract_sequences.py` before `samtools faidx` call. Simple bounds check: if `start > end`, swap them.
```python
if start > end:
    start, end = end, start
```
Impact: sequence_error failures collapsed **1,658 → 1**, extraction improved **8,745 → 10,402** (+1,657 sequences).

**Note:** This is a band-aid. The root cause (which resolver emits backwards coords?) should be fixed in Stage 2 for correctness, but the fix is warranted defensively here.

### 3. Sources that resolve without coordinates — resolved 🟢
`gramene` was updated to request `fl=*` and now emits coordinates. `noncode_2016` was reclassified to `matched_not_found` because it is existence-only by design.

### 4. `missing_coordinates` — resolved 🟢
The merge guard now prevents rows with empty `chrom` or missing `start`/`end` from counting as resolved. There are no remaining `missing_coordinates` failures in the latest extraction summary.

---

## Overcomplicated Patterns (delete / simplify)

### A. URL-propagation refactor status — now load-bearing ✅
_Updated 2026-07-12:_ this is no longer dead plumbing.
- `merge_resolved.py` now fills plant URLs via `assembly_name` and fallback `organism` join.
- URL-backed rows with missing accession now get deterministic URL-slug cache keys.
- `prepare_accession_list` now emits a manifest (`cache_key`, `fasta_url`) instead of bare accession list.
- `download_assembly.py` now prefers manifest `fasta_url` and falls back to NCBI lookup only when needed.
- Fan-out wildcard now accepts URL-slug cache keys.

Net effect: extraction improved **12,163 → 13,199** and `assembly_not_cached` dropped **1,322 → 314**.

### B. Two download scripts, one deprecated 🟠
`download_assemblies.py` (338 lines) is marked **DEPRECATED**, superseded by `download_assembly.py` + `aggregate_downloads.py`, and kept only "for test compatibility" (`test_download_assemblies_phase4.py`). The helper funcs (`ncbi_fasta_url`, `is_ncbi_assembly_accession`, …) are **duplicated** in `download_assembly.py`.

**Fix:** point the test at `download_assembly.py`, delete `download_assemblies.py`. One source of truth for the FTP URL logic.

### C. NONCODE has 4 rules for one database 🟠
`resolve_noncode` + `resolve_noncode_v4` + `resolve_noncode_2016` + `resolve_noncode_assembly_accessions`, and `noncode_2016` was existence-only (now reclassified away from resolved). The v4/2016 fallbacks add noise; collapse: keep v5 + v4 (both yield coords), drop 2016 to the failure report unless it demonstrably adds extractable hits.

### D. Scale check (not urgent, but note it) 🟡
~14,300 LOC across 30+ scripts / 31 rule files for "parse ID → look up coords → slice FASTA." Much of it is irreducible (heterogeneous DBs), but the half-built refactor (A), the duplicate downloader (B), and redundant NONCODE tiers (C) are removable now. `resolve_abandoned_accessions.py` (1,142 lines, ~8h runtime) earns its 2,002 hits but is the obvious next target if runtime becomes a problem — it downloads full NCBI GTFs.

---

## Legitimate Data Issues (inherent, not pipeline bugs)

These are **not** fixable by better code — they stem from the source data itself. Distinguish them from the fixable bottlenecks above so we don't chase unrecoverable rows. Numbers are current-run failure counts.

### L1. Assembly-version / coordinate provenance mismatch (~1+ rows) 🟡
NONCODE v4 (and likely other legacy DBs) emit coordinates against an **older assembly** than the one we download. Chromosome names translate fine, but the coordinate can fall **past the end of the newer assembly's sequence** → `samtools faidx` returns an empty region. Confirmed on `NONDRET005521.2` (zebrafish, chr2:60.2Mb vs chr2 length 59.5Mb). Only 1 surfaced as `sequence_error` today, but the class is real. **Genuine fix would require per-DB assembly-version pinning + coordinate liftover — not worth it for the volume.** Optional: tag as `coord_out_of_bounds` for clearer reporting.

### L2. Suppressed / withdrawn NCBI assemblies 🟡
_Updated 2026-07-12 after the FTP-URL refactor + retry:_ the FTP-directory-listing resolver (`ftp_assembly_folder()`) **recovered most of the presumed-dead list** — `GCF_000001895.5` (`Rnor_6.0`), `GCF_034140825.1`, `GCF_000002775.4`, `GCF_000695525.1`, `GCF_020379485.1`, `GCF_028885655.2`, `GCF_000511025.2`, `GCF_001660625.3`, `GCF_964237555.1`, `GCF_000002295.2`, `GCF_000001635.26`, `GCA_000188115.2` all now download `ok`. They were never truly withdrawn — the old datasets-API path just couldn't build their URLs.

**Only 2 remain genuinely ungettable at the pinned version:**
- `GCF_000001215.3` (*D. melanogaster*, `Release_6_plus_MT`) — **HTTP 404**, superseded by `GCF_000001215.4`.
- `GCF_000001545.5` (*P. abelii*, `ponAbe2`) — **entire FTP folder removed**.

Recoverable **only** by bumping to a current assembly, which reintroduces L1-style coordinate risk (the resolved coords were computed against the old assembly report). Low ROI — accept as lost or bump case-by-case.

### L3. Cross-database chromosome-naming conventions (~roman-numeral subset of 592) 🟡
Worm/yeast/NONCODE use roman numerals (`chrIV`, `chrII`); NCBI assembly reports use arabic (`4`, `2`) in Assigned-Molecule. Arabic cases (`chr5`→`5`) are fixable (actionable #9); **roman↔arabic needs an explicit mapping table** because there is no authoritative column linking them. Partly fixable, partly a data-convention gap.

### L4. Upstream unresolved / pattern-unmatched IDs (~8,844) ⚪
Of 22,249 classified IDs, ~8,844 never resolve to coordinates at all (obsolete IDs, unsupported databases, malformed input). These never reach extraction. Out of scope for extraction fixes — a separate resolution-coverage question, tracked at the parse/resolve stages, not here.

---

## Resolution-coverage roadmap (matched_not_found — 5,993 rows)

_Added 2026-07-12. Source: [reports/reason_resolution_brainstorm/MASTER_ROI_RANKING.md](reports/reason_resolution_brainstorm/MASTER_ROI_RANKING.md) + per-reason files. This extends **L4** — it's the resolve-stage counterpart to the extract-stage audit above. Distinct from extraction failures (397); these are IDs that classify but never reach coordinates._

Current `matched_not_found` breakdown (from `results/matched_not_found.tsv`):

| Reason | Rows | ROI strategy | Effort / expected recovery |
|---|---|---|---|
| `not_found_in_gramene` | 1,935 | Legacy rice/maize crosswalk + versioned plant fallback (#2) | M-H / 40-75% |
| `missing_coordinates` | 1,483 | Organism-aware GTF fallback after NCBI coord failure (#3) | M / 20-40% |
| `phytozome_gff_no_match_oryza_sativa` | 1,025 | RAP/MSU legacy rice ID crosswalk (#2) — **not** the JGI-auth fix; these are ID-version mismatches vs the GFF3 we have | M-H / high |
| `matched_noncode2016_no_coordinates` | 674 | NONCODE normalized/fuzzy match + species fallback (#4) | M / 70-90% for the good-coverage subset |
| `phytozome_gff_no_match_zea_mays` | 644 | Legacy maize (Zm00001d vs B73 v5) crosswalk (#2) | M-H / high |
| `not_found_in_any_noncode` | 210 | NONCODE base-ID normalization (#4) | M / 40-90% |
| `assembly_mapping_failed:NC_008405.2` / `NC_008394.4` | 15 | ✅ **DONE** — static NC→GCF exception map (#5): rice Build 4.0 = `GCF_000005425.2` | **L / ~100%** |
| `worm_gtf_not_resolved` | 5 | Multi-release worm fallback (#7) | L-M / 60-80% |
| `sgd_gtf_not_resolved` | 1 | ✅ **DONE** — `Source:SGD;Acc:*` canonicalization + `dbxref` index key (#6) | L / ~100% |
| `phytozome_gff_no_match_citrus_sinensis` | 1 | ID-version mismatch (not a parser gap) — needs crosswalk, out of #6 scope | — |

**Highest ROI = Strategy #1 (Source Attribution Backbone), cross-cutting.** Persist transcript provenance (source tool/paper/DB/release) from Stage 0 through merge. Immediate payoff with *current* data — no new resolver logic: 38.5% of `matched_not_found` rows already carry a tool tag (`not_found_in_gramene` 767, `missing_coordinates` 697, phytozome zea 563, phytozome oryza 258). It routes every large bucket below with one architectural change.

**Recommended execution order:** #1 backbone → #2 legacy rice/maize + versioned plant fallback (biggest bucket volume) → #3 organism-aware GTF fallback → #4 NONCODE normalization → #5/#6 low-effort near-100% tails → #7 worm → #8 optional alignment deep-recovery (gated by config).

**Note on phytozome buckets:** `phytozome_gff_no_match_*` (1,669) is an **ID-version-mismatch** problem (legacy RAP/MSU/Zm IDs not present in the GFF3 build), *not* the JGI-auth gap just fixed. The auth fix adds *new species*; these buckets need *crosswalk tables* for species we already have.

**Stage 0 coverage gaps to close (prereq for #1):** CNIT + FEELnc datasets absent locally (`n_sequences_loaded = 0`); PreLnc now wired as a Stage 0 source (train/test FASTA for human/mouse/cow/Arabidopsis/rice/maize — broadens plant provenance).

**Highest-value external input:** a source mapping table from the paper (transcript ID → originating DB/tool/paper/release). Unlocks #1 and #2 routing precision immediately.

---

## Per-problem deep dives (subagent outputs)

Each concern was investigated in isolation; full findings in `audit/`:

| # | File | Verdict |
|---|---|---|
| 1 | [audit/01_chrom_naming.md](audit/01_chrom_naming.md) | Confirmed. Fetch `*_assembly_report.txt` at download; ~20-line translate helper in extract, **gated to GCF_/GCA_ only** (worm/yeast/fly/plant seqids already match). Report URL derivable from existing `ncbi_fasta_url()` folder logic. |
| 2 | [audit/02_extract_wiring.md](audit/02_extract_wiring.md) | Confirmed. **One-line fix**: `extract_sequences.smk:30` → `ncbi_chromosome_resolved.tsv`. Schemas identical (20 cols); the 5-row delta is unmappable `NT_479536.1` correctly diverted. Only other consumers are read-only report stats. |
| 3 | [audit/03_coordinateless_resolvers.md](audit/03_coordinateless_resolvers.md) | gramene: coords **are** cheaply fetchable — add `fl=*` to the existing API call (verified live, 25/25 return coords, zero extra requests) → **fetch, don't drop**. noncode_2016: existence-only by design, 674 unique but 0 extractable → **reclassify to matched_not_found**. Add a generic guard so "resolved" means "has coordinates". |
| 4 | [audit/04_url_refactor.md](audit/04_url_refactor.md) | **REVERT** the `fasta_url`/`gtf_url`/`assembly_name` plumbing + `fill_urls_from_table` — dead (0 consumers, `assembly_name` 0/15549). **Keep** the config `assembly_accession` enrichment. Revert **selectively**: commit `6683a76` also stripped the 8h NCBI GTF download — don't reintroduce it. |
| 5 | [audit/05_dedup_downloaders_noncode.md](audit/05_dedup_downloaders_noncode.md) | `download_assemblies.py` helpers are byte-identical to `download_assembly.py` and no rule uses it; **delete it** + delete/repoint 2 stale integration tests (no test actually imports it). NONCODE 4 rules → 2: `resolve_noncode_assembly_accessions` is a fake checkpoint duplicating merge's inline mapping (**fold into merge**); drop `resolve_noncode_2016`. |

## Actionables (in order)

1. ~~**[blocker] Fix extract input file**~~ ✅ **DONE (2026-07-12 11:34)** — `extract_sequences.smk:30` now reads `ncbi_chromosome_resolved.tsv`.
2. ~~**[blocker] Add chromosome-name translation**~~ ✅ **DONE (2026-07-12 11:34)** — report fetched at download (`download_assembly.py::fetch_assembly_report`), translated in `extract_sequences.py` via new `chrom_translation.py`; `tests/test_chrom_translation.py` (6 cases) passing. Improved extraction 795 → 8,745.
3. ~~**[blocker] Fix backwards coordinates**~~ ✅ **DONE (2026-07-12 11:57)** — added defensive swap in `extract_sequences.py` line 117. Improved extraction 8,745 → **10,402** (1,657 sequences gained). sequence_error 1,658 → 1.
4. ~~**Run one clean end-to-end pass**~~ ✅ **DONE (2026-07-12 ~12:00)** — full pipeline completed with all three blockers fixed. Results regenerated and trusted.
5. ~~**Fix coordinate-less resolvers**~~ ✅ **DONE (2026-07-12 ~12:xx)** — gramene now returns coordinates, noncode_2016 is routed to `matched_not_found`, and merge filters coordinate-less resolved rows. missing_coordinates is now gone from the latest failure report.
6. ~~**Revert the URL refactor selectively**~~ ⛔ **SUPERSEDED — see [Legitimate bottleneck analysis](#legitimate-bottleneck-analysis-2026-07-12).** Evidence shows the `fasta_url` column is the fix for ~1,800 extraction failures. **Finish the refactor (actionable #10), don't revert it.** Only genuinely-dead pieces (`assembly_name` 0/15549, `gtf_url`/`gtf_format` with no consumer) may still go; `fasta_url` stays.

### Bottleneck fixes (new — from 2026-07-12 analysis)

9. ~~**[bottleneck] Add Assigned-Molecule to chrom translation**~~ ✅ **DONE (2026-07-12)** — col `2` added to the alias map (gated to `assembled-molecule` rows to avoid scaffold clobber), case-insensitive lookups. Plus 2nd pass: `chr`-toggle now feeds candidates through the report map (pure `resolve_chrom_key()`). **Recovered chrom_not_found 1,633 → 90** across two reruns. Tests: `tests/test_chrom_translation.py` (15 cases) + fixture `tests/data/assigned_molecule_report.txt`. See [chrom_not_found analysis](#chrom_not_found-1633--mostly-a-one-line-translation-bug).
9b. ~~**[bottleneck] Decouple assembly-report fetch from the datasets API**~~ ✅ **DONE (2026-07-12)** — report URL now resolved from the FTP directory listing (`ftp_assembly_folder()` in [download_assembly.py](workflow/scripts/download_assembly.py)), not the rate-limited API that 429'd and starved already-cached genomes of their report. Pure path-math extracted to `ncbi_ftp_species_dir()` in [ncbi_assembly_utils.py](workflow/scripts/ncbi_assembly_utils.py) (`tests/test_ncbi_ftp_species_dir.py`). Live-verified HTTP 200. **Recovers the 76 `GCF_036512215.1` rows on rerun → chrom_not_found ~90 → ~14 (genuine data gaps only).**

    **▶ NEXT STEP:** rerun only the failed sentinels, not the whole download rule. Use [workflow/scripts/list_failed_download_targets.py](workflow/scripts/list_failed_download_targets.py) to turn `results/unresolved_assemblies.tsv` into a newline-delimited list of `resources/cache/<cache_key>/.download_done` targets, then rerun those files with `--rerun-triggers=mtime`. Important: either force the rerun or delete the listed `.download_done` files first, or Snakemake may consider them already satisfied. After that, `cut -f5 results/extraction_failed.tsv | sort | uniq -c` should confirm only the genuine residual gaps remain. The leftover organelle / legacy / sgd cases are genuine data gaps, not bugs — accept them as lost.
10. ~~**[bottleneck] Make `download_assembly.py` consume `fasta_url`**~~ ✅ **DONE (2026-07-12)** — implemented manifest-driven downloads (`cache_key`,`fasta_url`), URL-slug fan-out, and downloader URL-first behavior with NCBI fallback. Measured impact: **extracted 12,163 → 13,199**, `assembly_not_cached 1,322 → 314`.
11. ~~**[bottleneck] Populate `fasta_url` for URL-less non-NCBI rows**~~ 🟡 **PARTIALLY DONE (2026-07-12)** — plant stream now fills by `organism` fallback in merge (541/541). `plant_gtf` and most phytozome rows now carry URL metadata in merged output; however phytozome genome FASTA remains unavailable as a reliable source in current config workflow. **Decision:** defer phytozome (203 rows).
12. **[bottleneck, optional] Handle roman-numeral + dead assemblies** — roman↔arabic chrom mapping for worm/noncode `chrIV`; the "~48 dead assemblies" collapsed to **2** after the FTP-URL refactor + retry (2026-07-12): `GCF_000001215.3` (404, use `.4`) and `GCF_000001545.5` (`ponAbe2`, folder removed). Bump to current versions or accept as lost. Lowest ROI. See **L2**.

### Cleanup (unchanged)

13. **Delete `download_assemblies.py`**; delete/repoint the 2 stale integration-test assertions (no test imports it — `test_phase4` has its own copies).
14. **Collapse NONCODE 4 rules → 2** — fold `resolve_noncode_assembly_accessions` into merge, drop `resolve_noncode_2016`.
15. ~~**[reporting] Split conflated unresolved metric**~~ ✅ **DONE (2026-07-12)** — reports now distinguish:
  - `Unclassified` = `pattern_unmatched.tsv`
  - `Classified but unresolved` = `matched_not_found.tsv`
  - Legacy combined unresolved remains visible as a compatibility view.
    - Added explicit report equations/consistency checks:
      - `Input = Classified + Unclassified`
      - `Unresolved = Unclassified + Classified but unresolved`

  Implemented in [workflow/rules/report.smk](workflow/rules/report.smk) (new inputs), [workflow/scripts/generate_resolution_report.py](workflow/scripts/generate_resolution_report.py), [workflow/scripts/generate_report.py](workflow/scripts/generate_report.py), and [workflow/scripts/report_utils.py](workflow/scripts/report_utils.py). Cards, funnel rows, and next-action guidance now use the split counts.
  Reports regenerated: `results/report.html` and `results/resolution_report.html`.
  Current split snapshot: `pattern_unmatched=5008`, `matched_not_found=5993`, combined `unresolved=11001`.

**Impact of fixes 1–5:** extraction ~795 → **10,177**. Blockers eliminated.
**Observed impact of fixes 9–11:** extraction now **13,199**, with remaining failures concentrated in deferred phytozome + small non-plant URL gaps.

---

## Verification checklist (after fixes 1–3) ✅ VERIFIED

```bash
# Run completed 2026-07-12 ~12:00 and rerun after cleanup
snakemake --cores 8 extract_sequences

# Results (latest rerun):
grep -c '^>' results/output.fasta                     # 13,199 ✅
# NOTE: use FNR>1 (per-file), not NR>1 — otherwise the 3 extra file headers get counted as "fail_reason"
awk -F'\t' 'FNR>1{print $NF}' results/sequences/*.failed.tsv | sort | uniq -c
# 314 assembly_not_cached    → phytozome 203 + ensembl 53 + flybase 33 + noncode 25
#  82 chrom_not_found        → residual naming/data gaps
#    1 sequence_error        ✅ debugged — out-of-bounds coord (assembly-version mismatch), not a bug
# (the old "3 fail_reason" were miscounted header rows — gone with FNR>1)
```

**Verdict:** All blockers and bottleneck #10 are confirmed fixed. Remaining **397** failures are:
- Legitimate data issues (non-cached genomes, unmapped chromosomes)
- The 1 `sequence_error` is an out-of-bounds coordinate (assembly-version mismatch), debugged — not a bug
- The old "3 fail_reason" were miscounted header rows, not failures
- Not data bugs or pipeline wiring errors
- Mostly residual data gaps; the only large deferred block is phytozome.

---

## Next steps (bottlenecks first, then cleanup)

**Blockers resolved, #9 and #10 landed. Extraction at 13,199 (~85%). Remaining work is targeted cleanup + deferred phytozome policy.**

### Priority 0: Land the phytozome/JGI unblock (NEW 2026-07-12)
- JGI auth is now implemented + 5 new species wired. **Re-run parse → phytozome resolution → merge** to grow the resolved set with the new species (vitis/potato/amborella/chlamydomonas/physcomitrella). See [Phytozome/JGI unblock](#phytozomejgi-unblock-2026-07-12).
- To also *extract* phytozome coordinates, add genome-FASTA entries to the manifest (same token) — otherwise the new rows resolve but stay `assembly_not_cached`.

**⚠️ Unblock to run *now* while JGI restores are pending (2026-07-12):** `resolve_phytozome_gtf` needs **all** configured GFF3s present and sits **upstream of merge→extract→report**, so any single missing phytozome file blocks the whole pipeline (not just the phytozome stream). Two problems found + fixed so a full run can proceed today:
1. **Config path mismatch (fixed).** `citrus/sorghum/ricinus/oryza` pointed at `*.gene_exons.gff3.gz` filenames not on disk (disk has `*.gene.gff3.gz`). Because the `download_phytozome_gtf` script keys by *species name* while the resolve rule requests files by their `gtf:` *basename*, the mismatched paths triggered JGI downloads that fail with "manifest entry not found". Repointed the 4 `gtf:` paths to the existing `.gene.gff3.gz` files (verified they carry `mRNA` features the resolver reads). `touch`ed them so they postdate `manifest.json` (else mtime forces the same failing rebuild).
2. **JGI-only species deferred.** `solanum_tuberosum` (no local file), `vitis_vinifera` + `physcomitrella_patens` (**PURGED**, restore reqs pending, ≤24h) are commented out in `config/phytozome_gtf_sources.yaml` with dated re-enable notes. amborella + chlamydomonas are RESTORED and already downloaded, so they stay in.

Result: dry-run drops from 37 → 29 jobs with **0 phytozome download jobs** and no missing inputs. Run `snakemake --profile profiles/default --rerun-triggers=mtime`. When JGI restores vitis/physcomitrella (poll `files.jgi.doe.gov/request_archived_files/requests/652368`), uncomment them (and add solanum_tuberosum's GFF3), then rerun `resolve_phytozome_gtf` + downstream.

### Priority 1: Residual non-NCBI URL gaps (post-#10)
- `assembly_not_cached` residual is 314: phytozome 203 (extraction still needs genome FASTA — see Priority 0), ensembl 53, flybase 33, noncode(+v4) 25.
- For ensembl/flybase/noncode residuals, continue same manifest mechanism by filling reliable `fasta_url` in their resolver/config paths.
- When retrying downloads, target only the relevant `resources/cache/<cache_key>/.download_done` files derived from failed rows; do not force the already-ok assemblies.
- Prefer [workflow/scripts/list_failed_download_targets.py](workflow/scripts/list_failed_download_targets.py) over one-off shell parsing when generating the target list.

### ~~Chrom-translation fix (actionable #9)~~ ✅ DONE — chrom_not_found 1,633 → ~110.

### ~~NCBI-429 / dead-assembly retry~~ ✅ DONE (2026-07-12) — assembly_not_cached 1,857 → 1,322; 12/14 sentinels recovered, 2 genuine dead ends (L2).

### ~~Manifest-based URL download fan-out (actionable #10)~~ ✅ DONE (2026-07-12)
- `download_assembly.py` now consumes manifest `fasta_url` directly and supports URL-slug cache keys.
- `prepare_accession_list` now emits `results/assembly_download_manifest.tsv` (`cache_key`, `fasta_url`).
- `download_assembly` wildcard broadened for non-NCBI cache keys.
- Measured gain: extracted 12,163 → **13,199**.

### Priority 3: Cleanup (actionables #13–14)
- Delete deprecated `download_assemblies.py`; repoint `test_download_assemblies_phase4.py`.
- Collapse NONCODE 4 rules → 2 (fold `resolve_noncode_assembly_accessions` into merge, drop `resolve_noncode_2016`).
- **Do NOT** blanket-revert the URL refactor — `fasta_url` is now load-bearing (Priority 2).

### Priority 4: Report split snapshot
- Capture and paste the updated top-line split counts into the Evidence/Verification sections:
  - `unclassified = pattern_unmatched`
  - `classified_but_unresolved = matched_not_found`
  - `unresolved_total = unclassified + classified_but_unresolved`

**Files with blockers 1–3 fixes** (already committed):
- `workflow/rules/extract_sequences.smk` — input repoint (blocker #2)
- `workflow/scripts/download_assembly.py` — `fetch_assembly_report()` + 3 call sites (blocker #1)
- `workflow/scripts/chrom_translation.py` — **new** `load_chrom_translation()` helper (blocker #1)
- `workflow/scripts/extract_sequences.py` — translate `chrom` + defensive swap for backwards coords (blockers #1, #3)
- `tests/test_chrom_translation.py` + `tests/data/GCF_000001405.40_assembly_report.txt` — **new** (6 cases, passing)
- `workflow/scripts/gramene_resolver.py` — now requests and emits coordinates
- `workflow/scripts/resolve_noncode_2016.py` — now routes existence-only hits to unresolved
- `workflow/scripts/merge_resolved.py` — coordinate guard prevents invalid rows from counting as resolved

---

## Maintenance protocol — keep this audit the source of truth

**Before you end any investigation or fix, update this file.** The audit only stays useful if it matches reality.

1. **File findings under the right heading:** a bug/wiring problem → *Things That Break*; removable complexity → *Overcomplicated Patterns*; an inherent data limitation → *Legitimate Data Issues*.
2. **Update Actionables:** tick items you completed (`~~strikethrough~~ ✅ DONE (date)` with the measured impact), add new ones, and mark reversed conclusions **superseded** rather than deleting them.
3. **Refresh Next Steps** so the priority order reflects the current state.
4. **Re-run the counts** in the Verification checklist and paste the new numbers so the Evidence table stays honest.
5. **Leave a dated note** (`_Fixed YYYY-MM-DD_` / `_Investigated YYYY-MM-DD_`) so history is traceable.

If you dispatched subagents, fold their insights back here before finishing — the audit, not the chat log, is what the next agent reads.
