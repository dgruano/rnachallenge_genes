# Pipeline Audit — RNA Flanking Sequence Pipeline

_Date: 2026-07-14 · Branch: `main` · HEAD: `7ab78ce` (phytozome genome-FASTA extraction + assembly-collision inventory) · post-rerun refresh_

> **📌 This document is the source of truth for pipeline health.**
> Every agent (or person) who investigates or fixes anything **must finish by updating this file** — add findings to the right section (Things That Break / Overcomplicated Patterns / Legitimate Data Issues), tick or renumber the Actionables, and refresh Next Steps. Leave the audit consistent with reality before ending your turn. If a prior conclusion is wrong, mark it superseded rather than deleting the history.

## Plant resolution — namespace-keyed v44/Phytozome-v11 (2026-07-14, branch `fix/rnaplonc-ensembl-plants-v44`)

⏳ **CODE-DONE, URLs pending, NOT yet rerun.** Reworks how plant transcripts pick a
reference assembly, driven by the corrected source model. **Supersedes L6's AGPv3/v44
guesswork** with a *verified* namespace→assembly mapping. Spec:
[docs/superpowers/specs/2026-07-14-tool-aware-plant-resolution-design.md](docs/superpowers/specs/2026-07-14-tool-aware-plant-resolution-design.md).

**Corrected source model** (from `dani_notes.md` + user):
- **RNAPlonc → Phytozome v11** (its GreeNC source is Phytozome-based; v10 had no annotation
  files → v11 substitute). *No separate GreeNC resolver* — folded into the Phytozome path.
- **PreLnc coding → Ensembl Plants v44**; PreLnc non-coding → GreeNC (≈ Phytozome v11).

**Key correctness finding — resolve by namespace, not by (species, tool).** The worry that one
species needs different assemblies per tool is real, but the *namespace* already encodes it,
and it is **verified by ID presence**, not assumed:
- Maize: **0 of 100** input `GRMZM…` cores exist in `Zea_maysb73v4.AGPv4.gff3.gz`; `Zm00001d…`
  is native AGPv4. → `GRMZM` (AGPv3, RNAPlonc) and `Zm00001d` (AGPv4, PreLnc-coding) are
  **disjoint namespaces on different assemblies**. The tool split *is* the namespace split.
- Rice `Os…t` shared by PreLnc(v44) + LGC(r30): same IRGSP-1.0 assembly → coordinates agree;
  tool only picks an annotation *release*. No namespace maps to two assemblies in this input.

**Differentiating assembly versions** (the method, so this isn't re-litigated): (1) read the
header's inline tag when present — `Zm00001d…_T001 … chromosome:B73_RefGen_v4:…` = AGPv4;
(2) else ID prefix encodes the generation — `GRMZM`=RefGen_v3/AGPv3, `Zm00001d`=RefGen_v4/AGPv4,
`Zm00001e`=NAM5.0; (3) **always confirm by `zgrep -c -Ff <ids> <gff>.gz`** — near-zero ⇒ wrong
assembly. This check is now automated (below).

**Landed in code** (`90fe4d0`):
- `parse_ids.py` — splits maize `Zm00001d…` (→ v44/AGPv4) from `GRMZM…` (→ Phytozome v11/AGPv3);
  the 512 `Zm00001d` IDs were previously *unclassified*.
- `plant_gtf_sources.yaml` — PreLnc-coding (`arabidopsis`, `oryza_sativa`, maize `Zm00001d`)
  repointed to **release-44**, URLs are `FILL_ME__…` placeholders (grep `FILL_ME`).
- `phytozome_gtf_sources.yaml` — maize repointed **AGPv4→AGPv3 v11** (for `GRMZM`); **added
  `brachypodium_distachyon` + `manihot_esculenta`** (RNAPlonc `Bradi`/`Manes.`, previously
  `missing_coordinates`). `genome_id`/GFF filename are `FILL_ME__…` (fill from JGI portal).
- `resolution_guard.py` + wiring — each plant resolver computes a per-species match rate and
  **exits non-zero** if a species matches < `plant_resolution_min_match_rate` (default 0.02) of
  its IDs against a configured source (the wrong-assembly signature). Override with
  `plant_resolution_strict: false`. Unit test `tests/test_resolution_guard.py` (5 cases).

**Maize split works via the two resolvers, each keyed on species:** plant_gtf `zea_mays` = v44
(matches `Zm00001d`), phytozome `zea_mays` = AGPv3 v11 (matches `GRMZM`); each ID matches exactly
one file, the other routes to that resolver's *unresolved* output, and `merge_resolved` unions
the resolved rows. DAG validated (`snakemake -n --forcerun`): new species fan out into
download+resolve jobs; execution will fail only on the unfilled `FILL_ME` URLs/genome_ids.

**⏭ Blocked on @dgruano:** fill the `FILL_ME__` v44 URLs (`plant_gtf_sources.yaml`) and the
Phytozome v11 `genome_id`s (`phytozome_gtf_sources.yaml`, brachypodium/manihot/maize-AGPv3),
then rerun (command in Next Steps / commit). Expected recovery: 561 `GRMZM` + 512 `Zm00001d`
maize + 142 `Bradi` + 127 `Manes.` + rice/arabidopsis PreLnc-coding.

## Run refresh — 2026-07-14 (phytozome fix landed)

Latest `results/` (`grep -c '^>' output.fasta` + `FNR>1` failure tally). Two reruns on 2026-07-14: the first (post-`d546a0e`) only half-applied the phytozome fix; the **phytozome resolver regeneration (Priority 0) is now DONE** and the 503 rows extract.

| Signal | post-#10 | half-applied | **fix landed** | Δ vs prev |
|---|---|---|---|---|
| **Sequences extracted** | 13,199 | 15,319 | **15,820** | +501 ✅ |
| assembly_not_cached | 314 | 614 | **111** | −503 ✅ |
| chrom_not_found | 82 | 85 | **87** | +2 |
| sequence_error | 1 | 1 | 1 | — |
| **Real failure total** | 397 | 700 | **199** | −501 |

✅ **Phytozome cleared.** `resolved_ids.tsv` now carries per-species keys (`phytozome_vitis_vinifera` 140, `phytozome_sorghum_bicolor` 131, `phytozome_physcomitrella_patens` 122, `phytozome_oryza_sativa` 30, `phytozome_amborella_trichopoda` 34, `phytozome_citrus_sinensis` 26, `phytozome_ricinus_communis` 14, `phytozome_chlamydomonas_reinhardtii` 4, `phytozome_zea_mays` 2 = 503), matching the cached genomes. All extract except 2 `phytozome_zea_mays` rows (now a chrom-naming miss, not assembly).

**Remaining `assembly_not_cached` (111) — all previously-known, none phytozome:**

| bucket | accession | count | status |
|---|---|---|---|
| ensembl | `ARS-UCD2.0` (cattle) | 53 | Class 2 bare-name — needs `fasta_url` populated |
| flybase | `BDGP6` (drosophila) | 33 | Class 2 bare-name — needs `fasta_url` populated |
| noncode | `GCF_000001545.5` (ponAbe2) | 24 | L2 genuinely dead (FTP folder removed) |
| noncode_v4 | `GCF_000001215.3` (dmel) | 1 | L2 genuinely dead (404, use `.4`) |

So **86 of 111 are recoverable** (cattle 53 + fly 33 non-NCBI genomes); 25 are genuine dead ends.

**`chrom_not_found` (87) — dominated by plant_gtf tomato 70:** all `GCA_000188115.2` (tomato SL3.0), bare `1`–`12`. FASTA seqids are GenBank accessions `CM001064.2`–`CM001075.2` (from the `.fai`); this GCA **has** an NCBI `assembly_report.txt` that would map them, but the genome was cached via the Ensembl-Plants URL path which never fetched the report. Residual: noncode 9, ncbi 4, phytozome_zea_mays 2, plant_gtf potato 1, noncode_v4 1 — small data/naming gaps.

---

## TL;DR — UPDATED 2026-07-12 POST-FIX

**Status: Extraction working, four blockers/refinements fixed.**
- **Blocker #1 (chrom translation)** ✅ DONE — assembly reports fetched, chromosome names translated
- **Blocker #2 (extract file wiring)** ✅ DONE — extract now reads remapped `ncbi_chromosome_resolved.tsv`
- **Blocker #3 (backwards coordinates)** ✅ DONE — defensive swap added, sequence_error collapsed from 1,658 → 1
- **Cleanup #5 (coordinate-less resolvers)** ✅ DONE — gramene now returns coords, NONCODE2016 reclassified to matched_not_found, merge guards invalid resolved rows

**Current extraction:** **13,199 sequences** from 15,544 resolved IDs (~85% extraction rate) after landing bottleneck fix #10 (manifest-driven URL downloads) and rerunning. Up from 12,163.

Remaining failures were mostly legitimate data gaps. The only sizable unresolved download bucket was **phytozome (203)**, previously deferred due missing genome FASTA URLs/auth constraints — **JGI auth was implemented (2026-07-12), see [Phytozome/JGI unblock](#phytozomejgi-unblock-2026-07-12); the genome-FASTA extraction gap is now closed in code (2026-07-14), see [Phytozome genome-FASTA extraction](#phytozome-genome-fasta-extraction-2026-07-14).** Counts below predate that landing and need a post-rerun refresh.

---

## Phytozome genome-FASTA extraction (2026-07-14)

✅ **CODE-DONE + VERIFIED (2026-07-14 rerun).** The resolver now emits per-species `phytozome_<species>` keys into the resolved table, matching the cached genomes; **all 503 phytozome rows extract** (`assembly_not_cached` 614 → 111, none phytozome; extraction 15,319 → 15,820). See [Run refresh 2026-07-14](#run-refresh--2026-07-14-phytozome-fix-landed). Closes the last open half of the phytozome deferral — the 203 `assembly_not_cached` rows that resolved to coordinates but had no genome to slice against. Landed in `d546a0e` (`feat(phytozome): cache genome FASTAs so resolved rows extract`):

- **Per-species genome cache.** New `download_phytozome_fasta` rule ([download_phytozome_fasta.smk](workflow/rules/download_phytozome_fasta.smk) + [download_phytozome_fasta.py](workflow/scripts/download_phytozome_fasta.py)) fetches each species' JGI genome FASTA to `resources/cache/<assembly_accession>/genome.fasta` — exactly where `extract_sequences` looks. Fan-out is driven by the resolved phytozome species and feeds `.assemblies_ready`.
- **Collision fix.** The resolver now emits `assembly_accession = phytozome_<species>` instead of a shared literal `Phytozome`, so the 9 phytozome species no longer collapse into one cache dir. Version-matched JGI lookup keeps annotation↔assembly consistent.
- **Assembly/provenance collisions inventoried.** This work surfaced a set of accession→cache collisions across resolver streams — see [ASSEMBLY_COLLISIONS.md](ASSEMBLY_COLLISIONS.md). Class 1 (shared `Phytozome` key, 503 rows) is fixed by this commit; **Class 2 (~330 recoverable ensembl/flybase/wormbase/ncbi bare-name rows) is the highest-ROI next mapping fix.**
- **Deploy caveat:** delete any stale `resources/cache/Phytozome/` before the next run so extraction can't read an orphaned wrong genome.

**⏭ Rerun pending:** the fix is in code but not yet reflected in `results/` — the 13,199 / 397 counts and the phytozome-203 bucket throughout this doc are from a pre-`d546a0e` run. Re-run download→extract→report, then refresh the Evidence and Verification counts.

---

## Phytozome/JGI unblock (2026-07-12)

The phytozome deferral had two distinct blockers; the auth one is now closed:

1. **Annotation (GFF3) download — auth. ✅ UNBLOCKED.** `download_phytozome_gtf` ([resolve_phytozome_gtf.smk](workflow/rules/resolve_phytozome_gtf.smk)) now loads a JGI bearer token from `.env` (`JGI_SESSION_TOKEN` or `PHYTOZOME_BEARER`) and fails loudly if absent. Five new species wired into [config/phytozome_gtf_sources.yaml](config/phytozome_gtf_sources.yaml) + [resources/phytozome/manifest.json](resources/phytozome/manifest.json): amborella, chlamydomonas, physcomitrella, **vitis_vinifera** (`VIT_`/`GSVIVT`/`GTVIVG`), **solanum_tuberosum** (`PGSC`). `parse_ids.py` + `PREFIX_TO_SPECIES` route the new prefixes. This expands **transcript→coordinate resolution**, not extraction.
2. **Genome FASTA for extraction — ✅ CLOSED (2026-07-14).** The 203 `assembly_not_cached` phytozome rows resolved to coordinates but couldn't be sliced (the manifest downloaded GFF3 annotations, not genome FASTAs). Now fixed by the new `download_phytozome_fasta` rule + per-species cache keys — see [Phytozome genome-FASTA extraction](#phytozome-genome-fasta-extraction-2026-07-14). Rerun pending to reflect in `results/`.

**Net:** the token fix grows the *resolved* set (new plant/algae/moss/grape/potato transcripts get coordinates); closing the extraction gap for phytozome was a separate follow-up — **now landed (2026-07-14), see [Phytozome genome-FASTA extraction](#phytozome-genome-fasta-extraction-2026-07-14).**

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

## Run failures — 2026-07-13 (`fix/phytozome-jgi-download-id`)

A forced full run (`.snakemake/log/2026-07-13T004801…log`) completed 51/71 steps
then exited on three independent failures. Root-caused below. Two are real bugs
(NCBI GenBank, phytozome key mismatch); one is transient BioMart flakiness.

### TB-1. `resolve_ncbi_genbank` crashes the whole batch on 4 junk IDs ✅ DONE (2026-07-13)
- **Symptom:** `RuntimeError: Some IDs have invalid value and were omitted. Maximum ID value 18446744073709551615` in `_epost` ([ncbi_genbank_fetcher.py:206](workflow/scripts/ncbi_genbank_fetcher.py#L206)); the rule exits non-zero and takes down the run.
- **Root cause:** [resolve_ncbi_genbank.py:98](workflow/scripts/resolve_ncbi_genbank.py#L98) feeds **every** `db_source == ncbi` row straight to `Entrez.epost(db="nuccore", id=",".join(accessions))`. Of the 1,619 IDs, four are malformed — `xm_003`, `np_205`, `np_206`, `nc_201` (lowercase prefix, truncated numeric; not real accessions). NCBI rejects them, returns an `<ERROR>` element, and Biopython's `Entrez.read` **raises on it**, discarding the ~1,615 valid `XM_`/`NM_` accessions in the same batch. One bad ID poisons the whole epost.
- **Fix:** quarantine non-accession IDs before eposting. Filter `accessions` to a real accession shape (e.g. `^[A-Z]{2}_?\d{6,}(\.\d+)?$` after `.strip().upper()`) and write the rejects to the unresolved output instead of passing them to `fetch()`. Ideally also fix upstream in `parse_ids` so `xm_003`/`nc_201` classify as `unknown`, not `ncbi`.
- **Note:** these 4 junk IDs are the *only* reason the batch fails; the other 1,615 resolve fine once they're removed.
- **✅ Fix landed (2026-07-13):** [resolve_ncbi_genbank.py](workflow/scripts/resolve_ncbi_genbank.py) now filters `df_ncbi["transcript_id"]` through `ACCESSION_RE` (`^[A-Z]{1,2}_?\d{5,}(\.\d+)?$`, matched after `.strip().upper()`) **before** `fetcher.fetch()`. Non-matching IDs are routed to the unresolved output with `reason=invalid_accession` and never reach `Entrez.epost`, so one malformed ID can no longer poison the batch. The `fetch([])`/all-junk case is already guarded in [ncbi_genbank_fetcher.py:175](workflow/scripts/ncbi_genbank_fetcher.py#L175). Verified: the 4 junk IDs are rejected and real `XM_`/`NM_`/`NR_`/`XR_`/`AB` accessions pass. Upstream `parse_ids` classification left as-is (defensive filter at the epost boundary is sufficient).

### TB-2. `download_phytozome_gtf` — config key ≠ `gtf:` basename (empty logs) ✅ DONE (2026-07-13)
- **Symptom:** 4 SLURM jobs FAILED with **empty rule logs**: `Csinensis_154_v1.1.gene`, `Sbicolor_313_v3.1.gene`, `Rcommunis_119_v0.1.gene`, `Osativa_204_v7.0.gene`. The other 6 phytozome species succeeded.
- **Root cause:** the `{species}` wildcard is derived from the `gtf:` **filename** (`resolve_phytozome_gtf` requests `phytozome_sources[s]["gtf"]`, so Snakemake matches the download rule's `resources/phytozome/{species}.gff3.gz` output → `{species}=Csinensis_154_v1.1.gene`). But the config entries are keyed by **organism name** (`citrus_sinensis`, `sorghum_bicolor`, `ricinus_communis`, `oryza_sativa`). [download_phytozome_gtf.py:92](workflow/scripts/download_phytozome_gtf.py#L92) does `sources.get(species)` → `None` → no `genome_id`, no manifest entry, no URL → `WorkflowError` raised **before** any `log_path.write_text` (hence the empty logs). The 6 that succeeded have `gtf` basename == config key.
- **This reopens Priority 0 §"Config path mismatch (fixed)".** The 2026-07-12 fix repointed the 4 `gtf:` paths to on-disk `.gene.gff3.gz` files and relied on their mtime; the branch's genome_id-based rewrite (commit `55a8e86`) now re-resolves via `genome_id` keyed on the filename-wildcard, so the mismatch is live again — and a forced run rebuilds the `protected()` outputs regardless of the local files.
- **✅ Fix landed (2026-07-13) — folder-per-species layout.** Chosen over the two options above because it makes the `{species}` wildcard *genuinely* the species everywhere (config key, manifest key, wildcard all align), instead of papering over the mismatch.
  - **Layout.** Every Phytozome GFF3 now lives at `resources/phytozome/<species>/<source_file_name>.gff3.gz`. The `<species>` folder equals the config key, so it becomes the `{species}` wildcard and `sources.get(species)` resolves the `genome_id`. The inner file keeps the JGI source name for traceability. Rule output is `resources/phytozome/{species}/{gff}` with `wildcard_constraints: species=r"[^/]+"` ([resolve_phytozome_gtf.smk](workflow/rules/resolve_phytozome_gtf.smk)); log/benchmark carry both wildcards. The 6 already-working files were `mv`'d into their folders (kept basename, `touch`'d to postdate `manifest.json`); the 4 broken ones download fresh.
  - **`.gene`, not `gene_exons`.** The 4 broken `gtf:` paths stay on the `.gene.gff3.gz` variant. Pinned via new `portal_file_name` entries in [manifest.json](resources/phytozome/manifest.json) (the JGI heuristic otherwise prefers `gene_exons`). Bonus: `Sbicolor_313_v3.1.gene.gff3.gz` is `RESTORED` where its `gene_exons` is `PURGED`, so sorghum no longer needs a restore. Version verified against real input IDs: RNAChallenge Citrus IDs are `orange1.1g...m` (v1.1 namespace) → `Csinensis_154_v1.1.gene.gff3.gz`.
  - **No more silent failures.** `download_phytozome_gtf.py` wraps its body so any exception writes to the rule log before re-raising.
  - **Live-verified:** `prefer_name` resolves the `.gene` file for all 4 (citrus/sorghum/ricinus `RESTORED`; oryza `Osativa_204_v7.0.gene.gff3.gz` `PURGED` → first run fires a JGI restore, rerun ≤24 h). Dry-run: 4 download + 1 resolve job, no `ProtectedOutput`/`MissingInput` errors, 6 moved files up-to-date under `--rerun-triggers=mtime`.

### L5. `download_metadata_table` / solanum_lycopersicum — renamed mart dataset, NOT transient 🟡 ✅ DONE (2026-07-13)
- **Symptom:** `Empty response from BioMart for slycopersicum_eg_gene` on all 3 attempts. Only tomato failed; the other 4 species succeeded.
- **Cause (verified on live mart 2026-07-13):** **not** flakiness. Ensembl Plants moved tomato to a new assembly (**SL4.0**) and a new mart dataset, `slgca000188115v5cm_eg_gene`. The old `slycopersicum_eg_gene` slug (hard-coded in `download_metadata_table.py:34` **and** `biomart_plant_batch.py:52`) no longer exists → empty result. The other 4 species kept their `*_eg_gene` slugs, so they still work.
- **Why "fix the slug" is the wrong fix:** the SL4.0 dataset also changed the **ID namespace** — transcripts are `mRNA-Solyc…4.1`, genes `gene-Solyc…`, with bumped versions (our gene `Solyc06g068790.2.1` → mart `mRNA-Solyc06g068790.4.1`). BioMart filters on exact `ensembl_transcript_id`, so our SL3.0-era `Solyc…2.1` input IDs cannot match SL4.0 regardless of slug, and SL4.0 coordinates would *disagree* with the assembly the rest of the tomato path uses.
- **Resolution:** tomato is already fully resolved via the **plant_gtf** path (70 rows) from `release-30 / GCA_000188115.2` (SL3.0), whose `Solyc…2` namespace matches the input. The metadata table for tomato is therefore redundant *and* unsatisfiable from the current mart. See **L6** for the same root cause hitting the batch rule at scale.
- **✅ Fix landed (2026-07-13):** `solanum_lycopersicum` **commented out** (not deleted) of `external_metadata_tables` in [config/config.yaml](config/config.yaml), with an inline note pointing here. Kept as a comment so the intent is visible and it's trivially restorable if Ensembl Plants ever re-serves an SL3.0-compatible tomato mart. Dry-run confirms no `download_metadata_table` job for tomato; the other 4 species still run.

### L6. `biomart_plant_batch` resolves **0 rows** — release-iteration finds nothing 🔴 ⏸ DETACHED (2026-07-13)

> **⤷ Superseded in part (2026-07-14):** the AGPv3/v44 namespace→assembly reasoning below is
> now **verified and implemented** — see [Plant resolution — namespace-keyed v44/Phytozome-v11](#plant-resolution--namespace-keyed-v44phytozome-v11-2026-07-14-branch-fixrnaplonc-ensembl-plants-v44).
> Maize `GRMZM`→Phytozome v11 (AGPv3), `Zm00001d`→Ensembl Plants v44 (AGPv4), plus new
> Phytozome brachypodium/manihot. Rice-suffix tie-break (open decision #2) is moot: PreLnc(v44)
> and LGC(r30) share IRGSP-1.0, so coordinates agree.

**⏸ Detached as a quick fix (2026-07-13) — answers open decision #3 below.** Because the rule contributes **0 unique rows** (the 543 covered rows all resolve via plant_gtf/phytozome) while the flaky live-mart release loop makes runs slow and non-deterministic, [biomart_plant_batch.smk](workflow/rules/biomart_plant_batch.smk) no longer calls `biomart_plant_batch.py`. The rule now runs a trivial inline `run:` block that emits a **header-only** `biomart_resolved.tsv` (keeps `merge_resolved` satisfied) and passes every input ID through to `biomart_unresolved.tsv` with `reason=biomart_detached`. No network, `runtime` dropped 30→5 min / `mem_mb` 4096→512. Nothing else changed: `merge_resolved` is the only consumer of `biomart_resolved`, nothing consumes `biomart_unresolved`, and the plant GTF/Phytozome/Gramene fallbacks run off `external_*`, not BioMart. Dry-run parses clean. **Re-enable** by restoring `script: ../scripts/biomart_plant_batch.py` (script kept intact). This is a stop-gap, not the durable fix — the 1603-row gap still needs the pinned AGPv3/IRGSP GTFs (Tier 1 below).

The rule (`workflow/scripts/biomart_plant_batch.py`) tries a list of mart endpoints in order — `current, release-60, 59, 58, 57` (`RELEASE_ENDPOINTS`) — per species, retrying until one "works". On the last run it produced **1 header line and 0 data rows** in `biomart_resolved.tsv`; all **2906** input rows fell through to `biomart_unresolved.tsv` with reasons like `biomart_no_match_<species>_all_releases_exhausted`. The rule is currently dead weight.

**1 — Which transcripts depend on this rule.** Input is `results/external_unresolved.tsv` (plant IDs the local metadata-table fast path could not resolve). Of 2906 rows, **2146** are genuine plant-ID candidates (the other 760 are `not_plant_id`). By namespace:

| Namespace | Count | Species | Status |
|---|---|---|---|
| `Os…t…_NN` (RAP-DB) | 958 | rice | ❌ unresolved by any source |
| `GRMZM2G…_T0N` (AGPv3) | 640 | maize | ❌ unresolved by any source |
| `PGSC…` | 467 | potato | ✅ already via plant_gtf |
| `Solyc…`, `LOC_Os…`, `orange…` | ~81 | tomato/rice-MSU/citrus | ✅ already via plant_gtf / phytozome |

Net: **1603 of 2146 are resolved by *nothing*** (958 rice RAP + 640 maize GRMZM + a few stragglers); the remaining 543 are already covered elsewhere (485 plant_gtf, 58 phytozome). So BioMart's *unique* contribution today is **zero**.

**Root cause = same as L5:** the input IDs are **old-assembly namespaces that current Ensembl Plants no longer serves**. Maize `GRMZM2G…` is **AGPv3** (retired; current mart is `Zm-B73-NAM-5.0` → `Zm00001eb…` IDs). Rice `Os01t0102850_00` is **RAP-DB** native format; the live mart carries the hyphen form `Os01t0873800-01`. No live release matches → "all releases exhausted".

**2 — Origin tools for these 2146 candidates** (`results/tool_source_map.tsv`): PreLnc 1382, RNAplonc 1023, LGC 210, CPPred 6, PLEK 4, CPC2 3 (plus singletons). These are the **plant lncRNA classifiers** (RNAplonc and PreLnc are plant-specific), i.e. exactly the tool datasets this challenge is benchmarking — so the 1603 gap is worth closing.

**3 — Better source (no live BioMart, mirror the tomato precedent).** Both missing namespaces are available as **pinned Ensembl Plants FTP GTFs**, the same mechanism `plant_gtf_sources.yaml` already uses for tomato/potato. Assembly/release now **confirmed from the tool articles** (see `article_notes.md`, `dani_notes.md`):
- **Maize GRMZM →** `release-31/gtf/zea_mays/Zea_mays.AGPv3.31.gtf.gz`. Confirmed it carries `transcript_id "GRMZM2G059865_T01"` — drop-in `plant_gtf` source. Provenance: PreLnc/RNAplonc; GRMZM (AGPv3) is the **GreeNC / non-coding** lineage — PreLnc *coding* maize is v44 `Zm00001d…`, a **different namespace** that needs the v4 build, not this one.
- **Rice `Os…t…_0N` →** do **not** normalize to modern `-01`. The articles show LGC/PreLnc used **release-26** (coding, `IRGSP-1.0.26`) and **release-30** (ncrna, `IRGSP-1.0.30`), whose GTFs natively use the `Os…t…_0N` underscore format that matches our input. Pin release-26 + release-30 rice GTFs (URLs in the notes). This supersedes my earlier release-60 + normalization guess.

**Recommendation:** stop depending on the flaky live-mart release loop for plant coordinates. Add the AGPv3 maize GTF/FASTA to `plant_gtf_sources.yaml`, add RAP-ID normalization for rice, and treat `biomart_plant_batch` as a best-effort last resort (or retire it). This removes the 1603-row gap and the non-deterministic release iteration in one move.

#### Resolution strategy — primary + article fallback (for decision)

Two tiers, cheapest first. **Tier 1 handles the common case; Tier 2 is the fallback when the namespace does not uniquely pin an assembly.**

**Tier 1 — namespace → assembly (deterministic, no article needed).**
The ID prefix version-locks the build for these families, so we can pin the matching Ensembl Plants FTP GTF directly:
| Namespace | Assembly (implied by prefix) | Pinned source | Verified |
|---|---|---|---|
| `GRMZM2G…_T0N` | maize AGPv3 | `release-31/…/Zea_mays.AGPv3.31.gtf.gz` | ✅ carries `GRMZM2G059865_T01` |
| `Os…t…_NN` (RAP) | rice IRGSP-1.0 / RAP-DB | existing `release-60` rice GTF (+ `_NN`→`-NN` normalization) | ✅ base IDs present |
| `Solyc…2.1` | tomato SL3.0 / ITAG3 | `release-30` tomato GTF (already configured) | ✅ resolves 70 rows |
| `PGSC…` | potato | already configured | ✅ resolves 467 |

**Tier 2 — origin tool → source article → assembly (fallback).** When a namespace is ambiguous, retired, or a normalization guess is unverifiable, the tool's **publication** is often the only ground truth ("genome downloaded from &lt;DB&gt; release &lt;N&gt;"). We do not have the articles in-repo, but **@dgruano does**. The 1603 currently-unresolved rows trace back to only **three papers**, so this is a bounded manual step:

| Origin tool | Owns (unresolved rows) | Article should confirm |
|---|---|---|
| **PreLnc** | 958 rice (100%) + 400 maize | rice IRGSP/RAP release; maize AGPv build |
| **RNAplonc** | 561 maize (plant-specific tool) | maize AGPv build (expected AGPv3) |
| **LGC** | 191 rice (overlaps PreLnc) | rice release (cross-check PreLnc) |

**How the fallback plugs in:** the article-confirmed assembly just selects/confirms the `plant_gtf_sources.yaml` release to pin — it does **not** add a new resolver path. Tier 1 already proposes AGPv3 (maize) and IRGSP-1.0 (rice); Tier 2 is there to (a) confirm those two guesses against what PreLnc/RNAplonc actually used, and (b) resolve the rice `_00`↔`-01` suffix question if the FTP spot-check is inconclusive.

**Open decisions for @dgruano:**
1. Pin AGPv3 maize + rice normalization now (Tier 1), and use the articles only to *confirm* — or read the articles first, then pin?
2. Rice suffix `_00`: does PreLnc's source treat it as the primary transcript (`-01`) or a gene-level record? (Article or FTP spot-check.)
3. Retire `biomart_plant_batch` entirely, or keep it as a degraded last resort after the pinned GTFs?

### TB-3. Ensembl headers mislabelled `ncbi` via embedded-RefSeq substring ✅ FIXED — ⏸ DEFERRED MERGE (branch `fix/parse-ids-embedded-ncbi-substring`)
- **Symptom:** the first 11 rows of `results/ncbi_genbank_unresolved.tsv` are `invalid_accession` IDs like `NC_010`, `NT_003`, `NP_201`, `xm_003`, `np_205` — surfaced by the TB-1 quarantine filter. They are **not** NCBI IDs at all.
- **Root cause:** these are **Ensembl** transcripts whose IDs were extracted faultily. `extract_transcript_id` returns the whole underscore-joined header (`ENST00000570013.5_ENSG..._GANC_010_...`), which the anchored Ensembl `DB_PATTERN` can't match, so it falls to `find_embedded_accession`. There, the NCBI RefSeq `EMBEDDED_PATTERN` is listed **first** and searched unanchored with `\d+`, so it matches `NC_010` inside `GA``NC_010` (a gene-symbol token) before the Ensembl embedded pattern (parse_ids.py:520) ever runs. Mapping: `GANC_010`→`NC_010`, `PRNT_003`→`NT_003`, `PCNP_201`→`NP_201`, `Lexm_003`→`xm_003`, `Pianp_205`→`np_205`, `GMNC_002`→`NC_002`.
- **Fix:** tightened the NCBI embedded regex in [parse_ids.py](workflow/scripts/parse_ids.py) with a left word-boundary `(?<![A-Za-z0-9])` and `\d{6,}` (real RefSeq accessions never start mid-word and always carry ≥6 digits). The Ensembl embedded pattern now wins → these 11 route to the `ensembl` resolver. Verified: junk substrings return `None`, real embedded RefSeq (`gi|…|ref|NM_000014.6|`, `lcl|NC_000001.11`) still match.
- **⏸ Why deferred:** `parse_ids.py` is **Stage 1** — changing it invalidates the whole DAG and forces a full pipeline rerun. Per decision (2026-07-13) the fix is committed to isolated branch **`fix/parse-ids-embedded-ncbi-substring`** (commit `208b215`), to be **merged when a full rerun is acceptable**. Not on `fix/phytozome-jgi-download-id`.

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

## Refactors (deferred, non-blocking)

### R1. Rename `ncbi_chromosome_resolved.tsv` — misleading name 🟡
The file produced by `resolve_ncbi_chromosome_accessions` is **not** NCBI-only:
that rule passes *every* resolved row through and only *patches* NCBI chromosome
names, so the file is the full merged resolved table (carries all phytozome,
ensembl, noncode, … rows). The whole download+extract stage keys off it, which
makes the name actively misleading. Rename to something honest, e.g.
`resolved_chrom_patched.tsv`. Touches `resolve_ncbi_chromosome_accessions.smk`
(output), `download_assemblies.smk` + `extract_sequences.smk` (inputs),
`prepare_accession_list` + `aggregate_downloads.py` docstrings. Deferred because
it's a pure rename with no behavior change; do it alongside a run that already
regenerates these outputs. _Noted 2026-07-14 during phytozome-FASTA work (which
deliberately reads the honestly-named `resolved_ids.tsv` for its fan-out)._


### R2. Change formatted strings for snakemake pathvars
This is the appropriate coding pattern.

### R3. Double check hardcoded sources and mapping dicts
- PREFIX_TO_SPECIES
- others

### R4. Group output files in folders to avoid clutter

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
11. ~~**[bottleneck] Populate `fasta_url` for URL-less non-NCBI rows**~~ ✅ **DONE (2026-07-14)** — plant stream fills by `organism` fallback in merge (541/541); `plant_gtf` and phytozome rows carry URL metadata. The phytozome genome-FASTA gap (203 rows, previously deferred 2026-07-12) is now closed by `download_phytozome_fasta` — see [Phytozome genome-FASTA extraction](#phytozome-genome-fasta-extraction-2026-07-14). Rerun pending to confirm the bucket clears.
12. **[bottleneck, optional] Handle roman-numeral + dead assemblies** — roman↔arabic chrom mapping for worm/noncode `chrIV`; the "~48 dead assemblies" collapsed to **2** after the FTP-URL refactor + retry (2026-07-12): `GCF_000001215.3` (404, use `.4`) and `GCF_000001545.5` (`ponAbe2`, folder removed). Bump to current versions or accept as lost. Lowest ROI. See **L2**.

### Run-failure fixes (new — from 2026-07-13 run)

16. ~~**[blocker] Filter junk IDs before NCBI epost**~~ ✅ **DONE (2026-07-13)** — `resolve_ncbi_genbank.py` now quarantines non-accession IDs (`xm_003`, `np_205`, `np_206`, `nc_201`) via `ACCESSION_RE` before `fetcher.fetch()`/`Entrez.epost`, routing them to unresolved as `invalid_accession`. Unblocks the ~1,615 valid IDs that the crash was discarding. See **TB-1**.
17. ~~**[blocker] Fix phytozome config-key/`gtf:`-basename mismatch**~~ ✅ **DONE (2026-07-13)** — moved to a `resources/phytozome/<species>/<source_name>.gff3.gz` folder layout so the `{species}` wildcard equals the config key and `sources.get(species)` resolves `genome_id`; kept the `.gene` (not `gene_exons`) variant, pinned via `manifest.json` `portal_file_name`; wrapped `download_phytozome_gtf.py` so no failure is silent. citrus/sorghum/ricinus download immediately; oryza is PURGED (restore workflow). See **TB-2**.
18. **[not transient] Tomato mart dataset renamed** — `download_metadata_table` for solanum_lycopersicum fails because Plants moved tomato to SL4.0 (`slgca000188115v5cm_eg_gene`, new `mRNA-Solyc…` namespace); the old slug is gone. Tomato is already resolved via plant_gtf (SL3.0). Drop it from `external_metadata_tables`. Same root cause hits `biomart_plant_batch` at scale (1603 rice/maize rows) — see **L5**/**L6**.
19. **[deferred-merge] Merge Ensembl-mislabel parse fix** — ✅ fixed on branch **`fix/parse-ids-embedded-ncbi-substring`** (commit `208b215`); **not merged** because it's a Stage-1 change that forces a full rerun. Merge when a full rerun is scheduled. See **TB-3**.

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

### Priority -1: Fix 2026-07-13 run blockers (NEW — actionables #16–18)
- **#16 — NCBI epost crash (TB-1):** ✅ **DONE (2026-07-13)** — junk IDs filtered before epost in `resolve_ncbi_genbank.py`.
- **#17 — phytozome key mismatch (TB-2):** ✅ DONE (2026-07-13) — folder-per-species layout (`resources/phytozome/<species>/<source_name>.gff3.gz`), `.gene` pinned via manifest, failures now logged. oryza still needs a JGI restore (PURGED).
- **#18 — plant BioMart (L5/L6):** 🔴 not transient — Plants retired old assemblies/namespaces. Tomato covered by plant_gtf; 1603 rice(RAP)+maize(GRMZM) rows need pinned AGPv3 GTF + RAP normalization, not the mart.
- **#19 — Ensembl mislabel (TB-3):** ✅ fixed on branch `fix/parse-ids-embedded-ncbi-substring` (`208b215`), ⏸ **deferred merge** — Stage-1 change forces a full rerun; merge when one is scheduled.

### Priority 0: Regenerate the phytozome resolved table ✅ DONE (2026-07-14)
Resolver regenerated with per-species `phytozome_<species>` keys; **all 503 phytozome rows now extract** (`assembly_not_cached` 614 → 111, extraction 15,319 → **15,820**). Verified in [Run refresh 2026-07-14](#run-refresh--2026-07-14-phytozome-fix-landed). The phytozome deferral (open since Apr 27) is fully closed.

### Priority 0b: plant_gtf tomato chrom naming (✅ DONE — was 70 rows, the largest bucket)
`chrom_not_found` plant_gtf 71, all `GCA_000188115.2` (tomato SL3.0) bare `1`–`12`. The FASTA seqids are GenBank accessions `CM001064.2`–`CM001075.2` (verified from the `.fai`), so `1`→`CM001064.2` … `12`→`CM001075.2`. **Root cause:** the NCBI `assembly_report.txt` maps `1`→RefSeq `NC_015438.2`, but the cached Ensembl GCA-toplevel FASTA uses GenBank seqids (`CM001064.2`–`CM001075.2`). `resolve_chrom_key` didn't try the GenBank alias, causing `chrom_not_found`. **Fixed** by offering both RefSeq and GenBank accessions as candidates in `load_chrom_translation`. Verified: 0 remaining `chrom_not_found` rows for `GCA_000188115.2` in `results/extraction_failed.tsv`.

**⚠️ Unblock to run *now* while JGI restores are pending (2026-07-12):** `resolve_phytozome_gtf` needs **all** configured GFF3s present and sits **upstream of merge→extract→report**, so any single missing phytozome file blocks the whole pipeline (not just the phytozome stream). Two problems found + fixed so a full run can proceed today:
1. **Config path mismatch (fixed).** `citrus/sorghum/ricinus/oryza` pointed at `*.gene_exons.gff3.gz` filenames not on disk (disk has `*.gene.gff3.gz`). Because the `download_phytozome_gtf` script keys by *species name* while the resolve rule requests files by their `gtf:` *basename*, the mismatched paths triggered JGI downloads that fail with "manifest entry not found". Repointed the 4 `gtf:` paths to the existing `.gene.gff3.gz` files (verified they carry `mRNA` features the resolver reads). `touch`ed them so they postdate `manifest.json` (else mtime forces the same failing rebuild).
2. **JGI-only species deferred.** `solanum_tuberosum` (no local file), `vitis_vinifera` + `physcomitrella_patens` (**PURGED**, restore reqs pending, ≤24h) are commented out in `config/phytozome_gtf_sources.yaml` with dated re-enable notes. amborella + chlamydomonas are RESTORED and already downloaded, so they stay in.

Result: dry-run drops from 37 → 29 jobs with **0 phytozome download jobs** and no missing inputs. Run `snakemake --profile profiles/default --rerun-triggers=mtime`. When JGI restores vitis/physcomitrella (poll `files.jgi.doe.gov/request_archived_files/requests/652368`), uncomment them (and add solanum_tuberosum's GFF3), then rerun `resolve_phytozome_gtf` + downstream.

### Priority 1: Residual non-NCBI URL gaps 🟠 NEXT BIGGEST — 86 recoverable rows
- `assembly_not_cached` residual is **111**: phytozome cleared (Priority 0 ✅), leaving **ensembl `ARS-UCD2.0` (cattle, 53) + flybase `BDGP6` (fly, 33) = 86 recoverable**, plus 25 genuine L2 dead ends (ponAbe2 24 + dmel 1).
- These 86 are non-NCBI genomes with **bare assembly names** (`ARS-UCD2.0`, `BDGP6`) that don't carry a `fasta_url`, so the manifest fan-out has nothing to download. Fill reliable `fasta_url` for the cattle/fly assemblies in their resolver/config paths (Ensembl / FlyBase FTP), same manifest mechanism as plant_gtf. This is the largest remaining *recoverable* extraction bucket after the tomato-chrom fix (Priority 0b).
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
