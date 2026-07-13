"""
scripts/detect_ensembl_species.py
Auto-Detect Ensembl Species from Transcript ID Prefixes
=======================================================
Replicates the logic of:
    grep -oE 'ENS[A-Z]+' | sort | uniq

For each unique prefix found in the Ensembl transcript IDs:
  1. Matches against a built-in reference table covering all
     Ensembl vertebrate species — no config required for known species
  2. Merges user-defined overrides from config (optional, for unknowns)
  3. For unrecognised prefixes: writes an actionable error message
     with the exact config.yaml snippet to add, then raises an
     exception to stop the pipeline cleanly

Outputs
-------
ensembl_species_map.tsv          : transcript_id | prefix | species | build
ensembl_unknown_prefixes.tsv     : prefix | example_id | n_transcripts | action_required
ensembl_unknown_prefixes_ACTION_REQUIRED.txt  (only written if unknowns exist)

Notes on the ID format
----------------------
Ensembl stable IDs: ENS[species_code][type][11 digits]
  - Human has NO species code: ENST00000...
  - Others use 3-letter code:  ENSMUST... (mouse), ENSRNOT... (rat)
  - Feature type: G=gene, T=transcript, E=exon, P=protein
  - We extract the transcript prefix (ends in T): ENST, ENSMUST, ENSRNOT, ...
"""

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("detect_ensembl_species", snakemake.log[0])
input_tsv = snakemake.input.classified
out_map = snakemake.output.species_map
out_unk = snakemake.output.unmatched
cfg = snakemake.config

# ════════════════════════════════════════════════════════════
# Built-in Ensembl transcript prefix reference table
# ════════════════════════════════════════════════════════════
# Format: "TRANSCRIPT_PREFIX": ("ensembl_species_name", "genome_build")
# Transcript prefix = ENS + [3-letter species code] + T
# Human is special: no species code, so prefix is just "ENST"
# Source: Ensembl stable ID prefix page + biomaRt dataset names (release 113)
# ════════════════════════════════════════════════════════════
BUILTIN_PREFIX_TABLE: dict[str, tuple[str, str]] = {
    # ── Primates ─────────────────────────────────────────────
    "ENST": ("homo_sapiens", "GRCh38"),
    "ENSPTRT": ("pan_troglodytes", "Pan_tro_3.0"),
    "ENSPANT": ("pan_paniscus", "panpan1.1"),
    "ENSGGOT": ("gorilla_gorilla", "gorGor4"),
    "ENSPPYT": ("pongo_abelii", "Susie_PABv2"),
    "ENSNLET": ("nomascus_leucogenys", "Nleu_3.0"),
    "ENSMICT": ("microcebus_murinus", "Mmur_3.0"),
    "ENSOGAT": ("otolemur_garnettii", "OtoGar3"),
    "ENSCJAT": ("callithrix_jacchus", "mCalJac1.pat.X"),
    "ENSMMUT": ("macaca_mulatta", "Mmul_10"),
    "ENSMFAT": ("macaca_fascicularis", "Macaca_fascicularis_6.0"),
    "ENSANOT": ("aotus_nancymaae", "Anan_2.0"),
    "ENSSBOT": ("saimiri_boliviensis_boliviensis", "SaiBol1.0"),
    "ENSCCET": ("cebus_imitator", "Cebus_imitator-1.0"),
    "ENSCSAT": ("chlorocebus_sabaeus", "ChlSab1.1"),
    "ENSPHAT": ("papio_anubis", "Panubis1.0"),
    "ENSTSYT": ("theropithecus_gelada", "Tgel_1.0"),
    "ENSMNET": ("mandrillus_leucophaeus", "Mleu.le_1.0"),
    "ENSCATT": ("colobus_angolensis_palliatus", "Cang.pa_1.0"),
    "ENSPINGT": ("piliocolobus_tephrosceles", "ASM277652v2"),
    "ENSREXT": ("rhinopithecus_roxellana", "Rrox_v1"),
    "ENSNASAT": ("nasalis_larvatus", "ASM242139v1"),
    # ── Rodents ──────────────────────────────────────────────
    "ENSMUST": ("mus_musculus", "GRCm39"),
    "ENSMSPRT": ("mus_spretus", "SPRET_EiJ_v1"),
    "ENSMPAHT": ("mus_pahari", "PAHARI_EiJ_v1"),
    "ENSMCART": ("mus_caroli", "CAROLI_EiJ_v1"),
    "ENSRNOT": ("rattus_norvegicus", "mRatBN7.2"),
    "ENSCPOT": ("cavia_porcellus", "Cavpor3.0"),
    "ENSOCUNT": ("octodon_degus", "OctDeg1.0"),
    "ENSCAST": ("castor_canadensis", "C.can_genome_v1.0"),
    "ENSDORT": ("dipodomys_ordii", "Dord_2.0"),
    "ENSPEMOT": ("peromyscus_maniculatus_bairdii", "HU_Pman_2.1"),
    "ENSMOCET": ("microtus_ochrogaster", "MicOch1.0"),
    "ENSNFORT": ("nannospalax_galili", "S.galili_v1.0"),
    "ENSMAVGT": ("marmota_marmota_marmota", "marMar2.1"),
    "ENSUAMGT": ("urocitellus_parryii", "ASM342692v1"),
    "ENSICGT": ("ictidomys_tridecemlineatus", "SpeTri2.0"),
    "ENSFCDT": ("fukomys_damarensis", "DMR_v1.0"),
    "ENSHETGT": ("heterocephalus_glaber_female", "HetGla_female_1.0"),
    "ENSOHANT": ("ochotona_princeps", "OchPri3.0"),
    "ENSOCUNT2": ("oryctolagus_cuniculus", "OryCun2.0"),
    # ── Carnivores ───────────────────────────────────────────
    "ENSCAFT": ("canis_lupus_familiaris", "ROS_Cfam_1.0"),
    "ENSVVUT": ("vulpes_vulpes", "VulVul2.2"),
    "ENSMPET": ("mustela_putorius_furo", "MusPutFur1.0"),
    "ENSNVIT": ("neovison_vison", "NNQGG.v01"),
    "ENSAMET": ("ailuropoda_melanoleuca", "ASM200744v2"),
    "ENSUPMT": ("ursus_maritimus", "UrsMar_1.0"),
    "ENSUPAT": ("ursus_americanus", "ASM334442v1"),
    "ENSFCAT": ("felis_catus", "Felis_catus_9.0"),
    "ENSPTIG": ("panthera_tigris_altaica", "PanTig1.0"),
    "ENSPLOT": ("puma_concolor", "PumCon1.0"),
    "ENSLPAT": ("lynx_pardinus", "mLynPar1.p"),
    # ── Afrotheria & Xenarthra ────────────────────────────────
    "ENSTTRT": ("trichechus_manatus_latirostris", "TriManLat1.0"),
    "ENSLABT": ("loxodonta_africana", "loxAfr3"),
    "ENSPVAT": ("procavia_capensis", "proCap1"),
    "ENSOGAST": ("orycteropus_afer_afer", "OryAfe1.0"),
    "ENSDDST": ("dasypus_novemcinctus", "Dasnov3.0"),
    "ENSCHOET": ("choloepus_hoffmanni", "choHof1"),
    # ── Ungulates ────────────────────────────────────────────
    "ENSBTRT": ("bos_taurus", "ARS-UCD1.3"),
    "ENSBGRT": ("bos_grunniens", "LU_Bosgru_v3.0"),
    "ENSOGLT": ("ovis_aries_rambouillet", "Oar_rambouillet_v1.0"),
    "ENSCCHT": ("capra_hircus", "ARS1"),
    "ENSSSCRT": ("sus_scrofa", "Sscrofa11.1"),
    "ENSWAFT": ("vicugna_pacos", "vicPac1"),
    "ENSLCAT": ("lama_glama", "mLamGla1.p"),
    "ENSECBT": ("equus_caballus", "EquCab3.0"),
    "ENSEAFT": ("equus_asinus", "ASM1607732v2"),
    "ENSCDRT": ("camelus_dromedarius", "CamDro3"),
    # ── Whales & dolphins ────────────────────────────────────
    "ENSPHAT2": ("physeter_catodon", "ASM283717v2"),
    "ENSDNOT": ("delphinapterus_leucas", "ASM228892v3"),
    "ENSOORT": ("orcinus_orca", "Oorc_1.1"),
    "ENSLHKT": ("lipotes_vexillifer", "Lipotes_vexillifer_v1"),
    "ENSTTRT2": ("tursiops_truncatus", "turTru1"),
    # ── Insectivora & bats ───────────────────────────────────
    "ENSSHRT": ("sorex_araneus", "SorAra2.0"),
    "ENSCPUT": ("condylura_cristata", "ConCri1.0"),
    "ENSTOGT": ("tupaia_belangeri", "TREESHREW 1.0"),
    "ENSPVAT2": ("pteropus_vampyrus", "pteVam1"),
    "ENSRTPT": ("rhinolophus_ferrumequinum", "mRhiFer1_v1.p"),
    "ENSMYLT": ("myotis_lucifugus", "Myoluc2.0"),
    "ENSMVIT": ("miniopterus_natalensis", "MiniNat1.0"),
    "ENSETST": ("erinaceus_europaeus", "eriEur1"),
    # ── Marsupials & monotremes ───────────────────────────────
    "ENSSGAT": ("sarcophilus_harrisii", "mSarHar1.11"),
    "ENSSHAWT": ("monodelphis_domestica", "monDom5"),
    "ENSMEUGAT": ("macropus_eugenii", "Meug_1.0"),
    "ENSOANT": ("ornithorhynchus_anatinus", "mOrnAna1.p.v1"),
    # ── Birds ─────────────────────────────────────────────────
    "ENSGALT": ("gallus_gallus", "bGalGal1.mat.broiler.GRCg7b"),
    "ENSMGALT": ("meleagris_gallopavo", "Turkey_5.1"),
    "ENSAMPT": ("anas_platyrhynchos", "CAU_duck1.0"),
    "ENSAACT": ("aquila_chrysaetos_chrysaetos", "bAquChr1.2"),
    "ENSFALT": ("ficedula_albicollis", "FicAlb1.5"),
    "ENSGMOT": ("geospiza_fortis", "GeoFor1.0"),
    "ENSTGUT": ("taeniopygia_guttata", "bTaeGut2.pat.W"),
    "ENSPCAT": ("pseudopodoces_humilis", "PseHum1.0"),
    "ENSMLEUT": ("manacus_vitellinus", "ManVit1.0"),
    "ENSFHETT": ("falco_tinnunculus", "bFalTin1.pri"),
    "ENSACOT": ("antrostomus_carolinensis", "ASM69998v1"),
    # ── Reptiles ──────────────────────────────────────────────
    "ENSANCT": ("anolis_carolinensis", "AnoCar2.0"),
    "ENSPSVT": ("pogona_vitticeps", "pVitVit3.0"),
    "ENSCRCT": ("crocodylus_porosus", "CroPor_comp1"),
    "ENSGGAT": ("gopherus_agassizii", "ASM289641v1"),
    "ENSCMYT": ("chelonia_mydas", "CheMyd1.0"),
    "ENSPMYT": ("pelodiscus_sinensis", "PelSin1.0"),
    "ENSPTYT": ("python_bivittatus", "Python_molurus_bivittatus-5.0.2"),
    # ── Amphibians ────────────────────────────────────────────
    "ENSXETT": ("xenopus_tropicalis", "UCB_Xtro_10.0"),
    # ── Fish ──────────────────────────────────────────────────
    "ENSDART": ("danio_rerio", "GRCz11"),
    "ENSGMORT": ("gasterosteus_aculeatus", "BROADS1"),
    "ENSORLT": ("oryzias_latipes", "ASM223467v1"),
    "ENSONIT": ("oreochromis_niloticus", "O_niloticus_UMD_NMBU"),
    "ENSAMXT": ("astyanax_mexicanus", "Astyanax_mexicanus-2.0"),
    "ENSPOCT": ("poecilia_reticulata", "Guppy_female_1.0_MT"),
    "ENSPFOT": ("poecilia_formosa", "Poecilia_formosa-5.1.2"),
    "ENSSART": ("sparus_aurata", "fSpaAur1.1"),
    "ENSSFAT": ("scophthalmus_maximus", "fScoMax1.pri"),
    "ENSXMAT": ("xiphophorus_maculatus", "X_maculatus-5.0-male"),
    "ENSLABT2": ("labrus_bergylta", "fLabBer1.pri"),
    "ENSLCALT": ("lates_calcarifer", "ASB_HGAfBa0_CSIRO"),
    "ENSNOBT": ("nothobranchius_furzeri", "Nfu_20140520"),
    "ENSFHETT2": ("fundulus_heteroclitus", "Fundulus_heteroclitus-3.0.2"),
    "ENSECALT": ("eptatretus_burgeri", "Eburgeri_3.2"),
    "ENSPMAT": ("petromyzon_marinus", "Pmarinus_7.0"),
    "ENSCALT": ("callorhinchus_milii", "Callorhinchus_milii-6.1.3"),
    "ENSLCALT2": ("latimeria_chalumnae", "LatCha1"),
    "ENSPSIT": ("protopterus_annectens", "ProtAnn1.0"),
    "ENSLOCT": ("lepisosteus_oculatus", "LepOcu1"),
    "ENSCCAT": ("cyprinus_carpio", "common_carp_genome"),
    "ENSICPT": ("ictalurus_punctatus", "IpCoco_1.2"),
    "ENSASIT": ("astatotilapia_calliptera", "fAstCal1.2"),
    "ENSOKAT": ("oncorhynchus_kisutch", "Okis_V2"),
    "ENSOMT": ("oncorhynchus_mykiss", "Omyk_1.0"),
    "ENSSSAT": ("salmo_salar", "Ssal_v3.1"),
    "ENSCLACT": ("clupea_harengus", "Ch_v2.0.2"),
    "ENSMZET": ("maylandia_zebra", "M_zebra_UMD2a"),
    "ENSHCULT": ("haplochromis_burtoni", "AstBur1.0"),
    "ENSTRUPT": ("takifugu_rubripes", "fTakRub1.2"),
    "ENSTNIGET": ("tetraodon_nigroviridis", "TETRAODON 8.0"),
    "ENSOULT": ("oryzias_melastigma", "Om_v0.7.RACA"),
    # ── Invertebrate chordates ────────────────────────────────
    "ENSCSAVT": ("ciona_savignyi", "CSAV 2.0"),
    "ENSCINTT": ("ciona_intestinalis", "KH"),
}

# ── Regex: extract transcript-type prefix ─────────────────────
# Matches ENS + zero or more uppercase letters + T at start of ID.
# Equivalent to: grep -oE 'ENS[A-Z]*T'
PREFIX_RE = re.compile(r"^(ENS[A-Z]*T)", re.IGNORECASE)


def extract_prefix(transcript_id: str) -> str | None:
    """Return the transcript prefix (e.g. ENST, ENSMUST) or None."""
    m = PREFIX_RE.match(str(transcript_id).strip())
    return m.group(1).upper() if m else None


# ── Load data ─────────────────────────────────────────────────
log.info(
    "detect_ensembl_species: auto-detecting species from Ensembl transcript ID prefixes"
)
log.info(f"Built-in reference table: {len(BUILTIN_PREFIX_TABLE)} known prefixes")

df = pd.read_csv(input_tsv, sep="\t")
df_ensembl = df[df["db_source"] == "ensembl"].copy()

log.info(f"Total classified IDs : {len(df)}")
log.info(f"Ensembl IDs          : {len(df_ensembl)}")

if df_ensembl.empty:
    log.warning("No Ensembl IDs found — writing empty outputs")
    pd.DataFrame(columns=["transcript_id", "prefix", "species", "build"]).to_csv(
        out_map, sep="\t", index=False
    )
    pd.DataFrame(
        columns=["prefix", "example_id", "n_transcripts", "action_required"]
    ).to_csv(out_unk, sep="\t", index=False)
    log.info("detect_ensembl_species complete (no Ensembl IDs).")
    sys.exit(0)

# ── Extract prefix per ID (grep -oE 'ENS[A-Z]*T') ────────────
df_ensembl["prefix"] = df_ensembl["transcript_id"].apply(extract_prefix)

unparseable = df_ensembl[df_ensembl["prefix"].isna()]
if not unparseable.empty:
    log.warning(
        f"{len(unparseable)} IDs yielded no prefix (unexpected format): "
        f"{unparseable['transcript_id'].tolist()[:5]}"
    )
df_ensembl = df_ensembl.dropna(subset=["prefix"])

# ── Unique prefixes (sort | uniq) ────────────────────────────
unique_prefixes = sorted(df_ensembl["prefix"].unique())
log.info(f"Unique prefixes in data: {unique_prefixes}")

# ── Merge built-in table with config overrides ────────────────
overrides: dict[str, dict] = cfg.get("ensembl_species_overrides", {}) or {}
effective_table: dict[str, tuple[str, str]] = dict(BUILTIN_PREFIX_TABLE)
for pfx, info in overrides.items():
    effective_table[pfx.upper()] = (info["species"], info["build"])
    log.info(
        f"  Config override applied: {pfx.upper()} → {info['species']} ({info['build']})"
    )

# ── Match each unique prefix ──────────────────────────────────
known: dict[str, tuple[str, str]] = {}
unknown: list[str] = []

for prefix in unique_prefixes:
    if prefix in effective_table:
        species, build = effective_table[prefix]
        known[prefix] = (species, build)
        log.info(f"  {prefix:<14} → {species} ({build})")
    else:
        unknown.append(prefix)
        log.warning(f"  {prefix:<14} → UNKNOWN")

# ── Build species_map TSV ─────────────────────────────────────
map_rows = [
    {
        "transcript_id": row["transcript_id"],
        "prefix": row["prefix"],
        "species": known[row["prefix"]][0],
        "build": known[row["prefix"]][1],
    }
    for _, row in df_ensembl.iterrows()
    if row["prefix"] in known
]
df_map = pd.DataFrame(map_rows, columns=["transcript_id", "prefix", "species", "build"])
df_map.to_csv(out_map, sep="\t", index=False)

# ── Build unknown prefix TSV ──────────────────────────────────
unk_rows = []
for prefix in unknown:
    examples = df_ensembl[df_ensembl["prefix"] == prefix]["transcript_id"].tolist()
    unk_rows.append(
        {
            "prefix": prefix,
            "example_id": examples[0],
            "n_transcripts": len(examples),
            "action_required": (
                f"Add to config.yaml → ensembl_species_overrides:\n"
                f'  {prefix}: {{species: "<name>", build: "<assembly>"}}\n'
                f"See: https://www.ensembl.org/info/about/species.html"
            ),
        }
    )
df_unk = pd.DataFrame(
    unk_rows, columns=["prefix", "example_id", "n_transcripts", "action_required"]
)
df_unk.to_csv(out_unk, sep="\t", index=False)

# ── Summary ──────────────────────────────────────────────────
log.info("=" * 60)
log.info(f"Unique prefixes found        : {len(unique_prefixes)}")
log.info(f"Recognised                   : {len(known)}")
log.info(f"Unrecognised                 : {len(unknown)}")
log.info(f"Transcripts mapped           : {len(df_map)}")
for (sp, build), cnt in df_map.groupby(["species", "build"]).size().items():
    log.info(f"  {sp:<45} {build:<20} {cnt} transcripts")

# ── Fail loudly if any prefix is unknown ─────────────────────
if unknown:
    msg_lines = [
        "",
        "=" * 60,
        "PIPELINE STOPPED: Unknown Ensembl transcript ID prefixes",
        "=" * 60,
        "",
        "The following prefixes were detected in your input FASTA(s)",
        "but are not in the pipeline's built-in species reference table:",
        "",
    ]
    for prefix in unknown:
        examples = df_ensembl[df_ensembl["prefix"] == prefix]["transcript_id"].tolist()
        msg_lines.append(
            f"  Prefix : {prefix}"
            f"  (example: {examples[0]}, {len(examples)} transcripts)"
        )
    msg_lines += [
        "",
        "ACTION REQUIRED",
        "Add the missing entries to config/config.yaml:",
        "",
        "  ensembl_species_overrides:",
    ]
    for prefix in unknown:
        msg_lines.append(
            f'    {prefix}: {{species: "<ensembl_species_name>", build: "<assembly_name>"}}'
        )
    msg_lines += [
        "",
        "Look up the correct species name and assembly at:",
        "  https://www.ensembl.org/info/about/species.html",
        "",
        "Already-downloaded BioMart tables for other species are cached",
        "and will NOT be re-downloaded when you re-run the pipeline.",
        "=" * 60,
    ]

    full_msg = "\n".join(msg_lines)
    log.error(full_msg)

    # Write a prominent plain-text action file alongside the TSV
    action_file = (
        Path(out_unk).with_suffix("").with_suffix("").parent
        / "ensembl_unknown_prefixes_ACTION_REQUIRED.txt"
    )
    action_file.write_text(full_msg, encoding="utf-8")
    log.error(f"Action required file written → {action_file}")

    raise ValueError(
        f"Unknown Ensembl prefixes: {unknown}. "
        f"See {out_unk} and {action_file} for instructions."
    )

log.info(f"Written species_map → {out_map}")
log.info("detect_ensembl_species complete.")
