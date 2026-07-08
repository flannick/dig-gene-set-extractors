from __future__ import annotations

import csv
import urllib.request
from pathlib import Path

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

GTEX_V8_MEDIAN_TPM_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
)

# ── N-glycosylation: Set A ────────────────────────────────────────────────────
# Initiation enzymes: dolichol-linked oligosaccharide assembly (ER), OST transfer,
# and early trimming → high-mannose / immature glycoforms.
# Sources: Stanley P et al. N-Glycans. In: Varki A et al. (eds.) Essentials of Glycobiology
# 4th ed. CSHL Press 2022 (Ch. 9); UniProtKB pathway annotations.
GLYCO_A_INITIATION = [
    "DPAGT1", "ALG1", "ALG2", "ALG3", "ALG5", "ALG6", "ALG8", "ALG9", "ALG10", "ALG10B",
    "ALG11", "ALG12", "ALG13", "ALG14", "DPM1", "DPM2", "DPM3", "MPDU1", "DOLK", "SRD5A3",
    "RFT1", "STT3A", "STT3B", "RPN1", "RPN2", "DDOST", "DAD1", "OSTC", "TUSC3", "MAGT1",
    "OST4", "KRTCAP2", "MOGS", "GANAB", "PRKCSH", "MAN1B1",
]

# ── N-glycosylation: Set B — mature (complex/hybrid) ─────────────────────────
# Golgi maturation enzymes that elaborate the high-mannose N-glycan precursor:
# mannosidases trim Man8/9→Man5; MGAT enzymes add GlcNAc branches; galactosyl-,
# sialyl-, and fucosyltransferases build the antennae. Also includes O-fucosylation
# (POFUT1/2) on EGF/TSP1-domain proteins (Notch pathway).
# Mucin-type O-GalNAc enzymes (GALNT family, C1GALT1, GCNT1/3) are in Set E.
# Sources:
#   Stanley P et al. N-Glycans. In: Varki A et al. (eds.) Essentials of Glycobiology
#     4th ed. CSHL Press 2022 (Ch. 9).
#   Haltiwanger RS & Lowe JB. Role of glycosylation in development.
#     Annu Rev Biochem. 2004;73:491-537. doi:10.1146/annurev.biochem.73.011303.074043
GLYCO_B_MATURATION = [
    "MAN1A1", "MAN1A2", "MAN1C1", "MAN2A1", "MAN2A2",
    "MGAT1", "MGAT2", "MGAT3", "MGAT4A", "MGAT4B", "MGAT4C", "MGAT5", "MGAT5B",
    "B4GALT1", "B4GALT2", "B4GALT3", "B4GALT4", "B4GALT5", "B4GALT6",
    "ST3GAL1", "ST3GAL2", "ST3GAL3", "ST3GAL4", "ST3GAL5", "ST3GAL6",
    "ST6GAL1", "ST6GAL2",
    "ST8SIA1", "ST8SIA2", "ST8SIA4",
    "FUT1", "FUT2", "FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT8", "FUT9", "FUT10", "FUT11",
    "POFUT1", "POFUT2",
]

# ── O-glycosylation: Set C ────────────────────────────────────────────────────
# O-GlcNAc cycling enzymes: nucleocytoplasmic O-GlcNAc addition (OGT), removal (OGA/MGEA5),
# and EGF-domain-specific O-GlcNAc transferase (EOGT).
# This is a distinct, non-Golgi O-glycosylation pathway that dynamically modifies Ser/Thr
# in response to nutrient/metabolic state (UDP-GlcNAc availability).
# Sources:
#   Hart GW et al. Cross talk between O-GlcNAcylation and phosphorylation.
#     Annu Rev Biochem. 2011;80:825-58. doi:10.1146/annurev-biochem-060608-102511
#   OGlcNAc Atlas (https://oglcnac.mcw.edu/) — MCW/NIH curated database.
#   OGlcNAc Atlas (https://www.oglcnacatlas.com/) — Yang/Qian lab.
#   GlyGen O-GlcNAc data: https://data.glygen.org/ln2data/releases/data/current/reviewed/
#     human_protein_xref_oglcnac_atlas.csv, human_protein_xref_oglcnac_mcw.csv
GLYCO_C_OGLCNAC = [
    "OGT",    # O-GlcNAc transferase (adds GlcNAc to Ser/Thr; cytoplasmic/nuclear)
    "MGEA5",  # O-GlcNAcase / OGA (removes O-GlcNAc; also known as OGA)
    "EOGT",   # EGF-domain-specific O-GlcNAc transferase (ER-lumenal, extracellular proteins)
]

# ── O-glycosylation: Set D ────────────────────────────────────────────────────
# Proteoglycan O-xylosylation pathway: initiates heparan sulfate (HS) and
# chondroitin/dermatan sulfate (CS/DS) chains on proteoglycan core proteins.
# XYLT1/2 add the first xylose to Ser; subsequent enzymes build the tetrasaccharide
# linker (GalI-GalII-GlcA), then EXT/EXTL polymerases extend the HS chain;
# CSGALNACT1/2 initiate CS chains.
# Sources:
#   Esko JD, Kimata K, Lindahl U. Proteoglycans and Sulfated Glycosaminoglycans.
#     In: Varki A et al. Essentials of Glycobiology 4th ed. CSHL Press 2022 (Ch. 17).
#   Presto J et al. Heparan sulfate biosynthesis enzymes EXT1 and EXT2 affect NDST1
#     expression and heparan sulfate sulfation. J Biol Chem. 2008;283:16983-91.
#   Kitagawa H et al. Chondroitin sulfate biosynthesis.
#     Curr Opin Struct Biol. 2008;18:597-603.
GLYCO_D_PROTEOGLYCAN = [
    # Linker tetrasaccharide synthesis (Ser-Xyl-Gal-Gal-GlcA):
    "XYLT1", "XYLT2",       # xylosyltransferases I and II (Ser-Xyl; initiation)
    "B4GALT7",              # galactosyltransferase I (Xyl-GalI)
    "B3GALT6",              # galactosyltransferase II (GalI-GalII)
    "B4GAT1",               # glucuronyltransferase (GalII-GlcA; linker completion)
    # Heparan sulfate chain polymerization:
    "EXT1", "EXT2",         # HS co-polymerases (GlcA-GlcNAc alternating units)
    "EXTL1", "EXTL2", "EXTL3",  # EXT-like HS transferases
    # Chondroitin/dermatan sulfate initiation:
    "CSGALNACT1", "CSGALNACT2",  # CS GalNAc-transferases I and II (GlcA-GalNAc)
]

# ── O-glycosylation: Set E — mucin-type O-GalNAc ─────────────────────────────
# Full polypeptide GalNAc-transferase (ppGalNAcT / GALNT) family: add α-GalNAc
# to Ser/Thr in the Golgi, producing the Tn antigen (immature O-GalNAc).
# Core 1-4 elaboration enzymes convert Tn → T antigen → sialylated / extended
# mucin O-glycans. The complete GALNT family is included so GTEx t-stat reveals
# which initiators dominate in each tissue (many are tissue-restricted).
# Sources:
#   Schjoldager KT et al. Global view of human protein glycosylation pathways and
#     functions. Nat Struct Mol Biol. 2020;27(4):303-310. doi:10.1038/s41594-020-0371-0
#   Brockhausen I, Wandall HH, Ten Hagen KG. O-GalNAc glycans.
#     In: Varki A et al. Essentials of Glycobiology 4th ed. 2022. Chapter 10.
#   CAZy family GT27: https://www.cazy.org/GT27.html
GLYCO_E_MUCIN_INITIATION = [
    # Initiation (Tn antigen — immature O-GalNAc):
    "GALNT1",  "GALNT2",  "GALNT3",  "GALNT4",  "GALNT5",
    "GALNT6",  "GALNT7",  "GALNT8",  "GALNT9",  "GALNT10",
    "GALNT11", "GALNT12", "GALNT13", "GALNT14", "GALNT15",
    "GALNT16", "GALNT17", "GALNT18", "GALNT20",
    # Core 1 synthesis and elaboration (Tn → T antigen → sialylated/branched):
    "C1GALT1",    # T-synthase: GalNAc-α-Ser/Thr → Galβ1-3GalNAc (core 1 / T antigen)
    "C1GALT1C1",  # COSMC: ER chaperone required for C1GALT1 folding
    "GCNT1",      # Core 2 GlcNAc-T: branches core 1 → core 2 (lymphoid/secretory)
    "GCNT3",      # Core 2/4 GlcNAc-T: broader tissue expression
    "B3GNT6",     # Core 3 GlcNAc-T: GalNAc-α-Ser/Thr → GlcNAcβ1-3GalNAc (core 3)
    # Sialylation of O-GalNAc glycans:
    "ST3GAL1", "ST3GAL2",  # Siaα2-3Gal on core 1 and 2 (T antigen sialylation)
    "ST6GALNAC1", "ST6GALNAC2", "ST6GALNAC3", "ST6GALNAC4", "ST6GALNAC6",
]

# ── O-glycosylation: Set F — O-mannose / alpha-dystroglycan pathway ───────────
# Multi-step O-mannosylation of alpha-dystroglycan (α-DG), the receptor for
# laminin and other ECM ligands. POMT1/2 add O-Man in the ER; downstream enzymes
# build an unusual matriglycan chain via CDP-ribitol intermediates.
# Defects cause alpha-dystroglycanopathies: Walker-Warburg syndrome, Fukuyama
# congenital muscular dystrophy, LGMD2I, and related disorders.
# ISPD and GMPPB synthesize the activated substrates.
# Sources:
#   Yoshida-Moriguchi T & Campbell KP. Matriglycan: a polymeric glycan that links
#     dystroglycan to the basement membrane. Glycobiology. 2015;25(10):1039-52.
#     doi:10.1093/glycob/cwv066
#   Kanagawa M et al. Identification of a Post-translational Modification with
#     Ribitol-Phosphate and Its Defect in Muscular Dystrophy. Cell.
#     2016;167(5):1339-1353. doi:10.1016/j.cell.2016.09.003
GLYCO_F_O_MANNOSE = [
    "POMT1", "POMT2",   # O-mannosyltransferase 1/2 (ER; Ser/Thr-O-Man; initiation)
    "POMGNT1",          # POMGnT1: Man-α-Ser/Thr → GlcNAcβ1-2Man (core M1)
    "POMGNT2",          # POMGnT2: Man-α-Ser/Thr → GlcNAcβ1-4Man (core M3, α-DG specific)
    "B3GALNT2",         # β-1,3-GalNAc-T2: GlcNAcβ-Man → GalNAcβ1-3GlcNAc (core M3)
    "POMK",             # Phosphomannose kinase: phosphorylates 6-O-Man
    "FKTN",             # Fukutin: first CDP-ribitol phosphotransferase
    "FKRP",             # Fukutin-related protein: second CDP-ribitol phosphotransferase
    "TMEM5",            # RXYLT1: xylosyltransferase on the ribitol-P anchor
    "B4GAT1",           # β-1,4-GlcA-T: glucuronosyltransferase (shared with Set D linker)
    "LARGE1",           # Bifunctional Xyl-T/GlcA-T: extends the matriglycan repeats
    "LARGE2",           # GYLTL1B: homolog of LARGE1
    "ISPD",             # CDP-ribitol pyrophosphorylase A (activated substrate synthesis)
    "GMPPB",            # GDP-mannose pyrophosphorylase B (GDP-Man substrate)
]

# Catalog metadata for provenance records
CATALOG_SOURCES = {
    "A_initiation_immature": {
        "glycosylation_type": "N-linked",
        "pathway_stage": "initiation",
        "glycoform": "high-mannose / immature",
        "catalog_reference": (
            "Stanley P et al. N-Glycans. In: Varki A et al. (eds.) Essentials of Glycobiology "
            "4th ed. Cold Spring Harbor Laboratory Press, 2022. Chapter 9."
        ),
    },
    "B_maturation_mature": {
        "glycosylation_type": "N-linked (+ O-fucosylation)",
        "pathway_stage": "maturation (Golgi GlcNAc-branching, galactosylation, sialylation, fucosylation)",
        "glycoform": "complex / hybrid / mature; O-fucose on EGF/TSP1 domains",
        "catalog_reference": (
            "Stanley P et al. N-Glycans. In: Varki A et al. (eds.) Essentials of Glycobiology "
            "4th ed. Cold Spring Harbor Laboratory Press, 2022. Chapter 9. "
            "Haltiwanger RS & Lowe JB. Role of glycosylation in development. "
            "Annu Rev Biochem. 2004;73:491-537. doi:10.1146/annurev.biochem.73.011303.074043"
        ),
    },
    "C_O_GlcNAc": {
        "glycosylation_type": "O-linked (O-GlcNAc)",
        "pathway_stage": "nucleocytoplasmic O-GlcNAc cycling",
        "glycoform": "O-GlcNAc (Ser/Thr)",
        "catalog_reference": (
            "Hart GW et al. Cross talk between O-GlcNAcylation and phosphorylation: roles in "
            "signaling, transcription, and chronic disease. Annu Rev Biochem. "
            "2011;80:825-58. doi:10.1146/annurev-biochem-060608-102511. "
            "OGlcNAc Atlas (https://oglcnac.mcw.edu/); "
            "GlyGen O-GlcNAc data "
            "(https://data.glygen.org/ln2data/releases/data/current/reviewed/"
            "human_protein_xref_oglcnac_atlas.csv)."
        ),
    },
    "D_proteoglycan_Oxylosylation": {
        "glycosylation_type": "O-linked (O-xylosylation / proteoglycan GAG chains)",
        "pathway_stage": "proteoglycan GAG chain initiation (linker) and extension (HS/CS/DS)",
        "glycoform": "heparan sulfate / chondroitin sulfate / dermatan sulfate",
        "catalog_reference": (
            "Esko JD, Kimata K, Lindahl U. Proteoglycans and Sulfated Glycosaminoglycans. "
            "In: Varki A et al. (eds.) Essentials of Glycobiology 4th ed. "
            "Cold Spring Harbor Laboratory Press, 2022. Chapter 17. "
            "Presto J et al. J Biol Chem. 2008;283(28):16983-91. doi:10.1074/jbc.M801451200. "
            "Kitagawa H et al. Chondroitin sulfate biosynthesis. "
            "Curr Opin Struct Biol. 2008;18(5):597-603. doi:10.1016/j.sbi.2008.09.002"
        ),
    },
    "E_mucin_OGalNAc_initiation": {
        "glycosylation_type": "O-linked (mucin-type O-GalNAc)",
        "pathway_stage": (
            "O-GalNAc initiation (GALNT family; Tn antigen) through core 1-4 elaboration "
            "(C1GALT1/COSMC, GCNT1/3, B3GNT6) and mucin sialylation (ST3GAL1/2, ST6GALNAC1-4/6)"
        ),
        "glycoform": "Tn antigen → T antigen → core 1-4 / sialylated O-GalNAc (mucins)",
        "catalog_reference": (
            "Schjoldager KT et al. Global view of human protein glycosylation pathways and "
            "functions. Nat Struct Mol Biol. 2020;27(4):303-310. doi:10.1038/s41594-020-0371-0. "
            "Brockhausen I, Wandall HH, Ten Hagen KG. O-GalNAc glycans. "
            "In: Varki A et al. Essentials of Glycobiology 4th ed. 2022. Chapter 10. "
            "CAZy family GT27: https://www.cazy.org/GT27.html"
        ),
    },
    "F_O_mannose_dystroglycan": {
        "glycosylation_type": "O-linked (O-mannose / alpha-dystroglycan matriglycan)",
        "pathway_stage": (
            "O-mannosylation initiation (POMT1/2; ER) through matriglycan polymerization "
            "(FKTN/FKRP CDP-ribitol transfer, TMEM5, LARGE1/2 Xyl-GlcA extension) "
            "and activated-sugar synthesis (ISPD, GMPPB)"
        ),
        "glycoform": "O-mannose → matriglycan [(Xylα1-3GlcAβ1-3)n]; laminin-binding epitope on α-DG",
        "catalog_reference": (
            "Yoshida-Moriguchi T & Campbell KP. Matriglycan: a polymeric glycan that links "
            "dystroglycan to the basement membrane. Glycobiology. 2015;25(10):1039-52. "
            "doi:10.1093/glycob/cwv066. "
            "Kanagawa M et al. Identification of a Post-translational Modification with "
            "Ribitol-Phosphate and Its Defect in Muscular Dystrophy. Cell. "
            "2016;167(5):1339-1353. doi:10.1016/j.cell.2016.09.003"
        ),
    },
}


def _safe_name(tissue: str, max_len: int = 60) -> str:
    return "".join(c if c.isalnum() else "_" for c in tissue)[:max_len]


def _maybe_download(path_or_url: str, out_dir: Path) -> Path:
    if path_or_url.startswith(("http://", "https://")):
        dl_dir = out_dir / "downloads"
        dl_dir.mkdir(parents=True, exist_ok=True)
        dest = dl_dir / Path(path_or_url.rstrip("/").split("/")[-1])
        if not dest.exists():
            urllib.request.urlretrieve(path_or_url, dest)
        return dest
    return Path(path_or_url)


def _load_tstat_matrix(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        tissues = header[1:]
        gene_tstat: dict[str, list[float]] = {}
        for row in reader:
            if not row or not row[0]:
                continue
            try:
                gene_tstat[row[0]] = [float(x) for x in row[1:]]
            except ValueError:
                continue
    return tissues, gene_tstat


def _tissue_status(tstat: float) -> str:
    if tstat >= 2.0:
        return "enriched"
    if tstat >= -2.0:
        return "present"
    return "low_null"


def _write_catalog_set(
    set_dir: Path,
    set_name: str,
    description: str,
    genes: list[str],
    gene_tstat: dict[str, list[float]],
    tissue_idx: int,
    tstat_path: Path,
    set_type: str,
    tissue: str,
    threshold: float,
) -> int:
    set_dir.mkdir(parents=True, exist_ok=True)

    rows = [(g, gene_tstat[g][tissue_idx] if g in gene_tstat else None) for g in genes]

    with (set_dir / "geneset.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\tgtex_tstat\ttissue_status\n")
        for g, t in rows:
            status = _tissue_status(t) if t is not None else "no_data"
            tval = f"{t:.6f}" if t is not None else "NA"
            fh.write(f"{g}\t{tval}\t{status}\n")

    gene_symbols = [g for g, _ in rows]
    write_gmt([(set_name, gene_symbols)], set_dir / "genesets.gmt")

    catalog_meta = CATALOG_SOURCES.get(set_type, {})
    file_rec = input_file_record(str(tstat_path), "gtex_tstat_tsv")
    meta = make_metadata(
        converter_name="glycomaturity_gtex",
        parameters={
            "set_type": set_type,
            "tissue": tissue,
            "expression_threshold": threshold,
            "gtex_version": "v8",
            "gtex_source_url": GTEX_V8_MEDIAN_TPM_URL,
            "glycosylation_type": catalog_meta.get("glycosylation_type", ""),
            "pathway_stage": catalog_meta.get("pathway_stage", ""),
            "glycoform": catalog_meta.get("glycoform", ""),
            "catalog_reference": catalog_meta.get("catalog_reference", ""),
            "catalog": "curated_glyco_enzymes",
            "absence_interpretation": "null (low t-stat is not negative evidence)",
        },
        data_type="transcriptomics",
        assay="bulk",
        organism="human",
        genome_build="GRCh38",
        files=[file_rec],
        gene_annotation={
            "mode": "curated",
            "source": "manual_glyco_enzyme_catalog",
            "citation": catalog_meta.get("catalog_reference", ""),
        },
        weights={
            "weight_type": "nonnegative",
            "normalization": {"method": "none", "target_sum": None},
            "aggregation": "catalog_membership",
        },
        summary={
            "n_input_features": len(genes),
            "n_genes": len(gene_symbols),
            "n_features_assigned": len(gene_symbols),
            "fraction_features_assigned": 1.0,
            "n_sets_emitted": 1,
        },
        gene_set_description=description,
        output_files=[
            {"path": "genesets.gmt", "role": "gmt_library"},
            {"path": "geneset.tsv", "role": "selected_program"},
            {"path": "geneset.meta.json", "role": "metadata_json"},
        ],
    )
    write_metadata(set_dir / "geneset.meta.json", meta)
    return len(gene_symbols)


def run(args) -> dict[str, object]:
    activate_runtime_context("glycomaturity_gtex", getattr(args, "provenance_overlay_json", None))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tstat_path = _maybe_download(args.gtex_tstat_tsv, out_dir)
    tissues, gene_tstat = _load_tstat_matrix(tstat_path)
    threshold = float(getattr(args, "expression_threshold", 2.0))

    n_sets = 0
    gene_counts: list[int] = []

    sets = [
        (
            "GlycoA_initiation",
            "A_initiation_immature",
            GLYCO_A_INITIATION,
            lambda t, st: (
                f"N-glycosylation initiation enzymes (dolichol assembly → OST transfer → ER "
                f"trimming; unglyco→immature/high-mannose) expression-annotated in {t} (GTEx V8). "
                f"Curated reference catalog; absence=null."
            ),
        ),
        (
            "GlycoB_maturation",
            "B_maturation_mature",
            GLYCO_B_MATURATION,
            lambda t, st: (
                f"N-glycosylation maturation enzymes (Golgi mannosidases, GlcNAc-branching "
                f"MGAT1-5, galactosylation B4GALT1-6, sialylation ST3/6GAL, polysialylation "
                f"ST8SIA, fucosylation FUT1-11; elaborates immature→complex/mature N-glycoforms) "
                f"expression-annotated in {t} (GTEx V8 median TPM). "
                f"Includes O-fucosylation (POFUT1/2; EGF/TSP1-domain Notch pathway). "
                f"Mucin O-GalNAc enzymes are in Set E. Curated reference catalog; absence=null."
            ),
        ),
        (
            "GlycoC_OGlcNAc",
            "C_O_GlcNAc",
            GLYCO_C_OGLCNAC,
            lambda t, st: (
                f"O-GlcNAc cycling enzymes (OGT adds, OGA/MGEA5 removes nucleocytoplasmic "
                f"O-GlcNAc on Ser/Thr; EOGT for EGF-domain substrates) expression-annotated "
                f"in {t} (GTEx V8). "
                f"Source: Hart et al. Annu Rev Biochem 2011; OGlcNAc Atlas (MCW/NIH); "
                f"GlyGen human_protein_xref_oglcnac_atlas.csv."
            ),
        ),
        (
            "GlycoD_proteoglycan",
            "D_proteoglycan_Oxylosylation",
            GLYCO_D_PROTEOGLYCAN,
            lambda t, st: (
                f"Proteoglycan O-xylosylation pathway enzymes (Ser-Xyl linker synthesis via "
                f"XYLT1/2; tetrasaccharide completion via B4GALT7/B3GALT6/B4GAT1; heparan "
                f"sulfate extension via EXT1/2/EXTL1-3; chondroitin sulfate initiation via "
                f"CSGALNACT1/2) expression-annotated in {t} (GTEx V8). "
                f"Source: Esko et al. Essentials of Glycobiology 4th ed. Ch. 17."
            ),
        ),
        (
            "GlycoE_mucin_OGalNAc",
            "E_mucin_OGalNAc_initiation",
            GLYCO_E_MUCIN_INITIATION,
            lambda t, st: (
                f"Mucin-type O-GalNAc glycosylation enzymes: full GALNT family (GALNT1-20; "
                f"ppGalNAcTs add the immature Tn antigen α-GalNAc to Ser/Thr) and core "
                f"elaboration enzymes (C1GALT1/COSMC → T antigen; GCNT1/3 → core 2/4; "
                f"B3GNT6 → core 3; ST3GAL1/2, ST6GALNAC1-4/6 → sialylated mucin O-glycans) "
                f"expression-annotated in {t} (GTEx V8 median TPM). "
                f"Source: Schjoldager et al. Nat Struct Mol Biol 2020; "
                f"Essentials of Glycobiology 4th ed. Ch. 10; CAZy GT27."
            ),
        ),
        (
            "GlycoF_Omannose_dystroglycan",
            "F_O_mannose_dystroglycan",
            GLYCO_F_O_MANNOSE,
            lambda t, st: (
                f"O-mannose / alpha-dystroglycan (α-DG) matriglycan pathway: POMT1/2 (ER "
                f"O-Man initiation), POMGNT1/2, B3GALNT2, POMK, FKTN/FKRP (CDP-ribitol "
                f"transfer), TMEM5/RXYLT1, B4GAT1, LARGE1/2 (Xyl-GlcA matriglycan extension), "
                f"ISPD/GMPPB (substrate synthesis) expression-annotated in {t} (GTEx V8 "
                f"median TPM). Defects cause alpha-dystroglycanopathies (Walker-Warburg "
                f"syndrome, Fukuyama CMD, LGMD2I). "
                f"Source: Yoshida-Moriguchi & Campbell, Glycobiology 2015; "
                f"Kanagawa et al. Cell 2016."
            ),
        ),
    ]

    for ti, tissue in enumerate(tissues):
        safe_tissue = _safe_name(tissue)
        for prefix, set_type, genes, desc_fn in sets:
            set_name = f"{prefix}_{safe_tissue}"
            n = _write_catalog_set(
                set_dir=out_dir / set_name,
                set_name=set_name,
                description=desc_fn(tissue, set_type),
                genes=genes,
                gene_tstat=gene_tstat,
                tissue_idx=ti,
                tstat_path=tstat_path,
                set_type=set_type,
                tissue=tissue,
                threshold=threshold,
            )
            gene_counts.append(n)
            n_sets += 1

    return {
        "n_peaks": sum(gene_counts),
        "n_genes": sum(gene_counts),
        "n_groups": n_sets,
        "n_genes_min": min(gene_counts) if gene_counts else 0,
        "n_genes_max": max(gene_counts) if gene_counts else 0,
        "out_dir": str(out_dir),
    }
