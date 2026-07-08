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

# Curated N-glycosylation enzyme catalogs (manually curated reference sets).
# A = initiation enzymes: dolichol assembly, OST transfer, ER trimming → immature/high-mannose.
# B = maturation enzymes: Golgi processing and terminal capping → mature glycoforms.
GLYCO_A_INITIATION = [
    "DPAGT1", "ALG1", "ALG2", "ALG3", "ALG5", "ALG6", "ALG8", "ALG9", "ALG10", "ALG10B",
    "ALG11", "ALG12", "ALG13", "ALG14", "DPM1", "DPM2", "DPM3", "MPDU1", "DOLK", "SRD5A3",
    "RFT1", "STT3A", "STT3B", "RPN1", "RPN2", "DDOST", "DAD1", "OSTC", "TUSC3", "MAGT1",
    "OST4", "KRTCAP2", "MOGS", "GANAB", "PRKCSH", "MAN1B1",
]

GLYCO_B_MATURATION = [
    "MAN1A1", "MAN1A2", "MAN1C1", "MAN2A1", "MAN2A2",
    "MGAT1", "MGAT2", "MGAT3", "MGAT4A", "MGAT4B", "MGAT4C", "MGAT5", "MGAT5B",
    "B4GALT1", "B4GALT2", "B4GALT3", "B4GALT4", "B4GALT5", "B4GALT6", "B4GALT7",
    "ST3GAL1", "ST3GAL2", "ST3GAL3", "ST3GAL4", "ST3GAL5", "ST3GAL6",
    "ST6GAL1", "ST6GAL2",
    "ST6GALNAC1", "ST6GALNAC2", "ST6GALNAC3", "ST6GALNAC4", "ST6GALNAC6",
    "ST8SIA1", "ST8SIA2", "ST8SIA4",
    "FUT1", "FUT2", "FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT8", "FUT9", "FUT10", "FUT11",
    "POFUT1", "POFUT2",
    "GALNT1", "GALNT2", "GALNT3", "GALNT4", "GALNT6", "GALNT7", "GALNT10", "GALNT12",
    "C1GALT1", "C1GALT1C1", "GCNT1", "GCNT3",
]


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
    converter_name: str,
    threshold: float,
    set_type: str,
    tissue: str,
) -> int:
    set_dir.mkdir(parents=True, exist_ok=True)

    rows = [
        (g, gene_tstat[g][tissue_idx] if g in gene_tstat else None)
        for g in genes
    ]

    with (set_dir / "geneset.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\tgtex_tstat\ttissue_status\n")
        for g, t in rows:
            status = _tissue_status(t) if t is not None else "no_data"
            tval = f"{t:.6f}" if t is not None else "NA"
            fh.write(f"{g}\t{tval}\t{status}\n")

    gene_symbols = [g for g, _ in rows]
    write_gmt([(set_name, gene_symbols)], set_dir / "genesets.gmt")

    file_rec = input_file_record(str(tstat_path), "gtex_tstat_tsv")
    meta = make_metadata(
        converter_name=converter_name,
        parameters={
            "set_type": set_type,
            "tissue": tissue,
            "expression_threshold": threshold,
            "gtex_version": "v8",
            "gtex_source_url": GTEX_V8_MEDIAN_TPM_URL,
            "catalog": "curated_nglyco_enzymes",
        },
        data_type="transcriptomics",
        assay="bulk",
        organism="human",
        genome_build="GRCh38",
        files=[file_rec],
        gene_annotation={"mode": "curated", "source": "manual_nglyco_catalog"},
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
    converter_name = "glycomaturity_gtex"

    for ti, tissue in enumerate(tissues):
        safe_tissue = _safe_name(tissue)

        n = _write_catalog_set(
            set_dir=out_dir / f"GlycoA_initiation_{safe_tissue}",
            set_name=f"GlycoA_initiation_{safe_tissue}",
            description=(
                f"N-glycosylation initiation enzymes (dolichol assembly → OST transfer → ER trimming; "
                f"unglyco→immature/high-mannose) expression-annotated in {tissue} (GTEx V8). "
                f"Curated reference catalog; absence=null (low t-stat is not negative evidence)."
            ),
            genes=GLYCO_A_INITIATION,
            gene_tstat=gene_tstat,
            tissue_idx=ti,
            tstat_path=tstat_path,
            converter_name=converter_name,
            threshold=threshold,
            set_type="A_initiation_immature",
            tissue=tissue,
        )
        gene_counts.append(n)
        n_sets += 1

        n = _write_catalog_set(
            set_dir=out_dir / f"GlycoB_maturation_{safe_tissue}",
            set_name=f"GlycoB_maturation_{safe_tissue}",
            description=(
                f"N-glycosylation maturation enzymes (Golgi processing + terminal capping; "
                f"immature→mature glycoforms) expression-annotated in {tissue} (GTEx V8). "
                f"Curated reference catalog; absence=null."
            ),
            genes=GLYCO_B_MATURATION,
            gene_tstat=gene_tstat,
            tissue_idx=ti,
            tstat_path=tstat_path,
            converter_name=converter_name,
            threshold=threshold,
            set_type="B_maturation_mature",
            tissue=tissue,
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
