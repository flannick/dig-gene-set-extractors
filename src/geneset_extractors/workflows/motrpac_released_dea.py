from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from geneset_extractors.workflows.gtex_runtime_common import write_tsv, write_workflow_provenance_graph


LEGACY_TISSUE_TERMS: dict[str, str] = {
    "BLOOD": "T30-Blood-RNA",
    "BLOOD-RNA": "T30-Blood-RNA",
    "HIPPOC": "T52-Hippocampus",
    "HIPPOCAMPUS": "T52-Hippocampus",
    "CORTEX": "T53-Cortex",
    "HYPOTH": "T54-Hypothalamus",
    "HYPOTHALAMUS": "T54-Hypothalamus",
    "SKM-GN": "T55-Gastrocnemius",
    "SKMGN": "T55-Gastrocnemius",
    "GASTROCNEMIUS": "T55-Gastrocnemius",
    "SKM-VL": "T56-Vastus-Lateralis",
    "SKMVL": "T56-Vastus-Lateralis",
    "VASTUS-LATERALIS": "T56-Vastus-Lateralis",
    "HEART": "T58-Heart",
    "KIDNEY": "T59-Kidney",
    "ADRNL": "T60-Adrenal",
    "ADRENAL": "T60-Adrenal",
    "COLON": "T61-Colon",
    "SPLEEN": "T62-Spleen",
    "TESTES": "T63-Testes",
    "OVARY": "T64-Ovaries",
    "OVARIES": "T64-Ovaries",
    "LUNG": "T66-Lung",
    "SMLINT": "T67-Small-Intestine",
    "SMALL-INTESTINE": "T67-Small-Intestine",
    "SMALL_INTESTINE": "T67-Small-Intestine",
    "LIVER": "T68-Liver",
    "BAT": "T69-Brown-Adipose",
    "BROWN-ADIPOSE": "T69-Brown-Adipose",
    "WAT-SC": "T70-White-Adipose",
    "WATSC": "T70-White-Adipose",
    "WHITE-ADIPOSE": "T70-White-Adipose",
    "VENACV": "T99-Vena-Cava",
    "VENA-CAVA": "T99-Vena-Cava",
    "VENA_CAVA": "T99-Vena-Cava",
}


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _require_dir(path: Path, label: str) -> None:
    if not path.is_dir():
        raise NotADirectoryError(f"Missing {label}: {path}")


def _threshold(score: float) -> int | None:
    if score > 0:
        return 1
    if score < 0:
        return -1
    return None


def _load_feature_annotation(feature_annot: Path) -> pd.DataFrame:
    _require_file(feature_annot, "TRNSCRPT_FEATURE_ANNOT.txt")
    return pd.read_csv(feature_annot, sep="\t")


def _build_gene_mapper(feature_df: pd.DataFrame) -> dict[str, str]:
    required = {"gene_id", "gene_name"}
    missing = required - set(feature_df.columns)
    if missing:
        raise ValueError(f"Feature annotation is missing required columns: {sorted(missing)}")

    gene_mapper: dict[str, str] = {}
    for _, row in feature_df.iterrows():
        gene_mapper[row["gene_id"]] = row["gene_name"]
    return gene_mapper


def _load_symbol_map(mapping_file: Path) -> dict[str, str]:
    _require_file(mapping_file, "Harmonizome mapping file")
    df = pd.read_csv(mapping_file, sep="\t", header=None, dtype=str, keep_default_na=False)

    if df.shape[1] >= 3:
        input_col, output_col = 1, 2
    elif df.shape[1] == 2:
        input_col, output_col = 0, 1
    else:
        raise ValueError(
            "Mapping file must have either two tab-delimited columns "
            "(input_symbol, approved_symbol) or at least three columns "
            "where columns 2 and 3 are input_symbol and approved_symbol."
        )

    df[input_col] = df[input_col].astype(str).str.strip()
    df[output_col] = df[output_col].astype(str).str.strip()
    df = df[(df[input_col] != "") & (df[output_col] != "")]
    return dict(zip(df[input_col], df[output_col]))


def _load_timewise_dea(
    dea_dir: Path,
    gene_mapper: dict[str, str],
    symbol_map: dict[str, str],
    padj_max: float,
) -> tuple[pd.DataFrame, dict[str, int], list[Path]]:
    _require_dir(dea_dir, "MoTrPAC DEA directory")

    motrpac = pd.DataFrame([])
    dea_paths: list[Path] = []
    for rnaseqfile in sorted(os.listdir(dea_dir)):
        if "timewise" in rnaseqfile:
            rnaseq_path = dea_dir / rnaseqfile
            dea_paths.append(rnaseq_path)
            rnaseq = pd.read_csv(rnaseq_path, sep="\t").get(
                ["feature_ID", "tissue", "sex", "comparison_group", "adj_p_value", "logFC"]
            )
            if rnaseq is None:
                raise ValueError(f"{rnaseq_path} is missing one or more required timewise columns")
            rnaseq["term"] = rnaseq["tissue"] + "_" + rnaseq["sex"] + "_" + rnaseq["comparison_group"]
            motrpac = pd.concat([motrpac, rnaseq])

    if motrpac.empty:
        raise ValueError(f"No DEA files containing 'timewise' were found in {dea_dir}")

    before_gene_mapper = len(motrpac)
    motrpac["feature_ID"] = motrpac["feature_ID"].map(gene_mapper)
    motrpac = motrpac.dropna()

    before_padj = len(motrpac)
    motrpac = motrpac[["term", "feature_ID", "adj_p_value", "logFC"]]
    motrpac = motrpac[motrpac["adj_p_value"] < padj_max]

    before_symbolmap = len(motrpac)
    motrpac["feature_ID"] = motrpac["feature_ID"].apply(str.upper).map(symbol_map)
    motrpac = motrpac.dropna()
    motrpac["logFC"] = motrpac["logFC"].apply(_threshold)
    motrpac.columns = ["term", "gene", "adj_p_value", "threshold"]

    audit = {
        "timewise_rows_loaded": before_gene_mapper,
        "timewise_rows_after_gene_mapper": before_padj,
        "timewise_rows_after_padj": before_symbolmap,
        "timewise_rows_after_symbol_map": len(motrpac),
    }
    return motrpac, audit, dea_paths


def _load_training_dea(
    dea_dir: Path,
    gene_mapper: dict[str, str],
    symbol_map: dict[str, str],
    padj_max: float,
) -> tuple[pd.DataFrame, dict[str, int], list[Path]]:
    _require_dir(dea_dir, "MoTrPAC DEA directory")

    motrpac_training = pd.DataFrame([])
    dea_paths: list[Path] = []
    for rnaseqfile in sorted(os.listdir(dea_dir)):
        if "training" in rnaseqfile:
            rnaseq_path = dea_dir / rnaseqfile
            dea_paths.append(rnaseq_path)
            rnaseq = pd.read_csv(rnaseq_path, sep="\t").get(["feature_ID", "tissue", "adj_p_value"])
            if rnaseq is None:
                raise ValueError(f"{rnaseq_path} is missing one or more required training columns")
            motrpac_training = pd.concat([motrpac_training, rnaseq])

    if motrpac_training.empty:
        raise ValueError(f"No DEA files containing 'training' were found in {dea_dir}")

    before_gene_mapper = len(motrpac_training)
    motrpac_training["feature_ID"] = motrpac_training["feature_ID"].map(gene_mapper)
    motrpac_training = motrpac_training.dropna()

    before_padj = len(motrpac_training)
    motrpac_training = motrpac_training[["tissue", "feature_ID", "adj_p_value"]]
    motrpac_training = motrpac_training[motrpac_training["adj_p_value"] < padj_max]

    before_symbolmap = len(motrpac_training)
    motrpac_training["feature_ID"] = motrpac_training["feature_ID"].apply(str.upper).map(symbol_map)
    motrpac_training = motrpac_training.dropna()
    motrpac_training["threshold"] = 1

    motrpac_training["tissue"] = motrpac_training["tissue"] + "_consensus"
    motrpac_training.columns = ["term", "gene", "adj_p_value", "threshold"]

    audit = {
        "training_rows_loaded": before_gene_mapper,
        "training_rows_after_gene_mapper": before_padj,
        "training_rows_after_padj": before_symbolmap,
        "training_rows_after_symbol_map": len(motrpac_training),
    }
    return motrpac_training, audit, dea_paths


def _combine_and_standardize(motrpac: pd.DataFrame, motrpac_training: pd.DataFrame) -> pd.DataFrame:
    combined = pd.concat([motrpac, motrpac_training]).reset_index(drop=True)
    combined["adj_p_value"] = combined["adj_p_value"].apply(lambda x: np.log10(x) * -1)
    combined["adj_p_value"] = combined["adj_p_value"].mul(combined["threshold"])
    return combined


def _legacy_base_term(term: str) -> str:
    term = str(term).strip()
    parts = term.split("_")
    tissue = parts[0]
    tissue_key = tissue.upper().replace(" ", "-")
    base = LEGACY_TISSUE_TERMS.get(tissue_key, tissue)
    if len(parts) >= 2 and parts[1].lower() == "consensus":
        return f"{base}_Consensus"
    if len(parts) >= 3:
        sex = parts[1].capitalize()
        week = parts[2].upper()
        return f"{base}_{sex}_{week}"
    return base


def _build_signed_term_rows(processed_df: pd.DataFrame) -> list[dict[str, str]]:
    rows_out: list[dict[str, str]] = []
    for _, row in processed_df.iterrows():
        term = str(row.get("term", "")).strip()
        gene = str(row.get("gene", "")).strip()
        score_value = float(row.get("adj_p_value", 0.0) or 0.0)
        sign_value = float(row.get("threshold", 0.0) or 0.0)
        if not term or not gene or sign_value == 0.0:
            continue
        rows_out.append(
            {
                "term": _legacy_base_term(term),
                "gene_id": gene,
                "gene_symbol": gene,
                "score": str(score_value),
                "sign": str(sign_value),
            }
        )
    return rows_out


def run(args) -> dict[str, object]:
    feature_annot = Path(args.feature_annot).resolve()
    dea_dir = Path(args.dea_dir).resolve()
    mapping_file = Path(args.mapping_file).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    _require_file(feature_annot, "feature annotation")
    _require_dir(dea_dir, "DEA directory")
    _require_file(mapping_file, "mapping file")

    feature_df = _load_feature_annotation(feature_annot)
    gene_mapper = _build_gene_mapper(feature_df)
    symbol_map = _load_symbol_map(mapping_file)
    timewise_df, timewise_audit, timewise_paths = _load_timewise_dea(
        dea_dir=dea_dir,
        gene_mapper=gene_mapper,
        symbol_map=symbol_map,
        padj_max=float(args.padj_max),
    )
    training_df, training_audit, training_paths = _load_training_dea(
        dea_dir=dea_dir,
        gene_mapper=gene_mapper,
        symbol_map=symbol_map,
        padj_max=float(args.padj_max),
    )
    processed_df = _combine_and_standardize(timewise_df, training_df)

    processed_path = out_dir / "motrpac_processed.tsv"
    signed_term_path = out_dir / "motrpac_signed_term_gene.tsv"
    audit_path = out_dir / "motrpac_processing_audit.tsv"
    processed_df.to_csv(processed_path, sep="\t", index=False)
    audit_rows = (
        [{"metric": key, "value": value} for key, value in timewise_audit.items()]
        + [{"metric": key, "value": value} for key, value in training_audit.items()]
        + [
            {"metric": "combined_rows", "value": len(processed_df)},
            {"metric": "unique_terms", "value": int(processed_df["term"].nunique())},
            {"metric": "unique_genes", "value": int(processed_df["gene"].nunique())},
            {"metric": "padj_max", "value": float(args.padj_max)},
        ]
    )
    write_tsv(audit_path, audit_rows, ["metric", "value"])
    signed_term_rows = _build_signed_term_rows(processed_df)
    write_tsv(signed_term_path, signed_term_rows, ["term", "gene_id", "gene_symbol", "score", "sign"])

    input_paths: list[tuple[Path, str]] = [
        (feature_annot, "feature_annot"),
        (mapping_file, "mapping_file"),
    ]
    for dea_path in timewise_paths:
        input_paths.append((dea_path, "timewise_dea_tsv"))
    for dea_path in training_paths:
        input_paths.append((dea_path, "training_dea_tsv"))
    write_workflow_provenance_graph(
        workflow_name="motrpac_released_dea",
        module_name="geneset_extractors.workflows.motrpac_released_dea",
        output_dir=out_dir,
        focus_output_path=signed_term_path,
        output_paths=[
            (signed_term_path, "signed_term_gene_tsv"),
            (processed_path, "processed_tsv"),
            (audit_path, "processing_audit_tsv"),
        ],
        input_paths=input_paths,
        parameters={
            "padj_max": float(args.padj_max),
            "n_processed_rows": len(processed_df),
            "n_signed_term_rows": len(signed_term_rows),
            "n_terms": int(processed_df["term"].nunique()),
            "n_genes": int(processed_df["gene"].nunique()),
        },
    )
    return {
        "n_rows": len(processed_df),
        "n_terms": int(processed_df["term"].nunique()),
        "n_genes": int(processed_df["gene"].nunique()),
        "out_dir": str(out_dir),
    }
