from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

import numpy as np

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


def _build_gene_mapper(feature_annot: Path) -> dict[str, str]:
    with feature_annot.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"gene_id", "gene_name"}
    if not rows:
        raise ValueError(f"No feature annotation rows found in {feature_annot}")
    missing = required - set(rows[0].keys())
    if missing:
        raise ValueError(f"Feature annotation is missing required columns: {sorted(missing)}")
    return {
        str(row.get("gene_id", "")).strip(): str(row.get("gene_name", "")).strip()
        for row in rows
        if str(row.get("gene_id", "")).strip() and str(row.get("gene_name", "")).strip()
    }


def _load_symbol_map(mapping_file: Path) -> dict[str, str]:
    rows: list[list[str]] = []
    with mapping_file.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = [list(row) for row in reader]
    if not rows:
        raise ValueError(f"Mapping file was empty: {mapping_file}")
    width = max(len(row) for row in rows)
    if width >= 3:
        input_col, output_col = 1, 2
    elif width == 2:
        input_col, output_col = 0, 1
    else:
        raise ValueError("Mapping file must have at least two tab-delimited columns")
    symbol_map: dict[str, str] = {}
    for row in rows:
        if len(row) <= max(input_col, output_col):
            continue
        source = str(row[input_col]).strip()
        target = str(row[output_col]).strip()
        if source and target:
            symbol_map[source] = target
    return symbol_map


def _load_timewise_rows(
    *,
    dea_dir: Path,
    gene_mapper: dict[str, str],
    symbol_map: dict[str, str],
    padj_max: float,
) -> tuple[list[dict[str, Any]], dict[str, int], list[Path]]:
    rows_out: list[dict[str, Any]] = []
    dea_paths: list[Path] = []
    rows_loaded = 0
    rows_after_gene_mapper = 0
    rows_after_padj = 0
    rows_after_symbol_map = 0
    for dea_path in sorted(dea_dir.iterdir()):
        if "timewise" not in dea_path.name or not dea_path.is_file():
            continue
        dea_paths.append(dea_path)
        with dea_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                rows_loaded += 1
                feature_id = str(row.get("feature_ID", "")).strip()
                mapped_gene = gene_mapper.get(feature_id)
                if not mapped_gene:
                    continue
                rows_after_gene_mapper += 1
                try:
                    adj_p_value = float(str(row.get("adj_p_value", "")).strip())
                except ValueError:
                    continue
                if adj_p_value >= float(padj_max):
                    continue
                rows_after_padj += 1
                human_gene = symbol_map.get(mapped_gene.upper())
                if not human_gene:
                    continue
                rows_after_symbol_map += 1
                try:
                    logfc = float(str(row.get("logFC", "")).strip())
                except ValueError:
                    continue
                sign = _threshold(logfc)
                if sign is None:
                    continue
                term = "_".join(
                    [
                        str(row.get("tissue", "")).strip(),
                        str(row.get("sex", "")).strip(),
                        str(row.get("comparison_group", "")).strip(),
                    ]
                )
                rows_out.append(
                    {
                        "term": term,
                        "gene": human_gene,
                        "adj_p_value": adj_p_value,
                        "threshold": sign,
                    }
                )
    if not dea_paths:
        raise ValueError(f"No DEA files containing 'timewise' were found in {dea_dir}")
    audit = {
        "timewise_rows_loaded": rows_loaded,
        "timewise_rows_after_gene_mapper": rows_after_gene_mapper,
        "timewise_rows_after_padj": rows_after_padj,
        "timewise_rows_after_symbol_map": rows_after_symbol_map,
    }
    return rows_out, audit, dea_paths


def _load_training_rows(
    *,
    dea_dir: Path,
    gene_mapper: dict[str, str],
    symbol_map: dict[str, str],
    padj_max: float,
) -> tuple[list[dict[str, Any]], dict[str, int], list[Path]]:
    rows_out: list[dict[str, Any]] = []
    dea_paths: list[Path] = []
    rows_loaded = 0
    rows_after_gene_mapper = 0
    rows_after_padj = 0
    rows_after_symbol_map = 0
    for dea_path in sorted(dea_dir.iterdir()):
        if "training" not in dea_path.name or not dea_path.is_file():
            continue
        dea_paths.append(dea_path)
        with dea_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                rows_loaded += 1
                feature_id = str(row.get("feature_ID", "")).strip()
                mapped_gene = gene_mapper.get(feature_id)
                if not mapped_gene:
                    continue
                rows_after_gene_mapper += 1
                try:
                    adj_p_value = float(str(row.get("adj_p_value", "")).strip())
                except ValueError:
                    continue
                if adj_p_value >= float(padj_max):
                    continue
                rows_after_padj += 1
                human_gene = symbol_map.get(mapped_gene.upper())
                if not human_gene:
                    continue
                rows_after_symbol_map += 1
                term = f"{str(row.get('tissue', '')).strip()}_consensus"
                rows_out.append(
                    {
                        "term": term,
                        "gene": human_gene,
                        "adj_p_value": adj_p_value,
                        "threshold": 1,
                    }
                )
    if not dea_paths:
        raise ValueError(f"No DEA files containing 'training' were found in {dea_dir}")
    audit = {
        "training_rows_loaded": rows_loaded,
        "training_rows_after_gene_mapper": rows_after_gene_mapper,
        "training_rows_after_padj": rows_after_padj,
        "training_rows_after_symbol_map": rows_after_symbol_map,
    }
    return rows_out, audit, dea_paths


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


def _build_signed_term_rows(processed_rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    rows_out: list[dict[str, str]] = []
    for row in processed_rows:
        term = str(row.get("term", "")).strip()
        gene = str(row.get("gene", "")).strip()
        if not term or not gene:
            continue
        score_value = float(row.get("adj_p_value", 0.0) or 0.0)
        sign_value = float(row.get("threshold", 0.0) or 0.0)
        if sign_value == 0.0:
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

    gene_mapper = _build_gene_mapper(feature_annot)
    symbol_map = _load_symbol_map(mapping_file)
    timewise_rows, timewise_audit, timewise_paths = _load_timewise_rows(
        dea_dir=dea_dir,
        gene_mapper=gene_mapper,
        symbol_map=symbol_map,
        padj_max=float(args.padj_max),
    )
    training_rows, training_audit, training_paths = _load_training_rows(
        dea_dir=dea_dir,
        gene_mapper=gene_mapper,
        symbol_map=symbol_map,
        padj_max=float(args.padj_max),
    )
    processed_rows = list(timewise_rows) + list(training_rows)
    for row in processed_rows:
        standardized = np.log10(float(row["adj_p_value"])) * -1
        row["adj_p_value"] = standardized * float(row["threshold"])

    processed_path = out_dir / "motrpac_processed.tsv"
    signed_term_path = out_dir / "motrpac_signed_term_gene.tsv"
    audit_path = out_dir / "motrpac_processing_audit.tsv"
    write_tsv(processed_path, processed_rows, ["term", "gene", "adj_p_value", "threshold"])
    audit_rows = (
        [{"metric": key, "value": value} for key, value in timewise_audit.items()]
        + [{"metric": key, "value": value} for key, value in training_audit.items()]
        + [
            {"metric": "combined_rows", "value": len(processed_rows)},
            {"metric": "unique_terms", "value": len({str(row["term"]) for row in processed_rows})},
            {"metric": "unique_genes", "value": len({str(row["gene"]) for row in processed_rows})},
            {"metric": "padj_max", "value": float(args.padj_max)},
        ]
    )
    write_tsv(audit_path, audit_rows, ["metric", "value"])
    signed_term_rows = _build_signed_term_rows(processed_rows)
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
            "n_processed_rows": len(processed_rows),
            "n_signed_term_rows": len(signed_term_rows),
            "n_terms": len({str(row["term"]) for row in processed_rows}),
            "n_genes": len({str(row["gene"]) for row in processed_rows}),
        },
    )
    return {
        "n_rows": len(processed_rows),
        "n_terms": len({str(row["term"]) for row in processed_rows}),
        "n_genes": len({str(row["gene"]) for row in processed_rows}),
        "out_dir": str(out_dir),
    }
