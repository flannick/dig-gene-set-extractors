"""Convert released IGVF Perturb-seq differential-expression tables to signed GMTs.

This workflow preserves the legacy IGVF PS1 behavior: it applies the
per-analysis-set schema declared by the wrapper, retains significant rows when
a p-value threshold is supplied, assigns direction from a signed effect or a
fold-change ratio, ranks each term/direction by magnitude, and emits the
legacy-compatible signed term-gene intermediate plus the authoritative GMT.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from geneset_extractors.extractors.converters import signed_term_gene
from geneset_extractors.workflows.gtex_runtime_common import write_tsv, write_workflow_provenance_graph


TERM_PREFIX = "IGVF_Perturb_Seq"
SIGNED_TSV_NAME = "igvf_perturbseq_signed_term_gene.tsv"


def _setting(row: dict[str, str], name: str) -> str:
    value = str(row.get(name, "")).strip()
    return "" if value.upper() == "NA" else value


def _load_settings(path: Path, analysis_set_id: str) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    matches = [row for row in rows if str(row.get("analysis_set_id", "")).strip() == analysis_set_id]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one analysis_set_id={analysis_set_id!r} in {path}")
    return {str(key): str(value) for key, value in matches[0].items()}


def _read_table(path: Path, sep: str) -> pd.DataFrame:
    resolved_sep = None if sep in {"", "auto"} else sep
    return pd.read_csv(path, sep=resolved_sep, engine="python" if resolved_sep is None else None, dtype=str)


def _number(value: object, label: str) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Could not parse {label} value {value!r} as a number") from exc


def _signed_rows(table: pd.DataFrame, settings: dict[str, str]) -> list[dict[str, str]]:
    term_column = _setting(settings, "term_column")
    symbol_column = _setting(settings, "gene_symbol_column")
    gene_id_column = _setting(settings, "gene_id_column")
    effect_column = _setting(settings, "effect_column")
    ratio_column = _setting(settings, "ratio_column")
    score_column = _setting(settings, "score_column")
    pvalue_column = _setting(settings, "pvalue_column")
    pvalue_max = _setting(settings, "pvalue_max")
    required = [term_column, symbol_column]
    required.extend(column for column in (effect_column, ratio_column, score_column, pvalue_column) if column)
    missing = [column for column in required if column not in table.columns]
    if missing:
        raise ValueError(f"Input table is missing declared columns: {', '.join(missing)}")
    if bool(effect_column) == bool(ratio_column):
        raise ValueError("Declare exactly one of effect_column or ratio_column")
    out = pd.DataFrame()
    out["term"] = table[term_column].astype(str)
    out["gene"] = table[symbol_column].astype(str).str.upper()
    out["gene_id"] = table[gene_id_column].astype(str) if gene_id_column and gene_id_column in table.columns else out["gene"]
    if effect_column:
        effect = pd.to_numeric(table[effect_column], errors="coerce")
        out["sign"] = np.where(effect > 0, 1, -1)
        magnitude = effect.abs()
    else:
        ratio = pd.to_numeric(table[ratio_column], errors="coerce")
        out["sign"] = np.where(ratio > 1, 1, -1)
        magnitude = np.log2(ratio.where(ratio > 0)).abs()
    out["signed_score"] = pd.to_numeric(table[score_column], errors="coerce").abs() * out["sign"] if score_column else magnitude * out["sign"]
    out = out.dropna(subset=["signed_score"])
    if pvalue_column and pvalue_max:
        out = out[pd.to_numeric(table[pvalue_column], errors="coerce") <= float(pvalue_max)]
    score_threshold = _setting(settings, "score_threshold")
    if score_threshold:
        out = out[out["signed_score"].abs() >= float(score_threshold)]
    out = out[out["gene"].notna() & (out["gene"] != "NAN") & (out["gene"] != "")]
    out["_absz"] = out["signed_score"].abs()
    out = out.sort_values("_absz", ascending=False).drop_duplicates(["term", "gene"]).drop(columns="_absz")
    top_k = _setting(settings, "top_k_per_direction")
    if top_k:
        out = out.assign(_absz=out["signed_score"].abs()).sort_values("_absz", ascending=False).groupby(["term", "sign"], group_keys=False).head(int(top_k)).drop(columns="_absz")
    return [
        {"gene": str(row.gene), "gene_id": str(row.gene_id), "term": str(row.term), "score": str(abs(float(row.signed_score))), "signed_score": str(float(row.signed_score)), "sign": str(int(row.sign))}
        for row in out.reset_index(drop=True).itertuples(index=False)
    ]


def _write_workflow_gmt(rows: list[dict[str, str]], path: Path, min_size: int) -> None:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["term"], row["sign"])].append(row)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for (term, sign), group in sorted(grouped.items()):
            seen: set[str] = set()
            genes: list[str] = []
            for gene in sorted(row["gene"] for row in group):
                if gene not in seen:
                    seen.add(gene)
                    genes.append(gene)
            if len(genes) >= min_size:
                label = "Up" if sign == "1" else "Down"
                handle.write("\t".join([f"{term}_{label}", *genes]) + "\n")


def run(args) -> dict[str, object]:
    expression_tsv = Path(args.expression_tsv).resolve()
    manifest_tsv = Path(args.analysis_set_manifest).resolve()
    out_dir = Path(args.out_dir).resolve()
    if not expression_tsv.is_file():
        raise FileNotFoundError(f"Missing expression TSV: {expression_tsv}")
    if not manifest_tsv.is_file():
        raise FileNotFoundError(f"Missing analysis-set manifest: {manifest_tsv}")
    settings = _load_settings(manifest_tsv, args.analysis_set_id)
    rows = _signed_rows(_read_table(expression_tsv, _setting(settings, "sep") or "auto"), settings)
    workflow_dir = out_dir / "workflow"
    extractor_dir = out_dir / "extractor"
    workflow_dir.mkdir(parents=True, exist_ok=True)
    processed_path = workflow_dir / "igvf_perturbseq_processed.tsv"
    signed_path = workflow_dir / SIGNED_TSV_NAME
    write_tsv(processed_path, rows, ["gene", "gene_id", "term", "score", "signed_score", "sign"])
    signed_rows = sorted(({"term": row["term"], "gene_id": row["gene_id"], "gene_symbol": row["gene"], "score": row["score"], "sign": row["sign"]} for row in rows), key=lambda row: (row["term"], -int(row["sign"]), row["gene_symbol"]))
    write_tsv(signed_path, signed_rows, ["term", "gene_id", "gene_symbol", "score", "sign"])
    workflow_gmt = workflow_dir / str(args.gmt_name)
    _write_workflow_gmt(rows, workflow_gmt, int(args.min_gmt_size))
    write_workflow_provenance_graph(
        workflow_name="igvf_perturbseq", module_name=__name__, output_dir=workflow_dir, focus_output_path=signed_path,
        output_paths=[(signed_path, "signed_term_gene_tsv"), (processed_path, "processed_tsv"), (workflow_gmt, "workflow_gmt")],
        input_paths=[(expression_tsv, "released_igvf_differential_expression_tsv"), (manifest_tsv, "analysis_set_manifest")],
        parameters={"analysis_set_id": args.analysis_set_id, "min_gmt_size": int(args.min_gmt_size), "n_rows": len(rows)},
    )
    converter_args = type("SignedTermGeneArgs", (), {
        "table_tsv": str(signed_path), "out_dir": str(extractor_dir), "organism": args.organism, "genome_build": args.genome_build,
        "term_column": "term", "term_prefix": TERM_PREFIX, "gene_id_column": "gene_id", "gene_symbol_column": "gene_symbol",
        "score_column": "score", "sign_column": "sign", "emit_mode": "grouped_rows", "gmt_name_separator": "_",
        "gmt_signed_labels": "up_dn", "gmt_min_genes": int(args.min_gmt_size), "emit_small_gene_sets": False,
        "emit_gmt": True, "gmt_prefer_symbol": True, "gmt_require_symbol": True, "gmt_format": "classic", "provenance_overlay_json": None,
    })()
    signed_term_gene.run(converter_args)
    return {"n_rows": len(rows), "out_dir": str(out_dir)}
