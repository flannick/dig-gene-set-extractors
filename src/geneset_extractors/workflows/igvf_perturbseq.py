"""IGVF Perturb-seq long differential-expression to signed gene-set workflow."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from geneset_extractors.core.provenance import activate_runtime_context
from geneset_extractors.workflows.gtex_runtime_common import write_tsv, write_workflow_provenance_graph


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _load_gene_symbol_mapping(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    mapping = df.set_index(1)[2] if df.shape[1] >= 3 else df.set_index(0)[1]
    if df.shape[1] < 2:
        raise ValueError(f"Mapping file {path} has {df.shape[1]} column(s); expected at least 2.")
    return mapping.dropna().astype(str).to_dict()


def _zscore(row: pd.Series) -> pd.Series:
    std = row.std()
    return row * np.nan if pd.isna(std) or std == 0 else (row - row.mean()) / std


def _load_matrix(expression_tsv: Path, mapping_file: Path | None, orientation: str) -> pd.DataFrame:
    wide = pd.read_csv(expression_tsv, compression="gzip" if str(expression_tsv).endswith(".gz") else None, sep="\t")
    wide = wide.set_index("Unnamed: 0" if "Unnamed: 0" in wide.columns else wide.columns[0])
    if orientation == "perturbation_by_gene":
        wide = wide.T
    elif orientation != "gene_by_perturbation":
        raise ValueError(f"Unsupported orientation: {orientation}")
    wide.index = wide.index.astype(str).str.upper().map(_load_gene_symbol_mapping(mapping_file)) if mapping_file else wide.index.astype(str).str.upper()
    return wide[wide.index.notna() & ~wide.index.duplicated()]


def _read_table_any(path: Path, sep: str) -> pd.DataFrame:
    if sep == "auto":
        stem = str(path)[:-3] if str(path).endswith(".gz") else str(path)
        sep = "," if stem.lower().endswith(".csv") else "\t"
    return pd.read_csv(path, compression="gzip" if str(path).endswith(".gz") else None, sep=sep, dtype=str)


def _load_long_de(path: Path, args) -> pd.DataFrame:
    df = _read_table_any(path, str(getattr(args, "sep", "auto")))
    term_col, symbol_col = str(args.term_column), str(args.gene_symbol_column)
    for column, label in ((term_col, "term_column"), (symbol_col, "gene_symbol_column")):
        if column not in df.columns:
            raise ValueError(f"Column '{column}' ({label}) not found. Available: {list(df.columns)}")
    effect_col, ratio_col = getattr(args, "effect_column", None), getattr(args, "ratio_column", None)
    if not effect_col and not ratio_col:
        raise ValueError("long_de mode requires either --effect_column or --ratio_column for direction.")
    out = pd.DataFrame({"Perturbation": df[term_col].astype(str), "Gene": df[symbol_col].astype(str).str.upper()})
    gene_id_col = getattr(args, "gene_id_column", None)
    out["gene_id"] = df[gene_id_col].astype(str) if gene_id_col and gene_id_col in df.columns else out["Gene"]
    if effect_col:
        effect = pd.to_numeric(df[effect_col], errors="coerce")
        out["sign"], magnitude = np.where(effect > 0, 1, -1), effect.abs()
    else:
        ratio = pd.to_numeric(df[ratio_col], errors="coerce")
        out["sign"], magnitude = np.where(ratio > 1, 1, -1), np.log2(ratio.where(ratio > 0)).abs()
    score_col = getattr(args, "score_column", None)
    out["z"] = (pd.to_numeric(df[score_col], errors="coerce").abs() if score_col and score_col in df.columns else magnitude) * out["sign"]
    out = out.dropna(subset=["z"])
    pvalue_col, pvalue_max = getattr(args, "pvalue_column", None), getattr(args, "pvalue_max", None)
    if pvalue_col and pvalue_col in df.columns and pvalue_max is not None:
        out = out[pd.to_numeric(df[pvalue_col], errors="coerce") <= float(pvalue_max)]
    if getattr(args, "score_threshold", None) is not None:
        out = out[out["z"].abs() >= float(args.score_threshold)]
    out = out[out["Gene"].notna() & ~out["Gene"].isin(["", "NAN"])]
    out = out.assign(_absz=out["z"].abs()).sort_values("_absz", ascending=False).drop_duplicates(["Perturbation", "Gene"])
    top_k = getattr(args, "top_k_per_direction", None)
    if top_k:
        out = out.groupby(["Perturbation", "sign"], group_keys=False).head(int(top_k))
    # Pandas 3 preserves missing scalar values through ``astype(str)`` in a
    # StringArray; the legacy pandas behavior emitted their literal ``nan``
    # text.  Use Python conversion to preserve that historical contract.
    out["Perturbation"] = out["Perturbation"].map(str)
    out["Gene"] = out["Gene"].map(str)
    return out.drop(columns="_absz").reset_index(drop=True)


def _write_combined_gmt(long: pd.DataFrame, output_file: Path, min_gmt_size: int) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8", newline="\n") as handle:
        for direction, suffix in ((1, "_Up"), (-1, "_Down")):
            subset = long[long["sign"] == direction]
            for perturbation in sorted(subset["Perturbation"].unique()):
                genes = sorted(subset.loc[subset["Perturbation"] == perturbation, "Gene"].unique())
                if len(genes) >= min_gmt_size:
                    handle.write("\t".join([f"{perturbation}{suffix}", *genes]) + "\n")


def run(args) -> dict[str, object]:
    activate_runtime_context("igvf_perturbseq", getattr(args, "provenance_overlay_json", None))
    expression_tsv, out_dir = Path(args.expression_tsv).resolve(), Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    _require_file(expression_tsv, "expression TSV")
    mapping_file = Path(args.mapping_file).resolve() if getattr(args, "mapping_file", None) else None
    if mapping_file:
        _require_file(mapping_file, "mapping file")
    if str(getattr(args, "input_mode", "matrix")) == "long_de":
        long = _load_long_de(expression_tsv, args)
    else:
        matrix = _load_matrix(expression_tsv, mapping_file, str(getattr(args, "orientation", "perturbation_by_gene")))
        zmat = matrix.apply(_zscore, axis=1)
        long = zmat[abs(zmat) >= float(args.z_threshold)].stack().sort_values().to_frame().reset_index()
        long.columns = ["Gene", "Perturbation", "z"]
        long["sign"], long["gene_id"] = np.where(long["z"] > 0, 1, -1), long["Gene"]
    processed_rows = [{"gene": str(r.Gene), "gene_id": str(r.gene_id), "term": str(r.Perturbation), "score": str(abs(float(r.z))), "signed_score": str(float(r.z)), "sign": str(int(r.sign))} for r in long.itertuples()]
    processed_path = out_dir / "igvf_perturbseq_processed.tsv"
    write_tsv(processed_path, processed_rows, ["gene", "gene_id", "term", "score", "signed_score", "sign"])
    notebook_gmt = out_dir / str(args.gmt_name)
    _write_combined_gmt(long, notebook_gmt, int(args.min_gmt_size))
    signed_rows = sorted(({"term": r["term"], "gene_id": r["gene_id"], "gene_symbol": r["gene"], "score": r["score"], "sign": r["sign"]} for r in processed_rows), key=lambda r: (r["term"], -int(r["sign"]), r["gene_symbol"]))
    signed_path = out_dir / "igvf_perturbseq_signed_term_gene.tsv"
    write_tsv(signed_path, signed_rows, ["term", "gene_id", "gene_symbol", "score", "sign"])
    source_input_id = str(getattr(args, "source_input_id", "expression_tsv"))
    inputs = [(expression_tsv, source_input_id)] + ([(mapping_file, "mapping_file")] if mapping_file else [])
    write_workflow_provenance_graph(workflow_name="igvf_perturbseq", module_name="geneset_extractors.workflows.igvf_perturbseq", output_dir=out_dir, focus_output_path=signed_path, output_paths=[(signed_path, "table_tsv"), (processed_path, "processed_tsv"), (notebook_gmt, "notebook_combined_gmt")], input_paths=inputs, parameters={"input_mode": str(getattr(args, "input_mode", "matrix")), "source_input_id": source_input_id, "z_threshold": float(args.z_threshold), "min_gmt_size": int(args.min_gmt_size), "orientation": str(getattr(args, "orientation", "perturbation_by_gene")), "n_rows": len(processed_rows)})
    return {"n_rows": len(processed_rows), "out_dir": str(out_dir)}
