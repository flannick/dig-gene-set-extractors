from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from geneset_extractors.workflows.gtex_runtime_common import write_tsv, write_workflow_provenance_graph


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _load_gene_symbol_mapping(path: Path) -> dict[str, str]:
    """Mapping file: 2 columns (id, symbol) or >=3 columns (col1=id, col2=symbol)."""
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if df.shape[1] >= 3:
        mapping = df.set_index(1)[2]
    elif df.shape[1] == 2:
        mapping = df.set_index(0)[1]
    else:
        raise ValueError(
            f"Mapping file {path} has {df.shape[1]} column(s); expected 2 columns or at least 3 columns."
        )
    mapping = mapping.dropna()
    mapping.index = mapping.index.astype(str)
    return mapping.astype(str).to_dict()


def _zscore(row: pd.Series) -> pd.Series:
    std = row.std()
    if pd.isna(std) or std == 0:
        return row * np.nan
    return (row - row.mean()) / std


def _load_gene_by_perturbation(
    expression_tsv: Path,
    mapping_file: Path | None,
    orientation: str,
) -> pd.DataFrame:
    """Return a genes (rows) x perturbations (columns) matrix of effect scores.

    orientation:
      - "perturbation_by_gene": input rows are perturbations, columns are genes (matches the
        LINCS released-matrix layout). The matrix is transposed so genes become rows.
      - "gene_by_perturbation": input rows are genes, columns are perturbations (used as-is).
    """
    compression = "gzip" if str(expression_tsv).endswith(".gz") else None
    wide = pd.read_csv(expression_tsv, compression=compression, sep="\t")
    index_col = "Unnamed: 0" if "Unnamed: 0" in wide.columns else wide.columns[0]
    wide = wide.set_index(index_col)
    if orientation == "perturbation_by_gene":
        wide = wide.T
    elif orientation != "gene_by_perturbation":
        raise ValueError(f"Unsupported orientation: {orientation}")
    wide = wide.rename_axis("Gene", axis=0).rename_axis("Perturbation", axis=1)
    if mapping_file is not None:
        genemapping = _load_gene_symbol_mapping(mapping_file)
        wide.index = wide.index.astype(str).str.upper().map(genemapping)
    else:
        wide.index = wide.index.astype(str).str.upper()
    wide = wide[wide.index.notna()]
    wide = wide[~wide.index.duplicated()]
    return wide


def _threshold_long(gene_by_pert: pd.DataFrame, z_threshold: float) -> pd.DataFrame:
    zmat = gene_by_pert.apply(_zscore, axis=1)
    # pandas 3's stack implementation retains masked NaNs, unlike the prior
    # default. Drop them explicitly so only threshold-passing z-scores become
    # signed records; otherwise NaNs are incorrectly classified as down.
    long = zmat[abs(zmat) >= z_threshold].stack().dropna().sort_values().to_frame().reset_index()
    long.columns = ["Gene", "Perturbation", "z"]
    long = long[long["Gene"].notna()]
    long["Gene"] = long["Gene"].astype(str).str.upper()
    long["sign"] = long["z"].apply(lambda x: 1 if x > 0 else -1)
    return long


def _read_table_any(path: Path, sep: str) -> pd.DataFrame:
    compression = "gzip" if str(path).endswith(".gz") else None
    if sep == "auto":
        stem = str(path)[:-3] if str(path).endswith(".gz") else str(path)
        sep = "," if stem.lower().endswith(".csv") else "\t"
    return pd.read_csv(path, compression=compression, sep=sep, dtype=str)


def _load_long_de(path: Path, args) -> pd.DataFrame:
    """Build a long (term, Gene, score, sign) frame from a tidy per-perturbation DE table.

    Direction (up/down) is taken from a signed effect column (``--effect_column``) or from a
    fold-change *ratio* column (``--ratio_column``: >1 is up, <1 is down). Magnitude for
    ranking/thresholding comes from ``--score_column`` (absolute value); if omitted it falls
    back to the absolute signed effect. Optional significance filtering uses ``--pvalue_column``.
    """
    sep = str(getattr(args, "sep", "auto"))
    df = _read_table_any(path, sep)

    term_col = str(args.term_column)
    symbol_col = str(args.gene_symbol_column)
    gene_id_col = getattr(args, "gene_id_column", None)
    effect_col = getattr(args, "effect_column", None)
    ratio_col = getattr(args, "ratio_column", None)
    score_col = getattr(args, "score_column", None)
    pvalue_col = getattr(args, "pvalue_column", None)

    for required, label in [(term_col, "term_column"), (symbol_col, "gene_symbol_column")]:
        if required not in df.columns:
            raise ValueError(f"Column '{required}' ({label}) not found. Available: {list(df.columns)}")
    if not effect_col and not ratio_col:
        raise ValueError("long_de mode requires either --effect_column or --ratio_column for direction.")

    out = pd.DataFrame()
    out["Perturbation"] = df[term_col].astype(str)
    out["Gene"] = df[symbol_col].astype(str).str.upper()
    out["gene_id"] = df[gene_id_col].astype(str) if gene_id_col and gene_id_col in df.columns else out["Gene"]

    if effect_col:
        if effect_col not in df.columns:
            raise ValueError(f"Column '{effect_col}' (effect_column) not found. Available: {list(df.columns)}")
        effect = pd.to_numeric(df[effect_col], errors="coerce")
        out["sign"] = np.where(effect > 0, 1, -1)
        magnitude = effect.abs()
    else:
        if ratio_col not in df.columns:
            raise ValueError(f"Column '{ratio_col}' (ratio_column) not found. Available: {list(df.columns)}")
        ratio = pd.to_numeric(df[ratio_col], errors="coerce")
        # Fold-change ratio: >1 means up-regulated, <1 means down-regulated (log2 pivot at 1).
        out["sign"] = np.where(ratio > 1, 1, -1)
        magnitude = (np.log2(ratio.where(ratio > 0))).abs()

    if score_col and score_col in df.columns:
        out["z"] = pd.to_numeric(df[score_col], errors="coerce").abs() * out["sign"]
    else:
        out["z"] = magnitude * out["sign"]

    out = out.dropna(subset=["z"])
    if pvalue_col and pvalue_col in df.columns:
        pvals = pd.to_numeric(df[pvalue_col], errors="coerce")
        pmax = getattr(args, "pvalue_max", None)
        if pmax is not None:
            keep = pvals <= float(pmax)
            out = out[keep.reindex(out.index, fill_value=False)]

    score_threshold = getattr(args, "score_threshold", None)
    if score_threshold is not None:
        out = out[out["z"].abs() >= float(score_threshold)]

    out = out[out["Gene"].notna() & (out["Gene"] != "NAN") & (out["Gene"] != "")]
    # Collapse duplicate (term, gene) pairs to the strongest-magnitude observation.
    out["_absz"] = out["z"].abs()
    out = out.sort_values("_absz", ascending=False).drop_duplicates(["Perturbation", "Gene"]).drop(columns="_absz")

    top_k = getattr(args, "top_k_per_direction", None)
    if top_k:
        top_k = int(top_k)
        out = (
            out.assign(_absz=out["z"].abs())
            .sort_values("_absz", ascending=False)
            .groupby(["Perturbation", "sign"], group_keys=False)
            .head(top_k)
            .drop(columns="_absz")
        )
    return out.reset_index(drop=True)


def _write_combined_gmt(long: pd.DataFrame, output_file: Path, min_gmt_size: int) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8", newline="\n") as handle:
        for direction, suffix in [(1, "_Up"), (-1, "_Down")]:
            subset = long[long["sign"] == direction]
            for perturbation in sorted(subset["Perturbation"].unique()):
                genes = sorted(subset.loc[subset["Perturbation"] == perturbation, "Gene"].unique())
                if len(genes) >= min_gmt_size:
                    handle.write("\t".join([f"{perturbation}{suffix}", *genes]) + "\n")


def run(args) -> dict[str, object]:
    expression_tsv = Path(args.expression_tsv).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    _require_file(expression_tsv, "expression TSV")

    mapping_file = None
    if getattr(args, "mapping_file", None):
        mapping_file = Path(args.mapping_file).resolve()
        _require_file(mapping_file, "mapping file")

    z_threshold = float(args.z_threshold)
    min_gmt_size = int(args.min_gmt_size)
    gmt_name = str(args.gmt_name)
    orientation = str(getattr(args, "orientation", "perturbation_by_gene"))
    input_mode = str(getattr(args, "input_mode", "matrix"))

    if input_mode == "long_de":
        long = _load_long_de(expression_tsv, args)
    elif input_mode == "matrix":
        gene_by_pert = _load_gene_by_perturbation(expression_tsv, mapping_file, orientation)
        long = _threshold_long(gene_by_pert, z_threshold)
    else:
        raise ValueError(f"Unsupported input_mode: {input_mode}")

    has_gene_id = "gene_id" in long.columns
    processed_rows = [
        {
            "gene": str(row["Gene"]),
            "gene_id": str(row["gene_id"]) if has_gene_id else str(row["Gene"]),
            "term": str(row["Perturbation"]),
            "score": str(abs(float(row["z"]))),
            "signed_score": str(float(row["z"])),
            "sign": str(int(row["sign"])),
        }
        for _, row in long.iterrows()
    ]
    processed_path = out_dir / "igvf_perturbseq_processed.tsv"
    write_tsv(processed_path, processed_rows, ["gene", "gene_id", "term", "score", "signed_score", "sign"])

    notebook_gmt = out_dir / gmt_name
    _write_combined_gmt(long, notebook_gmt, min_gmt_size)

    signed_rows = sorted(
        [
            {
                "term": row["term"],
                "gene_id": row["gene_id"],
                "gene_symbol": row["gene"],
                "score": row["score"],
                "sign": row["sign"],
            }
            for row in processed_rows
        ],
        key=lambda row: (row["term"], -int(row["sign"]), row["gene_symbol"]),
    )
    signed_path = out_dir / "igvf_perturbseq_signed_term_gene.tsv"
    write_tsv(signed_path, signed_rows, ["term", "gene_id", "gene_symbol", "score", "sign"])

    input_paths = [(expression_tsv, "expression_tsv")]
    if mapping_file is not None:
        input_paths.append((mapping_file, "mapping_file"))

    write_workflow_provenance_graph(
        workflow_name="igvf_perturbseq",
        module_name="geneset_extractors.workflows.igvf_perturbseq",
        output_dir=out_dir,
        focus_output_path=signed_path,
        output_paths=[
            (signed_path, "signed_term_gene_tsv"),
            (processed_path, "processed_tsv"),
            (notebook_gmt, "notebook_combined_gmt"),
        ],
        input_paths=input_paths,
        parameters={
            "input_mode": input_mode,
            "z_threshold": z_threshold,
            "min_gmt_size": min_gmt_size,
            "orientation": orientation,
            "n_rows": len(processed_rows),
            "n_terms": int(long["Perturbation"].nunique()) if not long.empty else 0,
            "n_genes": int(long["Gene"].nunique()) if not long.empty else 0,
        },
    )
    return {"n_rows": len(processed_rows), "out_dir": str(out_dir)}
