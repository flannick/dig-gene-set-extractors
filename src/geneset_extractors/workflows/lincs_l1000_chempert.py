from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from geneset_extractors.workflows.gtex_runtime_common import write_tsv, write_workflow_provenance_graph


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _load_gene_symbol_mapping(path: Path) -> dict[str, str]:
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


def _zscore(gene: pd.Series) -> pd.Series:
    std = gene.std()
    if pd.isna(std) or std == 0:
        return gene * np.nan
    return (gene - gene.mean()) / std


def _preprocess_chempert(expression_tsv: Path, mapping_file: Path, z_threshold: float) -> pd.DataFrame:
    compression = "gzip" if str(expression_tsv).endswith(".gz") else None
    chempert_wide = pd.read_csv(expression_tsv, compression=compression, sep="\t")
    index_col = "Unnamed: 0" if "Unnamed: 0" in chempert_wide.columns else chempert_wide.columns[0]
    chempert_wide = chempert_wide.set_index(index_col).T.rename_axis("Gene", axis=0).rename_axis(
        "Chemical Perturbation", axis=1
    )
    genemapping = _load_gene_symbol_mapping(mapping_file)
    upper_index = chempert_wide.index.astype(str).str.upper()
    mapped_index = upper_index.map(genemapping)
    chempert_wide.index = mapped_index
    chempert_wide = chempert_wide[chempert_wide.index.duplicated() == False]
    chempert_wide = chempert_wide.apply(_zscore, axis=1)
    chempert = chempert_wide[abs(chempert_wide) >= z_threshold].stack().sort_values().to_frame().reset_index()
    chempert = chempert[chempert["Gene"].isna() == False]
    chempert.columns = ["Gene", "Chemical Perturbation", "z"]
    chempert["Gene"] = chempert["Gene"].astype(str).str.upper()
    chempert["threshold"] = chempert["z"].apply(lambda x: 1 if x > 0 else -1)
    return chempert


def _write_combined_gmt(chempert: pd.DataFrame, output_file: Path, min_gmt_size: int) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8", newline="\n") as handle:
        for direction, suffix in [(1, "_Up"), (-1, "_Down")]:
            subset = chempert[chempert["threshold"] == direction]
            for perturbation in sorted(subset["Chemical Perturbation"].unique()):
                genes = sorted(subset.loc[subset["Chemical Perturbation"] == perturbation, "Gene"].unique())
                if len(genes) >= min_gmt_size:
                    handle.write("\t".join([f"{perturbation}{suffix}", *genes]) + "\n")


def run(args) -> dict[str, object]:
    expression_tsv = Path(args.expression_tsv).resolve()
    mapping_file = Path(args.mapping_file).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    _require_file(expression_tsv, "expression TSV")
    _require_file(mapping_file, "mapping file")

    z_threshold = float(args.z_threshold)
    min_gmt_size = int(args.min_gmt_size)
    gmt_name = str(args.gmt_name)

    chempert = _preprocess_chempert(expression_tsv, mapping_file, z_threshold)
    processed_rows = [
        {
            "gene": str(row["Gene"]),
            "term": str(row["Chemical Perturbation"]),
            "score": str(abs(float(row["z"]))),
            "signed_score": str(float(row["z"])),
            "sign": str(int(row["threshold"])),
        }
        for _, row in chempert.iterrows()
    ]
    processed_path = out_dir / "lincs_l1000_processed.tsv"
    write_tsv(processed_path, processed_rows, ["gene", "term", "score", "signed_score", "sign"])

    notebook_gmt = out_dir / gmt_name
    _write_combined_gmt(chempert, notebook_gmt, min_gmt_size)

    signed_rows = sorted(
        [
            {
                "term": row["term"],
                "gene_id": row["gene"],
                "gene_symbol": row["gene"],
                "score": row["score"],
                "sign": row["sign"],
            }
            for row in processed_rows
        ],
        key=lambda row: (row["term"], -int(row["sign"]), row["gene_symbol"]),
    )
    signed_path = out_dir / "lincs_l1000_signed_term_gene.tsv"
    write_tsv(signed_path, signed_rows, ["term", "gene_id", "gene_symbol", "score", "sign"])

    write_workflow_provenance_graph(
        workflow_name="lincs_l1000_chempert",
        module_name="geneset_extractors.workflows.lincs_l1000_chempert",
        output_dir=out_dir,
        focus_output_path=signed_path,
        output_paths=[
            (signed_path, "signed_term_gene_tsv"),
            (processed_path, "processed_tsv"),
            (notebook_gmt, "notebook_combined_gmt"),
        ],
        input_paths=[
            (expression_tsv, "expression_tsv"),
            (mapping_file, "mapping_file"),
        ],
        parameters={
            "z_threshold": z_threshold,
            "min_gmt_size": min_gmt_size,
            "n_rows": len(processed_rows),
            "n_terms": int(chempert["Chemical Perturbation"].nunique()) if not chempert.empty else 0,
            "n_genes": int(chempert["Gene"].nunique()) if not chempert.empty else 0,
        },
    )
    return {"n_rows": len(processed_rows), "out_dir": str(out_dir)}
