from __future__ import annotations

from pathlib import Path

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


def _preprocess_crisprko(expression_tsv: Path, mapping_file: Path, top_n: int) -> pd.DataFrame:
    l1000 = pd.read_csv(expression_tsv, sep="\t", index_col=0)
    l1000 = l1000[l1000.index.map(lambda x: not pd.isna(x))]
    l1000 = l1000[l1000.index.map(lambda x: not str(x).startswith("BRDN"))]
    genemapping = _load_gene_symbol_mapping(mapping_file)
    l1000 = l1000[l1000.index.map(lambda x: x in genemapping)]
    l1000.index = l1000.index.map(lambda x: genemapping[x])
    l1000 = l1000.T
    l1000 = l1000[l1000.index.map(lambda x: x in genemapping)]
    l1000.index = l1000.index.map(lambda x: genemapping[x])
    l1000 = l1000.sort_index(axis=1)
    l1000 = l1000.rename_axis("Gene", axis=0).rename_axis("Gene KO", axis=1)
    l1000 = l1000.stack().reset_index()
    up = (
        l1000[l1000[0] > 0]
        .groupby("Gene KO")
        .apply(lambda x: x.sort_values(0, ascending=False).head(top_n))
        .reset_index(drop=True)
    )
    down = (
        l1000[l1000[0] < 0]
        .groupby("Gene KO")
        .apply(lambda x: x.sort_values(0).head(top_n))
        .reset_index(drop=True)
    )
    l1000 = pd.concat([up, down]).reset_index(drop=True)
    l1000["Threshold Value"] = l1000[0].apply(lambda x: 1 if x > 0 else -1)
    return l1000


def _write_combined_gmt(l1000: pd.DataFrame, output_file: Path, min_gmt_size: int) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    up_by_ko = l1000[l1000["Threshold Value"] == 1].groupby("Gene KO", sort=False)
    down_by_ko = l1000[l1000["Threshold Value"] == -1].groupby("Gene KO", sort=False)
    with output_file.open("w", encoding="utf-8", newline="\n") as handle:
        for attribute, group in up_by_ko:
            genes = group["Gene"].tolist()
            if len(genes) >= min_gmt_size:
                handle.write("\t".join([f"{attribute}_Up", *genes]) + "\n")
        for attribute, group in down_by_ko:
            genes = group["Gene"].tolist()
            if len(genes) >= min_gmt_size:
                handle.write("\t".join([f"{attribute}_Down", *genes]) + "\n")


def run(args) -> dict[str, object]:
    expression_tsv = Path(args.expression_tsv).resolve()
    mapping_file = Path(args.mapping_file).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    _require_file(expression_tsv, "expression TSV")
    _require_file(mapping_file, "mapping file")

    top_n = int(args.top_n)
    min_gmt_size = int(args.min_gmt_size)
    gmt_name = str(args.gmt_name)

    l1000 = _preprocess_crisprko(expression_tsv, mapping_file, top_n)
    processed_rows = [
        {
            "gene": str(row["Gene"]),
            "term": str(row["Gene KO"]),
            "score": str(abs(float(row[0]))),
            "signed_score": str(float(row[0])),
            "sign": str(int(row["Threshold Value"])),
        }
        for _, row in l1000.iterrows()
    ]
    processed_path = out_dir / "lincs_l1000_processed.tsv"
    write_tsv(processed_path, processed_rows, ["gene", "term", "score", "signed_score", "sign"])

    notebook_gmt = out_dir / gmt_name
    _write_combined_gmt(l1000, notebook_gmt, min_gmt_size)

    signed_rows = [
        {
            "term": row["term"],
            "gene_id": row["gene"],
            "gene_symbol": row["gene"],
            "score": row["score"],
            "sign": row["sign"],
        }
        for row in processed_rows
    ]
    signed_path = out_dir / "lincs_l1000_signed_term_gene.tsv"
    write_tsv(signed_path, signed_rows, ["term", "gene_id", "gene_symbol", "score", "sign"])

    write_workflow_provenance_graph(
        workflow_name="lincs_l1000_crisprko",
        module_name="geneset_extractors.workflows.lincs_l1000_crisprko",
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
            "top_n": top_n,
            "min_gmt_size": min_gmt_size,
            "n_rows": len(processed_rows),
            "n_terms": int(l1000["Gene KO"].nunique()) if not l1000.empty else 0,
            "n_genes": int(l1000["Gene"].nunique()) if not l1000.empty else 0,
        },
    )
    return {"n_rows": len(processed_rows), "out_dir": str(out_dir)}
