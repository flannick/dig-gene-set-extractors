from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

from geneset_extractors.workflows.gtex_runtime_common import write_tsv, write_workflow_provenance_graph


GENESHOT_URL_DEFAULT = "https://maayanlab.cloud/geneshot/api/associate"


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _load_human_gene_info(path: Path):
    df = pd.read_csv(path, sep="\t", dtype={"#tax_id": str, "tax_id": str})
    tax_col = "#tax_id" if "#tax_id" in df.columns else "tax_id"
    if tax_col in df.columns:
        df = df[df[tax_col].astype(str) == "9606"].copy()
    df["Symbol"] = df["Symbol"].map(str.upper)
    keep_cols = [c for c in ["Symbol", "GeneID", "description", "type_of_gene"] if c in df.columns]
    human_gene_info = df[keep_cols].set_index("Symbol")
    return df, human_gene_info


def _load_base_matrix(matrix_path: Path):
    matrix = pd.read_csv(matrix_path, sep="\t", compression="infer", index_col=0)
    if matrix.index.name is None:
        matrix.index.name = "Gene"
    asctb_full_gmt = (
        matrix.stack()
        .reset_index()
        .replace(0, pd.NA)
        .dropna()
        .groupby("level_1")["Gene"]
        .agg(list)
        .to_dict()
    )
    return matrix, asctb_full_gmt


def _query_geneshot_for_term(
    gene_list: list[str],
    *,
    url: str,
    timeout: int,
    retries: int,
    pause_seconds: float,
) -> dict[str, float]:
    payload = {"gene_list": list(gene_list), "similarity": "coexpression"}
    last_err: Optional[Exception] = None
    for attempt in range(retries + 1):
        try:
            response = requests.post(url, json=payload, timeout=timeout)
            response.raise_for_status()
            data = json.loads(response.text)
            sim_genes = pd.DataFrame(data["association"]).T.sort_values("simScore", ascending=False)
            if pause_seconds:
                time.sleep(pause_seconds)
            return sim_genes["simScore"].to_dict()
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            if attempt < retries:
                time.sleep(max(pause_seconds, 1.0) * (attempt + 1))
            else:
                raise RuntimeError("Geneshot query failed") from last_err
    raise AssertionError("unreachable")


def _build_or_load_geneshot_augmented(
    asctb_full_gmt: dict[str, list[str]],
    augmented_tsv: Path,
    *,
    geneshot_url: str,
    limit_terms: int | None,
    timeout: int,
    retries: int,
    pause_seconds: float,
) -> pd.DataFrame:
    if augmented_tsv.exists():
        return pd.read_csv(augmented_tsv, sep="\t", index_col="Gene")
    terms = list(asctb_full_gmt)
    if limit_terms is not None:
        terms = terms[:limit_terms]
    augment_gmt: dict[str, dict[str, float]] = {}
    for term in terms:
        augment_gmt[term] = _query_geneshot_for_term(
            list(asctb_full_gmt[term]),
            url=geneshot_url,
            timeout=timeout,
            retries=retries,
            pause_seconds=pause_seconds,
        )
    asctb_augmented = pd.DataFrame(index=augment_gmt.keys(), data=augment_gmt.values()).T.sort_index().rename_axis("Gene")
    asctb_augmented.to_csv(augmented_tsv, sep="\t")
    return asctb_augmented


def _create_augmented_gmt(
    asctb_full_gmt: dict[str, list[str]],
    asctb_augmented: pd.DataFrame,
    *,
    threshold: float,
    cap_multiplier: int,
) -> dict[str, set]:
    aug: dict[str, list[str]] = {}
    for cell in asctb_full_gmt:
        if cell not in asctb_augmented.columns:
            aug[cell] = []
            continue
        aug[cell] = (
            asctb_augmented[cell][asctb_augmented[cell] >= threshold]
            .dropna()
            .sort_values(ascending=False)
            .index[: cap_multiplier * len(asctb_full_gmt[cell])]
            .to_list()
        )
    augmented_gmt: dict[str, set] = {}
    for term in asctb_full_gmt:
        base = set(asctb_full_gmt[term])
        base.update(set(aug[term]))
        augmented_gmt[term] = base
    return augmented_gmt


def _create_asctbaug_dataframe(augmented_gmt: dict[str, set]) -> pd.DataFrame:
    frame = pd.Series(index=augmented_gmt.keys(), data=augmented_gmt.values()).explode().reset_index()
    frame.columns = ["Cell Type", "Gene"]
    return frame


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    input_matrix = Path(args.input_matrix).resolve()
    human_gene_info = Path(args.human_gene_info).resolve()
    _require_file(input_matrix, "input matrix")
    _require_file(human_gene_info, "human_gene_info")

    _raw_human_gene_info, human_gene_info_df = _load_human_gene_info(human_gene_info)
    _matrix, asctb_full_gmt = _load_base_matrix(input_matrix)

    augmented_tsv = out_dir / "asctb_geneshot_augmented_genes.tsv"
    asctb_augmented = _build_or_load_geneshot_augmented(
        asctb_full_gmt,
        augmented_tsv,
        geneshot_url=str(getattr(args, "geneshot_url", GENESHOT_URL_DEFAULT)),
        limit_terms=(int(args.limit_terms) if getattr(args, "limit_terms", None) else None),
        timeout=int(getattr(args, "request_timeout", 120)),
        retries=int(getattr(args, "request_retries", 2)),
        pause_seconds=float(getattr(args, "pause_seconds", 0.1)),
    )
    asctb_augmented = asctb_augmented[asctb_augmented.index.isin(human_gene_info_df.index)]
    augmented_gmt = _create_augmented_gmt(
        asctb_full_gmt,
        asctb_augmented,
        threshold=float(args.augmentation_threshold),
        cap_multiplier=int(args.cap_multiplier),
    )
    asctbaug = _create_asctbaug_dataframe(augmented_gmt)
    unsigned_rows = [
        {
            "term": str(row["Cell Type"]),
            "gene_id": str(row["Gene"]),
            "gene_symbol": str(row["Gene"]),
            "score": "1",
        }
        for _, row in asctbaug.iterrows()
    ]
    unsigned_path = out_dir / "hubmap_unsigned_term_gene.tsv"
    write_tsv(unsigned_path, unsigned_rows, ["term", "gene_id", "gene_symbol", "score"])

    write_workflow_provenance_graph(
        workflow_name="hubmap_asctb_augmented",
        module_name="geneset_extractors.workflows.hubmap_asctb_augmented",
        output_dir=out_dir,
        focus_output_path=unsigned_path,
        output_paths=[
            (unsigned_path, "unsigned_term_gene_tsv"),
            (augmented_tsv, "asctb_geneshot_augmented_genes_tsv"),
        ],
        input_paths=[
            (input_matrix, "input_matrix"),
            (human_gene_info, "human_gene_info"),
        ],
        parameters={
            "augmentation_threshold": float(args.augmentation_threshold),
            "cap_multiplier": int(args.cap_multiplier),
            "geneshot_url": str(getattr(args, "geneshot_url", GENESHOT_URL_DEFAULT)),
            "n_terms": len(augmented_gmt),
            "n_rows": len(unsigned_rows),
        },
    )
    return {"n_rows": len(unsigned_rows), "out_dir": str(out_dir)}
