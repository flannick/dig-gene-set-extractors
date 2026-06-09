from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

from geneset_extractors.workflows.gtex_runtime_common import write_tsv, write_workflow_provenance_graph


def _require_dir(path: Path, label: str) -> None:
    if not path.is_dir():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _detect_asctb_header_row(path: Path, max_lines: int = 50) -> int:
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for i, line in enumerate(handle):
            if i > max_lines:
                break
            if "CT/1" in line and ("AS/1" in line or "BGene" in line or "Gene" in line or "Protein" in line):
                return i
    return 10


def _prepare_raw_asctb_tables(raw_asctb_dir: Path, prepared_asctb_dir: Path) -> list[Path]:
    prepared_asctb_dir.mkdir(parents=True, exist_ok=True)
    prepared_paths: list[Path] = []
    for path in sorted(raw_asctb_dir.iterdir()):
        if not path.is_file() or path.suffix.lower() not in {".csv", ".tsv"}:
            continue
        header = _detect_asctb_header_row(path)
        sep = "\t" if path.suffix.lower() == ".tsv" else ","
        frame = pd.read_csv(path, header=header, sep=sep, dtype=str)
        out_path = prepared_asctb_dir / (f"{path.stem}.csv" if path.suffix.lower() == ".tsv" else path.name)
        frame.to_csv(out_path)
        prepared_paths.append(out_path)
    if not prepared_paths:
        raise RuntimeError(f"No raw ASCT+B tables found in {raw_asctb_dir}")
    return prepared_paths


def _load_asctb_tables(asctb_dir: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for table in tqdm(sorted(asctb_dir.iterdir()), desc="Combining ASCT+B tables"):
        if not table.is_file():
            continue
        frames.append(pd.read_csv(table, index_col=0))
    if not frames:
        raise RuntimeError(f"No prepared ASCT+B tables found in {asctb_dir}")
    return pd.concat(frames).reset_index(drop=True)


def _get_highest_resolution_cell_type(entry: pd.Series):
    cts = entry[["CT/1", "CT/2", "CT/3", "CT/4"]].dropna()
    return np.nan if len(cts) == 0 else (cts.index[-1], re.sub(r"(\w+)s$", r"\1", cts.iloc[-1]))


def _get_highest_resolution_cell_type_id(entry: pd.Series):
    ctids = entry[["CT/1/ID", "CT/2/ID", "CT/3/ID", "CT/4/ID"]].dropna()
    return np.nan if len(ctids) == 0 else ctids.iloc[-1]


def _add_labels(asctb: pd.DataFrame) -> pd.DataFrame:
    asctb = asctb.copy()
    asctb["Label"] = asctb.apply(_get_highest_resolution_cell_type, axis=1)
    asctb["CTID"] = asctb.apply(_get_highest_resolution_cell_type_id, axis=1)
    asctb = asctb.dropna(subset="Label")
    asctb["Label"] = (
        asctb["AS/1"].apply(str.capitalize)
        + "_"
        + asctb["Label"].apply(lambda x: x[0].replace("/", ""))
        + "_"
        + asctb["Label"].apply(lambda x: x[1].replace("_", " "))
    )
    asctb["Label"] = asctb["Label"].apply(lambda x: re.sub(r"\s\([^)]+\)", "", x))
    return asctb


def _get_marker_columns(asctb: pd.DataFrame):
    return asctb.columns[
        asctb.columns.map(
            lambda x: ("Gene" in x or "Protein" in x)
            and "LABEL" not in x
            and "ID" not in x
            and "ABBR" not in x
            and "NOTE" not in x
        )
    ]


def _get_all_genes_from_cols(entry: pd.Series, marker_cols) -> object:
    genes = set()
    entry_gene_cols = entry[marker_cols].dropna()
    for gene_col in entry_gene_cols:
        genes.update(set(str(gene_col).split(", ")))
    if len(genes) == 0:
        return np.nan
    return genes


def _add_raw_genes(asctb: pd.DataFrame) -> pd.DataFrame:
    asctb = asctb.copy()
    marker_cols = _get_marker_columns(asctb)
    asctb["Genes"] = asctb.apply(lambda row: _get_all_genes_from_cols(row, marker_cols), axis=1)
    asctb = asctb.dropna(subset="Genes")
    return asctb


def _load_geneinfo(human_gene_info: Path) -> pd.DataFrame:
    geneinfo = pd.read_csv(human_gene_info, sep="\t")
    geneinfo = geneinfo[geneinfo["#tax_id"] == 9606][geneinfo["type_of_gene"] == "protein-coding"].copy()
    geneinfo["Synonyms"] = geneinfo["Synonyms"].apply(str.split, sep="|")
    return geneinfo.explode("Synonyms")[["GeneID", "Symbol", "Synonyms", "description"]]


def _make_gene_mappers(geneinfo: pd.DataFrame):
    exact_map: dict[str, str] = {}
    synonym_map: dict[str, str] = {}
    for _, row in geneinfo.iterrows():
        symbol = str(row["Symbol"]).strip()
        if not symbol or symbol.lower() == "nan":
            continue
        exact_map[symbol.upper()] = symbol
        syn = str(row["Synonyms"]).strip()
        if not syn or syn == "-" or syn.lower() == "nan":
            continue
        synonym_map.setdefault(syn.upper(), symbol)
    return exact_map, synonym_map


def _clean_marker(gene_label: str) -> str:
    gene_label = str(gene_label).split(",")[0].strip()
    gene_label = re.sub(r"[+-]$", "", gene_label)
    gene_label = re.sub(r"\s\([^)]+\)", "", gene_label)
    return gene_label


def _resolve_gene_symbol(cleaned_marker: str, exact_map: dict[str, str], synonym_map: dict[str, str]):
    key = str(cleaned_marker).strip().upper()
    if not key:
        return None, "missing"
    if key in exact_map:
        return exact_map[key], "exact_symbol"
    if key in synonym_map:
        return synonym_map[key], "synonym"
    return None, "unmapped"


def _map_and_filter_genes(asctb: pd.DataFrame, exact_map: dict[str, str], synonym_map: dict[str, str]) -> pd.DataFrame:
    asctb = asctb.copy().explode("Genes")
    asctb["Raw Marker"] = asctb["Genes"]
    asctb["Clean Marker"] = asctb["Raw Marker"].map(_clean_marker)
    resolved = asctb["Clean Marker"].map(lambda x: _resolve_gene_symbol(x, exact_map, synonym_map))
    asctb["Genes"] = resolved.map(lambda x: x[0])
    asctb = asctb.dropna(subset="Genes").drop_duplicates(subset=["Label", "Genes"])
    return asctb[["Label", "CTID", "Genes"]].dropna().reset_index(drop=True)


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    human_gene_info = Path(args.human_gene_info).resolve()
    _require_file(human_gene_info, "human_gene_info")

    input_paths: list[tuple[Path, str]] = [(human_gene_info, "human_gene_info")]
    prepared_tables: list[Path] = []

    if getattr(args, "raw_asctb_dir", None):
        raw_asctb_dir = Path(args.raw_asctb_dir).resolve()
        _require_dir(raw_asctb_dir, "raw ASCT+B directory")
        input_paths.append((raw_asctb_dir, "raw_asctb_dir"))
        asctb_dir = out_dir / "ASCTB_Tables"
        prepared_tables = _prepare_raw_asctb_tables(raw_asctb_dir, asctb_dir)
    else:
        asctb_dir = Path(args.asctb_dir).resolve()
        _require_dir(asctb_dir, "ASCTB directory")
        input_paths.append((asctb_dir, "asctb_dir"))

    asctb = _load_asctb_tables(asctb_dir)
    asctb = _add_labels(asctb)
    asctb = _add_raw_genes(asctb)
    geneinfo = _load_geneinfo(human_gene_info)
    exact_map, synonym_map = _make_gene_mappers(geneinfo)
    asctb = _map_and_filter_genes(asctb, exact_map, synonym_map)

    edge_rows = [
        {
            "Gene": str(row["Genes"]),
            "Gene ID": "",
            "Cell Type": str(row["Label"]),
            "Cell Type ID": str(row["CTID"]),
            "Threshold": 1,
        }
        for _, row in asctb.iterrows()
    ]
    edge_path = out_dir / "gene_attribute_edges.txt.gz"
    pd.DataFrame(edge_rows).to_csv(edge_path, sep="\t", compression="gzip", index=False)

    unsigned_rows = [
        {
            "term": str(row["Label"]),
            "gene_id": str(row["Genes"]),
            "gene_symbol": str(row["Genes"]),
            "score": "1",
        }
        for _, row in asctb.iterrows()
    ]
    unsigned_path = out_dir / "hubmap_unsigned_term_gene.tsv"
    write_tsv(unsigned_path, unsigned_rows, ["term", "gene_id", "gene_symbol", "score"])

    write_workflow_provenance_graph(
        workflow_name="hubmap_asctb",
        module_name="geneset_extractors.workflows.hubmap_asctb",
        output_dir=out_dir,
        focus_output_path=unsigned_path,
        output_paths=[
            (unsigned_path, "unsigned_term_gene_tsv"),
            (edge_path, "gene_attribute_edges_tsv_gz"),
            *[(path, "prepared_asctb_table_csv") for path in prepared_tables],
        ],
        input_paths=input_paths,
        parameters={
            "n_terms": int(asctb["Label"].nunique()) if not asctb.empty else 0,
            "n_genes": int(asctb["Genes"].nunique()) if not asctb.empty else 0,
            "n_edges": len(unsigned_rows),
        },
    )
    return {"n_rows": len(unsigned_rows), "out_dir": str(out_dir)}
