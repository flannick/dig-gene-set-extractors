from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_ref_ubiquity_tsv(path: str | Path) -> list[dict[str, object]]:
    p = Path(path)
    rows: list[dict[str, object]] = []
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("ref ubiquity table is empty")
        fields = {f.strip() for f in reader.fieldnames if f}
        required_base = {"chrom", "start", "end"}
        if not required_base.issubset(fields):
            raise ValueError(
                "ref ubiquity table must include columns chrom,start,end and idf_ref or (df_ref,n_ref)"
            )
        has_idf = "idf_ref" in fields
        has_df = "df_ref" in fields and "n_ref" in fields
        if not has_idf and not has_df:
            raise ValueError(
                "ref ubiquity table must include idf_ref or both df_ref and n_ref columns"
            )
        for row in reader:
            chrom = str(row.get("chrom", "")).strip()
            if not chrom:
                continue
            start = int(row["start"])
            end = int(row["end"])
            if has_idf and str(row.get("idf_ref", "")).strip():
                idf = float(row["idf_ref"])
            else:
                df_ref = float(row.get("df_ref", 0.0))
                n_ref = float(row.get("n_ref", 0.0))
                idf = math.log((n_ref + 1.0) / (df_ref + 1.0)) + 1.0
            rows.append({"chrom": chrom, "start": start, "end": end, "idf_ref": float(idf)})
    return rows


def read_atlas_gene_stats_tsv(path: str | Path) -> dict[str, tuple[float, float]]:
    p = Path(path)
    out: dict[str, tuple[float, float]] = {}
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("atlas gene stats table is empty")
        fields = {f.strip() for f in reader.fieldnames if f}
        if "gene_id" not in fields:
            raise ValueError("atlas gene stats table must include gene_id")
        median_col = "median_score" if "median_score" in fields else ("median" if "median" in fields else None)
        mad_col = "mad_score" if "mad_score" in fields else ("mad" if "mad" in fields else None)
        if median_col is None or mad_col is None:
            raise ValueError("atlas gene stats table must include median_score and mad_score (or median and mad)")
        for row in reader:
            gene_id = str(row.get("gene_id", "")).strip()
            if not gene_id:
                continue
            out[gene_id] = (float(row[median_col]), float(row[mad_col]))
    return out
