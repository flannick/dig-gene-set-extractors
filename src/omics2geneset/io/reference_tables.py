from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path
from collections import defaultdict

from omics2geneset.core.atac_programs import ATLAS_SCORE_DEFINITION_DEFAULT


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _as_float(value: object) -> float | None:
    text = str(value).strip()
    if not text:
        return None
    lowered = text.lower()
    if lowered in {"na", "nan", "none", "null"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _median(values: list[float]) -> float:
    if not values:
        raise ValueError("cannot compute median of empty list")
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2 == 1:
        return float(ordered[mid])
    return float((ordered[mid - 1] + ordered[mid]) / 2.0)


def _infer_ccre_bed_path(ubiquity_path: Path) -> Path:
    stem = ubiquity_path.name.lower()
    candidates: list[Path] = []
    if "hg19" in stem or "grch37" in stem:
        candidates.extend(
            [
                ubiquity_path.with_name("encode_ccre_hg19.bed.gz"),
                ubiquity_path.with_name("encode_ccre_hg19.bed"),
            ]
        )
    if "hg38" in stem or "grch38" in stem:
        candidates.extend(
            [
                ubiquity_path.with_name("encode_ccre_hg38.bed.gz"),
                ubiquity_path.with_name("encode_ccre_hg38.bed"),
            ]
        )
    if "mm10" in stem or "grcm38" in stem:
        candidates.extend(
            [
                ubiquity_path.with_name("encode_ccre_mm10.bed.gz"),
                ubiquity_path.with_name("encode_ccre_mm10.bed"),
            ]
        )
    candidates.extend(
        [
            ubiquity_path.with_name("encode_ccre.bed.gz"),
            ubiquity_path.with_name("encode_ccre.bed"),
        ]
    )
    for path in candidates:
        if path.exists():
            return path
    raise ValueError(
        "ref ubiquity ccre_id schema requires a sibling cCRE BED file "
        "(expected encode_ccre_hg19.bed.gz, encode_ccre_hg38.bed.gz, or encode_ccre_mm10.bed.gz)"
    )


def _load_ccre_intervals_by_id(ccre_bed_path: Path) -> dict[str, tuple[str, int, int]]:
    out: dict[str, tuple[str, int, int]] = {}
    with _open_text(ccre_bed_path) as fh:
        for line in fh:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            cols = text.split("\t")
            if len(cols) < 5:
                continue
            chrom = cols[0].strip()
            ccre_id = cols[4].strip() or cols[3].strip()
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                continue
            if not chrom or not ccre_id:
                continue
            out[ccre_id] = (chrom, start, end)
    if not out:
        raise ValueError(f"cCRE BED had no usable rows: {ccre_bed_path}")
    return out


def _row_value(row: dict[str, str], fields: dict[str, str], key: str) -> str:
    col = fields.get(key)
    if not col:
        return ""
    return str(row.get(col, ""))


def read_ref_ubiquity_tsv(path: str | Path) -> list[dict[str, object]]:
    p = Path(path)
    rows: list[dict[str, object]] = []
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("ref ubiquity table is empty")
        fields = {f.strip().lower(): f.strip() for f in reader.fieldnames if f and f.strip()}
        has_coords = {"chrom", "start", "end"}.issubset(set(fields))
        has_ccre_id = "ccre_id" in fields
        has_idf = ("idf_ref" in fields) or ("idf" in fields)
        has_df_n = ("df_ref" in fields and "n_ref" in fields) or ("df" in fields and "n" in fields)

        if has_coords and (has_idf or has_df_n):
            for row in reader:
                chrom = _row_value(row, fields, "chrom").strip()
                if not chrom:
                    continue
                try:
                    start = int(_row_value(row, fields, "start"))
                    end = int(_row_value(row, fields, "end"))
                except ValueError:
                    continue
                idf = _as_float(_row_value(row, fields, "idf_ref"))
                if idf is None:
                    idf = _as_float(_row_value(row, fields, "idf"))
                if idf is None:
                    df_ref = _as_float(_row_value(row, fields, "df_ref"))
                    n_ref = _as_float(_row_value(row, fields, "n_ref"))
                    if df_ref is None or n_ref is None:
                        df_ref = _as_float(_row_value(row, fields, "df"))
                        n_ref = _as_float(_row_value(row, fields, "n"))
                    if df_ref is None or n_ref is None:
                        continue
                    idf = math.log((n_ref + 1.0) / (df_ref + 1.0)) + 1.0
                rows.append({"chrom": chrom, "start": start, "end": end, "idf_ref": float(idf)})
            if not rows:
                raise ValueError("ref ubiquity table had coordinate schema but no usable rows")
            return rows

        if has_ccre_id and (has_idf or has_df_n):
            ccre_bed = _infer_ccre_bed_path(p)
            ccre_intervals = _load_ccre_intervals_by_id(ccre_bed)
            for row in reader:
                ccre_id = _row_value(row, fields, "ccre_id").strip()
                if not ccre_id:
                    continue
                interval = ccre_intervals.get(ccre_id)
                if interval is None:
                    continue
                idf = _as_float(_row_value(row, fields, "idf_ref"))
                if idf is None:
                    idf = _as_float(_row_value(row, fields, "idf"))
                if idf is None:
                    df_ref = _as_float(_row_value(row, fields, "df_ref"))
                    n_ref = _as_float(_row_value(row, fields, "n_ref"))
                    if df_ref is None or n_ref is None:
                        df_ref = _as_float(_row_value(row, fields, "df"))
                        n_ref = _as_float(_row_value(row, fields, "n"))
                    if df_ref is None or n_ref is None:
                        continue
                    idf = math.log((n_ref + 1.0) / (df_ref + 1.0)) + 1.0
                chrom, start, end = interval
                rows.append({"chrom": chrom, "start": start, "end": end, "idf_ref": float(idf)})
            if not rows:
                raise ValueError(
                    "ref ubiquity table had ccre_id schema but no usable rows "
                    "(check placeholder NA values and cCRE BED compatibility)"
                )
            return rows

        raise ValueError(
            "ref ubiquity table must be either: "
            "(chrom,start,end with idf_ref or df_ref/n_ref) OR "
            "(ccre_id with idf/idf_ref or df/n plus sibling cCRE BED)"
        )


def read_atlas_gene_stats_tsv(path: str | Path) -> dict[str, tuple[float, float]]:
    by_score_definition = read_atlas_gene_stats_by_score_definition_tsv(path)
    if ATLAS_SCORE_DEFINITION_DEFAULT in by_score_definition:
        return by_score_definition[ATLAS_SCORE_DEFINITION_DEFAULT]
    if len(by_score_definition) == 1:
        only_key = next(iter(by_score_definition.keys()))
        return by_score_definition[only_key]
    raise ValueError(
        "atlas stats table contains multiple score_definition values; "
        "use read_atlas_gene_stats_by_score_definition_tsv()"
    )


def read_atlas_gene_stats_by_score_definition_tsv(path: str | Path) -> dict[str, dict[str, tuple[float, float]]]:
    p = Path(path)
    out: dict[str, dict[str, tuple[float, float]]] = {}
    per_key_gene_scores: dict[tuple[str, str], list[float]] = defaultdict(list)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("atlas gene stats table is empty")
        fields = {f.strip().lower(): f.strip() for f in reader.fieldnames if f and f.strip()}
        if "gene_id" not in fields:
            raise ValueError("atlas gene stats table must include gene_id")
        score_def_col = fields.get("score_definition") or fields.get("score_def")
        median_col = fields.get("median_score") or fields.get("median")
        mad_col = fields.get("mad_score") or fields.get("mad")
        if median_col and mad_col:
            for row in reader:
                gene_id = _row_value(row, fields, "gene_id").strip()
                if not gene_id:
                    continue
                median = _as_float(row.get(median_col, ""))
                mad = _as_float(row.get(mad_col, ""))
                if median is None or mad is None:
                    continue
                score_def = (
                    str(row.get(score_def_col, "")).strip()
                    if score_def_col
                    else ATLAS_SCORE_DEFINITION_DEFAULT
                )
                if not score_def:
                    score_def = ATLAS_SCORE_DEFINITION_DEFAULT
                out.setdefault(score_def, {})[gene_id] = (float(median), float(mad))
            if not out:
                raise ValueError("atlas gene stats table had median/mad schema but no usable rows")
            return out

        score_col = fields.get("score")
        if score_col:
            for row in reader:
                gene_id = _row_value(row, fields, "gene_id").strip()
                if not gene_id:
                    continue
                score = _as_float(row.get(score_col, ""))
                if score is None:
                    continue
                score_def = (
                    str(row.get(score_def_col, "")).strip()
                    if score_def_col
                    else ATLAS_SCORE_DEFINITION_DEFAULT
                )
                if not score_def:
                    score_def = ATLAS_SCORE_DEFINITION_DEFAULT
                per_key_gene_scores[(score_def, gene_id)].append(float(score))
            for (score_def, gene_id), values in per_key_gene_scores.items():
                median = _median(values)
                deviations = [abs(v - median) for v in values]
                mad = _median(deviations)
                out.setdefault(score_def, {})[gene_id] = (float(median), float(mad))
            if not out:
                raise ValueError("atlas reference profile table had score schema but no usable rows")
            return out

        raise ValueError(
            "atlas table must include either median_score/mad_score (or median/mad), "
            "or long-format columns gene_id and score"
        )
    return out
