from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class MatrixData:
    sample_ids: list[str]
    feature_ids: list[str]
    counts: np.ndarray  # samples x features
    gene_symbols: list[str]
    matrix_orientation: str
    sample_id_column: str
    feature_id_column: str
    source_path: str


@dataclass(frozen=True)
class TableData:
    fieldnames: list[str]
    rows: list[dict[str, str]]


class MetadataJoinError(ValueError):
    pass



def _clean(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()



def _parse_float(raw: Any) -> float:
    text = _clean(raw)
    if not text:
        return 0.0
    try:
        value = float(text)
    except ValueError:
        return 0.0
    if not np.isfinite(value):
        return 0.0
    return float(value)



def read_tsv_rows(path: str | Path, *, delimiter: str = "\t") -> TableData:
    rows: list[dict[str, str]] = []
    with Path(path).open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"TSV has no header: {path}")
        fieldnames = [str(x) for x in reader.fieldnames]
        for row in reader:
            rows.append({str(k): _clean(v) for k, v in row.items() if k is not None})
    return TableData(fieldnames=fieldnames, rows=rows)



def write_tsv(path: str | Path, rows: list[dict[str, Any]], *, fieldnames: list[str] | None = None) -> None:
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        keyset: set[str] = set()
        for row in rows:
            keyset.update(str(k) for k in row.keys())
        fieldnames = sorted(keyset)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})



def read_matrix_tsv(
    path: str | Path,
    *,
    orientation: str,
    sample_id_column: str,
    feature_id_column: str,
    gene_symbol_column: str | None = None,
    delimiter: str = "\t",
) -> MatrixData:
    matrix_path = Path(path)
    with matrix_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Matrix TSV is empty: {path}") from exc
        header = [str(x) for x in header]

        if orientation == "sample_by_gene":
            if sample_id_column not in header:
                raise ValueError(
                    f"sample_id_column '{sample_id_column}' not found in matrix header. "
                    f"Available columns: {', '.join(header[:20])}"
                )
            sample_idx = header.index(sample_id_column)
            feature_ids = [col for i, col in enumerate(header) if i != sample_idx]
            if not feature_ids:
                raise ValueError("No feature columns remained after removing sample_id_column.")
            gene_symbols = list(feature_ids)
            sample_ids: list[str] = []
            rows: list[list[float]] = []
            for row in reader:
                if not row:
                    continue
                sample_id = _clean(row[sample_idx] if sample_idx < len(row) else "")
                if not sample_id:
                    continue
                sample_ids.append(sample_id)
                values: list[float] = []
                for i in range(len(header)):
                    if i == sample_idx:
                        continue
                    raw = row[i] if i < len(row) else ""
                    values.append(_parse_float(raw))
                rows.append(values)
            if not sample_ids:
                raise ValueError(f"No samples were parsed from matrix: {path}")
            counts = np.asarray(rows, dtype=float)
            return MatrixData(
                sample_ids=sample_ids,
                feature_ids=feature_ids,
                counts=counts,
                gene_symbols=gene_symbols,
                matrix_orientation=orientation,
                sample_id_column=sample_id_column,
                feature_id_column=feature_id_column,
                source_path=str(matrix_path),
            )

        if orientation != "gene_by_sample":
            raise ValueError(f"Unsupported matrix orientation: {orientation}")

        if feature_id_column not in header:
            raise ValueError(
                f"feature_id_column '{feature_id_column}' not found in matrix header. "
                f"Available columns: {', '.join(header[:20])}"
            )
        feature_idx = header.index(feature_id_column)
        symbol_idx = header.index(gene_symbol_column) if gene_symbol_column and gene_symbol_column in header else None
        sample_ids = [
            col
            for i, col in enumerate(header)
            if i not in {feature_idx} and i != symbol_idx
        ]
        if not sample_ids:
            raise ValueError("No sample columns remained after removing feature metadata columns.")
        feature_ids: list[str] = []
        gene_symbols: list[str] = []
        gene_rows: list[list[float]] = []
        for row in reader:
            if not row:
                continue
            feature_id = _clean(row[feature_idx] if feature_idx < len(row) else "")
            if not feature_id:
                continue
            feature_ids.append(feature_id)
            gene_symbols.append(_clean(row[symbol_idx]) if symbol_idx is not None and symbol_idx < len(row) else feature_id)
            values: list[float] = []
            for i in range(len(header)):
                if i == feature_idx or i == symbol_idx:
                    continue
                raw = row[i] if i < len(row) else ""
                values.append(_parse_float(raw))
            gene_rows.append(values)
        if not feature_ids:
            raise ValueError(f"No features were parsed from matrix: {path}")
        counts = np.asarray(gene_rows, dtype=float).T
        return MatrixData(
            sample_ids=sample_ids,
            feature_ids=feature_ids,
            counts=counts,
            gene_symbols=gene_symbols,
            matrix_orientation=orientation,
            sample_id_column=sample_id_column,
            feature_id_column=feature_id_column,
            source_path=str(matrix_path),
        )



def table_index(rows: list[dict[str, str]], key_column: str) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        key = _clean(row.get(key_column))
        if not key:
            continue
        out[key] = row
    return out



def join_metadata(
    primary_rows: list[dict[str, str]],
    *,
    primary_key: str,
    secondary_rows: list[dict[str, str]] | None,
    secondary_key: str | None,
    suffix: str = "",
) -> tuple[list[dict[str, str]], dict[str, int]]:
    if not secondary_rows:
        return [dict(row) for row in primary_rows], {"n_joined": 0, "n_missing_join": 0}
    if not secondary_key:
        raise MetadataJoinError("secondary_key is required when secondary metadata is provided")
    idx = table_index(secondary_rows, secondary_key)
    joined: list[dict[str, str]] = []
    n_joined = 0
    n_missing = 0
    for row in primary_rows:
        out = dict(row)
        key = _clean(row.get(primary_key))
        extra = idx.get(key)
        if extra is None:
            n_missing += 1
            joined.append(out)
            continue
        n_joined += 1
        for col, value in extra.items():
            if col == secondary_key:
                continue
            target = col if col not in out else f"{col}{suffix}"
            out[target] = _clean(value)
        joined.append(out)
    return joined, {"n_joined": n_joined, "n_missing_join": n_missing}



def subset_matrix(matrix: MatrixData, sample_ids: list[str]) -> MatrixData:
    idx_by_sample = {sample_id: i for i, sample_id in enumerate(matrix.sample_ids)}
    keep_idx = [idx_by_sample[sid] for sid in sample_ids if sid in idx_by_sample]
    counts = matrix.counts[np.asarray(keep_idx, dtype=int), :]
    keep_samples = [matrix.sample_ids[i] for i in keep_idx]
    return MatrixData(
        sample_ids=keep_samples,
        feature_ids=list(matrix.feature_ids),
        counts=counts,
        gene_symbols=list(matrix.gene_symbols),
        matrix_orientation=matrix.matrix_orientation,
        sample_id_column=matrix.sample_id_column,
        feature_id_column=matrix.feature_id_column,
        source_path=matrix.source_path,
    )



def write_matrix_tsv(path: str | Path, matrix: MatrixData) -> None:
    rows: list[dict[str, Any]] = []
    for sample_id, values in zip(matrix.sample_ids, matrix.counts, strict=False):
        row: dict[str, Any] = {matrix.sample_id_column: sample_id}
        for feature_id, value in zip(matrix.feature_ids, values, strict=False):
            row[feature_id] = float(value)
        rows.append(row)
    fieldnames = [matrix.sample_id_column] + list(matrix.feature_ids)
    write_tsv(path, rows, fieldnames=fieldnames)
