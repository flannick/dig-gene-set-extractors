from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from geneset_extractors.preprocessing.rnaseq.de_io import MatrixData


@dataclass(frozen=True)
class PseudobulkResult:
    matrix: MatrixData
    metadata_rows: list[dict[str, str]]
    n_cells_input: int
    n_cells_used: int
    n_pseudobulks: int
    n_dropped_small_groups: int



def _clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()



def build_pseudobulk(
    matrix: MatrixData,
    *,
    cell_metadata_rows: list[dict[str, str]],
    cell_id_column: str,
    donor_column: str,
    group_columns: list[str],
    min_cells_per_pseudobulk: int,
    pseudobulk_id_column: str = "sample_id",
) -> PseudobulkResult:
    meta_by_cell = {str(row.get(cell_id_column, "")).strip(): row for row in cell_metadata_rows if str(row.get(cell_id_column, "")).strip()}
    sample_index = {sample_id: i for i, sample_id in enumerate(matrix.sample_ids)}
    grouped: dict[tuple[tuple[str, str], ...], list[str]] = {}
    for cell_id in matrix.sample_ids:
        meta = meta_by_cell.get(cell_id)
        if meta is None:
            continue
        donor = _clean(meta.get(donor_column))
        if not donor:
            continue
        key_parts = [(donor_column, donor)]
        for col in group_columns:
            key_parts.append((col, _clean(meta.get(col))))
        grouped.setdefault(tuple(key_parts), []).append(cell_id)

    kept_keys: list[tuple[tuple[str, str], ...]] = []
    metadata_rows: list[dict[str, str]] = []
    matrix_rows: list[np.ndarray] = []
    n_cells_used = 0
    n_dropped_small = 0

    for key in sorted(grouped):
        cell_ids = grouped[key]
        if len(cell_ids) < int(min_cells_per_pseudobulk):
            n_dropped_small += 1
            continue
        kept_keys.append(key)
        idx = np.asarray([sample_index[cell_id] for cell_id in cell_ids], dtype=int)
        matrix_rows.append(np.asarray(matrix.counts[idx, :], dtype=float).sum(axis=0))
        n_cells_used += len(cell_ids)
        meta_row = {pseudobulk_id_column: "__".join(f"{k}={v}" for k, v in key)}
        for col, value in key:
            meta_row[col] = value
        meta_row["n_cells"] = str(len(cell_ids))
        metadata_rows.append(meta_row)

    if not matrix_rows:
        raise ValueError(
            "No pseudobulks remained after aggregation. Lower --min_cells_per_pseudobulk or check donor/cell metadata."
        )

    counts = np.vstack(matrix_rows)
    pb_matrix = MatrixData(
        sample_ids=[row[pseudobulk_id_column] for row in metadata_rows],
        feature_ids=list(matrix.feature_ids),
        counts=counts,
        gene_symbols=list(matrix.gene_symbols),
        matrix_orientation="sample_by_gene",
        sample_id_column=pseudobulk_id_column,
        feature_id_column=matrix.feature_id_column,
        source_path=matrix.source_path,
    )
    return PseudobulkResult(
        matrix=pb_matrix,
        metadata_rows=metadata_rows,
        n_cells_input=len(matrix.sample_ids),
        n_cells_used=n_cells_used,
        n_pseudobulks=len(metadata_rows),
        n_dropped_small_groups=n_dropped_small,
    )
