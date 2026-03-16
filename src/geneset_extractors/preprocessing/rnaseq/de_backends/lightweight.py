from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np

from geneset_extractors.preprocessing.rnaseq.de_design import ComparisonSpec


@dataclass(frozen=True)
class LightweightContrastResult:
    rows: list[dict[str, Any]]
    summary: dict[str, Any]



def bh_adjust(pvalues: list[float | None]) -> list[float | None]:
    indexed = [(idx, float(p)) for idx, p in enumerate(pvalues) if p is not None and math.isfinite(float(p))]
    if not indexed:
        return [None for _ in pvalues]
    indexed.sort(key=lambda item: item[1])
    n = len(indexed)
    adjusted = [None for _ in pvalues]
    running = 1.0
    for rank in range(n - 1, -1, -1):
        idx, p = indexed[rank]
        value = min(running, (p * n) / float(rank + 1))
        running = value
        adjusted[idx] = min(1.0, max(0.0, value))
    return adjusted



def _normal_two_sided_pvalue(stat: float) -> float:
    return float(math.erfc(abs(float(stat)) / math.sqrt(2.0)))



def _safe_float(value: object) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        out = float(text)
    except ValueError:
        return None
    if not math.isfinite(out):
        return None
    return float(out)



def _is_numeric_column(values: list[str]) -> bool:
    seen = False
    for value in values:
        parsed = _safe_float(value)
        if parsed is None:
            continue
        seen = True
    return seen



def _dummy_encode(values: list[str]) -> tuple[np.ndarray, list[str]]:
    levels = [level for level in sorted({value for value in values if value})]
    if len(levels) <= 1:
        return np.zeros((len(values), 0), dtype=float), []
    base = levels[0]
    columns = [level for level in levels if level != base]
    matrix = np.zeros((len(values), len(columns)), dtype=float)
    for i, value in enumerate(values):
        for j, level in enumerate(columns):
            matrix[i, j] = 1.0 if value == level else 0.0
    return matrix, [f"{base}->{level}" for level in columns]



def _build_design_matrix(
    metadata_rows: list[dict[str, str]],
    *,
    sample_ids: list[str],
    spec: ComparisonSpec,
    covariates: list[str],
    batch_columns: list[str],
) -> tuple[np.ndarray, int, list[str]]:
    row_by_sample = {str(row.get("sample_id", "")).strip(): row for row in metadata_rows}
    rows = [row_by_sample[sid] for sid in sample_ids]
    columns = [np.ones(len(sample_ids), dtype=float)]
    column_names = ["Intercept"]
    group_indicator = np.asarray(
        [1.0 if str(row.get(spec.group_column, "")).strip() == spec.group_a else 0.0 for row in rows],
        dtype=float,
    )
    columns.append(group_indicator)
    column_names.append(f"{spec.group_column}:{spec.group_a}")

    for col in covariates + batch_columns:
        values = [str(row.get(col, "")).strip() for row in rows]
        numeric_values = [_safe_float(value) for value in values]
        if all(value is not None for value in numeric_values) and any(value is not None for value in numeric_values):
            columns.append(np.asarray([float(value or 0.0) for value in numeric_values], dtype=float))
            column_names.append(col)
            continue
        dummy_matrix, dummy_names = _dummy_encode(values)
        for j in range(dummy_matrix.shape[1]):
            columns.append(dummy_matrix[:, j])
            column_names.append(f"{col}:{dummy_names[j]}")

    design = np.column_stack(columns)
    return design, 1, column_names



def _filter_low_expression(counts: np.ndarray, *, group_mask: np.ndarray, other_mask: np.ndarray) -> np.ndarray:
    if counts.size == 0:
        return np.zeros(counts.shape[1], dtype=bool)
    lib_sizes = counts.sum(axis=1)
    lib_sizes = np.where(lib_sizes <= 0.0, 1.0, lib_sizes)
    cpm = (counts / lib_sizes[:, None]) * 1_000_000.0
    min_group = max(2, int(min(group_mask.sum(), other_mask.sum())))
    return ((cpm[group_mask, :] >= 1.0).sum(axis=0) >= min_group) | ((cpm[other_mask, :] >= 1.0).sum(axis=0) >= min_group)



def _log_cpm(counts: np.ndarray) -> np.ndarray:
    lib_sizes = counts.sum(axis=1)
    lib_sizes = np.where(lib_sizes <= 0.0, 1.0, lib_sizes)
    cpm = (counts / lib_sizes[:, None]) * 1_000_000.0
    return np.log2(cpm + 0.5)



def _ols_group_effect(y: np.ndarray, design: np.ndarray, effect_idx: int) -> tuple[float, float, float]:
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ beta
    resid = y - fitted
    n = design.shape[0]
    p = design.shape[1]
    df = max(n - p, 1)
    sigma2 = float(np.dot(resid, resid) / float(df))
    xtx_inv = np.linalg.pinv(design.T @ design)
    se = math.sqrt(max(float(xtx_inv[effect_idx, effect_idx]) * sigma2, 1e-12))
    effect = float(beta[effect_idx])
    stat = effect / se if se > 0 else 0.0
    pvalue = _normal_two_sided_pvalue(stat)
    return effect, stat, pvalue



def run_lightweight_contrast(
    *,
    counts: np.ndarray,
    feature_ids: list[str],
    gene_symbols: list[str],
    sample_ids: list[str],
    metadata_rows: list[dict[str, str]],
    spec: ComparisonSpec,
    covariates: list[str],
    batch_columns: list[str],
    low_expression_filter: bool,
) -> LightweightContrastResult:
    row_by_sample = {str(row.get("sample_id", "")).strip(): row for row in metadata_rows}
    rows = [row_by_sample[sid] for sid in sample_ids]
    group_values = np.asarray([str(row.get(spec.group_column, "")).strip() for row in rows], dtype=object)
    group_a_mask = group_values == spec.group_a
    if spec.comparison_kind == "group_vs_rest" or spec.group_b is None:
        group_b_mask = group_values != spec.group_a
        group_b_label = "rest"
    else:
        group_b_mask = group_values == spec.group_b
        group_b_label = spec.group_b
    if int(group_a_mask.sum()) < 2 or int(group_b_mask.sum()) < 2:
        raise ValueError(
            f"Comparison {spec.comparison_id} needs at least 2 samples per group after filtering; "
            f"got n_group_a={int(group_a_mask.sum())}, n_group_b={int(group_b_mask.sum())}."
        )
    keep_mask = np.asarray(group_a_mask | group_b_mask, dtype=bool)
    counts_sub = np.asarray(counts[keep_mask, :], dtype=float)
    sample_ids_sub = [sid for sid, keep in zip(sample_ids, keep_mask, strict=False) if keep]
    metadata_sub = [row for row, keep in zip(rows, keep_mask, strict=False) if keep]
    group_a_sub = np.asarray([str(row.get(spec.group_column, "")).strip() == spec.group_a for row in metadata_sub], dtype=bool)
    group_b_sub = ~group_a_sub if spec.group_b is None or spec.comparison_kind == "group_vs_rest" else np.asarray([str(row.get(spec.group_column, "")).strip() == spec.group_b for row in metadata_sub], dtype=bool)

    if low_expression_filter:
        feature_keep = _filter_low_expression(counts_sub, group_mask=group_a_sub, other_mask=group_b_sub)
    else:
        feature_keep = np.ones(counts_sub.shape[1], dtype=bool)
    if not feature_keep.any():
        raise ValueError(f"Comparison {spec.comparison_id} retained no genes after low-expression filtering")

    counts_kept = counts_sub[:, feature_keep]
    feature_ids_kept = [feature_id for feature_id, keep in zip(feature_ids, feature_keep, strict=False) if keep]
    gene_symbols_kept = [gene_symbol for gene_symbol, keep in zip(gene_symbols, feature_keep, strict=False) if keep]
    log_cpm = _log_cpm(counts_kept)
    design, effect_idx, column_names = _build_design_matrix(
        metadata_sub,
        sample_ids=sample_ids_sub,
        spec=spec,
        covariates=covariates,
        batch_columns=batch_columns,
    )

    raw_rows: list[dict[str, Any]] = []
    raw_pvalues: list[float | None] = []
    for j, gene_id in enumerate(feature_ids_kept):
        y = log_cpm[:, j]
        logfc, stat, pvalue = _ols_group_effect(y, design, effect_idx)
        mean_expr = float(y.mean())
        raw_rows.append(
            {
                "comparison_id": spec.comparison_id,
                "gene_id": gene_id,
                "gene_symbol": gene_symbols_kept[j] or gene_id,
                "logFC": logfc,
                "stat": stat,
                "pvalue": pvalue,
                "group_a": spec.group_a,
                "group_b": group_b_label,
                "stratum": "|".join(f"{k}={v}" for k, v in sorted(spec.strata.items())),
                "backend": "lightweight",
                "n_group_a": int(group_a_sub.sum()),
                "n_group_b": int(group_b_sub.sum()),
                "mean_expr": mean_expr,
                "model_formula": " + ".join(column_names),
            }
        )
        raw_pvalues.append(pvalue)

    padj_values = bh_adjust(raw_pvalues)
    rows_out: list[dict[str, Any]] = []
    for row, padj in zip(raw_rows, padj_values, strict=False):
        out = dict(row)
        out["padj"] = padj if padj is not None else ""
        rows_out.append(out)

    return LightweightContrastResult(
        rows=rows_out,
        summary={
            "comparison_id": spec.comparison_id,
            "backend": "lightweight",
            "n_group_a": int(group_a_sub.sum()),
            "n_group_b": int(group_b_sub.sum()),
            "n_features_input": counts.shape[1],
            "n_features_retained": len(feature_ids_kept),
            "low_expression_filter": bool(low_expression_filter),
        },
    )
