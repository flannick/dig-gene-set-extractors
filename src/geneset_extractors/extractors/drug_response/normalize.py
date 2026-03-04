from __future__ import annotations

from statistics import NormalDist
import statistics

from geneset_extractors.extractors.drug_response.io import ResponseRecord


_LOWER_IS_MORE_SENSITIVE = {"logfold_change", "auc", "ic50"}


def resolve_response_direction(
    *,
    response_metric: str,
    response_direction: str | None,
) -> str:
    explicit = "" if response_direction is None else str(response_direction).strip()
    if explicit:
        return explicit
    metric = str(response_metric).strip().lower()
    if metric in _LOWER_IS_MORE_SENSITIVE:
        return "lower_is_more_sensitive"
    return "higher_is_more_sensitive"


def aggregate_response_records(
    records: list[ResponseRecord],
) -> tuple[dict[str, dict[str, float]], dict[str, object]]:
    by_pair: dict[tuple[str, str], list[float]] = {}
    for rec in records:
        key = (str(rec.sample_id), str(rec.drug_id))
        by_pair.setdefault(key, []).append(float(rec.response))

    matrix: dict[str, dict[str, float]] = {}
    n_duplicates = 0
    for (sample_id, drug_id), values in by_pair.items():
        if len(values) > 1:
            n_duplicates += 1
        value = float(sum(values) / float(len(values)))
        matrix.setdefault(sample_id, {})
        matrix[sample_id][drug_id] = value
    summary = {
        "n_input_records": len(records),
        "n_sample_drug_pairs": len(by_pair),
        "n_duplicate_sample_drug_pairs": n_duplicates,
        "n_samples": len(matrix),
        "n_drugs": len({d for row in matrix.values() for d in row}),
    }
    return matrix, summary


def orient_matrix(
    matrix: dict[str, dict[str, float]],
    *,
    response_direction: str,
) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    sign = -1.0 if str(response_direction) == "lower_is_more_sensitive" else 1.0
    for sample_id, row in matrix.items():
        out[sample_id] = {drug_id: sign * float(value) for drug_id, value in row.items()}
    return out


def _rank_normal(values: list[tuple[str, float]]) -> dict[str, float]:
    if not values:
        return {}
    ranked = sorted(values, key=lambda kv: (float(kv[1]), str(kv[0])))
    n = len(ranked)
    normal = NormalDist()
    out: dict[str, float] = {}
    for idx, (sample_id, _value) in enumerate(ranked, start=1):
        p = (float(idx) - 0.5) / float(n)
        out[sample_id] = float(normal.inv_cdf(p))
    return out


def transform_matrix_by_drug(
    matrix: dict[str, dict[str, float]],
    *,
    response_transform: str,
    eps: float = 1e-6,
) -> tuple[dict[str, dict[str, float]], dict[str, object]]:
    samples = sorted(matrix)
    drugs = sorted({drug for row in matrix.values() for drug in row})
    out: dict[str, dict[str, float]] = {sample: {} for sample in samples}
    drug_stats: dict[str, dict[str, float]] = {}

    for drug_id in drugs:
        values = [(sample_id, float(matrix[sample_id][drug_id])) for sample_id in samples if drug_id in matrix[sample_id]]
        if not values:
            continue
        raw_vals = [v for _s, v in values]
        if response_transform == "none":
            for sample_id, value in values:
                out[sample_id][drug_id] = float(value)
            drug_stats[drug_id] = {
                "n_samples": len(values),
                "median": float(statistics.median(raw_vals)),
                "mad": float(statistics.median([abs(v - statistics.median(raw_vals)) for v in raw_vals])),
            }
            continue
        if response_transform == "rank_normal_score":
            transformed = _rank_normal(values)
            for sample_id, z in transformed.items():
                out[sample_id][drug_id] = float(z)
            drug_stats[drug_id] = {
                "n_samples": len(values),
                "median": float(statistics.median(raw_vals)),
                "mad": float(statistics.median([abs(v - statistics.median(raw_vals)) for v in raw_vals])),
            }
            continue
        if response_transform != "robust_z_mad":
            raise ValueError(f"Unsupported response_transform: {response_transform}")
        med = float(statistics.median(raw_vals))
        mad = float(statistics.median([abs(v - med) for v in raw_vals]))
        denom = mad + float(eps)
        for sample_id, value in values:
            out[sample_id][drug_id] = float((float(value) - med) / denom)
        drug_stats[drug_id] = {"n_samples": len(values), "median": med, "mad": mad}

    n_total = len(samples) * len(drugs) if samples and drugs else 0
    n_present = sum(len(row) for row in out.values())
    summary = {
        "response_transform": response_transform,
        "n_samples": len(samples),
        "n_drugs": len(drugs),
        "n_values_present": n_present,
        "n_values_missing": (n_total - n_present if n_total else 0),
        "fraction_values_present": (float(n_present) / float(n_total) if n_total else 0.0),
        "drug_stats": drug_stats,
    }
    return out, summary
