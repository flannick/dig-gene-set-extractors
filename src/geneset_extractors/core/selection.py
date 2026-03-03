from __future__ import annotations

import math


def _sorted_items(scores: dict[str, float]) -> list[tuple[str, float]]:
    return sorted(scores.items(), key=lambda kv: (-float(kv[1]), str(kv[0])))


def select_top_k(scores: dict[str, float], k: int) -> list[str]:
    if k <= 0:
        return []
    return [gene_id for gene_id, _ in _sorted_items(scores)[:k]]


def select_quantile(scores: dict[str, float], q: float) -> list[str]:
    if not (0.0 < q < 1.0):
        raise ValueError("q must be in (0,1)")
    if not scores:
        return []
    k = max(1, int(math.ceil(len(scores) * q)))
    return select_top_k(scores, k)


def select_threshold(scores: dict[str, float], min_score: float) -> list[str]:
    return [gene_id for gene_id, score in _sorted_items(scores) if float(score) >= float(min_score)]


def softmax_temperature(scores: dict[str, float], temperature: float) -> dict[str, float]:
    if temperature <= 0:
        raise ValueError("temperature must be > 0")
    if not scores:
        return {}
    max_score = max(float(v) for v in scores.values())
    exps = {
        g: math.exp((float(s) - max_score) / temperature)
        for g, s in scores.items()
    }
    z = sum(exps.values())
    if z <= 0:
        n = len(scores)
        return {g: 1.0 / n for g in scores}
    return {g: v / z for g, v in exps.items()}


def _effective_n(weights: dict[str, float]) -> float:
    s2 = sum(float(v) * float(v) for v in weights.values())
    if s2 <= 0:
        return 0.0
    return 1.0 / s2


def find_temperature_for_neff(scores: dict[str, float], target_neff: float, tol: float = 1e-3) -> float:
    if not scores:
        raise ValueError("scores is empty")
    if target_neff <= 0:
        raise ValueError("target_neff must be > 0")

    lo = 1e-6
    hi = 1e6
    for _ in range(100):
        mid = (lo + hi) / 2.0
        w = softmax_temperature(scores, mid)
        neff = _effective_n(w)
        if abs(neff - target_neff) <= tol:
            return mid
        if neff < target_neff:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def ranked_gene_ids(scores: dict[str, float]) -> list[str]:
    return [gene_id for gene_id, _ in _sorted_items(scores)]


def within_set_l1_weights(scores: dict[str, float], selected_gene_ids: list[str]) -> dict[str, float]:
    if not selected_gene_ids:
        return {}
    total = sum(float(scores.get(g, 0.0)) for g in selected_gene_ids)
    if total > 0:
        return {g: float(scores.get(g, 0.0)) / total for g in selected_gene_ids}
    uniform = 1.0 / len(selected_gene_ids)
    return {g: uniform for g in selected_gene_ids}


def global_l1_weights(scores: dict[str, float]) -> dict[str, float]:
    total = sum(float(v) for v in scores.values())
    if total <= 0:
        return dict(scores)
    return {g: float(v) / total for g, v in scores.items()}
