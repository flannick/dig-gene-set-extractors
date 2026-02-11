from __future__ import annotations


def transform_peak_weight(weight: float, mode: str) -> float:
    if mode == "signed":
        return weight
    if mode == "abs":
        return abs(weight)
    if mode == "positive":
        return max(weight, 0.0)
    if mode == "negative":
        return max(-weight, 0.0)
    raise ValueError(f"Unsupported peak weight transform: {mode}")


def score_genes(
    peak_weights: list[float],
    links: list[dict[str, object]],
    peak_weight_transform: str,
) -> dict[str, float]:
    gene_weights: dict[str, float] = {}
    transformed = [transform_peak_weight(float(w), peak_weight_transform) for w in peak_weights]
    for link in links:
        p = int(link["peak_index"])
        g = str(link["gene_id"])
        lw = float(link["link_weight"])
        gene_weights[g] = gene_weights.get(g, 0.0) + transformed[p] * lw
    return gene_weights
