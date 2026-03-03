from __future__ import annotations


def normalize(gene_weights: dict[str, float], method: str) -> dict[str, float]:
    if method == "none":
        return dict(gene_weights)
    if method == "l1" or method == "within_set_l1":
        total = sum(gene_weights.values())
        if total <= 0:
            return dict(gene_weights)
        return {k: v / total for k, v in gene_weights.items()}
    raise ValueError(f"Unsupported normalization method: {method}")
