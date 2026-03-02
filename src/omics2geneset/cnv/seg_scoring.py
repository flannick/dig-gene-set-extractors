from __future__ import annotations

import math


def focal_penalty(
    segment_length_bp: int,
    focal_length_scale_bp: float,
    focal_length_alpha: float,
) -> float:
    length = max(1.0, float(segment_length_bp))
    scale = max(1.0, float(focal_length_scale_bp))
    alpha = max(1e-9, float(focal_length_alpha))
    return 1.0 / (1.0 + (length / scale) ** alpha)


def gene_count_penalty(n_genes_in_segment: int, mode: str) -> float:
    n = max(1, int(n_genes_in_segment))
    if mode == "none":
        return 1.0
    if mode == "inv_sqrt":
        return 1.0 / math.sqrt(float(n))
    if mode == "inv":
        return 1.0 / float(n)
    raise ValueError(f"Unsupported gene_count_penalty mode: {mode}")


def purity_corrected_amplitude(
    amplitude: float,
    purity: float | None,
    *,
    use_purity_correction: bool,
    purity_floor: float,
    max_abs_amplitude: float,
) -> float:
    value = float(amplitude)
    if use_purity_correction and purity is not None:
        value = value / max(float(purity), float(purity_floor))
    cap = abs(float(max_abs_amplitude))
    if cap > 0:
        value = max(-cap, min(cap, value))
    return value

