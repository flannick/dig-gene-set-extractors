from __future__ import annotations

import math


PROGRAM_LINKED_ACTIVITY = "linked_activity"
PROGRAM_PROMOTER_ACTIVITY = "promoter_activity"
PROGRAM_DISTAL_ACTIVITY = "distal_activity"
PROGRAM_ENHANCER_BIAS = "enhancer_bias"
PROGRAM_TFIDF_DISTAL = "tfidf_distal"


_BULK_ALLOWED = (
    PROGRAM_LINKED_ACTIVITY,
    PROGRAM_PROMOTER_ACTIVITY,
    PROGRAM_DISTAL_ACTIVITY,
    PROGRAM_ENHANCER_BIAS,
)
_SC_ALLOWED = (
    PROGRAM_LINKED_ACTIVITY,
    PROGRAM_PROMOTER_ACTIVITY,
    PROGRAM_DISTAL_ACTIVITY,
    PROGRAM_ENHANCER_BIAS,
    PROGRAM_TFIDF_DISTAL,
)

_BULK_PRESETS: dict[str, tuple[str, ...]] = {
    "none": (PROGRAM_LINKED_ACTIVITY,),
    "default": (
        PROGRAM_LINKED_ACTIVITY,
        PROGRAM_PROMOTER_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
        PROGRAM_ENHANCER_BIAS,
    ),
    "all": (
        PROGRAM_LINKED_ACTIVITY,
        PROGRAM_PROMOTER_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
        PROGRAM_ENHANCER_BIAS,
    ),
}
_SC_PRESETS: dict[str, tuple[str, ...]] = {
    "none": (PROGRAM_LINKED_ACTIVITY,),
    "default": (
        PROGRAM_LINKED_ACTIVITY,
        PROGRAM_PROMOTER_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
        PROGRAM_ENHANCER_BIAS,
        PROGRAM_TFIDF_DISTAL,
    ),
    "all": (
        PROGRAM_LINKED_ACTIVITY,
        PROGRAM_PROMOTER_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
        PROGRAM_ENHANCER_BIAS,
        PROGRAM_TFIDF_DISTAL,
    ),
}


def parse_csv_tokens(value: str | None) -> list[str]:
    text = "" if value is None else str(value).strip()
    if not text:
        return []
    return [tok.strip() for tok in text.split(",") if tok.strip()]


def resolve_program_methods(
    assay: str,
    program_preset: str | None,
    program_methods_csv: str | None,
) -> list[str]:
    if assay == "bulk":
        allowed = _BULK_ALLOWED
        presets = _BULK_PRESETS
    elif assay == "single_cell":
        allowed = _SC_ALLOWED
        presets = _SC_PRESETS
    else:
        raise ValueError(f"Unsupported assay for ATAC program methods: {assay}")

    preset_name = "default" if program_preset is None else str(program_preset).strip()
    if not preset_name:
        preset_name = "default"
    if preset_name not in presets:
        raise ValueError(f"Unsupported program_preset: {preset_name}")

    requested = parse_csv_tokens(program_methods_csv)
    expanded: list[str] = []
    if requested:
        for tok in requested:
            if tok == "all":
                expanded.extend(allowed)
            else:
                expanded.append(tok)
    else:
        expanded.extend(presets[preset_name])

    out: list[str] = []
    seen: set[str] = set()
    for method in expanded:
        if method not in allowed:
            raise ValueError(f"Unsupported program method for assay={assay}: {method}")
        if method in seen:
            continue
        out.append(method)
        seen.add(method)

    if PROGRAM_LINKED_ACTIVITY not in seen:
        out.insert(0, PROGRAM_LINKED_ACTIVITY)
    return out


def promoter_peak_indices(promoter_links: list[dict[str, object]]) -> set[int]:
    return {int(link["peak_index"]) for link in promoter_links}


def mask_peak_weights(
    peak_weights: list[float],
    include_indices: set[int] | None = None,
    exclude_indices: set[int] | None = None,
) -> list[float]:
    out: list[float] = []
    has_include = include_indices is not None
    has_exclude = exclude_indices is not None
    for i, w in enumerate(peak_weights):
        keep = True
        if has_include and i not in include_indices:
            keep = False
        if has_exclude and i in exclude_indices:
            keep = False
        out.append(float(w) if keep else 0.0)
    return out


def enhancer_bias_scores(
    promoter_scores: dict[str, float],
    distal_scores: dict[str, float],
    pseudocount: float = 1e-6,
) -> dict[str, float]:
    genes = set(promoter_scores) | set(distal_scores)
    eps = float(pseudocount)
    out: dict[str, float] = {}
    for gene_id in genes:
        prom = float(promoter_scores.get(gene_id, 0.0))
        dist = float(distal_scores.get(gene_id, 0.0))
        out[gene_id] = math.log2((dist + eps) / (prom + eps))
    return out


def compute_peak_idf(
    matrix_entries,
    n_peaks: int,
    n_cells: int,
) -> list[float]:
    df = [0] * n_peaks
    for peak_i, _cell_i, val in matrix_entries:
        if float(val) != 0.0:
            df[int(peak_i)] += 1
    denom = float(n_cells) + 1.0
    return [math.log(denom / (float(d) + 1.0)) + 1.0 for d in df]
