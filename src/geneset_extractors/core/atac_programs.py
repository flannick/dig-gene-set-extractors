from __future__ import annotations

import math

from geneset_extractors.core.reference_calibration import apply_peak_idf


PROGRAM_LINKED_ACTIVITY = "linked_activity"
PROGRAM_PROMOTER_ACTIVITY = "promoter_activity"
PROGRAM_DISTAL_ACTIVITY = "distal_activity"
PROGRAM_ENHANCER_BIAS = "enhancer_bias"
PROGRAM_TFIDF_DISTAL = "tfidf_distal"
PROGRAM_REF_UBIQUITY_PENALTY = "ref_ubiquity_penalty"
PROGRAM_ATLAS_RESIDUAL = "atlas_residual"
PROGRAM_PRESET_NONE = "none"
PROGRAM_PRESET_DEFAULT = "default"
PROGRAM_PRESET_CONNECTABLE = "connectable"
PROGRAM_PRESET_QC = "qc"
PROGRAM_PRESET_EXPERIMENTAL = "experimental"
PROGRAM_PRESET_ALL = "all"
ATLAS_SCORE_DEFINITION_DEFAULT = "__default__"
REFERENCE_PROGRAM_METHODS = (
    PROGRAM_REF_UBIQUITY_PENALTY,
    PROGRAM_ATLAS_RESIDUAL,
)
CONTRAST_METHOD_NONE = "none"
CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY = "auto_prefer_ref_ubiquity_else_none"
ALLOWED_CONTRAST_METHODS = (
    CONTRAST_METHOD_NONE,
    PROGRAM_REF_UBIQUITY_PENALTY,
    PROGRAM_ATLAS_RESIDUAL,
    CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY,
)
CALIBRATION_METHOD_NONE = CONTRAST_METHOD_NONE
CALIBRATION_METHOD_AUTO_PREFER_REF_UBIQUITY = CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY
ALLOWED_CALIBRATION_METHODS = ALLOWED_CONTRAST_METHODS


_BULK_ALLOWED = (
    PROGRAM_LINKED_ACTIVITY,
    PROGRAM_PROMOTER_ACTIVITY,
    PROGRAM_DISTAL_ACTIVITY,
    PROGRAM_ENHANCER_BIAS,
    PROGRAM_REF_UBIQUITY_PENALTY,
    PROGRAM_ATLAS_RESIDUAL,
)
_SC_ALLOWED = (
    PROGRAM_LINKED_ACTIVITY,
    PROGRAM_PROMOTER_ACTIVITY,
    PROGRAM_DISTAL_ACTIVITY,
    PROGRAM_ENHANCER_BIAS,
    PROGRAM_TFIDF_DISTAL,
    PROGRAM_REF_UBIQUITY_PENALTY,
    PROGRAM_ATLAS_RESIDUAL,
)

_BULK_PRESETS: dict[str, tuple[str, ...]] = {
    PROGRAM_PRESET_NONE: (PROGRAM_LINKED_ACTIVITY,),
    PROGRAM_PRESET_DEFAULT: (
        PROGRAM_LINKED_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
    ),
    PROGRAM_PRESET_CONNECTABLE: (
        PROGRAM_LINKED_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
    ),
    PROGRAM_PRESET_QC: (
        PROGRAM_PROMOTER_ACTIVITY,
    ),
    PROGRAM_PRESET_EXPERIMENTAL: (
        PROGRAM_LINKED_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
        PROGRAM_ENHANCER_BIAS,
    ),
    PROGRAM_PRESET_ALL: (
        PROGRAM_LINKED_ACTIVITY,
        PROGRAM_PROMOTER_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
        PROGRAM_ENHANCER_BIAS,
        PROGRAM_REF_UBIQUITY_PENALTY,
        PROGRAM_ATLAS_RESIDUAL,
    ),
}
_SC_PRESETS: dict[str, tuple[str, ...]] = {
    PROGRAM_PRESET_NONE: (PROGRAM_LINKED_ACTIVITY,),
    PROGRAM_PRESET_DEFAULT: (
        PROGRAM_LINKED_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
    ),
    PROGRAM_PRESET_CONNECTABLE: (
        PROGRAM_LINKED_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
    ),
    PROGRAM_PRESET_QC: (
        PROGRAM_PROMOTER_ACTIVITY,
    ),
    PROGRAM_PRESET_EXPERIMENTAL: (
        PROGRAM_LINKED_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
        PROGRAM_ENHANCER_BIAS,
        PROGRAM_TFIDF_DISTAL,
    ),
    PROGRAM_PRESET_ALL: (
        PROGRAM_LINKED_ACTIVITY,
        PROGRAM_PROMOTER_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
        PROGRAM_ENHANCER_BIAS,
        PROGRAM_TFIDF_DISTAL,
        PROGRAM_REF_UBIQUITY_PENALTY,
        PROGRAM_ATLAS_RESIDUAL,
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
    return out


def ensure_reference_program_methods(program_methods: list[str]) -> list[str]:
    out = list(program_methods)
    seen = set(out)
    for method in REFERENCE_PROGRAM_METHODS:
        if method not in seen:
            out.append(method)
            seen.add(method)
    return out


def remove_reference_program_methods(program_methods: list[str]) -> list[str]:
    blocked = set(REFERENCE_PROGRAM_METHODS)
    return [method for method in program_methods if method not in blocked]


def resolve_calibration_methods(
    calibration_methods_csv: str | None,
    legacy_program_methods_csv: str | None,
    use_reference_bundle: bool,
) -> list[str]:
    requested = parse_csv_tokens(calibration_methods_csv)
    legacy = parse_csv_tokens(legacy_program_methods_csv)
    legacy_ref = [tok for tok in legacy if tok in REFERENCE_PROGRAM_METHODS]

    expanded: list[str] = []
    if requested:
        for tok in requested:
            if tok == "all":
                expanded.extend((CONTRAST_METHOD_NONE, PROGRAM_REF_UBIQUITY_PENALTY, PROGRAM_ATLAS_RESIDUAL))
            elif tok in {"auto", CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY}:
                expanded.append(CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY)
            else:
                expanded.append(tok)
    elif legacy_ref:
        expanded.append(CONTRAST_METHOD_NONE)
        expanded.extend(legacy_ref)
    else:
        expanded.append(CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY)

    out: list[str] = []
    seen: set[str] = set()
    for method in expanded:
        if method not in ALLOWED_CALIBRATION_METHODS:
            raise ValueError(f"Unsupported calibration method: {method}")
        if method in seen:
            continue
        out.append(method)
        seen.add(method)

    if not use_reference_bundle:
        out2: list[str] = []
        for method in out:
            if method in {PROGRAM_REF_UBIQUITY_PENALTY, PROGRAM_ATLAS_RESIDUAL}:
                continue
            if method == CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY:
                if CONTRAST_METHOD_NONE not in out2:
                    out2.append(CONTRAST_METHOD_NONE)
                continue
            if method not in out2:
                out2.append(method)
        out = out2

    if not out:
        return [CONTRAST_METHOD_NONE]
    return out


def resolve_contrast_methods(
    contrast_methods_csv: str | None,
    legacy_program_methods_csv: str | None,
    use_reference_bundle: bool,
) -> list[str]:
    return resolve_calibration_methods(
        contrast_methods_csv,
        legacy_program_methods_csv,
        use_reference_bundle=use_reference_bundle,
    )


def calibration_method_enablement_hint(method: str, genome_build: str) -> str | None:
    if method not in REFERENCE_PROGRAM_METHODS and method != CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY:
        return None
    gb = genome_build.strip().lower()
    if gb in {"hg19", "grch37"}:
        preset = "atac_default_optional_hg19"
    elif gb in {"hg38", "grch38"}:
        preset = "atac_default_optional_hg38"
    else:
        preset = ""
    if preset:
        return (
            "Enable by providing a reference bundle via --resources_dir <bundle_root> "
            "(see docs/atac_reference_bundle.md), or fetch resources with "
            f"`geneset-extractors resources fetch --preset {preset}`."
        )
    return (
        "Enable by providing a compatible reference bundle via --resources_dir <bundle_root> "
        "(see docs/atac_reference_bundle.md)."
    )


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


def atlas_residual_scores(
    scores: dict[str, float],
    gene_stats: dict[str, tuple[float, float]],
    metric: str,
    eps: float,
    min_raw_quantile: float = 0.0,
    use_log1p: bool = True,
) -> dict[str, float]:
    out: dict[str, float] = {}
    tiny = float(eps)
    q = float(min_raw_quantile)
    if q < 0.0 or q > 1.0:
        raise ValueError(f"atlas min_raw_quantile must be in [0,1], got {q}")
    floor = None
    if q > 0.0 and scores:
        raw_vals = sorted(float(v) for v in scores.values())
        if raw_vals:
            idx = max(0, min(len(raw_vals) - 1, int(math.floor(q * (len(raw_vals) - 1)))))
            floor = raw_vals[idx]

    def _signed_log1p(x: float) -> float:
        if x == 0.0:
            return 0.0
        return math.copysign(math.log1p(abs(x)), x)

    for gene_id, score in scores.items():
        if gene_id not in gene_stats:
            continue
        if floor is not None and float(score) < float(floor):
            continue
        median, mad = gene_stats[gene_id]
        s = float(score)
        if metric == "logratio":
            numerator = (math.log1p(max(s, 0.0)) if use_log1p else s) + tiny
            denominator = (math.log1p(max(float(median), 0.0)) if use_log1p else float(median)) + tiny
            if numerator <= 0.0 or denominator <= 0.0:
                continue
            out[gene_id] = math.log(numerator / denominator)
            continue
        if metric == "zscore":
            if use_log1p:
                s_t = _signed_log1p(s)
                m_t = _signed_log1p(float(median))
                mad_t = abs(_signed_log1p(float(median) + float(mad)) - m_t)
                out[gene_id] = (s_t - m_t) / (mad_t + tiny)
            else:
                out[gene_id] = (s - float(median)) / (float(mad) + tiny)
            continue
        raise ValueError(f"Unsupported atlas metric: {metric}")
    return out


def normalize_program_preset(preset: str | None) -> str:
    name = "default" if preset is None else str(preset).strip()
    if not name:
        return PROGRAM_PRESET_DEFAULT
    if name == PROGRAM_PRESET_DEFAULT:
        return PROGRAM_PRESET_CONNECTABLE
    return name


def default_link_methods_for_preset(preset: str | None) -> tuple[str, ...]:
    normalized = normalize_program_preset(preset)
    if normalized == PROGRAM_PRESET_CONNECTABLE:
        return ("nearest_tss", "distance_decay")
    if normalized == PROGRAM_PRESET_QC:
        return ("promoter_overlap",)
    if normalized == PROGRAM_PRESET_EXPERIMENTAL:
        return ("nearest_tss", "distance_decay")
    return ("promoter_overlap", "nearest_tss", "distance_decay")


def preferred_link_method_for_program(preset: str | None, program_method: str) -> str | None:
    normalized = normalize_program_preset(preset)
    if normalized == PROGRAM_PRESET_CONNECTABLE:
        if program_method == PROGRAM_LINKED_ACTIVITY:
            return "nearest_tss"
        if program_method == PROGRAM_DISTAL_ACTIVITY:
            return "distance_decay"
        return None
    if normalized == PROGRAM_PRESET_QC:
        if program_method in {PROGRAM_PROMOTER_ACTIVITY, PROGRAM_LINKED_ACTIVITY}:
            return "promoter_overlap"
        return None
    return None


def link_methods_for_program(
    preset: str | None,
    program_method: str,
    link_methods: list[str],
) -> list[str]:
    preferred = preferred_link_method_for_program(preset, program_method)
    if preferred is None:
        return list(link_methods)
    return [method for method in link_methods if method == preferred]


def contrast_method_enablement_hint(method: str, genome_build: str) -> str | None:
    return calibration_method_enablement_hint(method, genome_build)


def resolve_auto_calibration_methods(
    requested: list[str],
    ref_ubiquity_ready: bool,
) -> tuple[list[str], str | None]:
    out: list[str] = []
    fallback_reason: str | None = None
    for method in requested:
        if method == CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY:
            chosen = PROGRAM_REF_UBIQUITY_PENALTY if ref_ubiquity_ready else CONTRAST_METHOD_NONE
            if not ref_ubiquity_ready:
                fallback_reason = (
                    "calibration_method=auto_prefer_ref_ubiquity_else_none could not load ref_ubiquity resource; "
                    "falling back to calibration_method=none"
                )
        else:
            chosen = method
        if chosen not in out:
            out.append(chosen)
    if not out:
        out.append(CONTRAST_METHOD_NONE)
    return out, fallback_reason


def resolve_auto_contrast_methods(
    requested: list[str],
    ref_ubiquity_ready: bool,
) -> tuple[list[str], str | None]:
    return resolve_auto_calibration_methods(requested, ref_ubiquity_ready=ref_ubiquity_ready)


def peak_values_for_calibration(
    peak_values: list[float],
    calibration_method: str,
    ref_peak_idf: list[float] | None,
) -> list[float]:
    if calibration_method == PROGRAM_REF_UBIQUITY_PENALTY:
        if ref_peak_idf is None:
            raise ValueError("ref_ubiquity contrast requested but ref_peak_idf is unavailable")
        return apply_peak_idf(peak_values, ref_peak_idf)
    return [float(v) for v in peak_values]


def peak_values_for_contrast(
    peak_values: list[float],
    contrast_method: str,
    ref_peak_idf: list[float] | None,
) -> list[float]:
    return peak_values_for_calibration(
        peak_values,
        calibration_method=contrast_method,
        ref_peak_idf=ref_peak_idf,
    )


def score_definition_key(program_method: str, link_method: str) -> str:
    return f"program={program_method}__link_method={link_method}"


def atlas_stats_for_score_definition(
    atlas_stats_by_definition: dict[str, dict[str, tuple[float, float]]],
    score_definition: str,
) -> dict[str, tuple[float, float]] | None:
    if score_definition in atlas_stats_by_definition:
        return atlas_stats_by_definition[score_definition]
    if ATLAS_SCORE_DEFINITION_DEFAULT in atlas_stats_by_definition:
        return atlas_stats_by_definition[ATLAS_SCORE_DEFINITION_DEFAULT]
    return None
