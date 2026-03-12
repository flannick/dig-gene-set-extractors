from __future__ import annotations

from dataclasses import dataclass
from math import erf, sqrt
import sys

import numpy as np

from .features import MASS_DEPENDENT_VARS, SubjectSummary, feature_base_variable


MASS_COVARIATE_ORDER = ("lean_mass", "total_mass", "subject_mass")


@dataclass
class ContrastSignature:
    contrast_id: str
    contrast_label: str
    group_label: str
    feature_scores: dict[str, float]
    feature_models: dict[str, dict[str, object]]
    metadata: dict[str, object]


def _normal_pvalue_from_t(t_value: float) -> float:
    z = abs(float(t_value)) / sqrt(2.0)
    return float(max(0.0, min(1.0, 2.0 * (1.0 - 0.5 * (1.0 + erf(z))))))


def _mean(values: list[float]) -> float:
    return float(sum(values) / float(len(values) or 1))


def _pooled_sd(a: list[float], b: list[float]) -> float:
    if len(a) < 2 or len(b) < 2:
        return 1.0
    mean_a = _mean(a)
    mean_b = _mean(b)
    ss = sum((x - mean_a) ** 2 for x in a) + sum((x - mean_b) ** 2 for x in b)
    denom = max(1, len(a) + len(b) - 2)
    return float(max(1e-8, sqrt(ss / float(denom))))


def _welch_like(a: list[float], b: list[float]) -> float:
    if not a or not b:
        return 0.0
    pooled = _pooled_sd(a, b)
    return float((_mean(a) - _mean(b)) / pooled)


def _fit_linear(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    fitted = x @ beta
    resid = y - fitted
    dof = max(1, int(len(y) - x.shape[1]))
    sigma2 = float(np.sum(resid ** 2) / float(dof))
    xtx_inv = np.linalg.pinv(x.T @ x)
    se = np.sqrt(np.maximum(1e-12, np.diag(xtx_inv) * sigma2))
    return beta, se, sigma2


def _mass_effect(
    *,
    values: list[float],
    groups: list[int],
    covariate: list[float],
    min_group_size: int,
) -> tuple[float, dict[str, object]]:
    n_group = sum(groups)
    n_other = len(groups) - n_group
    if n_group < min_group_size or n_other < min_group_size:
        score = _welch_like([v for v, g in zip(values, groups) if g == 1], [v for v, g in zip(values, groups) if g == 0])
        return score, {
            "mass_adjusted": False,
            "interaction_tested": False,
            "interaction_retained": False,
            "model_kind": "fallback_small_design",
            "interaction_pvalue": None,
        }

    y = np.asarray(values, dtype=float)
    g = np.asarray(groups, dtype=float)
    c = np.asarray(covariate, dtype=float)
    c_centered = c - float(np.mean(c))
    x_full = np.column_stack([np.ones_like(g), g, c_centered, g * c_centered])
    beta_full, se_full, _sigma2 = _fit_linear(y, x_full)
    interaction_t = float(beta_full[3] / se_full[3]) if se_full[3] > 0 else 0.0
    interaction_p = _normal_pvalue_from_t(interaction_t)
    if interaction_p < 0.05:
        effect = float(beta_full[1])
        return effect, {
            "mass_adjusted": True,
            "interaction_tested": True,
            "interaction_retained": True,
            "interaction_pvalue": interaction_p,
            "model_kind": "glm_interaction",
            "group_effect": effect,
        }

    x_simple = np.column_stack([np.ones_like(g), g, c_centered])
    beta_simple, se_simple, _sigma2 = _fit_linear(y, x_simple)
    group_t = float(beta_simple[1] / se_simple[1]) if se_simple[1] > 0 else 0.0
    return float(group_t), {
        "mass_adjusted": True,
        "interaction_tested": True,
        "interaction_retained": False,
        "interaction_pvalue": interaction_p,
        "model_kind": "ancova",
        "group_effect": float(beta_simple[1]),
    }


def choose_mass_covariate(subjects: list[SubjectSummary], explicit_covariate: str | None) -> tuple[str | None, bool]:
    if explicit_covariate:
        key = str(explicit_covariate).strip().lower().replace(".", "_")
        if any(subject.metadata.get(key) is not None for subject in subjects):
            return key, True
    for key in MASS_COVARIATE_ORDER:
        present = [subject.metadata.get(key) for subject in subjects if subject.metadata.get(key) is not None]
        if len(present) >= max(2, len(subjects) // 2):
            return key, True
    return None, False


def detect_design_class(subjects: list[SubjectSummary]) -> str:
    if any(str(subject.metadata.get("period", "")).strip() for subject in subjects):
        return "unsupported_crossover"
    return "simple_grouped"


def build_group_contrasts(
    *,
    subjects_by_id: dict[str, SubjectSummary],
    explicit_mass_covariate: str | None,
    min_group_size: int,
) -> tuple[list[ContrastSignature], list[str]]:
    warnings: list[str] = []
    subjects = [subjects_by_id[key] for key in sorted(subjects_by_id)]
    design_class = detect_design_class(subjects)
    if design_class == "unsupported_crossover":
        raise ValueError("Detected crossover/period-style design markers; this design class is not supported in v1")
    groups = sorted({subject.group for subject in subjects if subject.group})
    if not groups:
        groups = ["all"]
        for subject in subjects:
            subject.group = "all"
    mass_covariate, mass_adjusted_possible = choose_mass_covariate(subjects, explicit_mass_covariate)
    if not mass_adjusted_possible:
        warning = "no suitable mass/body-composition covariate found; mass-dependent variables will run in exploratory fallback mode"
        warnings.append(warning)
        print(f"warning: {warning}", file=sys.stderr)

    feature_names = sorted({feature for subject in subjects for feature in subject.feature_values})
    contrasts: list[ContrastSignature] = []
    baseline_only = len(groups) == 1
    group_iter = groups if not baseline_only else [groups[0]]
    for group in group_iter:
        feature_scores: dict[str, float] = {}
        feature_models: dict[str, dict[str, object]] = {}
        included_subjects = [subject for subject in subjects if subject.group == group or not baseline_only]
        if baseline_only:
            subject_subset = subjects
        else:
            subject_subset = subjects
        for feature_name in feature_names:
            feature_values = []
            group_flags = []
            cov_values = []
            for subject in subject_subset:
                value = subject.feature_values.get(feature_name)
                if value is None:
                    continue
                feature_values.append(float(value))
                group_flags.append(1 if subject.group == group else 0)
                cov = subject.metadata.get(mass_covariate) if mass_covariate else None
                cov_values.append(float(cov) if cov is not None else float("nan"))
            if not feature_values:
                continue
            base_var = feature_base_variable(feature_name)
            if baseline_only:
                ctr = _mean(feature_values)
                scale = _pooled_sd(feature_values, feature_values)
                feature_scores[feature_name] = float(ctr / max(scale, 1e-8))
                feature_models[feature_name] = {
                    "interaction_tested": False,
                    "interaction_retained": False,
                    "mass_adjusted": False,
                    "model_kind": "baseline",
                }
                continue
            if base_var in MASS_DEPENDENT_VARS and mass_covariate is not None and all(np.isfinite(cov_values)):
                score, model_info = _mass_effect(
                    values=feature_values,
                    groups=group_flags,
                    covariate=cov_values,
                    min_group_size=min_group_size,
                )
                feature_scores[feature_name] = float(score)
                feature_models[feature_name] = model_info
            else:
                if base_var in MASS_DEPENDENT_VARS and mass_covariate is None:
                    warnings.append(f"mass-dependent feature {feature_name} used fallback contrast because no covariate was available")
                selected = [value for value, flag in zip(feature_values, group_flags) if flag == 1]
                other = [value for value, flag in zip(feature_values, group_flags) if flag == 0]
                feature_scores[feature_name] = _welch_like(selected, other)
                feature_models[feature_name] = {
                    "interaction_tested": False,
                    "interaction_retained": False,
                    "mass_adjusted": False,
                    "model_kind": "anova" if base_var not in MASS_DEPENDENT_VARS else "fallback_no_covariate",
                }
        interaction_tested = sum(1 for row in feature_models.values() if bool(row.get("interaction_tested")))
        interaction_retained = sum(1 for row in feature_models.values() if bool(row.get("interaction_retained")))
        contrast_id = group if not baseline_only else "baseline"
        contrasts.append(
            ContrastSignature(
                contrast_id=contrast_id,
                contrast_label=(f"{group}_vs_rest" if not baseline_only else "baseline"),
                group_label=group,
                feature_scores=feature_scores,
                feature_models=feature_models,
                metadata={
                    "mass_covariate": mass_covariate,
                    "mass_adjusted": bool(mass_covariate),
                    "interaction_tested": interaction_tested,
                    "interaction_retained": interaction_retained,
                    "design_class": design_class,
                    "run_ids": sorted({subject.run_id for subject in subjects if subject.run_id}),
                    "n_subjects_total": len(subjects),
                    "n_subjects_group": sum(1 for subject in subjects if subject.group == group),
                    "n_subjects_other": sum(1 for subject in subjects if subject.group != group),
                },
            )
        )
    return contrasts, warnings
