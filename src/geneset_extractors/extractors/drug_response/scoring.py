from __future__ import annotations

import math

import numpy as np


def compute_ubiquity_weights(
    *,
    z_matrix: dict[str, dict[str, float]],
    drugs: list[str],
    method: str,
    tau: float,
    epsilon: float,
) -> tuple[dict[str, float], dict[str, object]]:
    if method == "none":
        return {drug: 1.0 for drug in drugs}, {"method": "none", "tau": tau, "epsilon": epsilon}
    if method != "fraction_active":
        raise ValueError(f"Unsupported ubiquity_penalty method: {method}")

    samples = sorted(z_matrix)
    weights: dict[str, float] = {}
    activity_fraction: dict[str, float] = {}
    for drug in drugs:
        if not samples:
            frac = 0.0
        else:
            active = 0
            for sample_id in samples:
                value = float(z_matrix.get(sample_id, {}).get(drug, 0.0))
                if value > float(tau):
                    active += 1
            frac = float(active) / float(len(samples))
        activity_fraction[drug] = frac
        weights[drug] = 1.0 / (float(epsilon) + frac)
    return weights, {
        "method": method,
        "tau": tau,
        "epsilon": epsilon,
        "activity_fraction": activity_fraction,
    }


def compute_polypharm_weights(
    *,
    drug_targets: dict[str, dict[str, float]],
    enabled: bool,
    t0: int,
) -> tuple[dict[str, float], dict[str, object]]:
    out: dict[str, float] = {}
    target_counts: dict[str, int] = {}
    if not enabled:
        for drug in drug_targets:
            out[drug] = 1.0
            target_counts[drug] = len(drug_targets[drug])
        return out, {"enabled": False, "t0": t0, "target_counts": target_counts}

    threshold = max(1, int(t0))
    for drug, genes in drug_targets.items():
        n_targets = len(genes)
        target_counts[drug] = n_targets
        out[drug] = 1.0 / (1.0 + math.log(1.0 + (float(n_targets) / float(threshold))))
    return out, {"enabled": True, "t0": threshold, "target_counts": target_counts}


def compute_target_ubiquity_weights(
    *,
    drug_targets: dict[str, dict[str, float]],
    method: str,
) -> tuple[dict[str, float], dict[str, object]]:
    mode = str(method).strip().lower()
    gene_to_drug_count: dict[str, int] = {}
    n_drugs = len(drug_targets)
    for targets in drug_targets.values():
        for gene_id in targets:
            gene_to_drug_count[gene_id] = int(gene_to_drug_count.get(gene_id, 0)) + 1
    if mode == "none":
        return (
            {gene_id: 1.0 for gene_id in gene_to_drug_count},
            {
                "method": "none",
                "n_drugs": n_drugs,
                "n_genes": len(gene_to_drug_count),
                "n_genes_affected": 0,
                "top_ubiquitous": [],
            },
        )
    if mode != "idf":
        raise ValueError(f"Unsupported target_ubiquity_penalty method: {method}")

    gene_weights: dict[str, float] = {}
    for gene_id, f_g in gene_to_drug_count.items():
        gene_weights[gene_id] = float(math.log((float(n_drugs) + 1.0) / (float(f_g) + 1.0)))
    top_ubiquitous = sorted(gene_to_drug_count.items(), key=lambda kv: (-int(kv[1]), str(kv[0])))[:20]
    return (
        gene_weights,
        {
            "method": "idf",
            "n_drugs": n_drugs,
            "n_genes": len(gene_to_drug_count),
            "n_genes_affected": sum(1 for f_g in gene_to_drug_count.values() if int(f_g) > 1),
            "top_ubiquitous": [{"gene_id": gene_id, "f": int(f_g)} for gene_id, f_g in top_ubiquitous],
        },
    )


def score_genes_target_weighted(
    *,
    program_drug_scores: dict[str, float],
    drug_targets: dict[str, dict[str, float]],
    drug_weights: dict[str, float],
    target_gene_weights: dict[str, float] | None,
    signed: bool,
) -> tuple[dict[str, float], dict[str, object]]:
    gene_scores: dict[str, float] = {}
    n_drugs_with_mapping = 0
    n_drugs_used = 0
    for drug_id, score in program_drug_scores.items():
        targets = drug_targets.get(drug_id)
        if not targets:
            continue
        n_drugs_with_mapping += 1
        contribution = float(score)
        if not signed:
            contribution = max(contribution, 0.0)
        if contribution == 0.0:
            continue
        n_drugs_used += 1
        drug_weight = float(drug_weights.get(drug_id, 1.0))
        for gene_id, target_weight in targets.items():
            gene_weight = 1.0 if target_gene_weights is None else float(target_gene_weights.get(gene_id, 1.0))
            delta = contribution * drug_weight * gene_weight * float(target_weight)
            gene_scores[gene_id] = float(gene_scores.get(gene_id, 0.0)) + delta
    return gene_scores, {
        "n_drugs_with_mapping": n_drugs_with_mapping,
        "n_drugs_used_nonzero": n_drugs_used,
        "signed": bool(signed),
    }


def _score_genes_sparse_deconvolution(
    *,
    program_drug_scores: dict[str, float],
    drug_targets: dict[str, dict[str, float]],
    drug_weights: dict[str, float],
    target_gene_weights: dict[str, float] | None,
    signed: bool,
    sparse_alpha: float,
) -> tuple[dict[str, float], str | None]:
    try:
        from sklearn.linear_model import Lasso  # type: ignore
    except ModuleNotFoundError:
        return {}, (
            "scoring_model=sparse_deconvolution requested but scikit-learn is not installed; "
            "falling back to target_weighted_sum"
        )

    genes = sorted({gene for targets in drug_targets.values() for gene in targets})
    drugs = sorted(program_drug_scores)
    if not genes or not drugs:
        return {}, None

    gene_index = {gene: idx for idx, gene in enumerate(genes)}
    X = np.zeros((len(drugs), len(genes)), dtype=float)
    y = np.zeros((len(drugs),), dtype=float)

    for row_idx, drug_id in enumerate(drugs):
        y[row_idx] = float(program_drug_scores.get(drug_id, 0.0))
        weighted = float(drug_weights.get(drug_id, 1.0))
        for gene_id, target_weight in drug_targets.get(drug_id, {}).items():
            col_idx = gene_index[gene_id]
            gene_weight = 1.0 if target_gene_weights is None else float(target_gene_weights.get(gene_id, 1.0))
            X[row_idx, col_idx] = weighted * gene_weight * float(target_weight)

    if signed:
        y_pos = np.maximum(y, 0.0)
        y_neg = np.maximum(-y, 0.0)
        model_pos = Lasso(alpha=float(sparse_alpha), positive=True, fit_intercept=False, max_iter=10000)
        model_neg = Lasso(alpha=float(sparse_alpha), positive=True, fit_intercept=False, max_iter=10000)
        model_pos.fit(X, y_pos)
        model_neg.fit(X, y_neg)
        coefs = model_pos.coef_ - model_neg.coef_
    else:
        y = np.maximum(y, 0.0)
        model = Lasso(alpha=float(sparse_alpha), positive=True, fit_intercept=False, max_iter=10000)
        model.fit(X, y)
        coefs = model.coef_
    scores = {gene: float(coef) for gene, coef in zip(genes, coefs) if float(coef) != 0.0}
    return scores, None


def score_genes(
    *,
    scoring_model: str,
    program_drug_scores: dict[str, float],
    drug_targets: dict[str, dict[str, float]],
    drug_weights: dict[str, float],
    target_gene_weights: dict[str, float] | None,
    signed: bool,
    sparse_alpha: float = 0.01,
) -> tuple[dict[str, float], dict[str, object]]:
    if scoring_model == "target_weighted_sum":
        scores, info = score_genes_target_weighted(
            program_drug_scores=program_drug_scores,
            drug_targets=drug_targets,
            drug_weights=drug_weights,
            target_gene_weights=target_gene_weights,
            signed=signed,
        )
        return scores, {"model_used": "target_weighted_sum", **info}

    if scoring_model == "sparse_deconvolution":
        sparse_scores, warning = _score_genes_sparse_deconvolution(
            program_drug_scores=program_drug_scores,
            drug_targets=drug_targets,
            drug_weights=drug_weights,
            target_gene_weights=target_gene_weights,
            signed=signed,
            sparse_alpha=float(sparse_alpha),
        )
        if sparse_scores:
            return sparse_scores, {"model_used": "sparse_deconvolution"}
        fallback_scores, info = score_genes_target_weighted(
            program_drug_scores=program_drug_scores,
            drug_targets=drug_targets,
            drug_weights=drug_weights,
            target_gene_weights=target_gene_weights,
            signed=signed,
        )
        payload = {"model_used": "target_weighted_sum", "fallback_from": "sparse_deconvolution", **info}
        if warning:
            payload["warning"] = warning
        return fallback_scores, payload

    raise ValueError(f"Unsupported scoring_model: {scoring_model}")
