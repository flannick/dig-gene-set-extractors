from __future__ import annotations


def build_reference_gene_maps(
    *,
    reference_metadata: dict[str, dict[str, str]],
    compound_targets: dict[str, dict[str, float]],
) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]], dict[str, object]]:
    compound_map: dict[str, dict[str, float]] = {}
    genetic_map: dict[str, dict[str, float]] = {}
    n_compound = 0
    n_genetic = 0
    n_missing = 0
    for perturbation_id, row in reference_metadata.items():
        perturbation_type = str(row.get("perturbation_type", "")).strip().lower()
        if perturbation_type == "compound":
            compound_id = str(row.get("compound_id", "") or row.get("perturbation_id", perturbation_id)).strip()
            mapping = compound_targets.get(compound_id) or compound_targets.get(perturbation_id)
            if mapping:
                compound_map[perturbation_id] = dict(mapping)
                n_compound += 1
            else:
                n_missing += 1
        elif perturbation_type in {"orf", "crispr"}:
            gene_symbol = str(row.get("gene_symbol", "")).strip().upper()
            if gene_symbol:
                genetic_map[perturbation_id] = {gene_symbol: 1.0}
                n_genetic += 1
            else:
                n_missing += 1
        else:
            n_missing += 1
    return compound_map, genetic_map, {
        "n_compound_refs_with_targets": n_compound,
        "n_genetic_refs": n_genetic,
        "n_refs_missing_mapping": n_missing,
    }


def accumulate_gene_scores(
    *,
    evidence_by_ref: dict[str, float],
    ref_to_gene: dict[str, dict[str, float]],
) -> dict[str, float]:
    scores: dict[str, float] = {}
    for ref_id, evidence in evidence_by_ref.items():
        if evidence <= 0.0:
            continue
        mapping = ref_to_gene.get(ref_id)
        if not mapping:
            continue
        for gene_symbol, weight in mapping.items():
            scores[gene_symbol] = float(scores.get(gene_symbol, 0.0)) + float(evidence) * float(weight)
    return scores


def l1_normalize_scores(scores: dict[str, float]) -> dict[str, float]:
    total = sum(float(v) for v in scores.values() if float(v) > 0.0)
    if total <= 0.0:
        return {}
    return {gene: float(value) / total for gene, value in scores.items() if float(value) > 0.0}
