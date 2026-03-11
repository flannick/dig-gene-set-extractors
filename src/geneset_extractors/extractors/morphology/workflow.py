from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path
import re
import sys
from statistics import median

from geneset_extractors.core.gmt import build_gmt_sets_from_rows, parse_int_list_csv, parse_mass_list_csv, parse_str_list_csv, resolve_gmt_out_path, write_gmt
from geneset_extractors.core.metadata import make_metadata, write_metadata
from geneset_extractors.core.qc import write_run_summary_files
from geneset_extractors.core.selection import global_l1_weights, ranked_gene_ids, select_quantile, select_threshold, select_top_k, within_set_l1_weights
from geneset_extractors.extractors.morphology.mapping import accumulate_gene_scores, bounded_reference_gene_mapping, build_reference_gene_maps, l1_normalize_scores
from geneset_extractors.extractors.morphology.similarity import align_feature_space, pairwise_similarity, standardize_profiles


_SAFE_COMPONENT_RE = re.compile(r"[^A-Za-z0-9._=-]+")
_ARTIFACT_PATTERNS = [re.compile(r"^GPR[0-9A-Z]*$"), re.compile(r"^HTR[0-9A-Z]*$"), re.compile(r"^OR[0-9A-Z]+$")]


@dataclass
class MorphologyWorkflowConfig:
    converter_name: str
    out_dir: Path
    organism: str
    genome_build: str
    dataset_label: str
    signature_name: str
    mode: str
    similarity_metric: str
    similarity_power: float
    polarity: str
    query_modality_column: str
    same_modality_first: bool
    cross_modality_penalty: float
    max_reference_neighbors: int
    adaptive_neighbors: bool
    min_effective_neighbors: int
    neighbor_evidence_drop_ratio: float
    mutual_neighbor_filter: bool
    min_similarity: float
    control_calibration: str
    control_residual_components: int
    control_min_profiles_for_residualization: int
    hubness_penalty: str
    gene_recurrence_penalty: str
    compound_weight: float
    genetic_weight: float
    select: str
    top_k: int
    quantile: float
    min_score: float
    normalize: str
    emit_full: bool
    emit_gmt: bool
    gmt_out: str | None
    gmt_prefer_symbol: bool
    gmt_require_symbol: bool
    gmt_biotype_allowlist: str
    gmt_min_genes: int
    gmt_max_genes: int
    gmt_topk_list: str
    gmt_mass_list: str
    gmt_format: str
    emit_small_gene_sets: bool
    min_specificity_confidence_to_emit_opposite: str


@dataclass
class MorphologyBranch:
    name: str
    neighbors: list[tuple[str, float, float]]
    evidence_by_ref: dict[str, float]
    gene_scores: dict[str, float]


def _bundle_gene_universe(
    *,
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
) -> set[str]:
    genes: set[str] = set()
    for mapping in list(compound_map.values()) + list(genetic_map.values()):
        genes.update(mapping)
    return genes


def _safe_component(value: str, fallback: str) -> str:
    out = _SAFE_COMPONENT_RE.sub("_", str(value)).strip("_")
    return out or fallback


def _select_gene_ids(scores: dict[str, float], cfg: MorphologyWorkflowConfig) -> list[str]:
    if cfg.select == "none":
        return ranked_gene_ids(scores)
    if cfg.select == "top_k":
        return select_top_k(scores, int(cfg.top_k))
    if cfg.select == "quantile":
        return select_quantile(scores, float(cfg.quantile))
    if cfg.select == "threshold":
        return select_threshold(scores, float(cfg.min_score))
    raise ValueError(f"Unsupported selection method: {cfg.select}")


def _selected_weights(scores: dict[str, float], selected_gene_ids: list[str], normalize: str) -> dict[str, float]:
    if normalize == "none":
        return {gene_id: float(scores.get(gene_id, 0.0)) for gene_id in selected_gene_ids}
    if normalize == "l1":
        global_weights = global_l1_weights(scores)
        return {gene_id: float(global_weights.get(gene_id, 0.0)) for gene_id in selected_gene_ids}
    if normalize == "within_set_l1":
        return within_set_l1_weights(scores, selected_gene_ids)
    raise ValueError(f"Unsupported normalization method: {normalize}")


def _confidence_rank(level: str) -> int:
    return {"low": 0, "medium": 1, "high": 2}.get(str(level).strip().lower(), -1)


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["gene_id", "score", "rank"]
    if any("weight" in row for row in rows):
        fieldnames.append("weight")
    if any(str(row.get("gene_symbol", "")).strip() for row in rows):
        fieldnames.append("gene_symbol")
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _warn_gmt_diagnostics(gmt_diagnostics: list[dict[str, object]]) -> list[dict[str, object]]:
    warnings_payload: list[dict[str, object]] = []
    for diag in gmt_diagnostics:
        code = str(diag.get("code", ""))
        if code not in {"small_gene_set_skipped", "small_gene_set_emitted", "marginal_gene_count", "no_positive_genes", "require_symbol_heavy_drop"}:
            continue
        message = f"{code}: {diag.get('base_name','unknown')}"
        if code == "small_gene_set_skipped":
            message = (
                f"{code} program={diag.get('query_id', 'unknown')} polarity={diag.get('polarity', 'unknown')} "
                f"n_genes={int(diag.get('n_genes', 0) or 0)} min_required={int(diag.get('min_genes', 0) or 0)}"
            )
        print(f"warning: {message}", file=sys.stderr)
        warnings_payload.append(
            {
                "code": code,
                "message": message,
                "n_genes": int(diag.get("n_genes", 0) or 0),
                "min_genes": int(diag.get("min_genes", 0) or 0),
                "reason": str(diag.get("reason", "")),
            }
        )
    return warnings_payload


def _collect_skipped_gmt_outputs(gmt_diagnostics: list[dict[str, object]]) -> list[dict[str, object]]:
    return [
        {
            "reason": str(diag.get("reason", "")),
            "code": str(diag.get("code", "")),
            "n_genes": int(diag.get("n_genes", 0) or 0),
        }
        for diag in gmt_diagnostics
        if str(diag.get("code", "")) in {"small_gene_set_skipped", "no_positive_genes"}
    ]


def _artifact_family_fraction(rows: list[dict[str, object]], top_n: int = 20) -> float:
    top_rows = rows[: min(top_n, len(rows))]
    if not top_rows:
        return 0.0
    hits = 0
    for row in top_rows:
        symbol = str(row.get("gene_symbol", "") or row.get("gene_id", "")).strip().upper()
        if any(pat.match(symbol) for pat in _ARTIFACT_PATTERNS):
            hits += 1
    return float(hits) / float(len(top_rows))


def _candidate_neighbors(
    query_sim: dict[str, float],
    *,
    polarity: str,
    min_similarity: float,
) -> list[tuple[str, float]]:
    retained: list[tuple[str, float]] = []
    for ref_id, sim in query_sim.items():
        sim_value = float(sim)
        if polarity == "similar":
            if sim_value < float(min_similarity):
                continue
            score = sim_value
        else:
            score = -sim_value
            if score < float(min_similarity):
                continue
        if score <= 0.0:
            continue
        retained.append((ref_id, score))
    return retained


def _retrieval_confidence(
    n_positive_neighbors: int,
    neg_similarity_fraction: float,
) -> str:
    if n_positive_neighbors >= 10 and neg_similarity_fraction < 0.4:
        return "high"
    if n_positive_neighbors >= 5 and neg_similarity_fraction < 0.7:
        return "medium"
    return "low"


def _specificity_confidence(
    selected_rows: list[dict[str, object]],
    gene_scores: dict[str, float],
    *,
    neighbor_primary_target_agreement: float,
    neighbor_top3_target_agreement: float,
    high_hub_mass_fraction: float,
) -> tuple[str, float, float, float, float, float]:
    if not selected_rows or not gene_scores:
        return "low", 0.0, 0.0, 0.0, 0.0, 0.0
    top10_gene_mass = sum(float(row.get("weight", 0.0)) for row in selected_rows[:10])
    total_gene_mass = sum(float(v) for v in gene_scores.values())
    neighbor_target_concentration = 0.0
    if total_gene_mass > 0.0:
        neighbor_target_concentration = max(float(v) for v in gene_scores.values()) / total_gene_mass
    strong_gene_concentration = top10_gene_mass >= 0.65
    strong_target_agreement = (
        neighbor_target_concentration >= 0.25
        or neighbor_primary_target_agreement >= 0.35
        or neighbor_top3_target_agreement >= 0.5
    )
    low_hub_dominance = high_hub_mass_fraction <= 0.5
    if strong_gene_concentration and strong_target_agreement and low_hub_dominance:
        return "high", top10_gene_mass, neighbor_target_concentration, neighbor_primary_target_agreement, neighbor_top3_target_agreement, high_hub_mass_fraction
    if (strong_gene_concentration or strong_target_agreement) and high_hub_mass_fraction <= 0.75:
        return "medium", top10_gene_mass, neighbor_target_concentration, neighbor_primary_target_agreement, neighbor_top3_target_agreement, high_hub_mass_fraction
    return "low", top10_gene_mass, neighbor_target_concentration, neighbor_primary_target_agreement, neighbor_top3_target_agreement, high_hub_mass_fraction


def _label_vote_summary_from_neighbors(
    *,
    neighbors: list[tuple[str, float, float]],
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    reference_metadata: dict[str, dict[str, str]],
    target_annotations: dict[str, dict[str, str]] | None,
    field: str,
    query_modality: str | None,
) -> dict[str, object]:
    if not target_annotations:
        return {
            "top_label": None,
            "top_label_mass": 0.0,
            "votes": [],
            "top_label_support_count": 0,
            "top_label_same_modality_count": 0,
        }
    label_mass: dict[str, float] = {}
    label_refs: dict[str, set[str]] = {}
    label_same_modality_refs: dict[str, set[str]] = {}
    total_mass = 0.0
    for ref_id, _base, evidence in neighbors:
        mapping, meta = bounded_reference_gene_mapping(
            ref_id=ref_id,
            compound_map=compound_map,
            genetic_map=genetic_map,
        )
        if not mapping:
            continue
        ref_modality = _reference_modality_for_same_modality(
            ref_id=ref_id,
            reference_metadata=reference_metadata,
            mapping_modality=str(meta.get("modality") or ""),
        )
        ref_labels: set[str] = set()
        label_mass_by_ref: dict[str, float] = {}
        for gene, weight in mapping.items():
            label = str((target_annotations.get(gene, {}) or {}).get(field, "")).strip()
            if not label:
                continue
            contribution = float(evidence) * float(weight)
            if contribution <= 0.0:
                continue
            label_mass_by_ref[label] = float(label_mass_by_ref.get(label, 0.0)) + contribution
            ref_labels.add(label)
        for label, contribution in label_mass_by_ref.items():
            bounded_contribution = min(float(evidence), float(contribution))
            label_mass[label] = float(label_mass.get(label, 0.0)) + bounded_contribution
            total_mass += bounded_contribution
        for label in ref_labels:
            label_refs.setdefault(label, set()).add(ref_id)
            if _is_same_modality_pair(query_modality=query_modality, ref_modality=ref_modality):
                label_same_modality_refs.setdefault(label, set()).add(ref_id)
    if not label_mass or total_mass <= 0.0:
        return {
            "top_label": None,
            "top_label_mass": 0.0,
            "votes": [],
            "top_label_support_count": 0,
            "top_label_same_modality_count": 0,
        }
    ordered = sorted(label_mass.items(), key=lambda item: (-float(item[1]), str(item[0])))
    top_label, top_mass_raw = ordered[0]
    votes = [
        {
            "label": label,
            "support_mass": float(mass),
            "support_fraction": float(mass) / float(total_mass),
            "support_count": len(label_refs.get(label, set())),
            "same_modality_support_count": len(label_same_modality_refs.get(label, set())),
        }
        for label, mass in ordered[:5]
    ]
    return {
        "top_label": top_label,
        "top_label_mass": float(top_mass_raw) / float(total_mass),
        "votes": votes,
        "top_label_support_count": len(label_refs.get(top_label, set())),
        "top_label_same_modality_count": len(label_same_modality_refs.get(top_label, set())),
    }


def _normalize_modality_token(value: str | None) -> str | None:
    token = str(value or "").strip().lower()
    if not token:
        return None
    if token in {"compound", "cmpd", "small_molecule"}:
        return "compound"
    if token in {"orf", "overexpression"}:
        return "orf"
    if token in {"crispr", "crispra", "crispri", "knockout"}:
        return "crispr"
    if token == "genetic":
        return "genetic"
    return token


def _reference_modality_for_same_modality(
    *,
    ref_id: str,
    reference_metadata: dict[str, dict[str, str]],
    mapping_modality: str,
) -> str | None:
    ref_meta_modality = _normalize_modality_token(reference_metadata.get(ref_id, {}).get("perturbation_type", ""))
    if ref_meta_modality:
        return ref_meta_modality
    return _normalize_modality_token(mapping_modality)


def _is_same_modality_pair(
    *,
    query_modality: str | None,
    ref_modality: str | None,
) -> bool:
    query_token = _normalize_modality_token(query_modality)
    ref_token = _normalize_modality_token(ref_modality)
    if not query_token or not ref_token:
        return False
    if query_token == ref_token:
        return True
    if query_token == "compound" or ref_token == "compound":
        return False
    return False


def _reference_supports_nominal_target(
    *,
    ref_id: str,
    query_nominal_genes: set[str],
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
) -> bool:
    if not query_nominal_genes:
        return False
    mapping, _meta = bounded_reference_gene_mapping(
        ref_id=ref_id,
        compound_map=compound_map,
        genetic_map=genetic_map,
    )
    return bool(query_nominal_genes & set(mapping))


def _protected_same_modality_reference_ids(
    *,
    raw_candidate_neighbors: list[tuple[str, float]],
    reference_metadata: dict[str, dict[str, str]],
    query_modality: str | None,
    query_nominal_genes: set[str],
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    same_modality_top_n: int,
    raw_similarity_ratio: float,
) -> set[str]:
    if not raw_candidate_neighbors:
        return set()
    protected: set[str] = set()
    top_raw = float(raw_candidate_neighbors[0][1])
    same_modality_seen = 0
    for ref_id, base in raw_candidate_neighbors:
        ref_modality = _normalize_modality_token(reference_metadata.get(ref_id, {}).get("perturbation_type", ""))
        if query_modality and _is_same_modality_pair(query_modality=query_modality, ref_modality=ref_modality):
            if same_modality_seen < same_modality_top_n or float(base) >= top_raw * raw_similarity_ratio:
                protected.add(ref_id)
            same_modality_seen += 1
        if _reference_supports_nominal_target(
            ref_id=ref_id,
            query_nominal_genes=query_nominal_genes,
            compound_map=compound_map,
            genetic_map=genetic_map,
        ):
            protected.add(ref_id)
    return protected


def _reference_labels_for_level(
    *,
    ref_id: str,
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    target_annotations: dict[str, dict[str, str]] | None,
    field: str,
) -> set[str]:
    if not target_annotations:
        return set()
    mapping, _meta = bounded_reference_gene_mapping(
        ref_id=ref_id,
        compound_map=compound_map,
        genetic_map=genetic_map,
    )
    return {
        str((target_annotations.get(gene, {}) or {}).get(field, "")).strip()
        for gene in mapping
        if str((target_annotations.get(gene, {}) or {}).get(field, "")).strip()
    }


def _coherent_same_modality_compound_prelabel_ids(
    *,
    raw_candidate_neighbors: list[tuple[str, float]],
    reference_metadata: dict[str, dict[str, str]],
    query_modality: str | None,
    query_nominal_genes: set[str],
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    target_annotations: dict[str, dict[str, str]] | None,
    qc_weights: dict[str, float],
    same_modality_top_n: int,
    raw_similarity_ratio: float,
) -> set[str]:
    if query_modality != "compound" or not raw_candidate_neighbors or not target_annotations:
        return set()
    top_raw = float(raw_candidate_neighbors[0][1])
    eligible_neighbors: list[tuple[str, float, float]] = []
    same_modality_seen = 0
    for ref_id, base in raw_candidate_neighbors[: max(40, same_modality_top_n * 8)]:
        ref_modality = _normalize_modality_token(reference_metadata.get(ref_id, {}).get("perturbation_type", ""))
        if not _is_same_modality_pair(query_modality=query_modality, ref_modality=ref_modality):
            continue
        if same_modality_seen < same_modality_top_n or float(base) >= top_raw * raw_similarity_ratio:
            eligible_neighbors.append((ref_id, float(base), float(base) * float(qc_weights.get(ref_id, 1.0))))
        same_modality_seen += 1
    if len(eligible_neighbors) < 2:
        return set()
    levels = ["pathway_seed", "target_class", "mechanism_label", "target_family"]
    chosen_level = None
    chosen_label = None
    best_priority = (-1, -1.0, -1, "")
    for level in levels:
        payload = _score_label_from_neighbors(
            neighbors=eligible_neighbors,
            compound_map=compound_map,
            genetic_map=genetic_map,
            reference_metadata=reference_metadata,
            target_annotations=target_annotations,
            field=level,
            query_modality=query_modality,
        )
        label = str(payload.get("top_label") or "")
        if not label:
            continue
        support_count = int(payload.get("top_label_support_count", 0) or 0)
        same_modality_count = int(payload.get("top_label_same_modality_count", 0) or 0)
        label_mass = float(payload.get("top_label_mass", 0.0) or 0.0)
        query_labels = _query_labels_for_level(
            query_nominal_genes=query_nominal_genes,
            target_annotations=target_annotations,
            field=level,
        )
        query_consistent = label in query_labels if query_labels else False
        if query_labels and not query_consistent:
            continue
        if same_modality_count < 2:
            continue
        if not query_consistent and support_count < 3:
            continue
        if not query_consistent and label_mass < 0.32:
            continue
        if query_consistent and label_mass < 0.16:
            continue
        priority = (1 if query_consistent else 0, label_mass, support_count, label)
        if priority > best_priority:
            best_priority = priority
            chosen_level = level
            chosen_label = label
    if not chosen_level or not chosen_label:
        return set()
    return {
        ref_id
        for ref_id, _base, _evidence in eligible_neighbors
        if chosen_label
        in _reference_labels_for_level(
            ref_id=ref_id,
            compound_map=compound_map,
            genetic_map=genetic_map,
            target_annotations=target_annotations,
            field=chosen_level,
        )
    }


def _recurrence_weight(
    gene: str,
    *,
    gene_idf: dict[str, float],
    mode: str,
    stage: str,
) -> float:
    if mode == "none":
        return 1.0
    idf = float(gene_idf.get(gene, 1.0))
    if stage == "nomination":
        return idf
    return idf


def _combine_modal_gene_scores(
    *,
    evidence_by_ref: dict[str, float],
    query_modality: str | None,
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    cfg: MorphologyWorkflowConfig,
    gene_idf: dict[str, float],
) -> dict[str, float]:
    compound_scores: dict[str, float] = {}
    genetic_scores: dict[str, float] = {}
    for ref_id, evidence in evidence_by_ref.items():
        if float(evidence) <= 0.0:
            continue
        mapping, meta = bounded_reference_gene_mapping(
            ref_id=ref_id,
            compound_map=compound_map,
            genetic_map=genetic_map,
        )
        if not mapping:
            continue
        target_scores = compound_scores if str(meta.get("modality")) == "compound" else genetic_scores
        for gene, weight in mapping.items():
            target_scores[gene] = float(target_scores.get(gene, 0.0)) + float(evidence) * float(weight)
    compound_norm = l1_normalize_scores(compound_scores)
    genetic_norm = l1_normalize_scores(genetic_scores)
    compound_neighbors = sum(1 for ref_id in evidence_by_ref if ref_id in compound_map)
    genetic_neighbors = sum(1 for ref_id in evidence_by_ref if ref_id in genetic_map)
    if compound_norm and genetic_norm:
        cw = float(cfg.compound_weight)
        gw = float(cfg.genetic_weight)
        if query_modality == "compound":
            cw *= 1.2
            gw *= 0.8
        elif query_modality in {"orf", "crispr"}:
            gw *= 1.5
            cw *= 0.5
        if compound_neighbors < 5 <= genetic_neighbors:
            cw *= 0.25
        if genetic_neighbors < 5 <= compound_neighbors:
            gw *= 0.25
        total = cw + gw
        cw = cw / total if total > 0.0 else 0.5
        gw = gw / total if total > 0.0 else 0.5
    elif compound_norm:
        cw, gw = 1.0, 0.0
    elif genetic_norm:
        cw, gw = 0.0, 1.0
    else:
        cw, gw = 0.0, 0.0
    genes = sorted(set(compound_norm) | set(genetic_norm))
    scores = {
        gene: float(cw * compound_norm.get(gene, 0.0) + gw * genetic_norm.get(gene, 0.0))
        for gene in genes
        if float(cw * compound_norm.get(gene, 0.0) + gw * genetic_norm.get(gene, 0.0)) > 0.0
    }
    if cfg.gene_recurrence_penalty != "none":
        scores = {
            gene: float(score) * _recurrence_weight(
                gene,
                gene_idf=gene_idf,
                mode=cfg.gene_recurrence_penalty,
                stage="ranking",
            )
            for gene, score in scores.items()
            if float(score) > 0.0
        }
    return {gene: float(score) for gene, score in scores.items() if float(score) > 0.0}


def _build_core_branch(
    *,
    pooled_neighbors: list[tuple[str, float, float]],
    allowed_genes: set[str],
    target_support_weighted: dict[str, float],
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    cfg: MorphologyWorkflowConfig,
    gene_idf: dict[str, float],
) -> MorphologyBranch:
    core_neighbors = _refs_supporting_any_gene(
        penalized_neighbors=pooled_neighbors,
        compound_map=compound_map,
        genetic_map=genetic_map,
        allowed_genes=allowed_genes,
    )
    if int(cfg.max_reference_neighbors) > 0:
        core_neighbors = core_neighbors[: int(cfg.max_reference_neighbors)]
    evidence_by_ref = {ref_id: float(evidence) for ref_id, _base, evidence in core_neighbors}
    direct_scores = {
        gene: float(target_support_weighted.get(gene, 0.0))
        for gene in allowed_genes
        if float(target_support_weighted.get(gene, 0.0)) > 0.0
    }
    total_direct = sum(float(v) for v in direct_scores.values())
    gene_scores = {
        gene: float(score) / total_direct
        for gene, score in direct_scores.items()
        if total_direct > 0.0 and float(score) > 0.0
    }
    return MorphologyBranch(
        name="core",
        neighbors=core_neighbors,
        evidence_by_ref=evidence_by_ref,
        gene_scores=gene_scores,
    )


def _build_mechanism_branch(
    *,
    pooled_neighbors: list[tuple[str, float, float]],
    query_modality: str | None,
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    cfg: MorphologyWorkflowConfig,
    gene_idf: dict[str, float],
) -> MorphologyBranch:
    mechanism_neighbors = list(pooled_neighbors)
    if int(cfg.max_reference_neighbors) > 0:
        mechanism_neighbors = mechanism_neighbors[: int(cfg.max_reference_neighbors)]
    evidence_by_ref = {ref_id: float(evidence) for ref_id, _base, evidence in mechanism_neighbors}
    gene_scores = _combine_modal_gene_scores(
        evidence_by_ref=evidence_by_ref,
        query_modality=query_modality,
        compound_map=compound_map,
        genetic_map=genetic_map,
        cfg=cfg,
        gene_idf=gene_idf,
    )
    return MorphologyBranch(
        name="mechanism",
        neighbors=mechanism_neighbors,
        evidence_by_ref=evidence_by_ref,
        gene_scores=gene_scores,
    )


def _score_label_from_neighbors(
    *,
    neighbors: list[tuple[str, float, float]],
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    reference_metadata: dict[str, dict[str, str]],
    target_annotations: dict[str, dict[str, str]] | None,
    field: str,
    query_modality: str | None,
) -> dict[str, object]:
    summary = _label_vote_summary_from_neighbors(
        neighbors=neighbors,
        compound_map=compound_map,
        genetic_map=genetic_map,
        reference_metadata=reference_metadata,
        target_annotations=target_annotations,
        field=field,
        query_modality=query_modality,
    )
    votes = list(summary.get("votes", []))
    if not votes:
        return {
            "level": field,
            "top_label": None,
            "top_label_mass": 0.0,
            "top_label_support_count": 0,
            "top_label_same_modality_count": 0,
            "summary": summary,
        }
    specificity_bonus = {
        "pathway_seed": 1.00,
        "target_class": 0.92,
        "mechanism_label": 0.85,
        "target_family": 0.75,
    }.get(field, 0.75)
    top_label = str(summary.get("top_label") or "")
    top_mass = float(summary.get("top_label_mass", 0.0) or 0.0)
    support_count = int(summary.get("top_label_support_count", 0) or 0)
    same_modality_count = int(summary.get("top_label_same_modality_count", 0) or 0)
    modality_factor = 1.0
    if query_modality in {"orf", "crispr"} and same_modality_count <= 0:
        modality_factor = 0.7
    disagreement_penalty = 0.0
    return {
        "level": field,
        "top_label": top_label or None,
        "top_label_mass": top_mass,
        "top_label_support_count": support_count,
        "top_label_same_modality_count": same_modality_count,
        "specificity_bonus": specificity_bonus,
        "modality_factor": modality_factor,
        "disagreement_penalty": disagreement_penalty,
        "summary": summary,
    }


def _query_labels_for_level(
    *,
    query_nominal_genes: set[str],
    target_annotations: dict[str, dict[str, str]] | None,
    field: str,
) -> set[str]:
    if not query_nominal_genes or not target_annotations:
        return set()
    return {
        str((target_annotations.get(gene, {}) or {}).get(field, "")).strip()
        for gene in query_nominal_genes
        if str((target_annotations.get(gene, {}) or {}).get(field, "")).strip()
    }


def _choose_expansion_label(
    *,
    mode: str,
    query_modality: str | None,
    raw_neighbors: list[tuple[str, float, float]],
    prelabel_neighbors: list[tuple[str, float, float]],
    retained_neighbors: list[tuple[str, float, float]],
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    reference_metadata: dict[str, dict[str, str]],
    target_annotations: dict[str, dict[str, str]] | None,
    query_nominal_genes: set[str],
    exact_target_support_present: bool,
) -> tuple[
    dict[str, object],
    dict[str, dict[str, object]],
    dict[str, dict[str, object]],
    dict[str, dict[str, object]],
]:
    decision = {
        "allow_expansion": False,
        "allow_family": False,
        "allow_mechanism": False,
        "reason": "mode_not_expandable",
        "chosen_level": None,
        "chosen_label": None,
        "expansion_confidence": 0.0,
        "raw_retained_mismatch": False,
        "raw_prelabel_mismatch": False,
        "require_same_modality_support": query_modality in {"orf", "crispr"},
    }
    if mode not in {"mechanism", "hybrid"} or not target_annotations:
        return decision, {}, {}, {}
    levels = ["pathway_seed", "target_class", "mechanism_label", "target_family"]
    raw_scores: dict[str, dict[str, object]] = {}
    prelabel_scores: dict[str, dict[str, object]] = {}
    retained_scores: dict[str, dict[str, object]] = {}
    same_modality_raw_neighbors = [
        item
        for item in raw_neighbors
        if _is_same_modality_pair(
            query_modality=query_modality,
            ref_modality=_reference_modality_for_same_modality(
                ref_id=str(item[0]),
                reference_metadata=reference_metadata,
                mapping_modality=str(
                    bounded_reference_gene_mapping(
                        ref_id=str(item[0]),
                        compound_map=compound_map,
                        genetic_map=genetic_map,
                    )[1].get("modality")
                    or ""
                ),
            ),
        )
    ]
    same_modality_raw_scores: dict[str, dict[str, object]] = {}
    compound_like = query_modality == "compound"
    threshold = 0.18 if compound_like else 0.24
    for level in levels:
        raw_payload = _score_label_from_neighbors(
            neighbors=raw_neighbors,
            compound_map=compound_map,
            genetic_map=genetic_map,
            reference_metadata=reference_metadata,
            target_annotations=target_annotations,
            field=level,
            query_modality=query_modality,
        )
        prelabel_payload = _score_label_from_neighbors(
            neighbors=prelabel_neighbors,
            compound_map=compound_map,
            genetic_map=genetic_map,
            reference_metadata=reference_metadata,
            target_annotations=target_annotations,
            field=level,
            query_modality=query_modality,
        )
        retained_payload = _score_label_from_neighbors(
            neighbors=retained_neighbors,
            compound_map=compound_map,
            genetic_map=genetic_map,
            reference_metadata=reference_metadata,
            target_annotations=target_annotations,
            field=level,
            query_modality=query_modality,
        )
        raw_scores[level] = raw_payload
        prelabel_scores[level] = prelabel_payload
        retained_scores[level] = retained_payload
        same_modality_raw_scores[level] = _score_label_from_neighbors(
            neighbors=same_modality_raw_neighbors,
            compound_map=compound_map,
            genetic_map=genetic_map,
            reference_metadata=reference_metadata,
            target_annotations=target_annotations,
            field=level,
            query_modality=query_modality,
        )
    for level in levels:
        raw_payload = raw_scores[level]
        same_modality_raw_payload = same_modality_raw_scores[level]
        prelabel_payload = prelabel_scores[level]
        retained_payload = retained_scores[level]
        label = prelabel_payload.get("top_label")
        if not label:
            continue
        retained_mass = float(prelabel_payload.get("top_label_mass", 0.0) or 0.0)
        support_count = int(prelabel_payload.get("top_label_support_count", 0) or 0)
        same_modality_count = int(prelabel_payload.get("top_label_same_modality_count", 0) or 0)
        raw_label = raw_payload.get("top_label")
        retained_label = retained_payload.get("top_label")
        raw_mass = float(raw_payload.get("top_label_mass", 0.0) or 0.0)
        mismatch = bool(raw_label and raw_label != label)
        retained_mismatch = bool(raw_label and retained_label and raw_label != retained_label)
        if mismatch:
            decision["raw_prelabel_mismatch"] = True
        if retained_mismatch:
            decision["raw_retained_mismatch"] = True
        support_factor = min(1.0, float(support_count) / 3.0)
        modality_factor = 1.0
        if query_modality in {"orf", "crispr"} and same_modality_count <= 0:
            modality_factor = 0.65
        disagreement_factor = 0.75 if mismatch else 1.0
        confidence = retained_mass * support_factor * float(prelabel_payload.get("specificity_bonus", 1.0)) * modality_factor * disagreement_factor
        if raw_mass > 0.0 and raw_label == label:
            confidence *= 1.0 + min(0.15, 0.5 * raw_mass)
        query_labels = _query_labels_for_level(
            query_nominal_genes=query_nominal_genes,
            target_annotations=target_annotations,
            field=level,
        )
        raw_support = int(raw_payload.get("top_label_support_count", 0) or 0)
        raw_same_modality = int(raw_payload.get("top_label_same_modality_count", 0) or 0)
        same_modality_raw_label = str(same_modality_raw_payload.get("top_label") or "")
        same_modality_raw_mass = float(same_modality_raw_payload.get("top_label_mass", 0.0) or 0.0)
        same_modality_raw_support = int(same_modality_raw_payload.get("top_label_support_count", 0) or 0)
        retained_support = int(retained_payload.get("top_label_support_count", 0) or 0)
        raw_query_consistent = bool(raw_label and raw_label in query_labels)
        same_modality_raw_query_consistent = bool(same_modality_raw_label and same_modality_raw_label in query_labels)
        prelabel_query_consistent = bool(label and label in query_labels)
        retained_query_consistent = bool(retained_label and retained_label in query_labels)
        if (
            compound_like
            and same_modality_raw_label
            and label
            and same_modality_raw_label != label
            and (same_modality_raw_query_consistent or not query_labels)
            and same_modality_raw_support >= max(2, support_count + 1)
            and same_modality_raw_mass >= max(0.18, retained_mass * 0.9)
            and support_count <= 2
            and (not prelabel_query_consistent or not retained_query_consistent)
        ):
            fallback_confidence = same_modality_raw_mass * min(1.0, float(same_modality_raw_support) / 4.0) * float(same_modality_raw_payload.get("specificity_bonus", 1.0))
            if fallback_confidence >= threshold:
                decision.update(
                    {
                        "allow_expansion": True,
                        "reason": "raw_same_modality_coherence_guard",
                        "chosen_level": level,
                        "chosen_label": same_modality_raw_label,
                        "expansion_confidence": fallback_confidence,
                        "allow_family": level == "target_family",
                        "allow_mechanism": level in {"mechanism_label", "target_class", "pathway_seed"},
                    }
                )
                return decision, raw_scores, prelabel_scores, retained_scores
        if (
            compound_like
            and not exact_target_support_present
            and same_modality_raw_query_consistent
            and not prelabel_query_consistent
            and same_modality_raw_mass >= 0.15
            and same_modality_raw_support >= 2
            and (retained_support <= 2 or not retained_query_consistent)
        ):
            fallback_confidence = same_modality_raw_mass * min(1.0, float(same_modality_raw_support) / 4.0) * float(same_modality_raw_payload.get("specificity_bonus", 1.0))
            if fallback_confidence >= threshold:
                decision.update(
                    {
                        "allow_expansion": True,
                        "reason": "raw_same_modality_query_consistent_label_fallback",
                        "chosen_level": level,
                        "chosen_label": same_modality_raw_label,
                        "expansion_confidence": fallback_confidence,
                        "allow_family": level == "target_family",
                        "allow_mechanism": level in {"mechanism_label", "target_class", "pathway_seed"},
                    }
                )
                return decision, raw_scores, prelabel_scores, retained_scores
        if (
            compound_like
            and not exact_target_support_present
            and raw_query_consistent
            and not prelabel_query_consistent
            and raw_mass >= 0.15
            and (raw_support >= 3 or raw_same_modality >= 2)
            and (retained_support <= 2 or not retained_query_consistent)
        ):
            fallback_confidence = raw_mass * min(1.0, float(max(raw_support, raw_same_modality)) / 4.0) * float(raw_payload.get("specificity_bonus", 1.0))
            if fallback_confidence >= threshold:
                decision.update(
                    {
                        "allow_expansion": True,
                        "reason": "raw_query_consistent_label_fallback",
                        "chosen_level": level,
                        "chosen_label": raw_label,
                        "expansion_confidence": fallback_confidence,
                        "allow_family": level == "target_family",
                        "allow_mechanism": level in {"mechanism_label", "target_class", "pathway_seed"},
                    }
                )
                return decision, raw_scores, prelabel_scores, retained_scores
        prelabel_payload["expansion_confidence"] = confidence
        if confidence >= threshold:
            decision.update(
                {
                    "allow_expansion": True,
                    "reason": "prelabel_coherent_label_support",
                    "chosen_level": level,
                    "chosen_label": label,
                    "expansion_confidence": confidence,
                    "allow_family": level == "target_family",
                    "allow_mechanism": level in {"mechanism_label", "target_class", "pathway_seed"},
                }
            )
            return decision, raw_scores, prelabel_scores, retained_scores
    decision["raw_retained_mismatch"] = any(
        str((raw_scores.get(level, {}) or {}).get("top_label") or "") != str((retained_scores.get(level, {}) or {}).get("top_label") or "")
        and bool((retained_scores.get(level, {}) or {}).get("top_label"))
        for level in levels
    )
    decision["raw_prelabel_mismatch"] = any(
        str((raw_scores.get(level, {}) or {}).get("top_label") or "") != str((prelabel_scores.get(level, {}) or {}).get("top_label") or "")
        and bool((prelabel_scores.get(level, {}) or {}).get("top_label"))
        for level in levels
    )
    if mode == "mechanism" and query_modality in {"orf", "crispr"}:
        for level in levels:
            query_labels = _query_labels_for_level(
                query_nominal_genes=query_nominal_genes,
                target_annotations=target_annotations,
                field=level,
            )
            if not query_labels:
                continue
            for payload in (prelabel_scores.get(level, {}), retained_scores.get(level, {})):
                label = str(payload.get("top_label") or "")
                if not label or label not in query_labels:
                    continue
                support_count = int(payload.get("top_label_support_count", 0) or 0)
                same_modality_count = int(payload.get("top_label_same_modality_count", 0) or 0)
                label_mass = float(payload.get("top_label_mass", 0.0) or 0.0)
                if same_modality_count <= 0 or support_count <= 0:
                    continue
                fallback_confidence = (
                    label_mass
                    * min(1.0, float(max(support_count, same_modality_count)) / 2.0)
                    * float(payload.get("specificity_bonus", 1.0) or 1.0)
                )
                if fallback_confidence < 0.14:
                    continue
                decision.update(
                    {
                        "allow_expansion": True,
                        "reason": "genetic_query_consistent_label_fallback",
                        "chosen_level": level,
                        "chosen_label": label,
                        "expansion_confidence": fallback_confidence,
                        "allow_family": level == "target_family",
                        "allow_mechanism": level in {"mechanism_label", "target_class", "pathway_seed"},
                    }
                )
                return decision, raw_scores, prelabel_scores, retained_scores
    decision["reason"] = "prelabel_label_support_too_weak"
    return decision, raw_scores, prelabel_scores, retained_scores


def _label_scores_by_modality(
    *,
    neighbors: list[tuple[str, float, float]],
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    reference_metadata: dict[str, dict[str, str]],
    target_annotations: dict[str, dict[str, str]] | None,
    query_modality: str | None,
) -> dict[str, dict[str, dict[str, object]]]:
    split_neighbors = {
        "all": neighbors,
        "compound": [
            item for item in neighbors
            if str(reference_metadata.get(item[0], {}).get("perturbation_type", "")).strip().lower() == "compound"
        ],
        "genetic": [
            item for item in neighbors
            if str(reference_metadata.get(item[0], {}).get("perturbation_type", "")).strip().lower() in {"orf", "crispr"}
        ],
    }
    out: dict[str, dict[str, dict[str, object]]] = {}
    for modality, modality_neighbors in split_neighbors.items():
        out[modality] = {}
        for level in ["pathway_seed", "target_class", "mechanism_label", "target_family"]:
            out[modality][level] = _score_label_from_neighbors(
                neighbors=modality_neighbors,
                compound_map=compound_map,
                genetic_map=genetic_map,
                reference_metadata=reference_metadata,
                target_annotations=target_annotations,
                field=level,
                query_modality=query_modality,
            )
    return out


def _candidate_genes_for_label(
    *,
    chosen_level: str,
    chosen_label: str,
    target_annotations: dict[str, dict[str, str]] | None,
    bundle_genes: set[str],
    local_neighbor_genes: set[str],
) -> tuple[list[str], str]:
    if not target_annotations or not chosen_level or not chosen_label:
        return [], "none"
    bundle_candidates = sorted(
        gene
        for gene in bundle_genes
        if str((target_annotations.get(gene, {}) or {}).get(chosen_level, "")).strip() == chosen_label
    )
    if bundle_candidates:
        return bundle_candidates, "bundle"
    local_candidates = sorted(
        gene
        for gene in local_neighbor_genes
        if str((target_annotations.get(gene, {}) or {}).get(chosen_level, "")).strip() == chosen_label
    )
    if local_candidates:
        return local_candidates, "local"
    global_candidates = sorted(
        gene
        for gene, row in target_annotations.items()
        if str((row or {}).get(chosen_level, "")).strip() == chosen_label
    )
    return global_candidates, "global"


def _apply_local_weighted_expansion(
    *,
    base_gene_scores: dict[str, float],
    branch_neighbors: list[tuple[str, float, float]],
    chosen_level: str | None,
    chosen_label: str | None,
    expansion_confidence: float,
    target_annotations: dict[str, dict[str, str]] | None,
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    full_bundle_compound_map: dict[str, dict[str, float]],
    full_bundle_genetic_map: dict[str, dict[str, float]],
    query_nominal_genes: set[str],
) -> tuple[dict[str, float], list[str], dict[str, object]]:
    if not chosen_level or not chosen_label or not target_annotations:
        return dict(base_gene_scores), [], {
            "candidate_universe_size": 0,
            "candidate_scope": "none",
            "local_supported_genes": 0,
            "bundle_supported_genes": 0,
            "bundle_candidate_genes": 0,
            "global_fallback_genes": 0,
            "chosen_level": None,
            "chosen_label": None,
            "expansion_mass_injected": 0.0,
        }
    bundle_genes = _bundle_gene_universe(compound_map=full_bundle_compound_map, genetic_map=full_bundle_genetic_map)
    local_neighbor_genes: set[str] = set()
    local_gene_support: dict[str, float] = {}
    for ref_id, _base, evidence in branch_neighbors:
        mapping = compound_map.get(ref_id) or genetic_map.get(ref_id) or {}
        for gene, weight in mapping.items():
            if str((target_annotations.get(gene, {}) or {}).get(chosen_level, "")).strip() != chosen_label:
                continue
            local_neighbor_genes.add(gene)
            local_gene_support[gene] = float(local_gene_support.get(gene, 0.0)) + float(evidence) * float(weight)
    candidate_genes, scope = _candidate_genes_for_label(
        chosen_level=chosen_level,
        chosen_label=chosen_label,
        target_annotations=target_annotations,
        bundle_genes=bundle_genes,
        local_neighbor_genes=local_neighbor_genes,
    )
    nominal_matching = sorted(
        gene
        for gene in query_nominal_genes
        if str((target_annotations.get(gene, {}) or {}).get(chosen_level, "")).strip() == chosen_label
    )
    global_only_candidates_added = 0
    if chosen_level in {"pathway_seed", "target_class", "mechanism_label"}:
        global_candidates = sorted(
            gene
            for gene, row in target_annotations.items()
            if str((row or {}).get(chosen_level, "")).strip() == chosen_label
        )
        extra = [gene for gene in global_candidates if gene not in candidate_genes]
        global_only_candidates_added = len(extra)
        candidate_genes = sorted(set(candidate_genes) | set(global_candidates))
    query_nominal_genes_dropped_by_scope = sorted(gene for gene in nominal_matching if gene not in candidate_genes)
    if nominal_matching:
        candidate_genes = sorted(set(candidate_genes) | set(nominal_matching))
        query_nominal_genes_dropped_by_scope = []
    bundle_candidate_genes = sum(
        1
        for gene in bundle_genes
        if str((target_annotations.get(gene, {}) or {}).get(chosen_level, "")).strip() == chosen_label
    )
    if not candidate_genes:
        return dict(base_gene_scores), [], {
            "candidate_universe_size": 0,
            "candidate_scope": "none",
            "local_supported_genes": len(local_neighbor_genes),
            "bundle_supported_genes": 0,
            "bundle_candidate_genes": bundle_candidate_genes,
            "global_fallback_genes": 0,
            "global_only_candidates_added": global_only_candidates_added,
            "query_nominal_genes_matching_label": nominal_matching,
            "query_nominal_genes_dropped_by_scope": query_nominal_genes_dropped_by_scope,
            "bundle_gene_universe_source": "full_bundle",
            "chosen_level": chosen_level,
            "chosen_label": chosen_label,
            "expansion_mass_injected": 0.0,
        }
    bundle_label_counts: dict[str, int] = {}
    for ref_id in set(compound_map) | set(genetic_map):
        mapping = compound_map.get(ref_id) or genetic_map.get(ref_id) or {}
        for gene in mapping:
            if gene in candidate_genes:
                bundle_label_counts[gene] = int(bundle_label_counts.get(gene, 0)) + 1
    weighted_candidates: dict[str, float] = {}
    for gene in candidate_genes:
        local = float(local_gene_support.get(gene, 0.0))
        bundle = float(bundle_label_counts.get(gene, 0))
        weighted_candidates[gene] = local + 0.25 * bundle
    for gene in nominal_matching:
        if gene in candidate_genes:
            weighted_candidates[gene] = max(float(weighted_candidates.get(gene, 0.0)), 0.5)
    total_weight = sum(float(v) for v in weighted_candidates.values())
    base_total = sum(float(v) for v in base_gene_scores.values())
    expansion_mass = max(1e-6, 0.35 * base_total * float(expansion_confidence))
    boosted = dict(base_gene_scores)
    expanded_genes: list[str] = []
    if total_weight <= 0.0:
        return boosted, [], {
            "candidate_universe_size": len(candidate_genes),
            "candidate_scope": scope,
            "local_supported_genes": len(local_neighbor_genes),
            "bundle_supported_genes": sum(1 for gene in candidate_genes if gene in bundle_genes),
            "bundle_candidate_genes": bundle_candidate_genes,
            "global_fallback_genes": len(candidate_genes) if scope == "global" else 0,
            "global_only_candidates_added": global_only_candidates_added,
            "query_nominal_genes_matching_label": nominal_matching,
            "query_nominal_genes_dropped_by_scope": query_nominal_genes_dropped_by_scope,
            "bundle_gene_universe_source": "full_bundle",
            "chosen_level": chosen_level,
            "chosen_label": chosen_label,
            "expansion_mass_injected": 0.0,
        }
    for gene, raw_weight in sorted(weighted_candidates.items(), key=lambda item: (-float(item[1]), str(item[0]))):
        if float(raw_weight) <= 0.0:
            continue
        boosted[gene] = float(boosted.get(gene, 0.0)) + expansion_mass * (float(raw_weight) / total_weight)
        expanded_genes.append(gene)
    return boosted, sorted(set(expanded_genes)), {
        "candidate_universe_size": len(candidate_genes),
        "candidate_scope": scope,
        "local_supported_genes": len(local_neighbor_genes),
        "bundle_supported_genes": sum(1 for gene in candidate_genes if gene in bundle_genes),
        "bundle_candidate_genes": bundle_candidate_genes,
        "global_fallback_genes": len(candidate_genes) if scope == "global" else 0,
        "global_only_candidates_added": global_only_candidates_added,
        "query_nominal_genes_matching_label": nominal_matching,
        "query_nominal_genes_dropped_by_scope": query_nominal_genes_dropped_by_scope,
        "bundle_gene_universe_source": "full_bundle",
        "chosen_level": chosen_level,
        "chosen_label": chosen_label,
        "expansion_mass_injected": expansion_mass,
    }


def _rows_from_scores(scores: dict[str, float], cfg: MorphologyWorkflowConfig) -> list[dict[str, object]]:
    selected_gene_ids = sorted(_select_gene_ids(scores, cfg), key=lambda g: (-float(scores.get(g, 0.0)), str(g)))
    weights = _selected_weights(scores, selected_gene_ids, cfg.normalize)
    return [
        {
            "gene_id": gene_id,
            "gene_symbol": gene_id,
            "score": float(scores.get(gene_id, 0.0)),
            "weight": float(weights.get(gene_id, 0.0)),
            "rank": rank,
        }
        for rank, gene_id in enumerate(selected_gene_ids, start=1)
    ]


def _choose_preferred_variant(
    *,
    mode: str,
    core_rows: list[dict[str, object]],
    expanded_rows: list[dict[str, object]],
    top_target_mass: float,
    expansion_confidence: float,
) -> tuple[str, str]:
    if mode == "direct_target":
        return "core", "exact_target_mode"
    if mode == "mechanism":
        return "expanded", "mechanism_mode"
    if expansion_confidence >= 0.22 and len(expanded_rows) > len(core_rows):
        return "expanded", "expansion_confidence"
    if top_target_mass >= 0.5 and core_rows:
        return "core", "exact_target_confidence"
    if expanded_rows:
        return "expanded", "fallback_expanded"
    return "core", "fallback_core"


def _merge_hybrid_scores(
    *,
    core_gene_scores: dict[str, float],
    expanded_gene_scores: dict[str, float],
) -> tuple[dict[str, float], dict[str, int]]:
    merged: dict[str, float] = {}
    for gene in sorted(set(core_gene_scores) | set(expanded_gene_scores)):
        core = float(core_gene_scores.get(gene, 0.0))
        expanded = float(expanded_gene_scores.get(gene, 0.0))
        if core > 0.0 and expanded > 0.0:
            merged[gene] = max(core, 0.75 * core + 0.25 * expanded)
        elif core > 0.0:
            merged[gene] = core
        elif expanded > 0.0:
            merged[gene] = expanded
    return (
        {gene: score for gene, score in merged.items() if score > 0.0},
        {
            "hybrid_core_gene_count": len(core_gene_scores),
            "hybrid_expanded_gene_count": len(expanded_gene_scores),
            "hybrid_added_only_gene_count": len([gene for gene in expanded_gene_scores if gene not in core_gene_scores]),
        },
    )


def _refs_supporting_any_gene(
    *,
    penalized_neighbors: list[tuple[str, float, float]],
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    allowed_genes: set[str],
) -> list[tuple[str, float, float]]:
    supported: list[tuple[str, float, float]] = []
    for ref_id, base, evidence in penalized_neighbors:
        mapping = compound_map.get(ref_id) or genetic_map.get(ref_id) or {}
        if allowed_genes & set(mapping):
            supported.append((ref_id, base, evidence))
    return supported


def _target_support_from_neighbors(
    *,
    penalized_neighbors: list[tuple[str, float, float]],
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
) -> tuple[dict[str, float], dict[str, int], dict[str, dict[str, int]], dict[str, float]]:
    target_support: dict[str, float] = {}
    target_count: dict[str, int] = {}
    target_by_modality: dict[str, dict[str, int]] = {}
    best_similarity: dict[str, float] = {}
    for ref_id, base, evidence in penalized_neighbors:
        mapping, meta = bounded_reference_gene_mapping(
            ref_id=ref_id,
            compound_map=compound_map,
            genetic_map=genetic_map,
        )
        if not mapping:
            continue
        modality = str(meta.get("modality") or "compound")
        for gene, weight in mapping.items():
            contribution = float(evidence) * float(weight)
            if contribution <= 0.0:
                continue
            target_support[gene] = float(target_support.get(gene, 0.0)) + contribution
            target_count[gene] = int(target_count.get(gene, 0)) + 1
            target_by_modality.setdefault(gene, {})
            target_by_modality[gene][modality] = int(target_by_modality[gene].get(modality, 0)) + 1
            best_similarity[gene] = max(float(best_similarity.get(gene, 0.0)), float(base))
    return target_support, target_count, target_by_modality, best_similarity


def _hubness_weight(
    *,
    ref_id: str,
    reference_metadata: dict[str, dict[str, str]],
    mode: str,
    hub_rank_penalties: dict[str, float],
) -> float:
    if mode == "none":
        return 1.0
    if mode == "inverse_rank":
        return float(hub_rank_penalties.get(ref_id, 1.0))
    raw = reference_metadata.get(ref_id, {}).get("hub_score")
    try:
        hub_score = float(raw) if raw not in {None, ""} else None
    except ValueError:
        hub_score = None
    if hub_score is None:
        return 1.0
    return 1.0 / max(1e-6, float(hub_score))


def _reference_gene_idf(
    *,
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
) -> dict[str, float]:
    ref_maps = [mapping for mapping in list(compound_map.values()) + list(genetic_map.values()) if mapping]
    n_refs = len(ref_maps)
    if n_refs <= 0:
        return {}
    df: dict[str, int] = {}
    for mapping in ref_maps:
        for gene in mapping:
            df[gene] = int(df.get(gene, 0)) + 1
    import math

    denom = math.log(float(n_refs) + 1.0) or 1.0
    return {
        gene: math.log((float(n_refs) + 1.0) / (float(count) + 1.0)) / denom
        for gene, count in df.items()
    }


def _routed_primary_target_agreement(
    *,
    retained_neighbors: list[tuple[str, float]],
    evidence_by_ref: dict[str, float],
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
) -> float:
    target_mass: dict[str, float] = {}
    total_mass = 0.0
    for ref_id, _base in retained_neighbors:
        evidence = float(evidence_by_ref.get(ref_id, 0.0))
        if evidence <= 0.0:
            continue
        mapping = compound_map.get(ref_id) or genetic_map.get(ref_id)
        if not mapping:
            continue
        top_weight = max(float(w) for w in mapping.values())
        top_genes = [gene for gene, weight in mapping.items() if float(weight) == top_weight]
        share = evidence / float(len(top_genes) or 1)
        for gene in top_genes:
            target_mass[gene] = float(target_mass.get(gene, 0.0)) + share
        total_mass += evidence
    if total_mass <= 0.0 or not target_mass:
        return 0.0
    return max(float(v) for v in target_mass.values()) / total_mass


def _routed_topk_target_agreement(
    *,
    retained_neighbors: list[tuple[str, float]],
    evidence_by_ref: dict[str, float],
    compound_map: dict[str, dict[str, float]],
    genetic_map: dict[str, dict[str, float]],
    top_k: int,
) -> float:
    target_mass: dict[str, float] = {}
    total_mass = 0.0
    for ref_id, _base in retained_neighbors:
        evidence = float(evidence_by_ref.get(ref_id, 0.0))
        if evidence <= 0.0:
            continue
        mapping = compound_map.get(ref_id) or genetic_map.get(ref_id)
        if not mapping:
            continue
        ordered = sorted(mapping.items(), key=lambda item: (-float(item[1]), str(item[0])))[: max(1, int(top_k))]
        denom = sum(float(weight) for _gene, weight in ordered) or float(len(ordered))
        for gene, weight in ordered:
            share = evidence * (float(weight) / float(denom))
            target_mass[gene] = float(target_mass.get(gene, 0.0)) + share
        total_mass += evidence
    if total_mass <= 0.0 or not target_mass:
        return 0.0
    top_mass = sum(
        float(weight)
        for _gene, weight in sorted(target_mass.items(), key=lambda item: (-float(item[1]), str(item[0])))[: max(1, int(top_k))]
    )
    return top_mass / total_mass


def _high_hub_mass_fraction(
    *,
    evidence_by_ref: dict[str, float],
    reference_metadata: dict[str, dict[str, str]],
) -> float:
    scored: list[tuple[str, float]] = []
    total_mass = 0.0
    for ref_id, evidence in evidence_by_ref.items():
        if float(evidence) <= 0.0:
            continue
        raw = reference_metadata.get(ref_id, {}).get("hub_score")
        try:
            hub_score = float(raw) if raw not in {None, ""} else None
        except ValueError:
            hub_score = None
        if hub_score is None:
            continue
        scored.append((ref_id, hub_score))
        total_mass += float(evidence)
    if total_mass <= 0.0 or len(scored) < 4:
        return 0.0
    scored.sort(key=lambda item: (float(item[1]), str(item[0])))
    cutoff_index = max(0, int(0.75 * len(scored)))
    high_hub_ids = {ref_id for ref_id, _hub_score in scored[cutoff_index:]}
    high_hub_mass = sum(float(evidence_by_ref.get(ref_id, 0.0)) for ref_id in high_hub_ids)
    return high_hub_mass / total_mass


def _pairwise_ref_similarity(
    ref_id_a: str,
    ref_id_b: str,
    *,
    ref_vectors: dict[str, list[float]],
    cache: dict[tuple[str, str], float],
) -> float:
    key = tuple(sorted((ref_id_a, ref_id_b)))
    if key in cache:
        return float(cache[key])
    vec_a = ref_vectors.get(ref_id_a)
    vec_b = ref_vectors.get(ref_id_b)
    if vec_a is None or vec_b is None:
        cache[key] = 0.0
        return 0.0
    dot = sum(float(a) * float(b) for a, b in zip(vec_a, vec_b, strict=False))
    norm_a = sum(float(a) * float(a) for a in vec_a) ** 0.5
    norm_b = sum(float(b) * float(b) for b in vec_b) ** 0.5
    sim = dot / max(1e-8, norm_a * norm_b)
    cache[key] = float(sim)
    return float(sim)


def _query_modality(
    *,
    query_id: str,
    query_membership: dict[str, list[str]],
    query_metadata_rows: dict[str, dict[str, str]] | None,
    modality_column: str,
) -> str | None:
    if not query_metadata_rows:
        return None
    modes = {
        str(query_metadata_rows.get(member_id, {}).get(modality_column, "")).strip().lower()
        for member_id in query_membership.get(query_id, [])
        if str(query_metadata_rows.get(member_id, {}).get(modality_column, "")).strip()
    }
    if len(modes) == 1:
        return next(iter(modes))
    return None


def _apply_control_mean_center(
    *,
    query_vectors: dict[str, list[float]],
    ref_vectors: dict[str, list[float]],
    reference_metadata: dict[str, dict[str, str]],
) -> tuple[dict[str, list[float]], dict[str, list[float]], dict[str, object] | None]:
    control_ids = [
        ref_id
        for ref_id in ref_vectors
        if str(reference_metadata.get(ref_id, {}).get("is_control", "")).strip().lower() in {"1", "true", "t", "yes"}
    ]
    if not control_ids:
        return query_vectors, ref_vectors, None
    n_features = len(next(iter(ref_vectors.values()))) if ref_vectors else 0
    if n_features <= 0:
        return query_vectors, ref_vectors, None
    control_mean = [0.0] * n_features
    for ref_id in control_ids:
        vec = ref_vectors.get(ref_id, [])
        for idx, value in enumerate(vec):
            control_mean[idx] += float(value)
    control_mean = [value / float(len(control_ids)) for value in control_mean]
    q_out = {
        query_id: [float(value) - float(control_mean[idx]) for idx, value in enumerate(vec)]
        for query_id, vec in query_vectors.items()
    }
    r_out = {
        ref_id: [float(value) - float(control_mean[idx]) for idx, value in enumerate(vec)]
        for ref_id, vec in ref_vectors.items()
    }
    return q_out, r_out, {"mode": "mean_center", "n_control_profiles": len(control_ids)}


def _dot(a: list[float], b: list[float]) -> float:
    return sum(float(x) * float(y) for x, y in zip(a, b, strict=False))


def _l2_norm(vec: list[float]) -> float:
    return math.sqrt(sum(float(x) * float(x) for x in vec))


def _normalize_axis(vec: list[float]) -> list[float] | None:
    norm = _l2_norm(vec)
    if norm <= 1e-8:
        return None
    return [float(x) / norm for x in vec]


def _orthogonalize(vec: list[float], axes: list[list[float]]) -> list[float]:
    out = [float(x) for x in vec]
    for axis in axes:
        proj = _dot(out, axis)
        if abs(proj) <= 1e-12:
            continue
        out = [float(value) - float(proj) * float(axis[idx]) for idx, value in enumerate(out)]
    return out


def _estimate_control_nuisance_axes(
    *,
    control_vectors: list[list[float]],
    n_components: int,
    n_iter: int = 25,
) -> tuple[list[list[float]], list[float]]:
    if not control_vectors or n_components <= 0:
        return [], []
    residual_controls = [[float(x) for x in vec] for vec in control_vectors]
    axes: list[list[float]] = []
    explained: list[float] = []
    total_energy = sum(_dot(vec, vec) for vec in residual_controls)
    if total_energy <= 1e-12:
        return [], []
    n_features = len(residual_controls[0])
    max_components = min(int(n_components), len(residual_controls), n_features)
    for _ in range(max_components):
        seed = [0.0] * n_features
        for vec in residual_controls:
            for idx, value in enumerate(vec):
                seed[idx] += float(value)
        axis = _normalize_axis(seed)
        if axis is None:
            seed = [float(residual_controls[0][idx]) for idx in range(n_features)]
            axis = _normalize_axis(seed)
        if axis is None:
            break
        for _iter in range(n_iter):
            updated = [0.0] * n_features
            for vec in residual_controls:
                proj = _dot(vec, axis)
                if abs(proj) <= 1e-12:
                    continue
                for idx, value in enumerate(vec):
                    updated[idx] += float(proj) * float(value)
            updated = _orthogonalize(updated, axes)
            normalized = _normalize_axis(updated)
            if normalized is None:
                break
            axis = normalized
        axis_energy = 0.0
        next_residuals: list[list[float]] = []
        for vec in residual_controls:
            proj = _dot(vec, axis)
            axis_energy += float(proj) * float(proj)
            next_residuals.append([float(value) - float(proj) * float(axis[idx]) for idx, value in enumerate(vec)])
        if axis_energy <= 1e-10:
            break
        axes.append(axis)
        explained.append(axis_energy / total_energy)
        residual_controls = next_residuals
    return axes, explained


def _apply_control_residualization(
    *,
    query_vectors: dict[str, list[float]],
    ref_vectors: dict[str, list[float]],
    reference_metadata: dict[str, dict[str, str]],
    n_components: int,
    min_profiles: int,
) -> tuple[dict[str, list[float]], dict[str, list[float]], dict[str, object] | None]:
    control_ids = [
        ref_id
        for ref_id in ref_vectors
        if str(reference_metadata.get(ref_id, {}).get("is_control", "")).strip().lower() in {"1", "true", "t", "yes"}
    ]
    if not control_ids:
        return query_vectors, ref_vectors, None
    n_features = len(next(iter(ref_vectors.values()))) if ref_vectors else 0
    if n_features <= 0:
        return query_vectors, ref_vectors, None
    if len(control_ids) < int(min_profiles):
        q_out, r_out, fallback = _apply_control_mean_center(
            query_vectors=query_vectors,
            ref_vectors=ref_vectors,
            reference_metadata=reference_metadata,
        )
        summary = dict(fallback or {})
        summary.update(
            {
                "mode": "mean_center",
                "requested_mode": "residualize_controls",
                "fallback_reason": "too_few_controls",
                "min_profiles_required": int(min_profiles),
            }
        )
        return q_out, r_out, summary
    control_mean = [0.0] * n_features
    for ref_id in control_ids:
        vec = ref_vectors.get(ref_id, [])
        for idx, value in enumerate(vec):
            control_mean[idx] += float(value)
    control_mean = [value / float(len(control_ids)) for value in control_mean]
    centered_controls = [
        [float(value) - float(control_mean[idx]) for idx, value in enumerate(ref_vectors.get(ref_id, []))]
        for ref_id in control_ids
    ]
    axes, explained = _estimate_control_nuisance_axes(
        control_vectors=centered_controls,
        n_components=max(1, int(n_components)),
    )
    if not axes:
        q_out, r_out, fallback = _apply_control_mean_center(
            query_vectors=query_vectors,
            ref_vectors=ref_vectors,
            reference_metadata=reference_metadata,
        )
        summary = dict(fallback or {})
        summary.update(
            {
                "mode": "mean_center",
                "requested_mode": "residualize_controls",
                "fallback_reason": "axis_estimation_failed",
                "requested_components": int(n_components),
            }
        )
        return q_out, r_out, summary

    def transform(vec: list[float]) -> list[float]:
        centered = [float(value) - float(control_mean[idx]) for idx, value in enumerate(vec)]
        return _orthogonalize(centered, axes)

    q_out = {query_id: transform(vec) for query_id, vec in query_vectors.items()}
    r_out = {ref_id: transform(vec) for ref_id, vec in ref_vectors.items()}
    return q_out, r_out, {
        "mode": "residualize_controls",
        "n_control_profiles": len(control_ids),
        "n_components": len(axes),
        "requested_components": int(n_components),
        "min_profiles_required": int(min_profiles),
        "explained_variance_fraction": explained,
    }


def run_morphology_workflow(
    *,
    cfg: MorphologyWorkflowConfig,
    query_profiles: dict[str, dict[str, float]],
    reference_profiles: dict[str, dict[str, float]],
    reference_metadata: dict[str, dict[str, str]],
    compound_targets: dict[str, dict[str, float]],
    target_annotations: dict[str, dict[str, str]] | None,
    feature_schema: list[str] | None,
    feature_stats: dict[str, tuple[float, float]] | None,
    query_membership: dict[str, list[str]],
    parse_summary: dict[str, object],
    input_files: list[dict[str, str]],
    resources_info: dict[str, object] | None,
    exclude_query_ids_from_reference: bool,
    query_metadata_rows: dict[str, dict[str, str]] | None,
) -> dict[str, object]:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    reference_profiles_all = dict(reference_profiles)
    reference_metadata_all = dict(reference_metadata)
    reference_profiles_effective = dict(reference_profiles_all)
    reference_metadata_effective = dict(reference_metadata_all)
    if exclude_query_ids_from_reference:
        overlap = set(query_profiles) & set(reference_profiles_effective)
        for ref_id in overlap:
            reference_profiles_effective.pop(ref_id, None)
            reference_metadata_effective.pop(ref_id, None)
    control_ids = {rid for rid, row in reference_metadata_all.items() if str(row.get("is_control", "")).strip().lower() in {"1", "true", "t", "yes"}}

    shared_features, feature_summary = align_feature_space(query_profiles, reference_profiles_effective, feature_schema=feature_schema)
    if not shared_features:
        raise ValueError("No shared numeric features between query and reference profiles")
    if float(feature_summary.get("fraction_reference_features_shared", 1.0)) < 0.8:
        print(
            f"warning: shared feature fraction is low ({feature_summary['fraction_reference_features_shared']:.3f}); results may be unstable",
            file=sys.stderr,
        )

    query_vectors, query_scale_summary = standardize_profiles(query_profiles, features=shared_features, feature_stats=feature_stats)
    ref_vectors_all, ref_scale_summary = standardize_profiles(reference_profiles_effective | {rid: reference_profiles_all[rid] for rid in control_ids if rid in reference_profiles_all}, features=shared_features, feature_stats=feature_stats)
    control_calibration_summary: dict[str, object] | None = None
    root_fallback_warning: dict[str, object] | None = None
    if cfg.control_calibration == "mean_center":
        query_vectors, ref_vectors_all, control_calibration_summary = _apply_control_mean_center(
            query_vectors=query_vectors,
            ref_vectors=ref_vectors_all,
            reference_metadata=reference_metadata_all,
        )
    elif cfg.control_calibration == "residualize_controls":
        query_vectors, ref_vectors_all, control_calibration_summary = _apply_control_residualization(
            query_vectors=query_vectors,
            ref_vectors=ref_vectors_all,
            reference_metadata=reference_metadata_all,
            n_components=int(cfg.control_residual_components),
            min_profiles=int(cfg.control_min_profiles_for_residualization),
        )
        if control_calibration_summary and str(control_calibration_summary.get("fallback_reason", "")).strip():
            print(
                "warning: control residualization fell back to mean-centering because "
                f"{control_calibration_summary['fallback_reason']}",
                file=sys.stderr,
            )
            root_fallback_warning = {
                "code": "control_residualization_fallback",
                "reason": str(control_calibration_summary["fallback_reason"]),
            }
        else:
            root_fallback_warning = None
    else:
        root_fallback_warning = None
    ref_vectors = dict(ref_vectors_all)
    if control_ids:
        print(f"warning: excluding {len(control_ids)} reference control profiles from similarity candidates", file=sys.stderr)
        for rid in control_ids:
            ref_vectors.pop(rid, None)
            reference_profiles_effective.pop(rid, None)
            reference_metadata_effective.pop(rid, None)
    similarities = pairwise_similarity(query_vectors, ref_vectors, metric=cfg.similarity_metric)

    qc_missing = 0
    qc_zero = 0
    qc_weights: dict[str, float] = {}
    for ref_id in ref_vectors:
        raw = reference_metadata_effective.get(ref_id, {}).get("qc_weight")
        try:
            val = float(raw) if raw not in {None, ""} else 1.0
        except ValueError:
            val = 1.0
            qc_missing += 1
        if val <= 0.0:
            qc_zero += 1
        qc_weights[ref_id] = max(float(val), 0.0)
    if qc_missing or qc_zero:
        print(f"warning: qc_weight missing_or_zero for {qc_missing + qc_zero} reference perturbations", file=sys.stderr)

    full_bundle_compound_map, full_bundle_genetic_map, full_mapping_summary = build_reference_gene_maps(
        reference_metadata=reference_metadata_all,
        compound_targets=compound_targets,
    )
    compound_map, genetic_map, mapping_summary = build_reference_gene_maps(
        reference_metadata=reference_metadata_effective,
        compound_targets=compound_targets,
    )
    gene_idf = _reference_gene_idf(compound_map=compound_map, genetic_map=genetic_map)
    compound_ref_weighting_stats = {
        "compound_target_weighting_mode": "bounded_per_reference",
        "n_compound_refs_renormalized": 0,
        "median_raw_target_count": 0.0,
    }
    raw_target_counts: list[float] = []
    for ref_id in compound_map:
        _mapping, meta = bounded_reference_gene_mapping(
            ref_id=ref_id,
            compound_map=compound_map,
            genetic_map=genetic_map,
        )
        if str(meta.get("modality")) != "compound":
            continue
        raw_target_counts.append(float(meta.get("raw_target_count", 0) or 0))
        if bool(meta.get("renormalized")):
            compound_ref_weighting_stats["n_compound_refs_renormalized"] = int(compound_ref_weighting_stats["n_compound_refs_renormalized"]) + 1
    if raw_target_counts:
        compound_ref_weighting_stats["median_raw_target_count"] = median(raw_target_counts)
    hub_scores_present = 0
    hub_rank_penalties: dict[str, float] = {}
    hub_order: list[tuple[str, float]] = []
    for ref_id, row in reference_metadata_effective.items():
        raw = row.get("hub_score")
        try:
            hub_score = float(raw) if raw not in {None, ""} else None
        except ValueError:
            hub_score = None
        if hub_score is None:
            continue
        hub_scores_present += 1
        hub_order.append((ref_id, hub_score))
    if hub_order:
        hub_order.sort(key=lambda item: (float(item[1]), str(item[0])))
        n_hubs = len(hub_order)
        for rank, (ref_id, _score) in enumerate(hub_order, start=1):
            hub_rank_penalties[ref_id] = 1.0 / (float(rank) / float(n_hubs))

    manifest_rows: list[dict[str, object]] = []
    combined_gmt_sets: list[tuple[str, list[str]]] = []
    root_warnings: list[dict[str, object]] = []
    skipped_programs: list[dict[str, object]] = []
    ref_ref_similarity_cache: dict[tuple[str, str], float] = {}
    gmt_topk_list = parse_int_list_csv(str(cfg.gmt_topk_list))
    gmt_mass_list = parse_mass_list_csv(str(cfg.gmt_mass_list))
    gmt_biotype_allowlist = parse_str_list_csv(str(cfg.gmt_biotype_allowlist))
    polarities = [cfg.polarity] if cfg.polarity in {"similar", "opposite"} else ["similar", "opposite"]
    if cfg.hubness_penalty != "none" and hub_scores_present == 0:
        print("warning: hub_score not present in reference metadata; proceeding without hubness penalty", file=sys.stderr)
        root_warnings.append({"code": "hub_score_missing", "mode": cfg.hubness_penalty})
    if cfg.control_calibration == "residualize_controls" and control_calibration_summary is None:
        print("warning: control residualization requested but no usable controls were available", file=sys.stderr)
        root_warnings.append({"code": "control_residualization_unavailable"})
    elif root_fallback_warning is not None:
        root_warnings.append(root_fallback_warning)
    if cfg.polarity in {"opposite", "both"}:
        print("warning: opposite-polarity morphology programs are experimental and should be interpreted cautiously", file=sys.stderr)
        root_warnings.append({"code": "opposite_polarity_experimental"})
    if cfg.mode in {"mechanism", "hybrid"} and not target_annotations:
        print(
            f"warning: mode={cfg.mode} requested without target annotations; mechanism expansion will be limited to broad target aggregation",
            file=sys.stderr,
        )
        root_warnings.append({"code": "target_annotations_missing_for_mode", "mode": cfg.mode})

    for query_id in sorted(query_vectors):
        query_sim = similarities.get(query_id, {})
        query_modality = _query_modality(
            query_id=query_id,
            query_membership=query_membership,
            query_metadata_rows=query_metadata_rows,
            modality_column=cfg.query_modality_column,
        )
        query_nominal_genes = {
            str(query_metadata_rows.get(member_id, {}).get("gene_symbol", "")).strip().upper()
            for member_id in query_membership.get(query_id, [])
            if query_metadata_rows
            and str(query_metadata_rows.get(member_id, {}).get("gene_symbol", "")).strip()
        }
        neg_frac = 0.0
        if query_sim:
            neg_frac = float(sum(1 for v in query_sim.values() if float(v) < 0.0)) / float(len(query_sim))
        if cfg.polarity == "similar" and neg_frac >= 0.3:
            print(
                f"warning: query={query_id} has many negative similarities (fraction={neg_frac:.2f}); opposite matches are being ignored",
                file=sys.stderr,
            )
            root_warnings.append(
                {
                    "code": "many_negative_similarities_ignored",
                    "query_id": query_id,
                    "fraction_negative_similarity": neg_frac,
                }
            )
        for polarity in polarities:
            program_warnings: list[dict[str, object]] = []
            candidate_neighbors = _candidate_neighbors(
                query_sim,
                polarity=polarity,
                min_similarity=float(cfg.min_similarity),
            )
            raw_candidate_neighbors = sorted(candidate_neighbors, key=lambda item: (-float(item[1]), str(item[0])))
            protected_direct_ref_ids = _protected_same_modality_reference_ids(
                raw_candidate_neighbors=raw_candidate_neighbors,
                reference_metadata=reference_metadata_effective,
                query_modality=query_modality,
                query_nominal_genes=query_nominal_genes,
                compound_map=compound_map,
                genetic_map=genetic_map,
                same_modality_top_n=5,
                raw_similarity_ratio=0.9,
            )
            protected_mechanism_ref_ids = _protected_same_modality_reference_ids(
                raw_candidate_neighbors=raw_candidate_neighbors,
                reference_metadata=reference_metadata_effective,
                query_modality=query_modality,
                query_nominal_genes=query_nominal_genes,
                compound_map=compound_map,
                genetic_map=genetic_map,
                same_modality_top_n=3,
                raw_similarity_ratio=0.95,
            )
            coherent_compound_prelabel_ref_ids = _coherent_same_modality_compound_prelabel_ids(
                raw_candidate_neighbors=raw_candidate_neighbors,
                reference_metadata=reference_metadata_effective,
                query_modality=query_modality,
                query_nominal_genes=query_nominal_genes,
                compound_map=compound_map,
                genetic_map=genetic_map,
                target_annotations=target_annotations,
                qc_weights=qc_weights,
                same_modality_top_n=6,
                raw_similarity_ratio=0.9,
            )
            same_modality_anchor_targets: set[str] = set()
            same_modality_anchor_pathway_seeds: set[str] = set()
            same_modality_anchor_target_classes: set[str] = set()
            same_modality_anchor_families: set[str] = set()
            same_modality_anchor_mechanisms: set[str] = set()
            if cfg.same_modality_first and query_modality:
                same_modality_candidates = [
                    (ref_id, base)
                    for ref_id, base in raw_candidate_neighbors
                    if str(reference_metadata_effective.get(ref_id, {}).get("perturbation_type", "")).strip().lower() == query_modality
                ]
                for ref_id, _base in same_modality_candidates[:5]:
                    mapping = compound_map.get(ref_id) or genetic_map.get(ref_id) or {}
                    ordered = sorted(mapping.items(), key=lambda item: (-float(item[1]), str(item[0])))[:3]
                    for gene, _weight in ordered:
                        same_modality_anchor_targets.add(gene)
                        if target_annotations:
                            pathway_seed = str((target_annotations.get(gene, {}) or {}).get("pathway_seed", "")).strip()
                            target_class = str((target_annotations.get(gene, {}) or {}).get("target_class", "")).strip()
                            family = str((target_annotations.get(gene, {}) or {}).get("target_family", "")).strip()
                            mechanism = str((target_annotations.get(gene, {}) or {}).get("mechanism_label", "")).strip()
                            if pathway_seed:
                                same_modality_anchor_pathway_seeds.add(pathway_seed)
                            if target_class:
                                same_modality_anchor_target_classes.add(target_class)
                            if family:
                                same_modality_anchor_families.add(family)
                            if mechanism:
                                same_modality_anchor_mechanisms.add(mechanism)
            penalized_neighbors: list[tuple[str, float, float]] = []
            for ref_id, base in candidate_neighbors:
                hub_penalty = _hubness_weight(
                    ref_id=ref_id,
                    reference_metadata=reference_metadata_effective,
                    mode=cfg.hubness_penalty,
                    hub_rank_penalties=hub_rank_penalties,
                )
                evidence = float(
                    (float(base) ** float(cfg.similarity_power))
                    * float(qc_weights.get(ref_id, 1.0))
                    * float(hub_penalty)
                )
                ref_modality = str(reference_metadata_effective.get(ref_id, {}).get("perturbation_type", "")).strip().lower()
                if cfg.same_modality_first and query_modality and ref_modality and ref_modality != query_modality:
                    mapping = compound_map.get(ref_id) or genetic_map.get(ref_id) or {}
                    overlap = bool(same_modality_anchor_targets & set(mapping))
                    pathway_overlap = False
                    class_overlap = False
                    family_overlap = False
                    mechanism_overlap = False
                    if target_annotations and mapping:
                        for gene in mapping:
                            pathway_seed = str((target_annotations.get(gene, {}) or {}).get("pathway_seed", "")).strip()
                            target_class = str((target_annotations.get(gene, {}) or {}).get("target_class", "")).strip()
                            family = str((target_annotations.get(gene, {}) or {}).get("target_family", "")).strip()
                            mechanism = str((target_annotations.get(gene, {}) or {}).get("mechanism_label", "")).strip()
                            if pathway_seed and pathway_seed in same_modality_anchor_pathway_seeds:
                                pathway_overlap = True
                            if target_class and target_class in same_modality_anchor_target_classes:
                                class_overlap = True
                            if family and family in same_modality_anchor_families:
                                family_overlap = True
                            if mechanism and mechanism in same_modality_anchor_mechanisms:
                                mechanism_overlap = True
                    cross_penalty = float(cfg.cross_modality_penalty)
                    if query_modality == "compound":
                        cross_penalty = max(cross_penalty, 0.8)
                    if overlap:
                        evidence *= 1.0
                    else:
                        recovery = cross_penalty
                        if pathway_overlap:
                            recovery = max(recovery, 0.95 if query_modality == "compound" else 0.8)
                        elif class_overlap:
                            recovery = max(recovery, 0.9 if query_modality == "compound" else 0.7)
                        elif mechanism_overlap:
                            recovery = max(recovery, 0.85 if query_modality == "compound" else 0.6)
                        elif family_overlap:
                            recovery = max(recovery, 0.8 if query_modality == "compound" else 0.5)
                        evidence *= recovery
                if evidence <= 0.0:
                    continue
                penalized_neighbors.append((ref_id, float(base), evidence))
            penalized_neighbors.sort(key=lambda item: (-float(item[2]), -float(item[1]), str(item[0])))
            all_penalized_neighbors = list(penalized_neighbors)
            if cfg.mutual_neighbor_filter and penalized_neighbors:
                anchors = [ref_id for ref_id, _base, _evidence in penalized_neighbors[: min(5, len(penalized_neighbors))]]
                filtered_neighbors: list[tuple[str, float, float]] = []
                for ref_id, base, evidence in penalized_neighbors:
                    if ref_id in anchors:
                        filtered_neighbors.append((ref_id, base, evidence))
                        continue
                    cohesion_vals = [
                        max(
                            0.0,
                            _pairwise_ref_similarity(ref_id, anchor_id, ref_vectors=ref_vectors, cache=ref_ref_similarity_cache),
                        )
                        for anchor_id in anchors
                        if anchor_id != ref_id
                    ]
                    cohesion = sum(cohesion_vals) / float(len(cohesion_vals) or 1)
                    adjusted = evidence * (0.5 if cohesion < 0.1 else 1.0)
                    if adjusted > 0.0:
                        filtered_neighbors.append((ref_id, base, adjusted))
                penalized_neighbors = sorted(filtered_neighbors, key=lambda item: (-float(item[2]), -float(item[1]), str(item[0])))
            if cfg.adaptive_neighbors and penalized_neighbors:
                adaptive_neighbors: list[tuple[str, float, float]] = []
                top_evidence = float(penalized_neighbors[0][2])
                for ref_id, base, evidence in penalized_neighbors:
                    if (
                        adaptive_neighbors
                        and len(adaptive_neighbors) >= int(cfg.min_effective_neighbors)
                        and float(evidence) < top_evidence * float(cfg.neighbor_evidence_drop_ratio)
                    ):
                        break
                    adaptive_neighbors.append((ref_id, base, evidence))
                penalized_neighbors = adaptive_neighbors
            penalized_by_ref = {ref_id: (ref_id, float(base), float(evidence)) for ref_id, base, evidence in all_penalized_neighbors}
            raw_by_ref = {ref_id: float(base) for ref_id, base in raw_candidate_neighbors}
            if cfg.mode == "direct_target":
                raw_pool_cap = max(int(cfg.max_reference_neighbors) * 15, int(cfg.min_effective_neighbors), 250) if int(cfg.max_reference_neighbors) > 0 else len(raw_candidate_neighbors)
                direct_pool_ids: list[str] = []
                top_raw = float(raw_candidate_neighbors[0][1]) if raw_candidate_neighbors else 0.0
                for ref_id, base in raw_candidate_neighbors[:raw_pool_cap]:
                    include = (
                        len(direct_pool_ids) < int(cfg.min_effective_neighbors)
                        or float(base) >= top_raw * float(cfg.neighbor_evidence_drop_ratio)
                        or ref_id in protected_direct_ref_ids
                    )
                    if include and ref_id in penalized_by_ref:
                        direct_pool_ids.append(ref_id)
                for ref_id in protected_direct_ref_ids:
                    if ref_id in penalized_by_ref and ref_id not in direct_pool_ids:
                        direct_pool_ids.append(ref_id)
                pooled_neighbors = sorted(
                    [penalized_by_ref[ref_id] for ref_id in direct_pool_ids if ref_id in penalized_by_ref],
                    key=lambda item: (-float(item[2]), -float(item[1]), str(item[0])),
                )
                raw_pooled_neighbors = [
                    (
                        ref_id,
                        float(raw_by_ref.get(ref_id, 0.0)),
                        float(raw_by_ref.get(ref_id, 0.0)) * float(qc_weights.get(ref_id, 1.0)),
                    )
                    for ref_id in direct_pool_ids
                    if ref_id in raw_by_ref
                ]
                prelabel_neighbors = list(pooled_neighbors)
                raw_prelabel_neighbors = list(raw_pooled_neighbors)
            else:
                prelabel_neighbors = list(all_penalized_neighbors)
                protected_prelabel_ref_ids = set(protected_mechanism_ref_ids) | set(coherent_compound_prelabel_ref_ids)
                if protected_prelabel_ref_ids:
                    for ref_id in protected_prelabel_ref_ids:
                        if ref_id in penalized_by_ref and all(existing_ref != ref_id for existing_ref, _b, _e in prelabel_neighbors):
                            prelabel_neighbors.append(penalized_by_ref[ref_id])
                    prelabel_neighbors.sort(key=lambda item: (-float(item[2]), -float(item[1]), str(item[0])))
                prelabel_cap = max(int(cfg.max_reference_neighbors) * 8, int(cfg.min_effective_neighbors), 50) if int(cfg.max_reference_neighbors) > 0 else len(prelabel_neighbors)
                prelabel_neighbors = prelabel_neighbors[:prelabel_cap]
                if protected_prelabel_ref_ids:
                    kept_prelabel_ids = {ref_id for ref_id, _base, _evidence in prelabel_neighbors}
                    for ref_id in sorted(protected_prelabel_ref_ids):
                        if ref_id in penalized_by_ref and ref_id not in kept_prelabel_ids:
                            prelabel_neighbors.append(penalized_by_ref[ref_id])
                    prelabel_neighbors = sorted(prelabel_neighbors, key=lambda item: (-float(item[2]), -float(item[1]), str(item[0])))
                prelabel_ref_ids = [ref_id for ref_id, _base, _evidence in prelabel_neighbors]
                raw_prelabel_neighbors = [
                    (
                        ref_id,
                        float(raw_by_ref.get(ref_id, 0.0)),
                        float(raw_by_ref.get(ref_id, 0.0)) * float(qc_weights.get(ref_id, 1.0)),
                    )
                    for ref_id in prelabel_ref_ids
                    if ref_id in raw_by_ref
                ]
                mechanism_neighbors = list(penalized_neighbors)
                protected_mechanism_pool_ids = set(protected_mechanism_ref_ids) | set(coherent_compound_prelabel_ref_ids)
                if protected_mechanism_pool_ids:
                    for ref_id in protected_mechanism_pool_ids:
                        if ref_id in penalized_by_ref and all(existing_ref != ref_id for existing_ref, _b, _e in mechanism_neighbors):
                            mechanism_neighbors.append(penalized_by_ref[ref_id])
                    mechanism_neighbors.sort(key=lambda item: (-float(item[2]), -float(item[1]), str(item[0])))
                pool_cap = max(int(cfg.max_reference_neighbors) * 5, int(cfg.min_effective_neighbors), 25) if int(cfg.max_reference_neighbors) > 0 else len(mechanism_neighbors)
                pooled_neighbors = mechanism_neighbors[:pool_cap]
                if protected_mechanism_pool_ids:
                    kept_pool_ids = {ref_id for ref_id, _base, _evidence in pooled_neighbors}
                    for ref_id in sorted(protected_mechanism_pool_ids):
                        if ref_id in penalized_by_ref and ref_id not in kept_pool_ids:
                            pooled_neighbors.append(penalized_by_ref[ref_id])
                    pooled_neighbors = sorted(pooled_neighbors, key=lambda item: (-float(item[2]), -float(item[1]), str(item[0])))
                pooled_ref_ids = [ref_id for ref_id, _base, _evidence in pooled_neighbors]
                raw_pooled_neighbors = [
                    (
                        ref_id,
                        float(raw_by_ref.get(ref_id, 0.0)),
                        float(raw_by_ref.get(ref_id, 0.0)) * float(qc_weights.get(ref_id, 1.0)),
                    )
                    for ref_id in pooled_ref_ids
                    if ref_id in raw_by_ref
                ]
            target_support, target_support_count, target_support_by_modality, target_best_similarity = _target_support_from_neighbors(
                penalized_neighbors=pooled_neighbors,
                compound_map=compound_map,
                genetic_map=genetic_map,
            )
            target_support_weighted = {
                gene: (
                    float(score)
                    * _recurrence_weight(
                        gene,
                        gene_idf=gene_idf,
                        mode=cfg.gene_recurrence_penalty,
                        stage="nomination",
                    )
                    if cfg.gene_recurrence_penalty != "none"
                    else float(score)
                )
                for gene, score in target_support.items()
            }
            top_target_candidates_raw = sorted(target_support.items(), key=lambda item: (-float(item[1]), str(item[0])))
            top_target_candidates_weighted = sorted(target_support_weighted.items(), key=lambda item: (-float(item[1]), str(item[0])))
            top_target_candidates = top_target_candidates_raw[:5]
            top_support = float(top_target_candidates_raw[0][1]) if top_target_candidates_raw else 0.0
            top_weighted_support = float(top_target_candidates_weighted[0][1]) if top_target_candidates_weighted else 0.0
            top_best_similarity = max((float(v) for v in target_best_similarity.values()), default=0.0)
            if query_nominal_genes and any(float(target_support.get(gene, 0.0)) > 0.0 for gene in query_nominal_genes):
                top_target_genes = {gene for gene in query_nominal_genes if float(target_support.get(gene, 0.0)) > 0.0}
            else:
                top_target_genes = {
                    gene
                    for gene, mass in top_target_candidates_raw[:5]
                    if top_support > 0.0
                    and (
                        float(mass) >= top_support * 0.8
                        or int(target_support_count.get(gene, 0)) >= 2
                        or float(target_best_similarity.get(gene, 0.0)) >= top_best_similarity * 0.95
                    )
                } or ({top_target_candidates_raw[0][0]} if top_target_candidates_raw else set())
                if top_weighted_support > 0.0:
                    top_target_genes.update(
                        gene
                        for gene, mass in top_target_candidates_weighted[:5]
                        if float(mass) >= top_weighted_support * 0.8
                    )
            core_branch = _build_core_branch(
                pooled_neighbors=pooled_neighbors,
                allowed_genes=top_target_genes,
                target_support_weighted=target_support_weighted,
                compound_map=compound_map,
                genetic_map=genetic_map,
                cfg=cfg,
                gene_idf=gene_idf,
            )
            mechanism_branch = _build_mechanism_branch(
                pooled_neighbors=pooled_neighbors,
                query_modality=query_modality,
                compound_map=compound_map,
                genetic_map=genetic_map,
                cfg=cfg,
                gene_idf=gene_idf,
            )
            active_branch = core_branch if cfg.mode == "direct_target" else mechanism_branch
            retained_neighbors = [(ref_id, base) for ref_id, base, _evidence in active_branch.neighbors]
            evidence_by_ref: dict[str, float] = dict(active_branch.evidence_by_ref)
            retained_by_modality: dict[str, int] = {}
            for ref_id, _score in retained_neighbors:
                modality = str(reference_metadata_effective.get(ref_id, {}).get("perturbation_type", "")).strip().lower() or "unknown"
                retained_by_modality[modality] = int(retained_by_modality.get(modality, 0)) + 1
            exact_target_support_present = bool(
                query_nominal_genes and any(float(target_support.get(gene, 0.0)) > 0.0 for gene in query_nominal_genes)
            )
            label_decision, raw_label_scores, prelabel_label_scores, retained_label_scores = _choose_expansion_label(
                mode=cfg.mode,
                query_modality=query_modality,
                raw_neighbors=raw_prelabel_neighbors,
                prelabel_neighbors=prelabel_neighbors,
                retained_neighbors=mechanism_branch.neighbors,
                compound_map=compound_map,
                genetic_map=genetic_map,
                reference_metadata=reference_metadata_effective,
                target_annotations=target_annotations,
                query_nominal_genes=query_nominal_genes,
                exact_target_support_present=exact_target_support_present,
            )
            raw_label_scores_by_modality = _label_scores_by_modality(
                neighbors=raw_prelabel_neighbors,
                compound_map=compound_map,
                genetic_map=genetic_map,
                reference_metadata=reference_metadata_effective,
                target_annotations=target_annotations,
                query_modality=query_modality,
            )
            prelabel_label_scores_by_modality = _label_scores_by_modality(
                neighbors=prelabel_neighbors,
                compound_map=compound_map,
                genetic_map=genetic_map,
                reference_metadata=reference_metadata_effective,
                target_annotations=target_annotations,
                query_modality=query_modality,
            )
            retained_label_scores_by_modality = _label_scores_by_modality(
                neighbors=mechanism_branch.neighbors,
                compound_map=compound_map,
                genetic_map=genetic_map,
                reference_metadata=reference_metadata_effective,
                target_annotations=target_annotations,
                query_modality=query_modality,
            )
            expansion_decision = dict(label_decision)
            expanded_gene_scores, expanded_family_genes, expansion_scope = _apply_local_weighted_expansion(
                base_gene_scores=mechanism_branch.gene_scores,
                branch_neighbors=mechanism_branch.neighbors,
                chosen_level=str(expansion_decision.get("chosen_level") or "") or None,
                chosen_label=str(expansion_decision.get("chosen_label") or "") or None,
                expansion_confidence=float(expansion_decision.get("expansion_confidence", 0.0) or 0.0),
                target_annotations=target_annotations,
                compound_map=compound_map,
                genetic_map=genetic_map,
                full_bundle_compound_map=full_bundle_compound_map,
                full_bundle_genetic_map=full_bundle_genetic_map,
                query_nominal_genes=query_nominal_genes,
            )
            expansion_decision.update(expansion_scope)
            top_pathway_seed = str((prelabel_label_scores.get("pathway_seed", {}) or {}).get("top_label") or "") or None
            top_pathway_seed_mass = float((prelabel_label_scores.get("pathway_seed", {}) or {}).get("top_label_mass", 0.0) or 0.0)
            top_target_class = str((prelabel_label_scores.get("target_class", {}) or {}).get("top_label") or "") or None
            top_target_class_mass = float((prelabel_label_scores.get("target_class", {}) or {}).get("top_label_mass", 0.0) or 0.0)
            top_family = str((prelabel_label_scores.get("target_family", {}) or {}).get("top_label") or "") or None
            top_family_mass = float((prelabel_label_scores.get("target_family", {}) or {}).get("top_label_mass", 0.0) or 0.0)
            top_mechanism = str((prelabel_label_scores.get("mechanism_label", {}) or {}).get("top_label") or "") or None
            top_mechanism_mass = float((prelabel_label_scores.get("mechanism_label", {}) or {}).get("top_label_mass", 0.0) or 0.0)
            chosen_level = str(expansion_decision.get("chosen_level") or "")
            chosen_label = str(expansion_decision.get("chosen_label") or "")
            chosen_confidence = float(expansion_decision.get("expansion_confidence", 0.0) or 0.0)
            if chosen_level == "pathway_seed" and chosen_label:
                top_pathway_seed = chosen_label
                top_pathway_seed_mass = chosen_confidence
            elif chosen_level == "target_class" and chosen_label:
                top_target_class = chosen_label
                top_target_class_mass = chosen_confidence
            elif chosen_level == "mechanism_label" and chosen_label:
                top_mechanism = chosen_label
                top_mechanism_mass = chosen_confidence
            elif chosen_level == "target_family" and chosen_label:
                top_family = chosen_label
                top_family_mass = chosen_confidence
            raw_family_summary = (raw_label_scores.get("target_family", {}) or {}).get("summary", {})
            prelabel_family_summary = (prelabel_label_scores.get("target_family", {}) or {}).get("summary", {})
            retained_family_summary = (retained_label_scores.get("target_family", {}) or {}).get("summary", {})
            raw_mechanism_summary = (raw_label_scores.get("mechanism_label", {}) or {}).get("summary", {})
            prelabel_mechanism_summary = (prelabel_label_scores.get("mechanism_label", {}) or {}).get("summary", {})
            retained_mechanism_summary = (retained_label_scores.get("mechanism_label", {}) or {}).get("summary", {})
            core_gene_scores = core_branch.gene_scores
            core_rows = _rows_from_scores(core_gene_scores, cfg)
            expanded_rows = _rows_from_scores(expanded_gene_scores, cfg)
            top_target_mass = float(top_target_candidates[0][1]) / float(sum(target_support.values()) or 1.0) if top_target_candidates else 0.0
            preferred_variant, preferred_variant_reason = _choose_preferred_variant(
                mode=cfg.mode,
                core_rows=core_rows,
                expanded_rows=expanded_rows,
                top_target_mass=top_target_mass,
                expansion_confidence=float(expansion_decision.get("expansion_confidence", 0.0) or 0.0),
            )
            hybrid_merge_stats = {"hybrid_merge_applied": False, "hybrid_core_gene_count": len(core_gene_scores), "hybrid_expanded_gene_count": len(expanded_gene_scores), "hybrid_added_only_gene_count": 0}
            if cfg.mode == "hybrid":
                gene_scores, hybrid_merge_stats = _merge_hybrid_scores(
                    core_gene_scores=core_gene_scores,
                    expanded_gene_scores=expanded_gene_scores,
                )
                hybrid_merge_stats["hybrid_merge_applied"] = True
                selected_rows = _rows_from_scores(gene_scores, cfg)
            else:
                selected_rows = expanded_rows if preferred_variant == "expanded" else core_rows
                gene_scores = expanded_gene_scores if preferred_variant == "expanded" else core_gene_scores
            full_gene_ids = sorted(gene_scores, key=lambda g: (-float(gene_scores[g]), str(g)))
            full_rows = [{"gene_id": gene_id, "gene_symbol": gene_id, "score": float(gene_scores[gene_id]), "rank": rank} for rank, gene_id in enumerate(full_gene_ids, start=1)]
            frac_artifact = _artifact_family_fraction(selected_rows)
            if frac_artifact >= 0.4:
                print(
                    f"warning: query={query_id} polarity={polarity} appears dominated by common receptor-family genes; interpret cautiously",
                    file=sys.stderr,
                )
                warning_payload = {
                    "code": "artifact_family_dominance",
                    "query_id": query_id,
                    "polarity": polarity,
                    "fraction": frac_artifact,
                }
                program_warnings.append(warning_payload)
                root_warnings.append(warning_payload)

            top_neighbor_ids = [ref_id for ref_id, _score in retained_neighbors[:5]]
            top_neighbor_sims = [float(score) for _ref_id, score in retained_neighbors[:5]]
            top_neighbor_modalities = [
                str(reference_metadata_effective.get(ref_id, {}).get("perturbation_type", "")).strip().lower()
                for ref_id in top_neighbor_ids
            ]
            retained_hub_scores: list[float] = []
            for ref_id, _score in retained_neighbors:
                raw_hub = reference_metadata_effective.get(ref_id, {}).get("hub_score")
                try:
                    parsed_hub = float(raw_hub) if raw_hub not in {None, ""} else None
                except ValueError:
                    parsed_hub = None
                if parsed_hub is not None:
                    retained_hub_scores.append(parsed_hub)
            total_evidence = sum(float(v) for v in evidence_by_ref.values())
            top5_mass = 0.0
            if total_evidence > 0.0:
                top5_mass = sum(float(evidence_by_ref.get(ref_id, 0.0)) for ref_id in top_neighbor_ids) / total_evidence
            n_positive_neighbors = sum(1 for value in query_sim.values() if float(value) > 0.0)
            n_negative_neighbors = sum(1 for value in query_sim.values() if float(value) < 0.0)
            retrieval_confidence = _retrieval_confidence(n_positive_neighbors, neg_frac)
            neighbor_primary_target_agreement = _routed_primary_target_agreement(
                retained_neighbors=retained_neighbors,
                evidence_by_ref=evidence_by_ref,
                compound_map=compound_map,
                genetic_map=genetic_map,
            )
            neighbor_top3_target_agreement = _routed_topk_target_agreement(
                retained_neighbors=retained_neighbors,
                evidence_by_ref=evidence_by_ref,
                compound_map=compound_map,
                genetic_map=genetic_map,
                top_k=3,
            )
            high_hub_mass_fraction = _high_hub_mass_fraction(
                evidence_by_ref=evidence_by_ref,
                reference_metadata=reference_metadata_effective,
            )
            specificity_confidence, top10_gene_mass, neighbor_target_concentration, neighbor_primary_target_agreement, neighbor_top3_target_agreement, high_hub_mass_fraction = _specificity_confidence(
                selected_rows,
                gene_scores,
                neighbor_primary_target_agreement=neighbor_primary_target_agreement,
                neighbor_top3_target_agreement=neighbor_top3_target_agreement,
                high_hub_mass_fraction=high_hub_mass_fraction,
            )
            total_target_support = float(sum(target_support.values()))
            top_target_support_mass = float(top_target_candidates[0][1]) / total_target_support if top_target_candidates and total_target_support > 0.0 else 0.0
            query_target_support_mass = 0.0
            query_target_best_rank = None
            if query_nominal_genes and target_support:
                ordered_targets = sorted(target_support.items(), key=lambda item: (-float(item[1]), str(item[0])))
                rank_map = {gene: rank for rank, (gene, _score) in enumerate(ordered_targets, start=1)}
                query_target_support_mass = max(float(target_support.get(gene, 0.0)) for gene in query_nominal_genes) / total_target_support if total_target_support > 0.0 else 0.0
                query_target_best_rank = min((rank_map.get(gene) for gene in query_nominal_genes if gene in rank_map), default=None)
            if retrieval_confidence == "high" and specificity_confidence == "low":
                warning_payload = {
                    "code": "diffuse_gene_evidence",
                    "query_id": query_id,
                    "polarity": polarity,
                    "suggestion": "neighbors are geometrically coherent but gene evidence is diffuse; tighten --max_reference_neighbors, inspect hubness/context, and check target routing",
                }
                print(
                    f"warning: query={query_id} polarity={polarity} neighbors are geometrically coherent but gene evidence is diffuse",
                    file=sys.stderr,
                )
                program_warnings.append(warning_payload)
                root_warnings.append(warning_payload)
            if retrieval_confidence == "low":
                warning_payload = {
                    "code": "low_retrieval_confidence",
                    "query_id": query_id,
                    "polarity": polarity,
                    "suggestion": "try --polarity both, verify cell line/timepoint coherence, and use a same-timepoint bundle",
                }
                print(
                    f"warning: query={query_id} polarity={polarity} has low retrieval confidence; check context matching and polarity",
                    file=sys.stderr,
                )
                program_warnings.append(warning_payload)
                root_warnings.append(warning_payload)
            if polarity == "opposite" and _confidence_rank(specificity_confidence) < _confidence_rank(cfg.min_specificity_confidence_to_emit_opposite):
                warning_payload = {
                    "code": "opposite_program_suppressed_low_specificity",
                    "query_id": query_id,
                    "polarity": polarity,
                    "specificity_confidence": specificity_confidence,
                    "min_required": cfg.min_specificity_confidence_to_emit_opposite,
                    "suggestion": "rerun with --min_specificity_confidence_to_emit_opposite low if you explicitly want exploratory opposite-polarity outputs",
                }
                print(
                    f"warning: query={query_id} polarity={polarity} suppressed because specificity_confidence={specificity_confidence} < {cfg.min_specificity_confidence_to_emit_opposite}",
                    file=sys.stderr,
                )
                program_warnings.append(warning_payload)
                root_warnings.append(warning_payload)
                skipped_programs.append(
                    {
                        "query_id": query_id,
                        "polarity": polarity,
                        "reason": "low_specificity_opposite_suppressed",
                        "specificity_confidence": specificity_confidence,
                        "min_required": cfg.min_specificity_confidence_to_emit_opposite,
                    }
                )
                continue

            safe_program_id = _safe_component(f"{query_id}__polarity={polarity}", "program")
            program_dir = out_dir / f"program={safe_program_id}"
            program_dir.mkdir(parents=True, exist_ok=True)
            _write_rows(program_dir / "geneset.tsv", selected_rows)
            if cfg.mode == "hybrid":
                _write_rows(program_dir / "geneset.core.tsv", core_rows)
                _write_rows(program_dir / "geneset.expanded.tsv", expanded_rows)
            if cfg.emit_full:
                _write_rows(program_dir / "geneset.full.tsv", full_rows)

            gmt_sets: list[tuple[str, list[str]]] = []
            gmt_plans: list[dict[str, object]] = []
            gmt_diagnostics: list[dict[str, object]] = []
            gmt_path = resolve_gmt_out_path(program_dir, cfg.gmt_out)
            if cfg.emit_gmt:
                base_name = (
                    f"{cfg.converter_name}__dataset={_safe_component(cfg.dataset_label, 'dataset')}"
                    f"__query={_safe_component(query_id, 'query')}__polarity={polarity}"
                    f"__similarity_metric={cfg.similarity_metric}"
                )
                gmt_sets, gmt_plans = build_gmt_sets_from_rows(
                    rows=full_rows,
                    base_name=base_name,
                    prefer_symbol=bool(cfg.gmt_prefer_symbol),
                    min_genes=int(cfg.gmt_min_genes),
                    max_genes=int(cfg.gmt_max_genes),
                    topk_list=gmt_topk_list,
                    mass_list=gmt_mass_list,
                    split_signed=False,
                    require_symbol=bool(cfg.gmt_require_symbol),
                    allowed_biotypes={x.lower() for x in gmt_biotype_allowlist} if gmt_biotype_allowlist else None,
                    emit_small_gene_sets=bool(cfg.emit_small_gene_sets),
                    diagnostics=gmt_diagnostics,
                    context={"converter": cfg.converter_name, "query_id": query_id, "polarity": polarity},
                )
                if gmt_sets:
                    write_gmt(gmt_sets, gmt_path, gmt_format=cfg.gmt_format)
                    combined_gmt_sets.extend(gmt_sets)
            gmt_warnings = _warn_gmt_diagnostics(gmt_diagnostics)
            program_warnings.extend(gmt_warnings)
            root_warnings.extend(gmt_warnings)

            output_files = [{"path": str(program_dir / "geneset.tsv"), "role": "selected_program"}]
            if cfg.mode == "hybrid":
                output_files.append({"path": str(program_dir / "geneset.core.tsv"), "role": "core_program"})
                output_files.append({"path": str(program_dir / "geneset.expanded.tsv"), "role": "expanded_program"})
            if cfg.emit_full:
                output_files.append({"path": str(program_dir / "geneset.full.tsv"), "role": "full_scores"})
            if cfg.emit_gmt and gmt_sets:
                output_files.append({"path": str(gmt_path), "role": "gmt"})

            run_summary_payload = {
                "converter": cfg.converter_name,
                "dataset_label": cfg.dataset_label,
                "program_id": query_id,
                "selected_direction": polarity,
                "n_query_profiles": len(query_membership.get(query_id, [])),
                "n_reference_profiles": len(ref_vectors),
                "n_positive_neighbors": n_positive_neighbors,
                "n_negative_neighbors": n_negative_neighbors,
                "n_retained_neighbors": len(retained_neighbors),
                "effective_neighbor_count": len(retained_neighbors),
                "mode": cfg.mode,
                "preferred_variant": preferred_variant,
                "query_modality": query_modality,
                "retained_neighbors_by_modality": retained_by_modality,
                "neg_similarity_fraction": neg_frac,
                "top_neighbor_ids": top_neighbor_ids,
                "top_neighbor_modalities": top_neighbor_modalities,
                "top_neighbor_similarities": top_neighbor_sims,
                "raw_candidate_neighbor_ids": [ref_id for ref_id, _score in raw_candidate_neighbors[:10]],
                "raw_candidate_neighbor_similarities": [float(score) for _ref_id, score in raw_candidate_neighbors[:10]],
                "raw_candidate_neighbors_detail": [
                    {
                        "ref_id": ref_id,
                        "similarity": float(score),
                        "modality": str(reference_metadata_effective.get(ref_id, {}).get("perturbation_type", "")).strip().lower(),
                        "top_targets": sorted(
                            (compound_map.get(ref_id) or genetic_map.get(ref_id) or {}).items(),
                            key=lambda item: (-float(item[1]), str(item[0])),
                        )[:3],
                    }
                    for ref_id, score in raw_candidate_neighbors[:10]
                ],
                "retained_neighbors_detail": [
                    {
                        "ref_id": ref_id,
                        "similarity": float(score),
                        "evidence": float(evidence_by_ref.get(ref_id, 0.0)),
                        "modality": str(reference_metadata_effective.get(ref_id, {}).get("perturbation_type", "")).strip().lower(),
                        "target_family": sorted(
                            {
                                str((target_annotations.get(gene, {}) or {}).get("target_family", "")).strip()
                                for gene in (compound_map.get(ref_id) or genetic_map.get(ref_id) or {})
                                if target_annotations and str((target_annotations.get(gene, {}) or {}).get("target_family", "")).strip()
                            }
                        ),
                        "mechanism_label": sorted(
                            {
                                str((target_annotations.get(gene, {}) or {}).get("mechanism_label", "")).strip()
                                for gene in (compound_map.get(ref_id) or genetic_map.get(ref_id) or {})
                                if target_annotations and str((target_annotations.get(gene, {}) or {}).get("mechanism_label", "")).strip()
                            }
                        ),
                    }
                    for ref_id, score in retained_neighbors[:10]
                ],
                "core_branch_neighbor_ids": [ref_id for ref_id, _base, _evidence in core_branch.neighbors[:10]],
                "mechanism_branch_neighbor_ids": [ref_id for ref_id, _base, _evidence in mechanism_branch.neighbors[:10]],
                "core_branch_neighbor_count": len(core_branch.neighbors),
                "mechanism_branch_neighbor_count": len(mechanism_branch.neighbors),
                "top5_neighbor_score_mass_fraction": top5_mass,
                "retrieval_confidence": retrieval_confidence,
                "specificity_confidence": specificity_confidence,
                "top10_gene_mass": top10_gene_mass,
                "neighbor_target_concentration": neighbor_target_concentration,
                "neighbor_primary_target_agreement": neighbor_primary_target_agreement,
                "neighbor_top3_target_agreement": neighbor_top3_target_agreement,
                "high_hub_mass_fraction": high_hub_mass_fraction,
                "top_target_candidates": [{"gene_symbol": gene, "support_mass": score, "weighted_support_mass": float(target_support_weighted.get(gene, score)), "support_count": int(target_support_count.get(gene, 0)), "best_similarity": float(target_best_similarity.get(gene, 0.0)), "support_by_modality": target_support_by_modality.get(gene, {})} for gene, score in top_target_candidates],
                "top_target_support_mass": top_target_support_mass,
                "query_nominal_targets": sorted(query_nominal_genes),
                "query_target_best_rank": query_target_best_rank,
                "query_target_support_mass": query_target_support_mass,
                "top_pathway_seed": top_pathway_seed,
                "top_pathway_seed_mass": top_pathway_seed_mass,
                "top_target_class": top_target_class,
                "top_target_class_mass": top_target_class_mass,
                "top_family": top_family,
                "top_family_mass": top_family_mass,
                "top_mechanism": top_mechanism,
                "top_mechanism_mass": top_mechanism_mass,
                "family_vote_summary_raw": raw_family_summary,
                "family_vote_summary_prelabel": prelabel_family_summary,
                "family_vote_summary_retained": retained_family_summary,
                "mechanism_vote_summary_raw": raw_mechanism_summary,
                "mechanism_vote_summary_prelabel": prelabel_mechanism_summary,
                "mechanism_vote_summary_retained": retained_mechanism_summary,
                "expansion_decision": expansion_decision,
                "label_scores_raw": raw_label_scores,
                "label_scores_prelabel": prelabel_label_scores,
                "label_scores_retained": retained_label_scores,
                "label_scores_raw_by_modality": raw_label_scores_by_modality,
                "label_scores_prelabel_by_modality": prelabel_label_scores_by_modality,
                "label_scores_retained_by_modality": retained_label_scores_by_modality,
                "expanded_added_genes": expanded_family_genes,
                "core_gene_count": len(core_rows),
                "expanded_gene_count": len(expanded_rows),
                "preferred_variant_reason": preferred_variant_reason,
                **hybrid_merge_stats,
                "compound_target_weighting_mode": compound_ref_weighting_stats["compound_target_weighting_mode"],
                "n_compound_refs_renormalized": compound_ref_weighting_stats["n_compound_refs_renormalized"],
                "median_raw_target_count": compound_ref_weighting_stats["median_raw_target_count"],
                "hubness_penalty": cfg.hubness_penalty,
                "hub_score_present": bool(hub_scores_present > 0),
                "hub_score_summary_retained": {
                    "min": min(retained_hub_scores) if retained_hub_scores else None,
                    "median": median(retained_hub_scores) if retained_hub_scores else None,
                    "max": max(retained_hub_scores) if retained_hub_scores else None,
                },
                "bundle_contexts": parse_summary.get("bundle_manifest", {}).get("contexts") if isinstance(parse_summary.get("bundle_manifest"), dict) else None,
                "bundle_summary": parse_summary.get("bundle_manifest", {}).get("summary") if isinstance(parse_summary.get("bundle_manifest"), dict) else None,
                "parse_summary": parse_summary,
                "feature_alignment": feature_summary,
                "query_standardization": query_scale_summary,
                "reference_standardization": ref_scale_summary,
                "control_calibration": control_calibration_summary,
                "mapping_summary": mapping_summary,
                "requested_gmt_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
                "skipped_gmt_outputs": _collect_skipped_gmt_outputs(gmt_diagnostics),
                "warnings": program_warnings,
            }
            run_json_path, run_txt_path = write_run_summary_files(program_dir, run_summary_payload)
            output_files.append({"path": str(run_json_path), "role": "run_summary_json"})
            output_files.append({"path": str(run_txt_path), "role": "run_summary_text"})

            meta = make_metadata(
                converter_name=cfg.converter_name,
                parameters={
                    "similarity_metric": cfg.similarity_metric,
                    "similarity_power": cfg.similarity_power,
                    "polarity": polarity,
                    "max_reference_neighbors": cfg.max_reference_neighbors,
                    "min_similarity": cfg.min_similarity,
                    "control_calibration": cfg.control_calibration,
                    "mode": cfg.mode,
                    "control_residual_components": cfg.control_residual_components,
                    "control_min_profiles_for_residualization": cfg.control_min_profiles_for_residualization,
                    "hubness_penalty": cfg.hubness_penalty,
                    "same_modality_first": cfg.same_modality_first,
                    "cross_modality_penalty": cfg.cross_modality_penalty,
                    "adaptive_neighbors": cfg.adaptive_neighbors,
                    "mutual_neighbor_filter": cfg.mutual_neighbor_filter,
                    "gene_recurrence_penalty": cfg.gene_recurrence_penalty,
                    "min_specificity_confidence_to_emit_opposite": cfg.min_specificity_confidence_to_emit_opposite,
                    "compound_weight": cfg.compound_weight,
                    "genetic_weight": cfg.genetic_weight,
                    "select": cfg.select,
                    "top_k": cfg.top_k,
                    "quantile": cfg.quantile,
                    "min_score": cfg.min_score,
                    "normalize": cfg.normalize,
                    "emit_gmt": cfg.emit_gmt,
                    "gmt_format": cfg.gmt_format,
                },
                data_type="morphology_profile",
                assay="morphology",
                organism=cfg.organism,
                genome_build=cfg.genome_build,
                files=input_files,
                gene_annotation={"mode": "none", "source": "reference_mapping", "gene_id_field": "gene_symbol"},
                weights={"weight_type": "nonnegative", "normalization": {"method": cfg.normalize, "target_sum": 1.0 if cfg.normalize in {"within_set_l1", "l1"} else None}, "aggregation": "reference_similarity_mapping"},
                summary={
                    "n_input_features": len(shared_features),
                    "n_genes": len(selected_rows),
                    "n_features_assigned": len(shared_features),
                    "fraction_features_assigned": 1.0,
                    "feature_alignment": feature_summary,
                    "parse_summary": parse_summary,
                    "mapping_summary": mapping_summary,
                    "resources": resources_info,
                    "warnings": program_warnings,
                    "retrieval_confidence": retrieval_confidence,
                    "specificity_confidence": specificity_confidence,
                    "mode": cfg.mode,
                    "preferred_variant": preferred_variant,
                    "n_positive_neighbors": n_positive_neighbors,
                    "n_negative_neighbors": n_negative_neighbors,
                    "effective_neighbor_count": len(retained_neighbors),
                    "query_modality": query_modality,
                    "neg_similarity_fraction": neg_frac,
                    "top_neighbor_ids": top_neighbor_ids,
                    "top_neighbor_modalities": top_neighbor_modalities,
                    "top_neighbor_similarities": top_neighbor_sims,
                    "raw_candidate_neighbor_ids": [ref_id for ref_id, _score in raw_candidate_neighbors[:10]],
                    "raw_candidate_neighbor_similarities": [float(score) for _ref_id, score in raw_candidate_neighbors[:10]],
                    "core_branch_neighbor_ids": [ref_id for ref_id, _base, _evidence in core_branch.neighbors[:10]],
                    "mechanism_branch_neighbor_ids": [ref_id for ref_id, _base, _evidence in mechanism_branch.neighbors[:10]],
                    "core_branch_neighbor_count": len(core_branch.neighbors),
                    "mechanism_branch_neighbor_count": len(mechanism_branch.neighbors),
                    "family_vote_summary_raw": raw_family_summary,
                    "family_vote_summary_prelabel": prelabel_family_summary,
                    "family_vote_summary_retained": retained_family_summary,
                    "mechanism_vote_summary_raw": raw_mechanism_summary,
                    "mechanism_vote_summary_prelabel": prelabel_mechanism_summary,
                    "mechanism_vote_summary_retained": retained_mechanism_summary,
                    "label_scores_raw": raw_label_scores,
                    "label_scores_prelabel": prelabel_label_scores,
                    "label_scores_retained": retained_label_scores,
                    "label_scores_raw_by_modality": raw_label_scores_by_modality,
                    "label_scores_prelabel_by_modality": prelabel_label_scores_by_modality,
                    "label_scores_retained_by_modality": retained_label_scores_by_modality,
                    "expansion_decision": expansion_decision,
                    "top10_gene_mass": top10_gene_mass,
                    "neighbor_target_concentration": neighbor_target_concentration,
                    "neighbor_primary_target_agreement": neighbor_primary_target_agreement,
                    "neighbor_top3_target_agreement": neighbor_top3_target_agreement,
                    "high_hub_mass_fraction": high_hub_mass_fraction,
                    "top_target_candidates": [{"gene_symbol": gene, "support_mass": score, "weighted_support_mass": float(target_support_weighted.get(gene, score)), "support_count": int(target_support_count.get(gene, 0)), "best_similarity": float(target_best_similarity.get(gene, 0.0)), "support_by_modality": target_support_by_modality.get(gene, {})} for gene, score in top_target_candidates],
                    "top_target_support_mass": top_target_support_mass,
                    "query_nominal_targets": sorted(query_nominal_genes),
                    "query_target_best_rank": query_target_best_rank,
                    "query_target_support_mass": query_target_support_mass,
                    "top_pathway_seed": top_pathway_seed,
                    "top_pathway_seed_mass": top_pathway_seed_mass,
                    "top_target_class": top_target_class,
                    "top_target_class_mass": top_target_class_mass,
                    "top_family": top_family,
                    "top_family_mass": top_family_mass,
                    "top_mechanism": top_mechanism,
                    "top_mechanism_mass": top_mechanism_mass,
                    "core_gene_count": len(core_rows),
                    "expanded_gene_count": len(expanded_rows),
                    "expanded_added_genes": expanded_family_genes,
                    "preferred_variant_reason": preferred_variant_reason,
                    **hybrid_merge_stats,
                    "compound_target_weighting_mode": compound_ref_weighting_stats["compound_target_weighting_mode"],
                    "n_compound_refs_renormalized": compound_ref_weighting_stats["n_compound_refs_renormalized"],
                    "median_raw_target_count": compound_ref_weighting_stats["median_raw_target_count"],
                    "hubness_penalty": cfg.hubness_penalty,
                    "same_modality_first": cfg.same_modality_first,
                    "cross_modality_penalty": cfg.cross_modality_penalty,
                    "adaptive_neighbors": cfg.adaptive_neighbors,
                    "mutual_neighbor_filter": cfg.mutual_neighbor_filter,
                    "gene_recurrence_penalty": cfg.gene_recurrence_penalty,
                    "hub_score_present": bool(hub_scores_present > 0),
                    "control_calibration": control_calibration_summary,
                    "hub_score_summary_retained": {
                        "min": min(retained_hub_scores) if retained_hub_scores else None,
                        "median": median(retained_hub_scores) if retained_hub_scores else None,
                        "max": max(retained_hub_scores) if retained_hub_scores else None,
                    },
                    "bundle_contexts": parse_summary.get("bundle_manifest", {}).get("contexts") if isinstance(parse_summary.get("bundle_manifest"), dict) else None,
                    "bundle_summary": parse_summary.get("bundle_manifest", {}).get("summary") if isinstance(parse_summary.get("bundle_manifest"), dict) else None,
                },
                program_extraction={
                    "selection_method": cfg.select,
                    "selection_params": {"k": cfg.top_k, "quantile": cfg.quantile, "min_score": cfg.min_score},
                    "normalize": cfg.normalize,
                    "n_selected_genes": len(selected_rows),
                    "preferred_variant": preferred_variant,
                    "core_n_selected_genes": len(core_rows),
                    "expanded_n_selected_genes": len(expanded_rows),
                    "score_definition": f"morphology reference similarity polarity={polarity}",
                    "retrieval_confidence": retrieval_confidence,
                    "specificity_confidence": specificity_confidence,
                    "neighbor_primary_target_agreement": neighbor_primary_target_agreement,
                    "neighbor_top3_target_agreement": neighbor_top3_target_agreement,
                    "high_hub_mass_fraction": high_hub_mass_fraction,
                    "effective_neighbor_count": len(retained_neighbors),
                },
                output_files=output_files,
                gmt={
                    "written": bool(cfg.emit_gmt and gmt_sets),
                    "path": str(gmt_path.relative_to(program_dir)) if cfg.emit_gmt and gmt_sets and gmt_path.is_relative_to(program_dir) else (str(gmt_path) if cfg.emit_gmt and gmt_sets else None),
                    "prefer_symbol": bool(cfg.gmt_prefer_symbol),
                    "require_symbol": bool(cfg.gmt_require_symbol),
                    "biotype_allowlist": gmt_biotype_allowlist,
                    "min_genes": int(cfg.gmt_min_genes),
                    "max_genes": int(cfg.gmt_max_genes),
                    "emit_small_gene_sets": bool(cfg.emit_small_gene_sets),
                    "requested_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
                    "emitted_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
                    "skipped_outputs": _collect_skipped_gmt_outputs(gmt_diagnostics),
                    "diagnostics": gmt_diagnostics,
                    "plans": gmt_plans,
                    "format": cfg.gmt_format,
                },
            )
            write_metadata(program_dir / "geneset.meta.json", meta)
            manifest_rows.append({
                "program_id": query_id,
                "query_id": query_id,
                "polarity": polarity,
                "mode": cfg.mode,
                "preferred_variant": preferred_variant,
                "n_query_profiles": len(query_membership.get(query_id, [])),
                "n_reference_profiles": len(ref_vectors),
                "n_retained_neighbors": len(retained_neighbors),
                "retrieval_confidence": retrieval_confidence,
                "specificity_confidence": specificity_confidence,
                "neighbor_primary_target_agreement": neighbor_primary_target_agreement,
                "neighbor_top3_target_agreement": neighbor_top3_target_agreement,
                "high_hub_mass_fraction": high_hub_mass_fraction,
                "top_target_support_mass": top_target_support_mass,
                "n_genes_selected": len(selected_rows),
                "path": str(program_dir.relative_to(out_dir)),
            })

    if not manifest_rows:
        raise ValueError("No morphology programs were emitted. Check feature overlap, mappings, and input tables.")
    with (out_dir / "manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "program_id",
                "query_id",
                "polarity",
                "mode",
                "preferred_variant",
                "n_query_profiles",
                "n_reference_profiles",
                "n_retained_neighbors",
                "retrieval_confidence",
                "specificity_confidence",
                "neighbor_primary_target_agreement",
                "neighbor_top3_target_agreement",
                "high_hub_mass_fraction",
                "top_target_support_mass",
                "n_genes_selected",
                "path",
            ],
        )
        writer.writeheader()
        for row in manifest_rows:
            writer.writerow(row)
    if cfg.emit_gmt and combined_gmt_sets:
        write_gmt(combined_gmt_sets, out_dir / "genesets.gmt", gmt_format=cfg.gmt_format)
    root_summary = {
        "converter": cfg.converter_name,
        "dataset_label": cfg.dataset_label,
        "mode": cfg.mode,
        "n_groups": len(manifest_rows),
        "parse_summary": parse_summary,
        "feature_alignment": feature_summary,
        "mapping_summary": mapping_summary,
        "bundle_contexts": parse_summary.get("bundle_manifest", {}).get("contexts") if isinstance(parse_summary.get("bundle_manifest"), dict) else None,
        "bundle_summary": parse_summary.get("bundle_manifest", {}).get("summary") if isinstance(parse_summary.get("bundle_manifest"), dict) else None,
        "resources": resources_info,
        "warnings": root_warnings,
        "skipped_programs": skipped_programs,
        "control_calibration": control_calibration_summary,
    }
    run_json_path, run_txt_path = write_run_summary_files(out_dir, root_summary)
    return {"n_groups": len(manifest_rows), "out_dir": str(out_dir), "run_summary_json": str(run_json_path), "run_summary_txt": str(run_txt_path)}
