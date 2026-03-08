from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re
import sys
from statistics import median

from geneset_extractors.core.gmt import build_gmt_sets_from_rows, parse_int_list_csv, parse_mass_list_csv, parse_str_list_csv, resolve_gmt_out_path, write_gmt
from geneset_extractors.core.metadata import make_metadata, write_metadata
from geneset_extractors.core.qc import write_run_summary_files
from geneset_extractors.core.selection import global_l1_weights, ranked_gene_ids, select_quantile, select_threshold, select_top_k, within_set_l1_weights
from geneset_extractors.extractors.morphology.mapping import accumulate_gene_scores, build_reference_gene_maps, l1_normalize_scores
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
    similarity_metric: str
    similarity_power: float
    polarity: str
    max_reference_neighbors: int
    min_similarity: float
    hubness_penalty: str
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


def _retained_neighbors(
    query_sim: dict[str, float],
    *,
    polarity: str,
    max_reference_neighbors: int,
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
    retained.sort(key=lambda item: (-float(item[1]), str(item[0])))
    if int(max_reference_neighbors) > 0:
        retained = retained[: int(max_reference_neighbors)]
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
) -> tuple[str, float, float]:
    if not selected_rows or not gene_scores:
        return "low", 0.0, 0.0
    top10_gene_mass = sum(float(row.get("weight", 0.0)) for row in selected_rows[:10])
    total_gene_mass = sum(float(v) for v in gene_scores.values())
    neighbor_target_concentration = 0.0
    if total_gene_mass > 0.0:
        neighbor_target_concentration = max(float(v) for v in gene_scores.values()) / total_gene_mass
    if top10_gene_mass >= 0.7 and neighbor_target_concentration >= 0.35:
        return "high", top10_gene_mass, neighbor_target_concentration
    if top10_gene_mass >= 0.7 or neighbor_target_concentration >= 0.35:
        return "medium", top10_gene_mass, neighbor_target_concentration
    return "low", top10_gene_mass, neighbor_target_concentration


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


def run_morphology_workflow(
    *,
    cfg: MorphologyWorkflowConfig,
    query_profiles: dict[str, dict[str, float]],
    reference_profiles: dict[str, dict[str, float]],
    reference_metadata: dict[str, dict[str, str]],
    compound_targets: dict[str, dict[str, float]],
    feature_schema: list[str] | None,
    feature_stats: dict[str, tuple[float, float]] | None,
    query_membership: dict[str, list[str]],
    parse_summary: dict[str, object],
    input_files: list[dict[str, str]],
    resources_info: dict[str, object] | None,
    exclude_query_ids_from_reference: bool,
) -> dict[str, object]:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    reference_profiles_effective = dict(reference_profiles)
    reference_metadata_effective = dict(reference_metadata)
    if exclude_query_ids_from_reference:
        overlap = set(query_profiles) & set(reference_profiles_effective)
        for ref_id in overlap:
            reference_profiles_effective.pop(ref_id, None)
            reference_metadata_effective.pop(ref_id, None)
    control_ids = {rid for rid, row in reference_metadata_effective.items() if str(row.get("is_control", "")).strip().lower() in {"1", "true", "t", "yes"}}
    if control_ids:
        print(f"warning: excluding {len(control_ids)} reference control profiles from similarity candidates", file=sys.stderr)
        for rid in control_ids:
            reference_profiles_effective.pop(rid, None)
            reference_metadata_effective.pop(rid, None)

    shared_features, feature_summary = align_feature_space(query_profiles, reference_profiles_effective, feature_schema=feature_schema)
    if not shared_features:
        raise ValueError("No shared numeric features between query and reference profiles")
    if float(feature_summary.get("fraction_reference_features_shared", 1.0)) < 0.8:
        print(
            f"warning: shared feature fraction is low ({feature_summary['fraction_reference_features_shared']:.3f}); results may be unstable",
            file=sys.stderr,
        )

    query_vectors, query_scale_summary = standardize_profiles(query_profiles, features=shared_features, feature_stats=feature_stats)
    ref_vectors, ref_scale_summary = standardize_profiles(reference_profiles_effective, features=shared_features, feature_stats=feature_stats)
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

    compound_map, genetic_map, mapping_summary = build_reference_gene_maps(
        reference_metadata=reference_metadata_effective,
        compound_targets=compound_targets,
    )
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
    gmt_topk_list = parse_int_list_csv(str(cfg.gmt_topk_list))
    gmt_mass_list = parse_mass_list_csv(str(cfg.gmt_mass_list))
    gmt_biotype_allowlist = parse_str_list_csv(str(cfg.gmt_biotype_allowlist))
    polarities = [cfg.polarity] if cfg.polarity in {"similar", "opposite"} else ["similar", "opposite"]
    if cfg.hubness_penalty != "none" and hub_scores_present == 0:
        print("warning: hub_score not present in reference metadata; proceeding without hubness penalty", file=sys.stderr)
        root_warnings.append({"code": "hub_score_missing", "mode": cfg.hubness_penalty})

    for query_id in sorted(query_vectors):
        query_sim = similarities.get(query_id, {})
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
            retained_neighbors = _retained_neighbors(
                query_sim,
                polarity=polarity,
                max_reference_neighbors=int(cfg.max_reference_neighbors),
                min_similarity=float(cfg.min_similarity),
            )
            evidence_by_ref: dict[str, float] = {}
            for ref_id, base in retained_neighbors:
                hub_penalty = _hubness_weight(
                    ref_id=ref_id,
                    reference_metadata=reference_metadata_effective,
                    mode=cfg.hubness_penalty,
                    hub_rank_penalties=hub_rank_penalties,
                )
                evidence_by_ref[ref_id] = float(
                    (float(base) ** float(cfg.similarity_power))
                    * float(qc_weights.get(ref_id, 1.0))
                    * float(hub_penalty)
                )
            retained_by_modality: dict[str, int] = {}
            for ref_id, _score in retained_neighbors:
                modality = str(reference_metadata_effective.get(ref_id, {}).get("perturbation_type", "")).strip().lower() or "unknown"
                retained_by_modality[modality] = int(retained_by_modality.get(modality, 0)) + 1
            compound_scores = accumulate_gene_scores(evidence_by_ref=evidence_by_ref, ref_to_gene=compound_map)
            genetic_scores = accumulate_gene_scores(evidence_by_ref=evidence_by_ref, ref_to_gene=genetic_map)
            compound_norm = l1_normalize_scores(compound_scores)
            genetic_norm = l1_normalize_scores(genetic_scores)
            compound_neighbors = int(retained_by_modality.get("compound", 0))
            genetic_neighbors = int(retained_by_modality.get("orf", 0) + retained_by_modality.get("crispr", 0))
            if compound_norm and genetic_norm:
                cw = float(cfg.compound_weight)
                gw = float(cfg.genetic_weight)
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
            gene_scores = {
                gene: float(cw * compound_norm.get(gene, 0.0) + gw * genetic_norm.get(gene, 0.0))
                for gene in genes
                if float(cw * compound_norm.get(gene, 0.0) + gw * genetic_norm.get(gene, 0.0)) > 0.0
            }
            selected_gene_ids = _select_gene_ids(gene_scores, cfg)
            selected_gene_ids = sorted(selected_gene_ids, key=lambda g: (-float(gene_scores.get(g, 0.0)), str(g)))
            weights = _selected_weights(gene_scores, selected_gene_ids, cfg.normalize)
            selected_rows = [
                {"gene_id": gene_id, "gene_symbol": gene_id, "score": float(gene_scores.get(gene_id, 0.0)), "weight": float(weights.get(gene_id, 0.0)), "rank": rank}
                for rank, gene_id in enumerate(selected_gene_ids, start=1)
            ]
            full_gene_ids = sorted(gene_scores, key=lambda g: (-float(gene_scores[g]), str(g)))
            full_rows = [
                {"gene_id": gene_id, "gene_symbol": gene_id, "score": float(gene_scores[gene_id]), "rank": rank}
                for rank, gene_id in enumerate(full_gene_ids, start=1)
            ]
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
            specificity_confidence, top10_gene_mass, neighbor_target_concentration = _specificity_confidence(
                selected_rows,
                gene_scores,
            )
            if retrieval_confidence == "high" and specificity_confidence == "low":
                warning_payload = {
                    "code": "diffuse_gene_evidence",
                    "query_id": query_id,
                    "polarity": polarity,
                    "suggestion": "neighbors are geometrically coherent but gene evidence is diffuse; inspect hubness/context and target routing",
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

            safe_program_id = _safe_component(f"{query_id}__polarity={polarity}", "program")
            program_dir = out_dir / f"program={safe_program_id}"
            program_dir.mkdir(parents=True, exist_ok=True)
            _write_rows(program_dir / "geneset.tsv", selected_rows)
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
                "retained_neighbors_by_modality": retained_by_modality,
                "neg_similarity_fraction": neg_frac,
                "top_neighbor_ids": top_neighbor_ids,
                "top_neighbor_modalities": top_neighbor_modalities,
                "top_neighbor_similarities": top_neighbor_sims,
                "top5_neighbor_score_mass_fraction": top5_mass,
                "retrieval_confidence": retrieval_confidence,
                "specificity_confidence": specificity_confidence,
                "top10_gene_mass": top10_gene_mass,
                "neighbor_target_concentration": neighbor_target_concentration,
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
                    "hubness_penalty": cfg.hubness_penalty,
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
                    "n_positive_neighbors": n_positive_neighbors,
                    "n_negative_neighbors": n_negative_neighbors,
                    "neg_similarity_fraction": neg_frac,
                    "top_neighbor_ids": top_neighbor_ids,
                    "top_neighbor_modalities": top_neighbor_modalities,
                    "top_neighbor_similarities": top_neighbor_sims,
                    "top10_gene_mass": top10_gene_mass,
                    "neighbor_target_concentration": neighbor_target_concentration,
                    "hubness_penalty": cfg.hubness_penalty,
                    "hub_score_present": bool(hub_scores_present > 0),
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
                    "score_definition": f"morphology reference similarity polarity={polarity}",
                    "retrieval_confidence": retrieval_confidence,
                    "specificity_confidence": specificity_confidence,
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
                "n_query_profiles": len(query_membership.get(query_id, [])),
                "n_reference_profiles": len(ref_vectors),
                "n_retained_neighbors": len(retained_neighbors),
                "retrieval_confidence": retrieval_confidence,
                "specificity_confidence": specificity_confidence,
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
                "n_query_profiles",
                "n_reference_profiles",
                "n_retained_neighbors",
                "retrieval_confidence",
                "specificity_confidence",
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
        "n_groups": len(manifest_rows),
        "parse_summary": parse_summary,
        "feature_alignment": feature_summary,
        "mapping_summary": mapping_summary,
        "bundle_contexts": parse_summary.get("bundle_manifest", {}).get("contexts") if isinstance(parse_summary.get("bundle_manifest"), dict) else None,
        "bundle_summary": parse_summary.get("bundle_manifest", {}).get("summary") if isinstance(parse_summary.get("bundle_manifest"), dict) else None,
        "resources": resources_info,
        "warnings": root_warnings,
    }
    run_json_path, run_txt_path = write_run_summary_files(out_dir, root_summary)
    return {"n_groups": len(manifest_rows), "out_dir": str(out_dir), "run_summary_json": str(run_json_path), "run_summary_txt": str(run_txt_path)}
