from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
import sys
from typing import Iterable

import numpy as np

from geneset_extractors.core.gmt import build_gmt_sets_from_rows, resolve_gmt_out_path, write_gmt
from geneset_extractors.core.metadata import enrich_manifest_row, input_file_record, make_metadata, write_metadata
from geneset_extractors.core.qc import write_run_summary_files
from geneset_extractors.core.selection import global_l1_weights, ranked_gene_ids, select_quantile, select_threshold, select_top_k, within_set_l1_weights
from geneset_extractors.extractors.calorimetry.bundle import bundle_resources_info
from geneset_extractors.extractors.calorimetry.contrasts import ContrastSignature, build_group_contrasts
from geneset_extractors.extractors.calorimetry.features import CalorimetryFeatureConfig, FeatureExtractionResult, PROGRAM_VARIABLES, extract_subject_features, feature_base_variable
from geneset_extractors.extractors.calorimetry.io import read_calr_data_csv, read_exclusions_tsv, read_session_csv
from geneset_extractors.extractors.calorimetry.orthology import (
    ORTHOLOG_RESOURCE_ID,
    map_gene_scores_to_human,
    ortholog_resource_record,
    read_mouse_human_orthologs_tsv,
    read_packaged_mouse_human_orthologs_tsv,
)
from geneset_extractors.extractors.calorimetry.ontology import (
    build_gene_edge_index,
    build_term_neighbors,
    build_term_template_index,
    read_phenotype_gene_edges_tsv,
    read_term_hierarchy_tsv,
    read_term_templates_tsv,
    route_terms_to_genes,
    score_terms,
    top_terms,
)


@dataclass
class CalRWorkflowConfig:
    converter_name: str
    out_dir: Path
    organism: str
    genome_build: str
    dataset_label: str
    signature_name: str
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
    gmt_split_signed: bool
    gmt_format: str
    emit_small_gene_sets: bool
    analysis_start_hour: float | None
    analysis_end_hour: float | None
    photoperiod_lights_on_hour: float | None
    photoperiod_hours_light: float
    exploratory_without_session: bool
    exclusions_tsv: str | None
    mass_covariate: str | None
    min_group_size: int
    output_gene_species: str = "human"
    ortholog_policy: str = "unique_only"
    similarity_metric: str = "cosine"
    similarity_floor: float = 0.0
    similarity_power: float = 1.0
    hubness_penalty: str = "inverse_linear"
    provenance_mismatch_penalty: float = 0.15
    reference_bundle_id: str | None = None


def _resources_info_object(resources_info: dict[str, object] | None) -> dict[str, object] | None:
    if resources_info is None:
        return {"manifest": None, "resources_dir": None, "used": [], "missing": [], "warnings": []}
    resources_info.setdefault("used", [])
    resources_info.setdefault("missing", [])
    resources_info.setdefault("warnings", [])
    return resources_info


def _append_resource_use(resources_info: dict[str, object] | None, record: dict[str, object]) -> None:
    payload = _resources_info_object(resources_info)
    if payload is None:
        return
    used = payload.setdefault("used", [])
    record_path = str(record.get("path", ""))
    record_id = str(record.get("id", ""))
    for existing in used:
        if str(existing.get("id", "")) == record_id and str(existing.get("path", "")) == record_path:
            return
    used.append(record)


def _load_mouse_human_orthologs(
    *,
    output_gene_species: str,
    mouse_human_orthologs_tsv: str | None,
    resources_info: dict[str, object] | None,
) -> tuple[dict[str, object] | None, dict[str, object]]:
    if str(output_gene_species).strip().lower() != "human":
        return None, {"output_gene_species": "source", "orthology_applied": False}
    if mouse_human_orthologs_tsv:
        orthologs = read_mouse_human_orthologs_tsv(mouse_human_orthologs_tsv)
        _append_resource_use(resources_info, ortholog_resource_record(path=str(mouse_human_orthologs_tsv), method="explicit_or_bundle"))
        return orthologs, {
            "output_gene_species": "human",
            "orthology_applied": True,
            "ortholog_resource": orthologs.get("summary", {}),
        }
    orthologs = read_packaged_mouse_human_orthologs_tsv()
    ortholog_summary = orthologs.get("summary", {})
    _append_resource_use(
        resources_info,
        ortholog_resource_record(path=str(ortholog_summary.get("path", "")), method="packaged_default"),
    )
    return orthologs, {
        "output_gene_species": "human",
        "orthology_applied": True,
        "ortholog_resource": ortholog_summary,
    }


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys()) if rows else ["gene_id", "gene_symbol", "score", "weight", "rank"]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _select_gene_ids(scores: dict[str, float], cfg: CalRWorkflowConfig) -> list[str]:
    method = str(cfg.select)
    if method == "none":
        return ranked_gene_ids(scores)
    if method == "top_k":
        return select_top_k(scores, int(cfg.top_k))
    if method == "quantile":
        return select_quantile(scores, float(cfg.quantile))
    if method == "threshold":
        return select_threshold(scores, float(cfg.min_score))
    raise ValueError(f"Unsupported selection method: {cfg.select}")


def _selected_weights(magnitude_scores: dict[str, float], selected_gene_ids: list[str], normalize: str) -> dict[str, float]:
    if normalize == "none":
        return {gene_id: float(magnitude_scores.get(gene_id, 0.0)) for gene_id in selected_gene_ids}
    if normalize == "l1":
        return global_l1_weights({gene_id: float(magnitude_scores.get(gene_id, 0.0)) for gene_id in selected_gene_ids})
    if normalize == "within_set_l1":
        return within_set_l1_weights(magnitude_scores, selected_gene_ids)
    raise ValueError(f"Unsupported normalization method: {normalize}")


def _program_feature_scores(feature_scores: dict[str, float], program_name: str) -> dict[str, float]:
    if program_name == "global":
        return dict(feature_scores)
    keep = PROGRAM_VARIABLES.get(program_name, set())
    return {
        feature_name: float(score)
        for feature_name, score in feature_scores.items()
        if feature_base_variable(feature_name) in keep
    }


def _read_feature_profiles(path: str | Path) -> tuple[dict[str, dict[str, float]], list[str], dict[str, object]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Reference feature table missing header: {path}")
        fieldnames = [str(name) for name in reader.fieldnames]
        id_col = fieldnames[0]
        feature_cols = fieldnames[1:]
        profiles: dict[str, dict[str, float]] = {}
        for row in reader:
            reference_id = str(row.get(id_col, "")).strip()
            if not reference_id:
                continue
            profiles[reference_id] = {
                feature: float(row.get(feature, 0.0) or 0.0)
                for feature in feature_cols
                if str(row.get(feature, "")).strip() not in {"", "NA", "NaN", "nan"}
            }
    return profiles, feature_cols, {"path": str(path), "n_profiles": len(profiles), "n_features": len(feature_cols)}


def _read_reference_metadata(path: str | Path) -> tuple[dict[str, dict[str, str]], dict[str, object]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Reference metadata table missing header: {path}")
        fieldnames = [str(name) for name in reader.fieldnames]
        id_col = fieldnames[0]
        rows = {}
        for row in reader:
            reference_id = str(row.get(id_col, "")).strip()
            if not reference_id:
                continue
            rows[reference_id] = {str(k): str(v or "") for k, v in row.items()}
    return rows, {"path": str(path), "n_rows": len(rows)}


def _read_feature_schema(path: str | Path | None) -> tuple[list[str] | None, dict[str, object] | None]:
    if not path:
        return None, None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    features = [str(row.get("feature_name", row.get("feature", ""))).strip() for row in rows if str(row.get("feature_name", row.get("feature", ""))).strip()]
    return features, {"path": str(path), "n_features": len(features)}


def _read_feature_stats(path: str | Path | None) -> tuple[dict[str, tuple[float, float]] | None, dict[str, object] | None]:
    if not path:
        return None, None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    stats: dict[str, tuple[float, float]] = {}
    for row in rows:
        feature = str(row.get("feature_name", row.get("feature", ""))).strip()
        if not feature:
            continue
        center = float(row.get("center", 0.0) or 0.0)
        scale = float(row.get("scale", 1.0) or 1.0)
        stats[feature] = (center, scale if abs(scale) > 1e-12 else 1.0)
    return stats, {"path": str(path), "n_features": len(stats)}


def _align_and_standardize(
    feature_scores: dict[str, float],
    *,
    feature_schema: list[str] | None,
    feature_stats: dict[str, tuple[float, float]] | None,
) -> tuple[dict[str, float], dict[str, object]]:
    feature_names = feature_schema or sorted(feature_scores)
    aligned: dict[str, float] = {}
    matched = 0
    for feature in feature_names:
        if feature not in feature_scores:
            aligned[feature] = 0.0
            continue
        matched += 1
        value = float(feature_scores[feature])
        if feature_stats and feature in feature_stats:
            center, scale = feature_stats[feature]
            aligned[feature] = float((value - center) / scale)
        else:
            aligned[feature] = value
    extra = sorted(set(feature_scores) - set(feature_names))
    return aligned, {
        "n_schema_features": len(feature_names),
        "n_query_features_matched": matched,
        "n_query_features_extra": len(extra),
        "extra_features": extra[:10],
    }


def _cosine_similarity(query: dict[str, float], ref: dict[str, float]) -> float:
    features = sorted(set(query) | set(ref))
    if not features:
        return 0.0
    q = np.asarray([float(query.get(feature, 0.0)) for feature in features], dtype=float)
    r = np.asarray([float(ref.get(feature, 0.0)) for feature in features], dtype=float)
    denom = float(np.linalg.norm(q) * np.linalg.norm(r))
    if denom <= 0.0:
        return 0.0
    return float(np.dot(q, r) / denom)


def _pearson_similarity(query: dict[str, float], ref: dict[str, float]) -> float:
    features = sorted(set(query) | set(ref))
    if not features:
        return 0.0
    q = np.asarray([float(query.get(feature, 0.0)) for feature in features], dtype=float)
    r = np.asarray([float(ref.get(feature, 0.0)) for feature in features], dtype=float)
    if q.size < 2:
        return 0.0
    q = q - float(np.mean(q))
    r = r - float(np.mean(r))
    denom = float(np.linalg.norm(q) * np.linalg.norm(r))
    if denom <= 0.0:
        return 0.0
    return float(np.dot(q, r) / denom)


def _similarity(query: dict[str, float], ref: dict[str, float], metric: str) -> float:
    if metric == "cosine":
        return _cosine_similarity(query, ref)
    if metric == "pearson":
        return _pearson_similarity(query, ref)
    raise ValueError(f"Unsupported similarity_metric: {metric}")


def _apply_hubness_penalty(score: float, hub_score: float, mode: str) -> float:
    if mode == "none":
        return float(score)
    if mode == "inverse_linear":
        return float(score) / (1.0 + max(0.0, float(hub_score)))
    if mode == "inverse_rank":
        return float(score) / (1.0 + max(0.0, float(hub_score)) * 2.0)
    raise ValueError(f"Unsupported hubness_penalty: {mode}")


def _provenance_penalty(reference_meta: dict[str, str], query_meta: dict[str, object], mismatch_penalty: float) -> tuple[float, list[str]]:
    penalty = 1.0
    reasons: list[str] = []
    query_cov = str(query_meta.get("mass_covariate", "")).strip()
    ref_cov = str(reference_meta.get("mass_covariate", "")).strip()
    if query_cov and ref_cov and query_cov != ref_cov:
        penalty *= max(0.0, 1.0 - mismatch_penalty)
        reasons.append("mass_covariate_mismatch")
    query_acclimation = str(query_meta.get("acclimation_state", "")).strip()
    ref_acclimation = str(reference_meta.get("acclimation_state", "")).strip()
    if query_acclimation and ref_acclimation and query_acclimation != ref_acclimation:
        penalty *= max(0.0, 1.0 - mismatch_penalty)
        reasons.append("acclimation_mismatch")
    q_temp = query_meta.get("ambient_temperature")
    ref_temp = reference_meta.get("ambient_temperature")
    try:
        if q_temp is not None and ref_temp not in {None, ""} and abs(float(q_temp) - float(ref_temp)) > 3.0:
            penalty *= max(0.0, 1.0 - mismatch_penalty)
            reasons.append("ambient_temperature_mismatch")
    except (TypeError, ValueError):
        pass
    return penalty, reasons


def _render_gene_rows(gene_scores: dict[str, float], gene_symbols: dict[str, str], cfg: CalRWorkflowConfig) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    magnitude = {gene_id: abs(float(score)) for gene_id, score in gene_scores.items()}
    selected_gene_ids = _select_gene_ids(magnitude, cfg)
    weights = _selected_weights(magnitude, selected_gene_ids, cfg.normalize)
    full_rows = []
    for rank, gene_id in enumerate(sorted(gene_scores, key=lambda gid: (-abs(float(gene_scores[gid])), str(gid))), start=1):
        full_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene_symbols.get(gene_id, gene_id),
                "score": float(gene_scores[gene_id]),
                "rank": rank,
            }
        )
    selected_rows = []
    for rank, gene_id in enumerate(selected_gene_ids, start=1):
        selected_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene_symbols.get(gene_id, gene_id),
                "score": float(gene_scores[gene_id]),
                "weight": float(weights.get(gene_id, 0.0)),
                "rank": rank,
            }
        )
    return selected_rows, full_rows


def _emit_program_output(
    *,
    cfg: CalRWorkflowConfig,
    out_dir: Path,
    contrast: ContrastSignature,
    program_name: str,
    mode_label: str,
    gene_scores: dict[str, float],
    gene_symbols: dict[str, str],
    input_files: list[dict[str, str]],
    summary_extra: dict[str, object],
    parse_summary: dict[str, object],
    resources_info: dict[str, object] | None,
) -> dict[str, object]:
    if not gene_scores:
        raise ValueError(f"No gene scores available for {contrast.contrast_id} {program_name} {mode_label}")
    selected_rows, full_rows = _render_gene_rows(gene_scores, gene_symbols, cfg)
    _write_rows(out_dir / "geneset.tsv", selected_rows)
    output_files = [{"path": str(out_dir / "geneset.tsv"), "role": "selected_gene_program"}]
    if cfg.emit_full:
        _write_rows(out_dir / "geneset.full.tsv", full_rows)
        output_files.append({"path": str(out_dir / "geneset.full.tsv"), "role": "full_scores"})

    gmt_sets: list[tuple[str, list[str]]] = []
    gmt_path: Path | None = None
    gmt_plans: list[dict[str, object]] = []
    gmt_diagnostics: list[dict[str, object]] = []
    if cfg.emit_gmt:
        gmt_sets, gmt_plans = build_gmt_sets_from_rows(
            full_rows,
            base_name=f"{cfg.converter_name}__contrast={contrast.contrast_id}__program={program_name}__mode={mode_label}",
            prefer_symbol=bool(cfg.gmt_prefer_symbol),
            min_genes=int(cfg.gmt_min_genes),
            max_genes=int(cfg.gmt_max_genes),
            topk_list=[int(value) for value in str(cfg.gmt_topk_list).split(",") if str(value).strip()],
            mass_list=[],
            split_signed=bool(cfg.gmt_split_signed),
            require_symbol=bool(cfg.gmt_require_symbol),
            emit_small_gene_sets=bool(cfg.emit_small_gene_sets),
            diagnostics=gmt_diagnostics,
            context={"program_method": program_name, "contrast_method": contrast.contrast_id, "link_method": mode_label},
        )
        if gmt_sets:
            gmt_path = resolve_gmt_out_path(out_dir, cfg.gmt_out)
            write_gmt(gmt_sets, gmt_path, gmt_format=cfg.gmt_format)
            output_files.append({"path": str(gmt_path), "role": "gmt"})

    run_payload = {
        "converter": cfg.converter_name,
        "dataset_label": cfg.dataset_label,
        "group": contrast.group_label,
        "primary_program_method": program_name,
        "selected_direction": mode_label,
        **contrast.metadata,
        **summary_extra,
        "resources": resources_info,
        "gmt_diagnostics": gmt_diagnostics,
    }
    run_json, run_txt = write_run_summary_files(out_dir, run_payload)
    output_files.append({"path": str(run_json), "role": "run_summary_json"})
    output_files.append({"path": str(run_txt), "role": "run_summary_text"})

    assigned_features = len(contrast.feature_scores)
    fraction_features_assigned = 1.0 if assigned_features > 0 else 0.0
    alignment_summary = summary_extra.get("alignment_summary")
    if isinstance(alignment_summary, dict):
        try:
            assigned_features = int(alignment_summary.get("n_query_features_matched", assigned_features))
            schema_features = int(alignment_summary.get("n_schema_features", len(contrast.feature_scores)))
            fraction_features_assigned = float(assigned_features) / float(schema_features) if schema_features > 0 else 0.0
        except (TypeError, ValueError):
            pass

    meta = make_metadata(
        converter_name=cfg.converter_name,
        parameters={
            "dataset_label": cfg.dataset_label,
            "signature_name": cfg.signature_name,
            "program": program_name,
            "mode": mode_label,
            "select": cfg.select,
            "top_k": cfg.top_k,
            "normalize": cfg.normalize,
            "reference_bundle_id": cfg.reference_bundle_id,
        },
        data_type="calorimetry",
        assay="indirect_calorimetry",
        organism=cfg.organism,
        genome_build=cfg.genome_build,
        files=input_files,
        gene_annotation={"mode": "none", "source": cfg.converter_name, "gene_id_field": "gene_id"},
        weights={
            "weight_type": "signed",
            "normalization": {"method": cfg.normalize, "target_sum": 1.0 if cfg.normalize in {"within_set_l1", "l1"} else None},
            "aggregation": mode_label,
        },
        summary={
            "n_input_features": len(contrast.feature_scores),
            "n_genes": len(selected_rows),
            "n_features_assigned": assigned_features,
            "fraction_features_assigned": fraction_features_assigned,
            **contrast.metadata,
            **summary_extra,
            "parse_summary": parse_summary,
            "resources": resources_info,
        },
        output_files=output_files,
        gmt={
            "written": bool(gmt_sets),
            "path": (
                str(gmt_path.relative_to(out_dir))
                if gmt_path is not None and gmt_path.is_relative_to(out_dir)
                else (str(gmt_path) if gmt_path is not None else None)
            ),
            "prefer_symbol": bool(cfg.gmt_prefer_symbol),
            "require_symbol": bool(cfg.gmt_require_symbol),
            "biotype_allowlist": cfg.gmt_biotype_allowlist,
            "min_genes": int(cfg.gmt_min_genes),
            "max_genes": int(cfg.gmt_max_genes),
            "emit_small_gene_sets": bool(cfg.emit_small_gene_sets),
            "plans": gmt_plans,
            "diagnostics": gmt_diagnostics,
        },
    )
    write_metadata(out_dir / "geneset.meta.json", meta)
    return {
        "path": str(out_dir),
        "n_genes": len(selected_rows),
        "meta": meta,
        "gmt_sets": gmt_sets,
    }


def _bundle_file(bundle_manifest_path: Path, payload: dict[str, object], key: str) -> str | None:
    files = payload.get("files", {}) if isinstance(payload, dict) else {}
    if not isinstance(files, dict):
        return None
    rel = str(files.get(key, "")).strip()
    if not rel:
        return None
    return str(bundle_manifest_path.parent / rel)


def run_calr_ontology_workflow(
    *,
    cfg: CalRWorkflowConfig,
    calr_data_csv: str,
    session_csv: str | None,
    exclusions_tsv: str | None,
    term_templates_tsv: str,
    phenotype_gene_edges_tsv: str,
    term_hierarchy_tsv: str | None,
    mouse_human_orthologs_tsv: str | None,
    input_files: list[dict[str, object]] | None,
    resources_info: dict[str, object] | None,
) -> dict[str, object]:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if cfg.output_gene_species == "human":
        resources_info = _resources_info_object(resources_info)
    data_fieldnames, data_rows = read_calr_data_csv(calr_data_csv)
    session_fieldnames, session_rows = read_session_csv(session_csv)
    exclusion_rows = read_exclusions_tsv(exclusions_tsv)
    extracted = extract_subject_features(
        data_rows=data_rows,
        data_fieldnames=data_fieldnames,
        session_rows=session_rows,
        session_fieldnames=session_fieldnames,
        exclusion_rows=exclusion_rows,
        cfg=CalorimetryFeatureConfig(
            analysis_start_hour=cfg.analysis_start_hour,
            analysis_end_hour=cfg.analysis_end_hour,
            photoperiod_lights_on_hour=cfg.photoperiod_lights_on_hour,
            photoperiod_hours_light=cfg.photoperiod_hours_light,
            exclusions_tsv=cfg.exclusions_tsv,
            exploratory_without_session=cfg.exploratory_without_session,
        ),
    )
    contrasts, contrast_warnings = build_group_contrasts(
        subjects_by_id=extracted.subjects,
        explicit_mass_covariate=cfg.mass_covariate,
        min_group_size=cfg.min_group_size,
    )
    template_rows, template_summary = read_term_templates_tsv(term_templates_tsv)
    edge_rows, edge_summary = read_phenotype_gene_edges_tsv(phenotype_gene_edges_tsv)
    hierarchy_rows, hierarchy_summary = read_term_hierarchy_tsv(term_hierarchy_tsv)
    templates = build_term_template_index(template_rows)
    gene_edges = build_gene_edge_index(edge_rows)
    term_neighbors = build_term_neighbors(hierarchy_rows)
    orthologs, orthology_base_summary = _load_mouse_human_orthologs(
        output_gene_species=cfg.output_gene_species,
        mouse_human_orthologs_tsv=mouse_human_orthologs_tsv,
        resources_info=resources_info,
    )

    if input_files is None:
        input_files = [input_file_record(calr_data_csv, "calr_data_csv")]
        if session_csv:
            input_files.append(input_file_record(session_csv, "session_csv"))
        if exclusions_tsv:
            input_files.append(input_file_record(exclusions_tsv, "exclusions_tsv"))
        input_files.append(input_file_record(term_templates_tsv, "term_templates_tsv"))
        input_files.append(input_file_record(phenotype_gene_edges_tsv, "phenotype_gene_edges_tsv"))
        if term_hierarchy_tsv:
            input_files.append(input_file_record(term_hierarchy_tsv, "term_hierarchy_tsv"))
        if mouse_human_orthologs_tsv:
            input_files.append(input_file_record(mouse_human_orthologs_tsv, "mouse_human_orthologs_tsv"))

    manifest_rows: list[dict[str, object]] = []
    combined_gmt: list[tuple[str, list[str]]] = []
    for contrast in contrasts:
        for program_name in ("global", "thermogenesis", "substrate_use", "intake_balance", "activity_circadian"):
            program_scores = _program_feature_scores(contrast.feature_scores, program_name)
            if not program_scores:
                continue
            scored_terms = score_terms(feature_scores=program_scores, templates=templates, program_name=program_name)
            if scored_terms and float(scored_terms[0].get("score", 0.0)) < 0.05:
                print(f"warning: sparse ontology support for contrast={contrast.contrast_id} program={program_name}", file=sys.stderr)
            for mode_label in ("core", "expanded"):
                gene_scores, gene_symbols, routing_summary = route_terms_to_genes(
                    scored_terms=scored_terms,
                    gene_edges=gene_edges,
                    term_neighbors=term_neighbors,
                    mode=mode_label,
                )
                orthology_summary = dict(orthology_base_summary)
                if cfg.output_gene_species == "human":
                    gene_scores, gene_symbols, mapping_summary = map_gene_scores_to_human(
                        gene_scores=gene_scores,
                        gene_symbols=gene_symbols,
                        orthologs=orthologs,
                        policy=cfg.ortholog_policy,
                        fallback_to_source=False,
                    )
                    orthology_summary.update(mapping_summary)
                if not gene_scores:
                    continue
                child_dir = out_dir / f"program={program_name}__mode={mode_label}__contrast={contrast.contrast_id}"
                result = _emit_program_output(
                    cfg=cfg,
                    out_dir=child_dir,
                    contrast=contrast,
                    program_name=program_name,
                    mode_label=mode_label,
                    gene_scores=gene_scores,
                    gene_symbols=gene_symbols,
                    input_files=input_files,
                    summary_extra={
                        "session_mode": extracted.metadata_summary["session_mode"],
                        "session_group_layout_mode": extracted.metadata_summary["session_group_layout_mode"],
                        "session_window_layout_mode": extracted.metadata_summary["session_window_layout_mode"],
                        "analysis_start": extracted.metadata_summary["analysis_start"],
                        "analysis_end": extracted.metadata_summary["analysis_end"],
                        "analysis_window_source": extracted.metadata_summary["analysis_window_source"],
                        "photoperiod_source": extracted.metadata_summary["photoperiod_source"],
                        "ee_source": extracted.metadata_summary["ee_source"],
                        "rer_source": extracted.metadata_summary["rer_source"],
                        "excluded_subjects": extracted.metadata_summary["excluded_subjects"],
                        "top_terms": top_terms(scored_terms),
                        "routing_summary": routing_summary,
                        "output_gene_species": cfg.output_gene_species,
                        "orthology_summary": orthology_summary,
                        "warnings": extracted.warnings + contrast_warnings,
                    },
                    parse_summary={
                        "feature_extraction": extracted.metadata_summary,
                        "term_templates": template_summary,
                        "phenotype_gene_edges": edge_summary,
                        "term_hierarchy": hierarchy_summary,
                    },
                    resources_info=resources_info,
                )
                manifest_rows.append(
                    {
                        "contrast_id": contrast.contrast_id,
                        "contrast_label": contrast.contrast_label,
                        "program": program_name,
                        "mode": mode_label,
                        "path": child_dir.relative_to(out_dir),
                        "n_genes": result["n_genes"],
                    }
                )
                combined_gmt.extend(result["gmt_sets"])
    with (out_dir / "manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
        rows = [enrich_manifest_row(out_dir, out_dir / str(row["path"]), dict(row)) for row in manifest_rows]
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=["contrast_id", "contrast_label", "program", "mode", "geneset_id", "label", "path", "meta_path", "provenance_path", "focus_node_id", "n_genes"],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "path": str(row["path"])})
    if combined_gmt and cfg.emit_gmt:
        write_gmt(combined_gmt, out_dir / "genesets.gmt", gmt_format=cfg.gmt_format)
    root_summary = {
        "converter": cfg.converter_name,
        "dataset_label": cfg.dataset_label,
        "session_mode": extracted.metadata_summary["session_mode"],
        "session_group_layout_mode": extracted.metadata_summary["session_group_layout_mode"],
        "session_window_layout_mode": extracted.metadata_summary["session_window_layout_mode"],
        "analysis_start": extracted.metadata_summary["analysis_start"],
        "analysis_end": extracted.metadata_summary["analysis_end"],
        "analysis_window_source": extracted.metadata_summary["analysis_window_source"],
        "photoperiod_source": extracted.metadata_summary["photoperiod_source"],
        "design_class": extracted.metadata_summary["design_class"],
        "run_ids": extracted.metadata_summary["run_ids"],
        "n_subjects": extracted.metadata_summary["n_subjects"],
        "n_outputs": len(manifest_rows),
        "output_gene_species": cfg.output_gene_species,
        "orthology_summary": orthology_base_summary,
        "warnings": extracted.warnings + contrast_warnings,
        "resources": resources_info,
    }
    write_run_summary_files(out_dir, root_summary)
    return {"out_dir": str(out_dir), "n_groups": len(manifest_rows)}


def run_calr_profile_workflow(
    *,
    cfg: CalRWorkflowConfig,
    calr_data_csv: str,
    session_csv: str | None,
    exclusions_tsv: str | None,
    reference_profiles_tsv: str,
    reference_metadata_tsv: str,
    feature_schema_tsv: str | None,
    feature_stats_tsv: str | None,
    mouse_human_orthologs_tsv: str | None,
    input_files: list[dict[str, object]] | None,
    resources_info: dict[str, object] | None,
) -> dict[str, object]:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if cfg.output_gene_species == "human":
        resources_info = _resources_info_object(resources_info)
    data_fieldnames, data_rows = read_calr_data_csv(calr_data_csv)
    session_fieldnames, session_rows = read_session_csv(session_csv)
    exclusion_rows = read_exclusions_tsv(exclusions_tsv)
    extracted = extract_subject_features(
        data_rows=data_rows,
        data_fieldnames=data_fieldnames,
        session_rows=session_rows,
        session_fieldnames=session_fieldnames,
        exclusion_rows=exclusion_rows,
        cfg=CalorimetryFeatureConfig(
            analysis_start_hour=cfg.analysis_start_hour,
            analysis_end_hour=cfg.analysis_end_hour,
            photoperiod_lights_on_hour=cfg.photoperiod_lights_on_hour,
            photoperiod_hours_light=cfg.photoperiod_hours_light,
            exclusions_tsv=cfg.exclusions_tsv,
            exploratory_without_session=cfg.exploratory_without_session,
        ),
    )
    contrasts, contrast_warnings = build_group_contrasts(
        subjects_by_id=extracted.subjects,
        explicit_mass_covariate=cfg.mass_covariate,
        min_group_size=cfg.min_group_size,
    )
    reference_profiles, _ref_features, ref_summary = _read_feature_profiles(reference_profiles_tsv)
    reference_metadata, ref_meta_summary = _read_reference_metadata(reference_metadata_tsv)
    feature_schema, schema_summary = _read_feature_schema(feature_schema_tsv)
    feature_stats, feature_stats_summary = _read_feature_stats(feature_stats_tsv)
    orthologs, orthology_base_summary = _load_mouse_human_orthologs(
        output_gene_species=cfg.output_gene_species,
        mouse_human_orthologs_tsv=mouse_human_orthologs_tsv,
        resources_info=resources_info,
    )
    if input_files is None:
        input_files = [input_file_record(calr_data_csv, "calr_data_csv")]
        if session_csv:
            input_files.append(input_file_record(session_csv, "session_csv"))
        if exclusions_tsv:
            input_files.append(input_file_record(exclusions_tsv, "exclusions_tsv"))
        input_files.append(input_file_record(reference_profiles_tsv, "reference_profiles_tsv"))
        input_files.append(input_file_record(reference_metadata_tsv, "reference_metadata_tsv"))
        if feature_schema_tsv:
            input_files.append(input_file_record(feature_schema_tsv, "feature_schema_tsv"))
        if feature_stats_tsv:
            input_files.append(input_file_record(feature_stats_tsv, "feature_stats_tsv"))
        if mouse_human_orthologs_tsv:
            input_files.append(input_file_record(mouse_human_orthologs_tsv, "mouse_human_orthologs_tsv"))

    manifest_rows: list[dict[str, object]] = []
    combined_gmt: list[tuple[str, list[str]]] = []
    for contrast in contrasts:
        for program_name in ("global", "thermogenesis", "substrate_use", "intake_balance", "activity_circadian"):
            raw_program_scores = _program_feature_scores(contrast.feature_scores, program_name)
            if not raw_program_scores:
                continue
            aligned_query, alignment_summary = _align_and_standardize(raw_program_scores, feature_schema=feature_schema, feature_stats=feature_stats)
            neighbors: list[dict[str, object]] = []
            for reference_id, profile in reference_profiles.items():
                aligned_ref, _ = _align_and_standardize(profile, feature_schema=feature_schema, feature_stats=feature_stats)
                sim = _similarity(aligned_query, aligned_ref, cfg.similarity_metric)
                if sim < float(cfg.similarity_floor):
                    continue
                ref_meta = reference_metadata.get(reference_id, {})
                qc_weight = float(ref_meta.get("qc_weight", 1.0) or 1.0)
                hub_score = float(ref_meta.get("hub_score", 0.0) or 0.0)
                prov_penalty, mismatch_reasons = _provenance_penalty(ref_meta, contrast.metadata, cfg.provenance_mismatch_penalty)
                evidence = _apply_hubness_penalty(max(0.0, sim) ** float(cfg.similarity_power), hub_score, cfg.hubness_penalty)
                evidence *= qc_weight * prov_penalty
                if evidence <= 0.0:
                    continue
                has_output_gene = bool(str(ref_meta.get("output_gene_id", "")).strip() or str(ref_meta.get("output_gene_symbol", "")).strip())
                neighbors.append(
                    {
                        "reference_id": reference_id,
                        "similarity": sim,
                        "evidence": evidence,
                        "gene_id": str(
                            ref_meta.get(
                                "output_gene_id",
                                ref_meta.get("gene_id", ref_meta.get("gene_symbol", reference_id)),
                            )
                        ).strip()
                        or reference_id,
                        "gene_symbol": str(
                            ref_meta.get(
                                "output_gene_symbol",
                                ref_meta.get("gene_symbol", ref_meta.get("gene_id", reference_id)),
                            )
                        ).strip()
                        or reference_id,
                        "gene_species": (
                            str(ref_meta.get("output_gene_species", "")).strip()
                            or ("human" if has_output_gene else "")
                            or str(ref_meta.get("gene_species", ref_meta.get("source_gene_species", cfg.organism))).strip()
                            or str(cfg.organism)
                        ),
                        "mismatch_reasons": mismatch_reasons,
                    }
                )
            neighbors.sort(key=lambda row: (-float(row["evidence"]), -float(row["similarity"]), str(row["reference_id"])))
            if neighbors and float(neighbors[0]["similarity"]) < 0.25:
                print(
                    f"warning: low retrieval confidence for contrast={contrast.contrast_id} program={program_name}; top_similarity={neighbors[0]['similarity']:.3f}",
                    file=sys.stderr,
                )
            gene_scores: dict[str, float] = {}
            gene_symbols: dict[str, str] = {}
            source_gene_scores: dict[str, float] = {}
            source_gene_symbols: dict[str, str] = {}
            for neighbor in neighbors[: max(1, min(25, len(neighbors)))]:
                gene_id = str(neighbor["gene_id"])
                gene_symbol = str(neighbor["gene_symbol"])
                gene_species = str(neighbor.get("gene_species", "")).strip().lower()
                if cfg.output_gene_species == "human" and gene_species not in {"human", ""}:
                    source_gene_scores[gene_id] = float(source_gene_scores.get(gene_id, 0.0)) + float(neighbor["evidence"])
                    source_gene_symbols[gene_id] = gene_symbol
                    continue
                gene_scores[gene_id] = float(gene_scores.get(gene_id, 0.0)) + float(neighbor["evidence"])
                gene_symbols[gene_id] = gene_symbol
            orthology_summary = dict(orthology_base_summary)
            if cfg.output_gene_species == "human" and source_gene_scores:
                mapped_scores, mapped_symbols, mapping_summary = map_gene_scores_to_human(
                    gene_scores=source_gene_scores,
                    gene_symbols=source_gene_symbols,
                    orthologs=orthologs,
                    policy=cfg.ortholog_policy,
                    fallback_to_source=False,
                )
                for gene_id, score in mapped_scores.items():
                    gene_scores[gene_id] = float(gene_scores.get(gene_id, 0.0)) + float(score)
                    if gene_id in mapped_symbols:
                        gene_symbols[gene_id] = str(mapped_symbols[gene_id])
                orthology_summary.update(mapping_summary)
            if not gene_scores:
                continue
            retrieval_confidence = "low"
            top_similarity = float(neighbors[0]["similarity"]) if neighbors else 0.0
            if top_similarity >= 0.55 and len(neighbors) >= 3:
                retrieval_confidence = "high"
            elif top_similarity >= 0.35 and len(neighbors) >= 2:
                retrieval_confidence = "medium"
            child_dir = out_dir / f"program={program_name}__mode=core__contrast={contrast.contrast_id}"
            result = _emit_program_output(
                cfg=cfg,
                out_dir=child_dir,
                contrast=contrast,
                program_name=program_name,
                mode_label="core",
                gene_scores=gene_scores,
                gene_symbols=gene_symbols,
                input_files=input_files,
                summary_extra={
                    "session_mode": extracted.metadata_summary["session_mode"],
                    "session_group_layout_mode": extracted.metadata_summary["session_group_layout_mode"],
                    "session_window_layout_mode": extracted.metadata_summary["session_window_layout_mode"],
                    "analysis_start": extracted.metadata_summary["analysis_start"],
                    "analysis_end": extracted.metadata_summary["analysis_end"],
                    "analysis_window_source": extracted.metadata_summary["analysis_window_source"],
                    "photoperiod_source": extracted.metadata_summary["photoperiod_source"],
                    "ee_source": extracted.metadata_summary["ee_source"],
                    "rer_source": extracted.metadata_summary["rer_source"],
                    "excluded_subjects": extracted.metadata_summary["excluded_subjects"],
                        "retrieval_confidence": retrieval_confidence,
                        "top_neighbors": neighbors[:10],
                        "alignment_summary": alignment_summary,
                        "output_gene_species": cfg.output_gene_species,
                        "orthology_summary": orthology_summary,
                        "warnings": extracted.warnings + contrast_warnings,
                    },
                parse_summary={
                    "feature_extraction": extracted.metadata_summary,
                    "reference_profiles": ref_summary,
                    "reference_metadata": ref_meta_summary,
                    "feature_schema": schema_summary,
                    "feature_stats": feature_stats_summary,
                },
                resources_info=resources_info,
            )
            manifest_rows.append(
                {
                    "contrast_id": contrast.contrast_id,
                    "contrast_label": contrast.contrast_label,
                    "program": program_name,
                    "mode": "core",
                    "path": child_dir.relative_to(out_dir),
                    "n_genes": result["n_genes"],
                    "retrieval_confidence": retrieval_confidence,
                }
            )
            combined_gmt.extend(result["gmt_sets"])
    with (out_dir / "manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
        rows = [enrich_manifest_row(out_dir, out_dir / str(row["path"]), dict(row)) for row in manifest_rows]
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "contrast_id",
                "contrast_label",
                "program",
                "mode",
                "geneset_id",
                "label",
                "path",
                "meta_path",
                "provenance_path",
                "focus_node_id",
                "n_genes",
                "retrieval_confidence",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "path": str(row["path"])})
    if combined_gmt and cfg.emit_gmt:
        write_gmt(combined_gmt, out_dir / "genesets.gmt", gmt_format=cfg.gmt_format)
    root_summary = {
        "converter": cfg.converter_name,
        "dataset_label": cfg.dataset_label,
        "session_mode": extracted.metadata_summary["session_mode"],
        "session_group_layout_mode": extracted.metadata_summary["session_group_layout_mode"],
        "session_window_layout_mode": extracted.metadata_summary["session_window_layout_mode"],
        "analysis_start": extracted.metadata_summary["analysis_start"],
        "analysis_end": extracted.metadata_summary["analysis_end"],
        "analysis_window_source": extracted.metadata_summary["analysis_window_source"],
        "photoperiod_source": extracted.metadata_summary["photoperiod_source"],
        "design_class": extracted.metadata_summary["design_class"],
        "run_ids": extracted.metadata_summary["run_ids"],
        "n_subjects": extracted.metadata_summary["n_subjects"],
        "n_outputs": len(manifest_rows),
        "reference_bundle_id": cfg.reference_bundle_id,
        "output_gene_species": cfg.output_gene_species,
        "orthology_summary": orthology_base_summary,
        "warnings": extracted.warnings + contrast_warnings,
        "resources": resources_info,
    }
    write_run_summary_files(out_dir, root_summary)
    return {"out_dir": str(out_dir), "n_groups": len(manifest_rows)}
