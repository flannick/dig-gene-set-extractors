from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re
import statistics
import sys

from geneset_extractors.core.gmt import (
    build_gmt_sets_from_rows,
    parse_int_list_csv,
    parse_mass_list_csv,
    parse_str_list_csv,
    resolve_gmt_out_path,
    write_gmt,
)
from geneset_extractors.core.metadata import make_metadata, write_metadata
from geneset_extractors.core.qc import write_run_summary_files
from geneset_extractors.core.selection import (
    global_l1_weights,
    ranked_gene_ids,
    select_quantile,
    select_threshold,
    select_top_k,
    within_set_l1_weights,
)
from geneset_extractors.extractors.drug_response.normalize import (
    aggregate_response_records,
    orient_matrix,
    transform_matrix_by_drug,
)
from geneset_extractors.extractors.drug_response.scoring import (
    compute_polypharm_weights,
    compute_target_ubiquity_weights,
    compute_ubiquity_weights,
    score_genes,
)
from geneset_extractors.extractors.drug_response.io import ResponseRecord
from geneset_extractors.io.gtf import read_genes_from_gtf


_SAFE_COMPONENT_RE = re.compile(r"[^A-Za-z0-9._-]+")
_LOW_TARGET_COVERAGE_WARN = 0.5
_BROAD_TARGET_PATTERNS = [
    re.compile(r"^GPR[0-9A-Z]*$"),
    re.compile(r"^HTR[0-9A-Z]*$"),
    re.compile(r"^DRD[0-9A-Z]*$"),
    re.compile(r"^ADRA[0-9A-Z]*$"),
    re.compile(r"^ADRB[0-9A-Z]*$"),
    re.compile(r"^CHRM[0-9A-Z]*$"),
    re.compile(r"^OR[0-9A-Z]+$"),
]


@dataclass
class DrugResponseWorkflowConfig:
    converter_name: str
    out_dir: Path
    organism: str
    genome_build: str
    dataset_label: str
    response_metric: str
    response_direction: str
    response_transform: str
    contrast_method: str
    case_control_within_group: bool
    min_group_size: int
    scoring_model: str
    sparse_alpha: float
    response_ubiquity_penalty: str
    target_ubiquity_penalty: str
    ubiquity_tau: float
    ubiquity_epsilon: float
    polypharm_downweight: bool
    polypharm_t0: int
    max_programs: int
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
    gmt_split_signed: bool | None
    gmt_format: str
    emit_small_gene_sets: bool
    gtf: str | None
    gtf_source: str | None
    gtf_gene_id_field: str


def _safe_component(value: str, fallback: str) -> str:
    out = _SAFE_COMPONENT_RE.sub("_", str(value)).strip("_")
    return out or fallback


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["gene_id", "score", "rank"]
    if any("weight" in row for row in rows):
        fieldnames.append("weight")
    if any(str(row.get("gene_symbol", "")).strip() for row in rows):
        fieldnames.append("gene_symbol")
    if any(str(row.get("gene_biotype", "")).strip() for row in rows):
        fieldnames.append("gene_biotype")
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _select_gene_ids(scores: dict[str, float], cfg: DrugResponseWorkflowConfig) -> list[str]:
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


def _warn_gmt_diagnostics(gmt_diagnostics: list[dict[str, object]]) -> list[dict[str, object]]:
    warnings_payload: list[dict[str, object]] = []
    for diag in gmt_diagnostics:
        code = str(diag.get("code", ""))
        base_name = str(diag.get("base_name", "unknown"))
        n_genes = int(diag.get("n_genes", 0) or 0)
        reason = str(diag.get("reason", "")).strip()
        suggestion = str(diag.get("suggestion", "")).strip()
        message = ""
        if code in {"small_gene_set_skipped", "no_positive_genes"}:
            message = f"skipped GMT output for {base_name}; n_genes={n_genes}. {reason}"
            print(f"warning: {message}", file=sys.stderr)
        elif code == "small_gene_set_emitted":
            message = f"emitted small GMT output for {base_name}; n_genes={n_genes}. {reason}"
            print(f"warning: {message}", file=sys.stderr)
        elif code == "marginal_gene_count":
            message = f"marginal GMT signal for {base_name}; n_genes={n_genes}. {reason}"
            print(f"warning: {message}", file=sys.stderr)
        elif code == "require_symbol_heavy_drop":
            message = f"GMT symbol requirement dropped many rows for {base_name}. {reason}"
            print(f"warning: {message}", file=sys.stderr)
        else:
            continue
        warnings_payload.append(
            {
                "code": code,
                "message": message,
                "base_name": base_name,
                "n_genes": n_genes,
            }
        )
        if suggestion:
            print(f"warning:   {suggestion}", file=sys.stderr)
    return warnings_payload


def _collect_skipped_gmt_outputs(gmt_diagnostics: list[dict[str, object]]) -> list[dict[str, object]]:
    skipped: list[dict[str, object]] = []
    for diag in gmt_diagnostics:
        code = str(diag.get("code", ""))
        if code not in {"small_gene_set_skipped", "no_positive_genes"}:
            continue
        skipped.append(
            {
                "reason": str(diag.get("reason", "")),
                "code": code,
                "n_genes": int(diag.get("n_genes", 0) or 0),
                "variant": str(diag.get("variant", "primary")),
            }
        )
    return skipped


def _mean(values: list[float]) -> float:
    if not values:
        return 0.0
    return float(sum(values) / float(len(values)))


def _is_broad_pharmacology_gene(gene_symbol: str) -> bool:
    token = str(gene_symbol).strip().upper()
    if not token:
        return False
    for pat in _BROAD_TARGET_PATTERNS:
        if pat.match(token):
            return True
    return False


def _build_group_index(
    *,
    samples: list[str],
    groups_by_sample: dict[str, str] | None,
) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    if groups_by_sample is None:
        return out
    for sample in samples:
        group = str(groups_by_sample.get(sample, "")).strip()
        if not group:
            continue
        out.setdefault(group, []).append(sample)
    return out


def _build_programs(
    *,
    contrast_method: str,
    z_matrix: dict[str, dict[str, float]],
    samples: list[str],
    drugs: list[str],
    groups_by_sample: dict[str, str] | None,
    case_control_by_sample: dict[str, bool] | None,
    case_control_within_group: bool,
    min_group_size: int,
) -> tuple[list[dict[str, object]], list[str], list[dict[str, object]], dict[str, dict[str, object]]]:
    warnings: list[str] = []
    programs: list[dict[str, object]] = []
    group_index = _build_group_index(samples=samples, groups_by_sample=groups_by_sample)
    group_qc_rows: list[dict[str, object]] = []
    groups_skipped: dict[str, dict[str, object]] = {}
    min_n = max(1, int(min_group_size))

    if contrast_method == "none":
        for sample_id in samples:
            row = z_matrix.get(sample_id, {})
            drug_scores = {drug: float(value) for drug, value in row.items()}
            programs.append(
                {
                    "program_id": f"sample={sample_id}",
                    "signed": False,
                    "contrast_method": "none",
                    "drug_scores": drug_scores,
                    "n_samples": 1,
                    "criterion": float(len(drug_scores)),
                    "group": "",
                }
            )
        return programs, warnings, group_qc_rows, groups_skipped

    if contrast_method == "group_mean":
        if not group_index:
            raise ValueError("contrast_method=group_mean requires groups (use --groups_tsv or group metadata)")
        for group in sorted(group_index):
            members = group_index[group]
            n_group = len(members)
            n_rest = len(samples) - n_group
            if n_group < min_n:
                msg = f"group_skipped group={group} n={n_group} reason=min_group_size"
                warnings.append(msg)
                groups_skipped[group] = {"n": n_group, "reason": "min_group_size"}
                group_qc_rows.append(
                    {
                        "group": group,
                        "n_group": n_group,
                        "n_rest": n_rest,
                        "n_case": 0,
                        "n_control": 0,
                        "sets_emitted": 0,
                        "reason_if_skipped": "min_group_size",
                    }
                )
                continue
            if n_group < 2:
                warnings.append(
                    f"group={group} has only {n_group} samples for group_mean; estimates may be unstable"
                )
            scores: dict[str, float] = {}
            for drug in drugs:
                values = [float(z_matrix[s].get(drug)) for s in members if drug in z_matrix[s]]
                if not values:
                    continue
                scores[drug] = _mean(values)
            criterion = _mean([abs(v) for v in scores.values()]) if scores else 0.0
            programs.append(
                {
                    "program_id": f"group={group}",
                    "signed": False,
                    "contrast_method": "group_mean",
                    "drug_scores": scores,
                    "n_samples": n_group,
                    "criterion": criterion,
                    "group": group,
                    "group_qc": {
                        "group": group,
                        "n_group": n_group,
                        "n_rest": n_rest,
                        "n_case": 0,
                        "n_control": 0,
                        "reason_if_skipped": "",
                    },
                }
            )
            group_qc_rows.append(
                {
                    "group": group,
                    "n_group": n_group,
                    "n_rest": n_rest,
                    "n_case": 0,
                    "n_control": 0,
                    "sets_emitted": 0,
                    "reason_if_skipped": "",
                }
            )
        return programs, warnings, group_qc_rows, groups_skipped

    if contrast_method == "group_vs_rest":
        if not group_index:
            raise ValueError("contrast_method=group_vs_rest requires groups (use --groups_tsv or group metadata)")
        for group in sorted(group_index):
            members = set(group_index[group])
            n_group = len(members)
            rest = [sample for sample in samples if sample not in members]
            n_rest = len(rest)
            if n_group < min_n:
                msg = f"group_skipped group={group} n={n_group} reason=min_group_size"
                warnings.append(msg)
                groups_skipped[group] = {"n": n_group, "reason": "min_group_size"}
                group_qc_rows.append(
                    {
                        "group": group,
                        "n_group": n_group,
                        "n_rest": n_rest,
                        "n_case": 0,
                        "n_control": 0,
                        "sets_emitted": 0,
                        "reason_if_skipped": "min_group_size",
                    }
                )
                continue
            if n_group < 2 or n_rest < 2:
                warnings.append(
                    f"group={group} has small contrast sizes (group={n_group}, rest={n_rest})"
                )
            scores: dict[str, float] = {}
            for drug in drugs:
                vals_group = [float(z_matrix[s].get(drug)) for s in members if drug in z_matrix[s]]
                vals_rest = [float(z_matrix[s].get(drug)) for s in rest if drug in z_matrix[s]]
                if not vals_group or not vals_rest:
                    continue
                scores[drug] = _mean(vals_group) - _mean(vals_rest)
            criterion = _mean([abs(v) for v in scores.values()]) if scores else 0.0
            programs.append(
                {
                    "program_id": f"group={group}__vs_rest",
                    "signed": True,
                    "contrast_method": "group_vs_rest",
                    "drug_scores": scores,
                    "n_samples": n_group,
                    "criterion": criterion,
                    "group": group,
                    "group_qc": {
                        "group": group,
                        "n_group": n_group,
                        "n_rest": n_rest,
                        "n_case": 0,
                        "n_control": 0,
                        "reason_if_skipped": "",
                    },
                }
            )
            group_qc_rows.append(
                {
                    "group": group,
                    "n_group": n_group,
                    "n_rest": n_rest,
                    "n_case": 0,
                    "n_control": 0,
                    "sets_emitted": 0,
                    "reason_if_skipped": "",
                }
            )
        return programs, warnings, group_qc_rows, groups_skipped

    if contrast_method == "case_control":
        if case_control_by_sample is None:
            raise ValueError("contrast_method=case_control requires --case_control_tsv")
        if case_control_within_group:
            if not group_index:
                raise ValueError(
                    "contrast_method=case_control with --case_control_within_group true requires groups"
                )
            for group in sorted(group_index):
                members = group_index[group]
                n_group = len(members)
                n_rest = len(samples) - n_group
                if n_group < min_n:
                    msg = f"group_skipped group={group} n={n_group} reason=min_group_size"
                    warnings.append(msg)
                    groups_skipped[group] = {"n": n_group, "reason": "min_group_size"}
                    group_qc_rows.append(
                        {
                            "group": group,
                            "n_group": n_group,
                            "n_rest": n_rest,
                            "n_case": 0,
                            "n_control": 0,
                            "sets_emitted": 0,
                            "reason_if_skipped": "min_group_size",
                        }
                    )
                    continue
                case_samples = [s for s in members if bool(case_control_by_sample.get(s, False))]
                control_samples = [s for s in members if s in case_control_by_sample and not case_control_by_sample[s]]
                if len(case_samples) < 1 or len(control_samples) < 1:
                    warnings.append(
                        f"group={group} skipped for case_control: case={len(case_samples)} control={len(control_samples)}"
                    )
                    groups_skipped[group] = {
                        "n": n_group,
                        "reason": "missing_case_or_control",
                    }
                    group_qc_rows.append(
                        {
                            "group": group,
                            "n_group": n_group,
                            "n_rest": n_rest,
                            "n_case": len(case_samples),
                            "n_control": len(control_samples),
                            "sets_emitted": 0,
                            "reason_if_skipped": "missing_case_or_control",
                        }
                    )
                    continue
                if len(case_samples) < 2 or len(control_samples) < 2:
                    warnings.append(
                        f"group={group} has small case/control sizes ({len(case_samples)}/{len(control_samples)})"
                    )
                scores: dict[str, float] = {}
                for drug in drugs:
                    vals_case = [float(z_matrix[s].get(drug)) for s in case_samples if drug in z_matrix[s]]
                    vals_ctrl = [float(z_matrix[s].get(drug)) for s in control_samples if drug in z_matrix[s]]
                    if not vals_case or not vals_ctrl:
                        continue
                    scores[drug] = _mean(vals_case) - _mean(vals_ctrl)
                criterion = _mean([abs(v) for v in scores.values()]) if scores else 0.0
                programs.append(
                    {
                    "program_id": f"group={group}__case_vs_control",
                    "signed": True,
                    "contrast_method": "case_control",
                    "drug_scores": scores,
                    "n_samples": len(case_samples) + len(control_samples),
                    "criterion": criterion,
                    "group": group,
                    "group_qc": {
                        "group": group,
                        "n_group": n_group,
                        "n_rest": n_rest,
                        "n_case": len(case_samples),
                        "n_control": len(control_samples),
                        "reason_if_skipped": "",
                    },
                }
            )
                group_qc_rows.append(
                    {
                        "group": group,
                        "n_group": n_group,
                        "n_rest": n_rest,
                        "n_case": len(case_samples),
                        "n_control": len(control_samples),
                        "sets_emitted": 0,
                        "reason_if_skipped": "",
                    }
                )
            return programs, warnings, group_qc_rows, groups_skipped

        case_samples = [s for s in samples if bool(case_control_by_sample.get(s, False))]
        control_samples = [s for s in samples if s in case_control_by_sample and not case_control_by_sample[s]]
        if len(case_samples) < 1 or len(control_samples) < 1:
            raise ValueError(
                "contrast_method=case_control requires at least one case and one control sample in --case_control_tsv"
            )
        if len(case_samples) < 2 or len(control_samples) < 2:
            warnings.append(
                f"global case_control has small sample counts (case={len(case_samples)}, control={len(control_samples)})"
            )
        scores = {}
        for drug in drugs:
            vals_case = [float(z_matrix[s].get(drug)) for s in case_samples if drug in z_matrix[s]]
            vals_ctrl = [float(z_matrix[s].get(drug)) for s in control_samples if drug in z_matrix[s]]
            if not vals_case or not vals_ctrl:
                continue
            scores[drug] = _mean(vals_case) - _mean(vals_ctrl)
        criterion = _mean([abs(v) for v in scores.values()]) if scores else 0.0
        programs.append(
            {
                "program_id": "case_vs_control",
                "signed": True,
                "contrast_method": "case_control",
                "drug_scores": scores,
                "n_samples": len(case_samples) + len(control_samples),
                "criterion": criterion,
                "group": "",
            }
        )
        return programs, warnings, group_qc_rows, groups_skipped

    raise ValueError(f"Unsupported contrast_method: {contrast_method}")


def _program_rows(
    *,
    scores: dict[str, float],
    selection_basis: dict[str, float],
    selected_gene_ids: list[str],
    weights: dict[str, float],
    gene_meta: dict[str, dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    selected_rows: list[dict[str, object]] = []
    for rank, gene_id in enumerate(selected_gene_ids, start=1):
        row: dict[str, object] = {
            "gene_id": gene_id,
            "score": float(scores.get(gene_id, 0.0)),
            "weight": float(weights.get(gene_id, 0.0)),
            "rank": rank,
        }
        symbol = str(gene_meta.get(gene_id, {}).get("gene_symbol", "")).strip()
        biotype = str(gene_meta.get(gene_id, {}).get("gene_biotype", "")).strip()
        if symbol:
            row["gene_symbol"] = symbol
        if biotype:
            row["gene_biotype"] = biotype
        selected_rows.append(row)

    full_gene_ids = sorted(
        [gene_id for gene_id, value in scores.items() if float(value) != 0.0],
        key=lambda gene_id: (-abs(float(scores[gene_id])), str(gene_id)),
    )
    full_rows: list[dict[str, object]] = []
    for rank, gene_id in enumerate(full_gene_ids, start=1):
        row = {
            "gene_id": gene_id,
            "score": float(scores.get(gene_id, 0.0)),
            "rank": rank,
        }
        symbol = str(gene_meta.get(gene_id, {}).get("gene_symbol", "")).strip()
        biotype = str(gene_meta.get(gene_id, {}).get("gene_biotype", "")).strip()
        if symbol:
            row["gene_symbol"] = symbol
        if biotype:
            row["gene_biotype"] = biotype
        full_rows.append(row)
    return selected_rows, full_rows


def run_drug_response_workflow(
    *,
    cfg: DrugResponseWorkflowConfig,
    response_records: list[ResponseRecord],
    drug_targets: dict[str, dict[str, float]],
    response_summary: dict[str, object],
    target_summary: dict[str, object],
    groups_by_sample: dict[str, str] | None,
    case_control_by_sample: dict[str, bool] | None,
    input_files: list[dict[str, str]],
) -> dict[str, object]:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    matrix_raw, aggregate_summary = aggregate_response_records(response_records)
    matrix_oriented = orient_matrix(matrix_raw, response_direction=cfg.response_direction)
    z_matrix, normalize_summary = transform_matrix_by_drug(
        matrix_oriented,
        response_transform=cfg.response_transform,
    )
    samples = sorted(z_matrix)
    drugs = sorted({drug for row in z_matrix.values() for drug in row})
    n_drugs = len(drugs)
    n_samples = len(samples)
    run_warnings: list[dict[str, object]] = []

    def _warn(code: str, message: str, **extra: object) -> None:
        print(f"warning: {message}", file=sys.stderr)
        payload: dict[str, object] = {"code": code, "message": message}
        payload.update(extra)
        run_warnings.append(payload)

    coverage_drugs = sum(1 for drug in drugs if drug in drug_targets and drug_targets[drug])
    coverage_fraction = float(coverage_drugs) / float(n_drugs) if n_drugs else 0.0
    if coverage_fraction < _LOW_TARGET_COVERAGE_WARN:
        _warn(
            "low_target_coverage",
            "low drug-to-target coverage: "
            f"{coverage_drugs}/{n_drugs} drugs have target mappings "
            f"(fraction={coverage_fraction:.3f}).",
            coverage_drugs=coverage_drugs,
            n_drugs=n_drugs,
            fraction=coverage_fraction,
        )

    ubiq_weights, ubiq_summary = compute_ubiquity_weights(
        z_matrix=z_matrix,
        drugs=drugs,
        method=cfg.response_ubiquity_penalty,
        tau=float(cfg.ubiquity_tau),
        epsilon=float(cfg.ubiquity_epsilon),
    )
    target_gene_weights, target_ubiquity_summary = compute_target_ubiquity_weights(
        drug_targets=drug_targets,
        method=cfg.target_ubiquity_penalty,
    )
    poly_weights, poly_summary = compute_polypharm_weights(
        drug_targets=drug_targets,
        enabled=bool(cfg.polypharm_downweight),
        t0=int(cfg.polypharm_t0),
    )
    combined_drug_weights = {
        drug: float(ubiq_weights.get(drug, 1.0)) * float(poly_weights.get(drug, 1.0))
        for drug in drugs
    }

    programs, program_warnings, group_qc_rows, groups_skipped = _build_programs(
        contrast_method=cfg.contrast_method,
        z_matrix=z_matrix,
        samples=samples,
        drugs=drugs,
        groups_by_sample=groups_by_sample,
        case_control_by_sample=case_control_by_sample,
        case_control_within_group=bool(cfg.case_control_within_group),
        min_group_size=int(cfg.min_group_size),
    )
    for warning in program_warnings:
        _warn("program_construction", warning)

    requested_program_count = len(programs)
    truncated = False
    if requested_program_count > int(cfg.max_programs):
        programs = sorted(
            programs,
            key=lambda p: (-float(p.get("criterion", 0.0)), str(p.get("program_id", ""))),
        )[: int(cfg.max_programs)]
        truncated = True
        _warn(
            "program_truncated",
            f"requested {requested_program_count} programs but max_programs={cfg.max_programs}; "
            "truncating to top programs by selection criterion. Increase --max_programs to keep more.",
            requested=requested_program_count,
            max_programs=int(cfg.max_programs),
        )

    gtf_symbol_to_biotype: dict[str, str] = {}
    if cfg.gtf:
        for gene in read_genes_from_gtf(cfg.gtf, gene_id_field=cfg.gtf_gene_id_field):
            symbol = str(gene.gene_symbol or "").strip().upper()
            if symbol and symbol not in gtf_symbol_to_biotype:
                gtf_symbol_to_biotype[symbol] = str(gene.gene_biotype or "")

    gmt_topk_list = parse_int_list_csv(str(cfg.gmt_topk_list))
    gmt_mass_list = parse_mass_list_csv(str(cfg.gmt_mass_list))
    gmt_biotype_allowlist = parse_str_list_csv(str(cfg.gmt_biotype_allowlist))

    manifest_rows: list[dict[str, object]] = []
    combined_gmt_sets: list[tuple[str, list[str]]] = []
    safe_ids_seen: set[str] = set()
    n_sparse_fallback = 0
    total_sets_emitted = 0
    set_size_values: list[int] = []
    n_pos_sets = 0
    n_neg_sets = 0
    broad_dominance_programs: list[dict[str, object]] = []

    for program in programs:
        program_id = str(program["program_id"])
        signed = bool(program.get("signed", False))
        program_scores = {str(d): float(v) for d, v in dict(program.get("drug_scores", {})).items()}
        gene_scores, scoring_info = score_genes(
            scoring_model=cfg.scoring_model,
            program_drug_scores=program_scores,
            drug_targets=drug_targets,
            drug_weights=combined_drug_weights,
            target_gene_weights=target_gene_weights,
            signed=signed,
            sparse_alpha=float(cfg.sparse_alpha),
        )
        warning = str(scoring_info.get("warning", "")).strip()
        if warning:
            _warn("scoring_warning", f"program={program_id}: {warning}", program_id=program_id)
            n_sparse_fallback += 1

        if signed:
            selection_basis = {gene_id: abs(float(score)) for gene_id, score in gene_scores.items() if float(score) != 0.0}
        else:
            selection_basis = {gene_id: float(score) for gene_id, score in gene_scores.items() if float(score) > 0.0}
        selected_gene_ids = _select_gene_ids(selection_basis, cfg)
        selected_gene_ids = sorted(
            [gene_id for gene_id in selected_gene_ids if gene_id in selection_basis],
            key=lambda gene_id: (-float(selection_basis.get(gene_id, 0.0)), str(gene_id)),
        )
        weights = _selected_weights(selection_basis, selected_gene_ids, cfg.normalize)

        gene_meta: dict[str, dict[str, object]] = {}
        for gene_id in set(gene_scores) | set(selected_gene_ids):
            symbol = str(gene_id).upper()
            gene_meta[gene_id] = {
                "gene_symbol": symbol,
                "gene_biotype": str(gtf_symbol_to_biotype.get(symbol, "")).strip(),
            }
        selected_rows, full_rows = _program_rows(
            scores=gene_scores,
            selection_basis=selection_basis,
            selected_gene_ids=selected_gene_ids,
            weights=weights,
            gene_meta=gene_meta,
        )
        top_n = min(20, len(selected_rows))
        if top_n > 0:
            top_rows = selected_rows[:top_n]
            broad_hits = 0
            for row in top_rows:
                symbol = str(row.get("gene_symbol", "") or row.get("gene_id", ""))
                if _is_broad_pharmacology_gene(symbol):
                    broad_hits += 1
            frac = float(broad_hits) / float(top_n)
            if frac >= 0.4:
                message = (
                    "top genes appear enriched for broad receptor-family targets "
                    f"(program={program_id}, broad_fraction_top{top_n}={frac:.2f}). "
                    "This can reflect library pharmacology; consider stronger target ubiquity "
                    "penalty, high-confidence target filtering, and balanced sampling."
                )
                _warn(
                    "broad_pharmacology_dominance",
                    message,
                    program_id=program_id,
                    top_n=top_n,
                    broad_hits=broad_hits,
                    broad_fraction=frac,
                )
                broad_dominance_programs.append(
                    {
                        "program_id": program_id,
                        "top_n": top_n,
                        "broad_hits": broad_hits,
                        "broad_fraction": frac,
                    }
                )

        safe_id = _safe_component(program_id, "program")
        dedup = safe_id
        suffix = 2
        while dedup in safe_ids_seen:
            dedup = f"{safe_id}_{suffix}"
            suffix += 1
        safe_ids_seen.add(dedup)
        program_dir = out_dir / f"program={dedup}"
        program_dir.mkdir(parents=True, exist_ok=True)
        _write_rows(program_dir / "geneset.tsv", selected_rows)
        if cfg.emit_full:
            _write_rows(program_dir / "geneset.full.tsv", full_rows)

        gmt_sets: list[tuple[str, list[str]]] = []
        gmt_plans: list[dict[str, object]] = []
        gmt_diagnostics: list[dict[str, object]] = []
        gmt_path = resolve_gmt_out_path(program_dir, cfg.gmt_out)
        split_signed = bool(cfg.gmt_split_signed) if cfg.gmt_split_signed is not None else signed
        if cfg.emit_gmt:
            base_name = (
                f"{cfg.converter_name}"
                f"__dataset={_safe_component(cfg.dataset_label, 'dataset')}"
                f"__program={_safe_component(program_id, 'program')}"
                f"__contrast_method={cfg.contrast_method}"
                f"__scoring_model={cfg.scoring_model}"
            )
            gmt_sets, gmt_plans = build_gmt_sets_from_rows(
                rows=full_rows,
                base_name=base_name,
                prefer_symbol=bool(cfg.gmt_prefer_symbol),
                min_genes=int(cfg.gmt_min_genes),
                max_genes=int(cfg.gmt_max_genes),
                topk_list=gmt_topk_list,
                mass_list=gmt_mass_list,
                split_signed=split_signed,
                require_symbol=bool(cfg.gmt_require_symbol),
                allowed_biotypes={x.lower() for x in gmt_biotype_allowlist} if gmt_biotype_allowlist else None,
                emit_small_gene_sets=bool(cfg.emit_small_gene_sets),
                diagnostics=gmt_diagnostics,
                context={
                    "converter": cfg.converter_name,
                    "program_id": program_id,
                    "contrast_method": cfg.contrast_method,
                },
            )
            if gmt_sets:
                write_gmt(gmt_sets, gmt_path, gmt_format=cfg.gmt_format)
                combined_gmt_sets.extend(gmt_sets)
                total_sets_emitted += len(gmt_sets)
                for set_name, genes in gmt_sets:
                    set_size_values.append(len(genes))
                    if "__pos__" in set_name:
                        n_pos_sets += 1
                    elif "__neg__" in set_name:
                        n_neg_sets += 1
        run_warnings.extend(_warn_gmt_diagnostics(gmt_diagnostics))

        output_files: list[dict[str, object]] = [{"path": str(program_dir / "geneset.tsv"), "role": "selected_program"}]
        if cfg.emit_full:
            output_files.append({"path": str(program_dir / "geneset.full.tsv"), "role": "full_scores"})
        if cfg.emit_gmt and gmt_sets:
            output_files.append({"path": str(gmt_path), "role": "gmt"})

        run_summary_payload = {
            "converter": cfg.converter_name,
            "dataset_label": cfg.dataset_label,
            "program_id": program_id,
            "contrast_method": cfg.contrast_method,
            "response_metric": cfg.response_metric,
            "response_direction": cfg.response_direction,
            "response_transform": cfg.response_transform,
            "scoring_model_requested": cfg.scoring_model,
            "scoring_model_used": scoring_info.get("model_used", cfg.scoring_model),
            "n_samples": n_samples,
            "n_cell_lines": n_samples,
            "n_drugs": n_drugs,
            "n_compounds": n_drugs,
            "n_response_rows": int(aggregate_summary.get("n_input_records", 0) or 0),
            "n_drugs_with_targets": coverage_drugs,
            "fraction_drugs_with_targets": coverage_fraction,
            "ubiquity_penalty": {
                "method": cfg.response_ubiquity_penalty,
                "tau": cfg.ubiquity_tau,
                "epsilon": cfg.ubiquity_epsilon,
            },
            "target_ubiquity_penalty": target_ubiquity_summary,
            "polypharm_downweight": {
                "enabled": cfg.polypharm_downweight,
                "t0": cfg.polypharm_t0,
            },
            "n_genes_selected": len(selected_rows),
            "requested_gmt_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
            "skipped_gmt_outputs": _collect_skipped_gmt_outputs(gmt_diagnostics),
            "warnings": run_warnings,
        }
        run_json_path, run_txt_path = write_run_summary_files(program_dir, run_summary_payload)
        output_files.append({"path": str(run_json_path), "role": "run_summary_json"})
        output_files.append({"path": str(run_txt_path), "role": "run_summary_text"})

        params = {
            "dataset_label": cfg.dataset_label,
            "response_metric": cfg.response_metric,
            "response_direction": cfg.response_direction,
            "response_transform": cfg.response_transform,
            "contrast_method": cfg.contrast_method,
            "case_control_within_group": cfg.case_control_within_group,
            "scoring_model": cfg.scoring_model,
            "sparse_alpha": cfg.sparse_alpha,
            "ubiquity_penalty": cfg.response_ubiquity_penalty,
            "response_ubiquity_penalty": cfg.response_ubiquity_penalty,
            "target_ubiquity_penalty": cfg.target_ubiquity_penalty,
            "ubiquity_tau": cfg.ubiquity_tau,
            "ubiquity_epsilon": cfg.ubiquity_epsilon,
            "polypharm_downweight": cfg.polypharm_downweight,
            "polypharm_t0": cfg.polypharm_t0,
            "min_group_size": cfg.min_group_size,
            "max_programs": cfg.max_programs,
            "program_id": program_id,
            "selection_method": cfg.select,
            "top_k": cfg.top_k,
            "quantile": cfg.quantile,
            "min_score": cfg.min_score,
            "normalize": cfg.normalize,
            "emit_full": cfg.emit_full,
            "emit_gmt": cfg.emit_gmt,
            "gmt_split_signed": split_signed,
            "gmt_format": cfg.gmt_format,
        }

        meta = make_metadata(
            converter_name=cfg.converter_name,
            parameters=params,
            data_type="drug_response",
            assay="drug_screen",
            organism=cfg.organism,
            genome_build=cfg.genome_build,
            files=input_files,
            gene_annotation=(
                {
                    "mode": "gtf",
                    "gtf_path": str(cfg.gtf),
                    "source": cfg.gtf_source or "user",
                    "gene_id_field": cfg.gtf_gene_id_field,
                }
                if cfg.gtf
                else {"mode": "none", "source": "none", "gene_id_field": "gene_symbol"}
            ),
            weights={
                "weight_type": ("signed" if signed else "nonnegative"),
                "normalization": {
                    "method": cfg.normalize,
                    "target_sum": 1.0 if cfg.normalize in {"within_set_l1", "l1"} else None,
                },
                "aggregation": str(scoring_info.get("model_used", cfg.scoring_model)),
            },
            summary={
                "n_input_features": int(n_drugs),
                "n_genes": int(len(selected_rows)),
                "n_features_assigned": int(coverage_drugs),
                "fraction_features_assigned": float(coverage_fraction),
                "n_samples": int(n_samples),
                "n_program_drugs": int(len(program_scores)),
                "criterion": float(program.get("criterion", 0.0) or 0.0),
                "target_summary": target_summary,
                "response_summary": response_summary,
                "aggregate_summary": aggregate_summary,
                "normalize_summary": {
                    "response_transform": normalize_summary.get("response_transform"),
                    "n_values_present": normalize_summary.get("n_values_present"),
                    "n_values_missing": normalize_summary.get("n_values_missing"),
                    "fraction_values_present": normalize_summary.get("fraction_values_present"),
                },
                "scoring_info": scoring_info,
                "target_ubiquity_summary": target_ubiquity_summary,
                "broad_dominance_programs": broad_dominance_programs,
                "groups_skipped": groups_skipped,
                "warnings": run_warnings,
                "program_truncation": {
                    "requested": requested_program_count,
                    "emitted": len(programs),
                    "truncated": truncated,
                },
            },
            program_extraction={
                "selection_method": cfg.select,
                "selection_params": {
                    "k": cfg.top_k,
                    "quantile": cfg.quantile,
                    "min_score": cfg.min_score,
                },
                "normalize": cfg.normalize,
                "n_selected_genes": len(selected_rows),
                "score_definition": (
                    "drug-response target aggregation with "
                    f"contrast_method={cfg.contrast_method}, scoring_model={cfg.scoring_model}"
                ),
            },
            output_files=output_files,
            gmt={
                "written": bool(cfg.emit_gmt and gmt_sets),
                "path": (
                    str(gmt_path.relative_to(program_dir))
                    if cfg.emit_gmt and gmt_sets and gmt_path.is_relative_to(program_dir)
                    else str(gmt_path)
                )
                if cfg.emit_gmt and gmt_sets
                else None,
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

        group_qc = dict(program.get("group_qc", {}))
        group_name = str(program.get("group", "")).strip()
        if group_name:
            for row in group_qc_rows:
                if str(row.get("group", "")) == group_name and not str(row.get("reason_if_skipped", "")):
                    row["sets_emitted"] = int(row.get("sets_emitted", 0) or 0) + int(len(gmt_sets))

        manifest_rows.append(
            {
                "program_id": program_id,
                "group": group_name,
                "contrast_method": cfg.contrast_method,
                "signed": str(signed).lower(),
                "n_samples": int(program.get("n_samples", 0) or 0),
                "n_group": int(group_qc.get("n_group", 0) or 0),
                "n_rest": int(group_qc.get("n_rest", 0) or 0),
                "n_case": int(group_qc.get("n_case", 0) or 0),
                "n_control": int(group_qc.get("n_control", 0) or 0),
                "criterion": float(program.get("criterion", 0.0) or 0.0),
                "sets_emitted": int(len(gmt_sets)),
                "path": str(program_dir.relative_to(out_dir)),
            }
        )

    manifest_path = out_dir / "manifest.tsv"
    if not manifest_rows:
        raise ValueError(
            "No programs were emitted after filtering/contrast setup. "
            "Check group sizes, contrast_method, target coverage, and min_group_size."
        )

    with manifest_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "program_id",
                "group",
                "contrast_method",
                "signed",
                "n_samples",
                "n_group",
                "n_rest",
                "n_case",
                "n_control",
                "criterion",
                "sets_emitted",
                "path",
            ],
        )
        writer.writeheader()
        for row in manifest_rows:
            writer.writerow(row)

    group_qc_path = out_dir / "group_qc.tsv"
    with group_qc_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=["group", "n_group", "n_rest", "n_case", "n_control", "sets_emitted", "reason_if_skipped"],
        )
        writer.writeheader()
        for row in sorted(group_qc_rows, key=lambda x: str(x.get("group", ""))):
            writer.writerow(
                {
                    "group": str(row.get("group", "")),
                    "n_group": int(row.get("n_group", 0) or 0),
                    "n_rest": int(row.get("n_rest", 0) or 0),
                    "n_case": int(row.get("n_case", 0) or 0),
                    "n_control": int(row.get("n_control", 0) or 0),
                    "sets_emitted": int(row.get("sets_emitted", 0) or 0),
                    "reason_if_skipped": str(row.get("reason_if_skipped", "")),
                }
            )

    if cfg.emit_gmt and combined_gmt_sets:
        write_gmt(combined_gmt_sets, out_dir / "genesets.gmt", gmt_format=cfg.gmt_format)

    target_stats = target_summary.get("targets_per_drug", {}) if isinstance(target_summary, dict) else {}
    promisc_stats = target_summary.get("promiscuity", {}) if isinstance(target_summary, dict) else {}
    groups_total = len(group_qc_rows)
    groups_emitted = sum(1 for row in group_qc_rows if not str(row.get("reason_if_skipped", "")).strip())
    groups_skipped_count = groups_total - groups_emitted
    emitted_rows = list(manifest_rows)
    group_rows_payload = [
        {
            "group": str(row.get("group", "")),
            "n_group": int(row.get("n_group", 0) or 0),
            "n_rest": int(row.get("n_rest", 0) or 0),
            "n_case": int(row.get("n_case", 0) or 0),
            "n_control": int(row.get("n_control", 0) or 0),
            "sets_emitted": int(row.get("sets_emitted", 0) or 0),
            "reason_if_skipped": str(row.get("reason_if_skipped", "")),
        }
        for row in sorted(group_qc_rows, key=lambda x: str(x.get("group", "")))
    ]
    set_size_summary = {
        "min": min(set_size_values) if set_size_values else 0,
        "median": statistics.median(set_size_values) if set_size_values else 0,
        "max": max(set_size_values) if set_size_values else 0,
    }
    root_summary_payload = {
        "converter": cfg.converter_name,
        "dataset_label": cfg.dataset_label,
        "response_metric": cfg.response_metric,
        "response_direction": cfg.response_direction,
        "response_transform": cfg.response_transform,
        "contrast_method": cfg.contrast_method,
        "scoring_model": cfg.scoring_model,
        "n_samples": n_samples,
        "n_cell_lines": n_samples,
        "n_drugs": n_drugs,
        "n_compounds": n_drugs,
        "n_response_rows": int(aggregate_summary.get("n_input_records", 0) or 0),
        "n_groups_total": groups_total,
        "n_groups_emitted": groups_emitted,
        "n_groups_skipped": groups_skipped_count,
        "groups": group_rows_payload,
        "target_table_stats": {
            "drugs_retained": coverage_drugs,
            "drugs_dropped_missing_targets": max(0, n_drugs - coverage_drugs),
            "drugs_dropped_promisc": int(promisc_stats.get("n_dropped", 0) or 0),
            "drugs_capped_promisc": int(promisc_stats.get("n_capped", 0) or 0),
            "p90_targets_per_drug": target_stats.get("p90", 0),
        },
        "target_summary": target_summary,
        "broad_dominance_programs": broad_dominance_programs,
        "gene_set_emission": {
            "total_sets_emitted": total_sets_emitted,
            "size_summary": set_size_summary,
            "n_pos_sets": n_pos_sets,
            "n_neg_sets": n_neg_sets,
        },
        "warnings": run_warnings,
    }
    run_json_path, run_txt_path = write_run_summary_files(out_dir, root_summary_payload)
    print(
        "group\tn_group\tn_rest\tsets_emitted\treason_if_skipped",
        file=sys.stdout,
    )
    for row in group_rows_payload:
        print(
            f"{row['group']}\t{row['n_group']}\t{row['n_rest']}\t{row['sets_emitted']}\t{row['reason_if_skipped']}",
            file=sys.stdout,
        )

    return {
        "n_groups": groups_emitted,
        "out_dir": str(out_dir),
        "n_sparse_fallback": n_sparse_fallback,
        "n_programs_requested": requested_program_count,
        "n_programs_emitted": len(emitted_rows),
        "programs_truncated": bool(truncated),
        "n_groups_skipped": groups_skipped_count,
        "run_summary_json": str(run_json_path),
        "run_summary_txt": str(run_txt_path),
    }
