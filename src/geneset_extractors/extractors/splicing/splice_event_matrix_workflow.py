from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path
import sys

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import enrich_manifest_row
from geneset_extractors.core.qc import write_run_summary_files
from geneset_extractors.extractors.rnaseq.deg_scoring import sanitize_name_component
from geneset_extractors.extractors.splicing.splice_event_diff_workflow import (
    SpliceRow,
    SplicingWorkflowConfig,
    read_tsv_rows,
    run_splice_event_diff_workflow,
)


@dataclass(frozen=True)
class SampleMetaRecord:
    line_no: int
    values: dict[str, str]


@dataclass(frozen=True)
class ContrastSpec:
    contrast_id: str
    contrast_label: str
    study_contrast: str
    sample_ids_a: tuple[str, ...]
    sample_ids_b: tuple[str, ...]
    condition_a: str | None
    condition_b: str | None
    group_label: str | None


@dataclass
class SpliceMatrixWorkflowConfig:
    converter_name: str
    out_dir: Path
    organism: str
    genome_build: str
    dataset_label: str
    signature_name: str
    psi_matrix_tsv: str
    sample_metadata_tsv: str
    event_metadata_tsv: str | None
    coverage_matrix_tsv: str | None
    matrix_format: str
    sample_id_column: str
    group_column: str | None
    condition_column: str | None
    study_contrast: str
    condition_a: str | None
    condition_b: str | None
    min_samples_per_condition: int
    effect_metric: str
    missing_value_policy: str
    min_present_per_condition: int
    tool_family: str
    event_id_column: str | None
    event_group_column: str | None
    event_type_column: str | None
    gene_id_column: str | None
    gene_symbol_column: str | None
    chrom_column: str | None
    start_column: str | None
    end_column: str | None
    strand_column: str | None
    score_column: str | None
    stat_column: str | None
    delta_psi_column: str | None
    psi_column: str | None
    padj_column: str | None
    pvalue_column: str | None
    probability_column: str | None
    read_support_column: str | None
    novel_flag_column: str | None
    annotation_status_column: str | None
    score_mode: str
    score_transform: str
    confidence_weight_mode: str
    min_probability: float
    min_read_support: float
    neglog10p_cap: float
    neglog10p_eps: float
    delta_psi_soft_floor: float
    delta_psi_soft_floor_mode: str
    event_dup_policy: str
    gene_aggregation: str
    gene_topk_events: int
    gene_burden_penalty_mode: str
    min_gene_burden_penalty: float
    gene_support_penalty_mode: str
    locus_density_penalty_mode: str
    locus_density_window_bp: int
    locus_density_top_n: int
    locus_density_penalty_mode_explicit: bool
    source_dataset: str | None
    bundle_same_dataset_policy: str
    bundle_prior_profile: str
    bundle_prior_profile_explicit: bool
    artifact_action: str
    artifact_action_explicit: bool
    min_nonself_source_datasets_for_event_prior: int
    max_dataset_fraction_for_event_prior: float
    ambiguous_gene_policy: str
    impact_mode: str
    impact_min: float
    impact_max: float
    select: str
    top_k: int
    quantile: float
    min_score: float
    normalize: str
    emit_full: bool
    emit_gmt: bool
    gmt_out: str | None
    gmt_format: str
    gmt_prefer_symbol: bool
    gmt_require_symbol: bool
    gmt_biotype_allowlist: str
    gmt_min_genes: int
    gmt_max_genes: int
    gmt_topk_list: str
    gmt_mass_list: str
    gmt_split_signed: bool
    emit_small_gene_sets: bool



def _clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()



def _parse_float(value: object) -> float | None:
    text = _clean(value)
    if not text or text.lower() in {"na", "nan", "none", "null"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None



def read_sample_metadata_tsv(path: str | Path) -> tuple[list[str], list[SampleMetaRecord]]:
    rows: list[SampleMetaRecord] = []
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Sample metadata TSV has no header: {path}")
        fieldnames = [str(x) for x in reader.fieldnames]
        for idx, row in enumerate(reader, start=2):
            rows.append(SampleMetaRecord(line_no=idx, values={k: str(v) for k, v in row.items() if k is not None}))
    return fieldnames, rows



def _sample_rows_by_id(rows: list[SampleMetaRecord], sample_id_column: str) -> dict[str, SampleMetaRecord]:
    out: dict[str, SampleMetaRecord] = {}
    for row in rows:
        sample_id = _clean(row.values.get(sample_id_column))
        if not sample_id:
            continue
        out[sample_id] = row
    if not out:
        raise ValueError(f"No non-empty sample ids found in sample metadata column '{sample_id_column}'.")
    return out


def _infer_source_dataset_from_sample_rows(rows: dict[str, SampleMetaRecord]) -> str | None:
    for candidate in ("study_id", "source_dataset"):
        values = sorted({_clean(row.values.get(candidate)) for row in rows.values() if _clean(row.values.get(candidate))})
        if len(values) == 1:
            return values[0]
    return None



def _mean(values: list[float]) -> float:
    return sum(values) / float(len(values))



def _sample_variance(values: list[float], mean_value: float) -> float:
    if len(values) < 2:
        return 0.0
    return sum((x - mean_value) ** 2 for x in values) / float(len(values) - 1)



def _welch_t(values_a: list[float], values_b: list[float]) -> float:
    mean_a = _mean(values_a)
    mean_b = _mean(values_b)
    var_a = _sample_variance(values_a, mean_a)
    var_b = _sample_variance(values_b, mean_b)
    denom = math.sqrt((var_a / float(len(values_a))) + (var_b / float(len(values_b))))
    if denom <= 0.0:
        return 0.0
    return (mean_a - mean_b) / denom



def _approx_two_sided_p_from_z(value: float) -> float:
    return math.erfc(abs(float(value)) / math.sqrt(2.0))


def _bh_adjust_pvalues(values: list[float | None]) -> list[float | None]:
    indexed = [(idx, float(val)) for idx, val in enumerate(values) if val is not None]
    if not indexed:
        return [None for _ in values]
    n = len(indexed)
    adjusted: list[float | None] = [None for _ in values]
    running = 1.0
    for rank, (idx, pvalue) in enumerate(sorted(indexed, key=lambda item: item[1]), start=1):
        candidate = min(1.0, pvalue * float(n) / float(rank))
        adjusted[idx] = candidate
    for idx, _pvalue in reversed(sorted(indexed, key=lambda item: item[1])):
        running = min(running, float(adjusted[idx] or 1.0))
        adjusted[idx] = running
    return adjusted



def _present_values(row: SpliceRow, sample_ids: tuple[str, ...], missing_value_policy: str) -> list[float] | None:
    parsed: list[float] = []
    missing = 0
    for sample_id in sample_ids:
        value = _parse_float(row.values.get(sample_id))
        if value is None:
            missing += 1
            continue
        parsed.append(float(value))
    if missing_value_policy == "drop" and missing > 0:
        return None
    return parsed



def _resolve_condition_levels(rows: dict[str, SampleMetaRecord], condition_column: str, condition_a: str | None, condition_b: str | None) -> tuple[str, str]:
    if condition_a and condition_b:
        return str(condition_a), str(condition_b)
    levels = sorted({_clean(r.values.get(condition_column)) for r in rows.values() if _clean(r.values.get(condition_column))})
    if len(levels) != 2:
        raise ValueError(
            f"Could not infer exactly two condition levels from column '{condition_column}'. Observed levels: {', '.join(levels) if levels else '(none)'}"
        )
    return levels[0], levels[1]



def _build_contrast_specs(cfg: SpliceMatrixWorkflowConfig, sample_rows: dict[str, SampleMetaRecord]) -> list[ContrastSpec]:
    sample_ids = tuple(sorted(sample_rows))
    specs: list[ContrastSpec] = []
    if cfg.study_contrast == "baseline":
        return [ContrastSpec("baseline", "baseline_all_samples", "baseline", sample_ids, (), None, None, None)]
    if cfg.study_contrast == "condition_a_vs_b":
        if not cfg.condition_column:
            raise ValueError("study_contrast=condition_a_vs_b requires --condition_column")
        condition_a, condition_b = _resolve_condition_levels(sample_rows, cfg.condition_column, cfg.condition_a, cfg.condition_b)
        sample_ids_a = tuple(sorted(sid for sid, row in sample_rows.items() if _clean(row.values.get(cfg.condition_column or "")) == condition_a))
        sample_ids_b = tuple(sorted(sid for sid, row in sample_rows.items() if _clean(row.values.get(cfg.condition_column or "")) == condition_b))
        return [ContrastSpec(f"{sanitize_name_component(condition_a)}_vs_{sanitize_name_component(condition_b)}", f"{condition_a}_vs_{condition_b}", "condition_a_vs_b", sample_ids_a, sample_ids_b, condition_a, condition_b, None)]
    if cfg.study_contrast == "group_vs_rest":
        if not cfg.group_column:
            raise ValueError("study_contrast=group_vs_rest requires --group_column")
        groups = sorted({_clean(r.values.get(cfg.group_column)) for r in sample_rows.values() if _clean(r.values.get(cfg.group_column))})
        for group in groups:
            sample_ids_a = tuple(sorted(sid for sid, row in sample_rows.items() if _clean(row.values.get(cfg.group_column or "")) == group))
            sample_ids_b = tuple(sorted(sid for sid, row in sample_rows.items() if _clean(row.values.get(cfg.group_column or "")) != group))
            specs.append(ContrastSpec(f"group_{sanitize_name_component(group)}_vs_rest", f"{group}_vs_rest", "group_vs_rest", sample_ids_a, sample_ids_b, group, "rest", group))
        return specs
    if cfg.study_contrast == "condition_within_group":
        if not cfg.group_column or not cfg.condition_column:
            raise ValueError("study_contrast=condition_within_group requires --group_column and --condition_column")
        condition_a, condition_b = _resolve_condition_levels(sample_rows, cfg.condition_column, cfg.condition_a, cfg.condition_b)
        groups = sorted({_clean(r.values.get(cfg.group_column)) for r in sample_rows.values() if _clean(r.values.get(cfg.group_column))})
        for group in groups:
            in_group = {sid: row for sid, row in sample_rows.items() if _clean(row.values.get(cfg.group_column or "")) == group}
            sample_ids_a = tuple(sorted(sid for sid, row in in_group.items() if _clean(row.values.get(cfg.condition_column or "")) == condition_a))
            sample_ids_b = tuple(sorted(sid for sid, row in in_group.items() if _clean(row.values.get(cfg.condition_column or "")) == condition_b))
            specs.append(ContrastSpec(f"group_{sanitize_name_component(group)}__{sanitize_name_component(condition_a)}_vs_{sanitize_name_component(condition_b)}", f"{group}:{condition_a}_vs_{condition_b}", "condition_within_group", sample_ids_a, sample_ids_b, condition_a, condition_b, group))
        return specs
    raise ValueError(f"Unsupported study_contrast: {cfg.study_contrast}")



def _merge_event_metadata(matrix_rows: list[SpliceRow], event_metadata_tsv: str | None, event_id_column: str | None) -> list[SpliceRow]:
    if not event_metadata_tsv:
        return matrix_rows
    with Path(event_metadata_tsv).open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"event_metadata_tsv has no header: {event_metadata_tsv}")
        id_col = event_id_column if event_id_column in reader.fieldnames else ("event_id" if "event_id" in reader.fieldnames else None)
        if id_col is None:
            raise ValueError("event_metadata_tsv must contain the configured event_id column or an event_id column")
        meta_map = {str(row.get(id_col, "")).strip(): {k: str(v) for k, v in row.items() if k is not None} for row in reader}
    merged: list[SpliceRow] = []
    id_col = event_id_column or "event_id"
    for row in matrix_rows:
        event_id = _clean(row.values.get(id_col))
        values = dict(row.values)
        if event_id and event_id in meta_map:
            for key, value in meta_map[event_id].items():
                values.setdefault(key, value)
        merged.append(SpliceRow(line_no=row.line_no, values=values))
    return merged



def _coverage_map(coverage_matrix_tsv: str | None) -> dict[tuple[int, str], float]:
    if not coverage_matrix_tsv:
        return {}
    _fieldnames, rows = read_tsv_rows(coverage_matrix_tsv)
    out: dict[tuple[int, str], float] = {}
    for row in rows:
        for key, value in row.values.items():
            parsed = _parse_float(value)
            if parsed is not None:
                out[(row.line_no, key)] = float(parsed)
    return out



def _row_statistic(row: SpliceRow, *, sample_ids_a: tuple[str, ...], sample_ids_b: tuple[str, ...], effect_metric: str, missing_value_policy: str, min_samples_per_condition: int, min_present_per_condition: int) -> dict[str, float] | None:
    values_a = _present_values(row, sample_ids_a, missing_value_policy)
    values_b = _present_values(row, sample_ids_b, missing_value_policy)
    if values_a is None or values_b is None:
        return None
    if len(values_a) < int(min_samples_per_condition) or len(values_b) < int(min_samples_per_condition):
        return None
    if len(values_a) < int(min_present_per_condition) or len(values_b) < int(min_present_per_condition):
        return None
    mean_a = _mean(values_a)
    mean_b = _mean(values_b)
    mean_diff = float(mean_a - mean_b)
    if effect_metric == "mean_diff":
        stat_value = mean_diff
    elif effect_metric == "welch_t":
        stat_value = _welch_t(values_a, values_b)
    else:
        raise ValueError(f"Unsupported effect_metric: {effect_metric}")
    return {
        "mean_a": float(mean_a),
        "mean_b": float(mean_b),
        "delta_psi": float(mean_diff),
        "stat": float(stat_value),
        "pvalue": float(_approx_two_sided_p_from_z(stat_value)),
        "n_a": float(len(values_a)),
        "n_b": float(len(values_b)),
    }



def _baseline_statistic(row: SpliceRow, *, sample_ids: tuple[str, ...], missing_value_policy: str, min_present_per_condition: int) -> dict[str, float] | None:
    values = _present_values(row, sample_ids, missing_value_policy)
    if values is None:
        return None
    if len(values) < int(min_present_per_condition):
        return None
    mean_value = _mean(values)
    return {
        "mean_a": float(mean_value),
        "mean_b": 0.0,
        "delta_psi": float(mean_value),
        "stat": float(mean_value),
        "pvalue": 1.0,
        "n_a": float(len(values)),
        "n_b": 0.0,
    }



def _contrast_rows_for_events(*, event_rows: list[SpliceRow], contrast: ContrastSpec, cfg: SpliceMatrixWorkflowConfig, coverage_lookup: dict[tuple[int, str], float]) -> tuple[list[SpliceRow], dict[str, int]]:
    contrast_rows: list[SpliceRow] = []
    raw_pvalues: list[float | None] = []
    summary = {"n_input_events": len(event_rows), "n_events_retained": 0, "n_events_dropped_missing": 0, "n_events_with_support": 0, "n_events_with_bh_padj": 0}
    for row in event_rows:
        if contrast.study_contrast == "baseline":
            stat = _baseline_statistic(row, sample_ids=contrast.sample_ids_a, missing_value_policy=cfg.missing_value_policy, min_present_per_condition=cfg.min_present_per_condition)
        else:
            stat = _row_statistic(row, sample_ids_a=contrast.sample_ids_a, sample_ids_b=contrast.sample_ids_b, effect_metric=cfg.effect_metric, missing_value_policy=cfg.missing_value_policy, min_samples_per_condition=cfg.min_samples_per_condition, min_present_per_condition=cfg.min_present_per_condition)
        if stat is None:
            summary["n_events_dropped_missing"] += 1
            continue
        values = {k: str(v) for k, v in row.values.items()}
        values["delta_psi"] = str(stat["delta_psi"])
        values["stat"] = str(stat["stat"])
        values["pvalue"] = str(stat["pvalue"])
        values["padj"] = ""
        raw_pvalues.append(float(stat["pvalue"]))
        support_vals: list[float] = []
        for sample_id in contrast.sample_ids_a + contrast.sample_ids_b:
            val = coverage_lookup.get((row.line_no, sample_id))
            if val is not None:
                support_vals.append(float(val))
        if support_vals:
            values["read_support"] = str(sum(support_vals) / float(len(support_vals)))
            summary["n_events_with_support"] += 1
        contrast_rows.append(SpliceRow(line_no=row.line_no, values=values))
        summary["n_events_retained"] += 1
    for row, padj in zip(contrast_rows, _bh_adjust_pvalues(raw_pvalues), strict=False):
        if padj is None:
            continue
        row.values["padj"] = str(padj)
        summary["n_events_with_bh_padj"] += 1
    return contrast_rows, summary



def _make_child_cfg(cfg: SpliceMatrixWorkflowConfig, contrast: ContrastSpec, out_dir: Path) -> SplicingWorkflowConfig:
    signature_bits = [cfg.signature_name, f"study_contrast={contrast.study_contrast}", f"contrast={contrast.contrast_id}"]
    if contrast.group_label:
        signature_bits.append(f"group={contrast.group_label}")
    signature_name = "__".join(signature_bits)
    dataset_label = f"{cfg.dataset_label}::{contrast.contrast_label}"
    return SplicingWorkflowConfig(
        converter_name=cfg.converter_name,
        out_dir=out_dir,
        organism=cfg.organism,
        genome_build=cfg.genome_build,
        dataset_label=dataset_label,
        signature_name=signature_name,
        tool_family=cfg.tool_family,
        event_id_column=cfg.event_id_column,
        event_group_column=cfg.event_group_column,
        event_type_column=cfg.event_type_column,
        gene_id_column=cfg.gene_id_column,
        gene_symbol_column=cfg.gene_symbol_column,
        chrom_column=cfg.chrom_column,
        start_column=cfg.start_column,
        end_column=cfg.end_column,
        strand_column=cfg.strand_column,
        score_column=cfg.score_column,
        stat_column=cfg.stat_column,
        delta_psi_column=cfg.delta_psi_column,
        psi_column=cfg.psi_column,
        padj_column=cfg.padj_column,
        pvalue_column=cfg.pvalue_column,
        probability_column=cfg.probability_column,
        read_support_column=cfg.read_support_column,
        novel_flag_column=cfg.novel_flag_column,
        annotation_status_column=cfg.annotation_status_column,
        score_mode=cfg.score_mode,
        score_transform=cfg.score_transform,
        confidence_weight_mode=cfg.confidence_weight_mode,
        min_probability=cfg.min_probability,
        min_read_support=cfg.min_read_support,
        neglog10p_cap=cfg.neglog10p_cap,
        neglog10p_eps=cfg.neglog10p_eps,
        delta_psi_soft_floor=cfg.delta_psi_soft_floor,
        delta_psi_soft_floor_mode=cfg.delta_psi_soft_floor_mode,
        event_dup_policy=cfg.event_dup_policy,
        gene_aggregation=cfg.gene_aggregation,
        gene_topk_events=cfg.gene_topk_events,
        gene_burden_penalty_mode=cfg.gene_burden_penalty_mode,
        min_gene_burden_penalty=cfg.min_gene_burden_penalty,
        gene_support_penalty_mode=cfg.gene_support_penalty_mode,
        locus_density_penalty_mode=cfg.locus_density_penalty_mode,
        locus_density_window_bp=cfg.locus_density_window_bp,
        locus_density_top_n=cfg.locus_density_top_n,
        locus_density_penalty_mode_explicit=cfg.locus_density_penalty_mode_explicit,
        source_dataset=cfg.source_dataset,
        bundle_same_dataset_policy=cfg.bundle_same_dataset_policy,
        bundle_prior_profile=cfg.bundle_prior_profile,
        bundle_prior_profile_explicit=cfg.bundle_prior_profile_explicit,
        artifact_action=cfg.artifact_action,
        artifact_action_explicit=cfg.artifact_action_explicit,
        min_nonself_source_datasets_for_event_prior=cfg.min_nonself_source_datasets_for_event_prior,
        max_dataset_fraction_for_event_prior=cfg.max_dataset_fraction_for_event_prior,
        ambiguous_gene_policy=cfg.ambiguous_gene_policy,
        impact_mode=cfg.impact_mode,
        impact_min=cfg.impact_min,
        impact_max=cfg.impact_max,
        select=cfg.select,
        top_k=cfg.top_k,
        quantile=cfg.quantile,
        min_score=cfg.min_score,
        normalize=cfg.normalize,
        emit_full=cfg.emit_full,
        emit_gmt=cfg.emit_gmt,
        gmt_out=cfg.gmt_out,
        gmt_format=cfg.gmt_format,
        gmt_prefer_symbol=cfg.gmt_prefer_symbol,
        gmt_require_symbol=cfg.gmt_require_symbol,
        gmt_biotype_allowlist=cfg.gmt_biotype_allowlist,
        gmt_min_genes=cfg.gmt_min_genes,
        gmt_max_genes=cfg.gmt_max_genes,
        gmt_topk_list=cfg.gmt_topk_list,
        gmt_mass_list=cfg.gmt_mass_list,
        gmt_split_signed=cfg.gmt_split_signed,
        emit_small_gene_sets=cfg.emit_small_gene_sets,
    )



def run_splice_event_matrix_workflow(
    *,
    cfg: SpliceMatrixWorkflowConfig,
    alias_map: dict[str, object] | None,
    ubiquity_map: dict[str, object] | None,
    event_ubiquity_by_dataset_map: dict[str, object] | None,
    impact_map: dict[str, object] | None,
    gene_burden_map: dict[str, object] | None,
    gene_burden_by_dataset_map: dict[str, object] | None,
    gene_locus_map: dict[str, object] | None,
    input_files: list[dict[str, str]],
    resources_info: dict[str, object] | None,
) -> dict[str, object]:
    if cfg.matrix_format != "wide_events_by_sample":
        raise ValueError(f"Unsupported matrix_format: {cfg.matrix_format}")
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    matrix_fieldnames, matrix_rows = read_tsv_rows(cfg.psi_matrix_tsv)
    matrix_rows = _merge_event_metadata(matrix_rows, cfg.event_metadata_tsv, cfg.event_id_column)
    _meta_fieldnames, meta_rows = read_sample_metadata_tsv(cfg.sample_metadata_tsv)
    sample_rows = _sample_rows_by_id(meta_rows, cfg.sample_id_column)
    available_sample_ids = [sid for sid in sample_rows if sid in matrix_fieldnames]
    if not available_sample_ids:
        raise ValueError("No sample IDs from sample_metadata_tsv were found as columns in psi_matrix_tsv.")
    sample_rows = {sid: sample_rows[sid] for sid in available_sample_ids}
    if not cfg.source_dataset:
        cfg.source_dataset = _infer_source_dataset_from_sample_rows(sample_rows)

    coverage_lookup = _coverage_map(cfg.coverage_matrix_tsv)
    contrasts = _build_contrast_specs(cfg, sample_rows)
    if not contrasts:
        raise ValueError("No splicing matrix contrasts could be constructed from the supplied metadata.")

    files = list(input_files)
    manifest_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    combined_gmt_sets: list[tuple[str, list[str]]] = []
    multiple = len(contrasts) > 1
    emitted = 0
    effective_bundle_profiles: set[str] = set()
    effective_artifact_actions: set[str] = set()
    effective_locus_modes: set[str] = set()
    quality_tiers: set[str] = set()
    n_tcga_public_mode_contrasts = 0
    n_outputs_withheld_for_artifact = 0

    for contrast in contrasts:
        contrast_rows, contrast_summary = _contrast_rows_for_events(event_rows=matrix_rows, contrast=contrast, cfg=cfg, coverage_lookup=coverage_lookup)
        if not contrast_rows:
            qc_rows.append({
                "contrast_id": contrast.contrast_id,
                "contrast_label": contrast.contrast_label,
                "study_contrast": contrast.study_contrast,
                "group_label": contrast.group_label or "",
                "condition_a": contrast.condition_a or "",
                "condition_b": contrast.condition_b or "",
                "n_samples_a": len(contrast.sample_ids_a),
                "n_samples_b": len(contrast.sample_ids_b),
                "n_input_events": contrast_summary["n_input_events"],
                "n_events_retained": contrast_summary["n_events_retained"],
                "n_events_dropped_missing": contrast_summary["n_events_dropped_missing"],
                "n_events_with_support": contrast_summary["n_events_with_support"],
                "n_events_with_bh_padj": contrast_summary["n_events_with_bh_padj"],
                "tcga_public_mode": False,
                "effective_bundle_prior_profile": "",
                "effective_artifact_action": "",
                "effective_locus_density_penalty_mode": "",
                "quality_tier": "skipped",
                "artifact_action_applied": "",
                "output_withheld_reason": "",
                "status": "skipped",
                "reason_if_skipped": "No events passed matrix contrast QC; likely too many missing values or too few samples per condition.",
                "path": "",
            })
            print(
                f"warning: skipped splicing contrast {contrast.contrast_label}; no events passed matrix contrast QC. Try lowering --min_samples_per_condition or --min_present_per_condition.",
                file=sys.stderr,
            )
            continue

        child_out_dir = out_dir if not multiple else out_dir / f"contrast={sanitize_name_component(contrast.contrast_id)}"
        child_cfg = _make_child_cfg(cfg, contrast, child_out_dir)
        _result = run_splice_event_diff_workflow(
            cfg=child_cfg,
            fieldnames=list(matrix_fieldnames) + ["delta_psi", "stat", "pvalue", "padj", "read_support"],
            rows=contrast_rows,
            alias_map=alias_map,
            ubiquity_map=ubiquity_map,
            event_ubiquity_by_dataset_map=event_ubiquity_by_dataset_map,
            impact_map=impact_map,
            gene_burden_map=gene_burden_map,
            gene_burden_by_dataset_map=gene_burden_by_dataset_map,
            gene_locus_map=gene_locus_map,
            input_files=files,
            resources_info=resources_info,
        )
        emitted += 1
        qc_rows.append({
            "contrast_id": contrast.contrast_id,
            "contrast_label": contrast.contrast_label,
            "study_contrast": contrast.study_contrast,
            "group_label": contrast.group_label or "",
            "condition_a": contrast.condition_a or "",
            "condition_b": contrast.condition_b or "",
            "n_samples_a": len(contrast.sample_ids_a),
            "n_samples_b": len(contrast.sample_ids_b),
            "n_input_events": contrast_summary["n_input_events"],
            "n_events_retained": contrast_summary["n_events_retained"],
            "n_events_dropped_missing": contrast_summary["n_events_dropped_missing"],
            "n_events_with_support": contrast_summary["n_events_with_support"],
            "n_events_with_bh_padj": contrast_summary["n_events_with_bh_padj"],
            "tcga_public_mode": bool(_result.get("tcga_public_mode", False)),
            "effective_bundle_prior_profile": str(_result.get("effective_bundle_prior_profile", "")),
            "effective_artifact_action": str(_result.get("effective_artifact_action", "")),
            "effective_locus_density_penalty_mode": str(_result.get("effective_locus_density_penalty_mode", "")),
            "quality_tier": str(_result.get("quality_tier", "")),
            "artifact_action_applied": str(_result.get("artifact_action_applied", "")),
            "output_withheld_reason": str(_result.get("output_withheld_reason", "")),
            "status": "emitted",
            "reason_if_skipped": "",
            "path": "." if not multiple else str(child_out_dir.relative_to(out_dir)),
        })
        if bool(_result.get("tcga_public_mode", False)):
            n_tcga_public_mode_contrasts += 1
        effective_bundle_profiles.add(str(_result.get("effective_bundle_prior_profile", "")))
        effective_artifact_actions.add(str(_result.get("effective_artifact_action", "")))
        effective_locus_modes.add(str(_result.get("effective_locus_density_penalty_mode", "")))
        quality_tiers.add(str(_result.get("quality_tier", "")))
        if str(_result.get("output_withheld_reason", "")):
            n_outputs_withheld_for_artifact += 1
        if multiple:
            manifest_rows.append({
                "contrast_id": contrast.contrast_id,
                "study_contrast": contrast.study_contrast,
                "group_label": contrast.group_label or "",
                "condition_a": contrast.condition_a or "",
                "condition_b": contrast.condition_b or "",
                "path": str(child_out_dir.relative_to(out_dir)),
            })
        child_gmt_path = child_out_dir / "genesets.gmt"
        if cfg.emit_gmt and child_gmt_path.exists():
            with child_gmt_path.open("r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        combined_gmt_sets.append((parts[0], parts[1:]))

    if multiple:
        with (out_dir / "manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
            rows = [enrich_manifest_row(out_dir, out_dir / str(row["path"]), row) for row in manifest_rows]
            writer = csv.DictWriter(
                fh,
                delimiter="\t",
                fieldnames=[
                    "contrast_id",
                    "study_contrast",
                    "group_label",
                    "condition_a",
                    "condition_b",
                    "geneset_id",
                    "label",
                    "path",
                    "meta_path",
                    "provenance_path",
                    "focus_node_id",
                ],
            )
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

    with (out_dir / "contrast_qc.tsv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "contrast_id",
                "contrast_label",
                "study_contrast",
                "group_label",
                "condition_a",
                "condition_b",
                "n_samples_a",
                "n_samples_b",
                "n_input_events",
                "n_events_retained",
                "n_events_dropped_missing",
                "n_events_with_support",
                "n_events_with_bh_padj",
                "tcga_public_mode",
                "effective_bundle_prior_profile",
                "effective_artifact_action",
                "effective_locus_density_penalty_mode",
                "quality_tier",
                "artifact_action_applied",
                "output_withheld_reason",
                "status",
                "reason_if_skipped",
                "path",
            ],
        )
        writer.writeheader()
        for row in qc_rows:
            writer.writerow(row)

    if multiple and cfg.emit_gmt and combined_gmt_sets:
        write_gmt(combined_gmt_sets, out_dir / "genesets.gmt", gmt_format=cfg.gmt_format)

    summary_payload = {
        "converter": cfg.converter_name,
        "dataset_label": cfg.dataset_label,
        "signature_name": cfg.signature_name,
        "study_contrast": cfg.study_contrast,
        "matrix_format": cfg.matrix_format,
        "effect_metric": cfg.effect_metric,
        "missing_value_policy": cfg.missing_value_policy,
        "delta_psi_soft_floor_mode": cfg.delta_psi_soft_floor_mode,
        "delta_psi_soft_floor": cfg.delta_psi_soft_floor,
        "gene_burden_penalty_mode": cfg.gene_burden_penalty_mode,
        "gene_support_penalty_mode": cfg.gene_support_penalty_mode,
        "bundle_prior_profile": cfg.bundle_prior_profile,
        "artifact_action": cfg.artifact_action,
        "effective_bundle_prior_profiles": sorted(value for value in effective_bundle_profiles if value),
        "effective_artifact_actions": sorted(value for value in effective_artifact_actions if value),
        "effective_locus_density_penalty_modes": sorted(value for value in effective_locus_modes if value),
        "quality_tiers": sorted(value for value in quality_tiers if value),
        "n_tcga_public_mode_contrasts": n_tcga_public_mode_contrasts,
        "n_outputs_withheld_for_artifact": n_outputs_withheld_for_artifact,
        "source_dataset": cfg.source_dataset,
        "bundle_same_dataset_policy": cfg.bundle_same_dataset_policy,
        "min_samples_per_condition": cfg.min_samples_per_condition,
        "min_present_per_condition": cfg.min_present_per_condition,
        "n_samples_total": len(sample_rows),
        "n_matrix_rows": len(matrix_rows),
        "n_events_with_bh_padj": sum(int(row["n_events_with_bh_padj"]) for row in qc_rows),
        "n_contrasts_total": len(contrasts),
        "n_contrasts_emitted": emitted,
        "n_contrasts_skipped": len(contrasts) - emitted,
        "resources": resources_info,
        "contrast_qc_path": "contrast_qc.tsv",
    }
    write_run_summary_files(out_dir, summary_payload)
    if emitted == 0:
        raise ValueError("No splicing matrix contrasts produced emit-ready outputs.")
    return {
        "out_dir": str(out_dir),
        "n_contrasts_total": len(contrasts),
        "n_contrasts_emitted": emitted,
        "grouped_output": multiple,
    }
