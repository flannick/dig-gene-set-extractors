from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path
import sys

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import input_file_record
from geneset_extractors.core.qc import write_run_summary_files
from geneset_extractors.extractors.proteomics.ptm_site_diff_workflow import (
    PTMRow,
    PTMWorkflowConfig,
    read_tsv_rows,
    run_ptm_site_diff_workflow,
)
from geneset_extractors.extractors.rnaseq.deg_scoring import sanitize_name_component


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
    output_subdir: str | None


@dataclass
class PTMMatrixWorkflowConfig:
    converter_name: str
    out_dir: Path
    organism: str
    genome_build: str
    dataset_label: str
    signature_name: str
    ptm_type: str
    ptm_matrix_tsv: str
    sample_metadata_tsv: str
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
    protein_matrix_tsv: str | None
    site_id_column: str | None
    site_group_column: str | None
    gene_id_column: str | None
    gene_symbol_column: str | None
    protein_accession_column: str | None
    residue_column: str | None
    position_column: str | None
    score_column: str | None
    stat_column: str | None
    logfc_column: str | None
    padj_column: str | None
    pvalue_column: str | None
    localization_prob_column: str | None
    peptide_count_column: str | None
    protein_logfc_column: str | None
    protein_stat_column: str | None
    score_mode: str
    score_transform: str
    protein_adjustment: str
    protein_adjustment_lambda: float
    confidence_weight_mode: str
    min_localization_prob: float
    site_dup_policy: str
    gene_aggregation: str
    gene_topk_sites: int
    ambiguous_gene_policy: str
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
    emit_small_gene_sets: bool
    neglog10p_cap: float
    neglog10p_eps: float


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


def _present_values(row: PTMRow, sample_ids: tuple[str, ...], missing_value_policy: str) -> list[float] | None:
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


def _site_statistic(
    row: PTMRow,
    *,
    sample_ids_a: tuple[str, ...],
    sample_ids_b: tuple[str, ...],
    effect_metric: str,
    missing_value_policy: str,
    min_samples_per_condition: int,
    min_present_per_condition: int,
) -> dict[str, float] | None:
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
        "log2fc": float(mean_diff),
        "stat": float(stat_value),
        "pvalue": float(_approx_two_sided_p_from_z(stat_value)),
        "n_a": float(len(values_a)),
        "n_b": float(len(values_b)),
    }


def _baseline_statistic(
    row: PTMRow,
    *,
    sample_ids: tuple[str, ...],
    missing_value_policy: str,
    min_present_per_condition: int,
) -> dict[str, float] | None:
    values = _present_values(row, sample_ids, missing_value_policy)
    if values is None:
        return None
    if len(values) < int(min_present_per_condition):
        return None
    mean_value = _mean(values)
    return {
        "mean_a": float(mean_value),
        "mean_b": 0.0,
        "log2fc": float(mean_value),
        "stat": float(mean_value),
        "pvalue": 1.0,
        "n_a": float(len(values)),
        "n_b": 0.0,
    }


def _resolve_condition_levels(
    rows: dict[str, SampleMetaRecord],
    condition_column: str,
    condition_a: str | None,
    condition_b: str | None,
) -> tuple[str, str]:
    if condition_a and condition_b:
        return str(condition_a), str(condition_b)
    levels = sorted({_clean(r.values.get(condition_column)) for r in rows.values() if _clean(r.values.get(condition_column))})
    if len(levels) != 2:
        raise ValueError(
            f"Could not infer exactly two condition levels from column '{condition_column}'. "
            f"Observed levels: {', '.join(levels) if levels else '(none)'}"
        )
    return levels[0], levels[1]


def _build_contrast_specs(cfg: PTMMatrixWorkflowConfig, sample_rows: dict[str, SampleMetaRecord]) -> tuple[list[ContrastSpec], list[dict[str, object]]]:
    qc_rows: list[dict[str, object]] = []
    sample_ids = tuple(sorted(sample_rows))
    specs: list[ContrastSpec] = []

    if cfg.study_contrast == "baseline":
        specs.append(
            ContrastSpec(
                contrast_id="baseline",
                contrast_label="baseline_all_samples",
                study_contrast="baseline",
                sample_ids_a=sample_ids,
                sample_ids_b=(),
                condition_a=None,
                condition_b=None,
                group_label=None,
                output_subdir=None,
            )
        )
        return specs, qc_rows

    if cfg.study_contrast == "condition_a_vs_b":
        if not cfg.condition_column:
            raise ValueError("study_contrast=condition_a_vs_b requires --condition_column")
        condition_a, condition_b = _resolve_condition_levels(sample_rows, cfg.condition_column, cfg.condition_a, cfg.condition_b)
        sample_ids_a = tuple(sorted(sid for sid, row in sample_rows.items() if _clean(row.values.get(cfg.condition_column or "")) == condition_a))
        sample_ids_b = tuple(sorted(sid for sid, row in sample_rows.items() if _clean(row.values.get(cfg.condition_column or "")) == condition_b))
        specs.append(
            ContrastSpec(
                contrast_id=f"{sanitize_name_component(condition_a)}_vs_{sanitize_name_component(condition_b)}",
                contrast_label=f"{condition_a}_vs_{condition_b}",
                study_contrast="condition_a_vs_b",
                sample_ids_a=sample_ids_a,
                sample_ids_b=sample_ids_b,
                condition_a=condition_a,
                condition_b=condition_b,
                group_label=None,
                output_subdir=None,
            )
        )
        return specs, qc_rows

    if cfg.study_contrast == "group_vs_rest":
        if not cfg.group_column:
            raise ValueError("study_contrast=group_vs_rest requires --group_column")
        groups = sorted({_clean(r.values.get(cfg.group_column)) for r in sample_rows.values() if _clean(r.values.get(cfg.group_column))})
        for group in groups:
            sample_ids_a = tuple(sorted(sid for sid, row in sample_rows.items() if _clean(row.values.get(cfg.group_column or "")) == group))
            sample_ids_b = tuple(sorted(sid for sid, row in sample_rows.items() if _clean(row.values.get(cfg.group_column or "")) != group))
            specs.append(
                ContrastSpec(
                    contrast_id=f"group_{sanitize_name_component(group)}_vs_rest",
                    contrast_label=f"{group}_vs_rest",
                    study_contrast="group_vs_rest",
                    sample_ids_a=sample_ids_a,
                    sample_ids_b=sample_ids_b,
                    condition_a=group,
                    condition_b="rest",
                    group_label=group,
                    output_subdir=f"group={sanitize_name_component(group)}",
                )
            )
        return specs, qc_rows

    if cfg.study_contrast == "condition_within_group":
        if not cfg.group_column or not cfg.condition_column:
            raise ValueError("study_contrast=condition_within_group requires --group_column and --condition_column")
        condition_a, condition_b = _resolve_condition_levels(sample_rows, cfg.condition_column, cfg.condition_a, cfg.condition_b)
        groups = sorted({_clean(r.values.get(cfg.group_column)) for r in sample_rows.values() if _clean(r.values.get(cfg.group_column))})
        for group in groups:
            in_group = {sid: row for sid, row in sample_rows.items() if _clean(row.values.get(cfg.group_column or "")) == group}
            sample_ids_a = tuple(sorted(sid for sid, row in in_group.items() if _clean(row.values.get(cfg.condition_column or "")) == condition_a))
            sample_ids_b = tuple(sorted(sid for sid, row in in_group.items() if _clean(row.values.get(cfg.condition_column or "")) == condition_b))
            specs.append(
                ContrastSpec(
                    contrast_id=(
                        f"group_{sanitize_name_component(group)}__"
                        f"{sanitize_name_component(condition_a)}_vs_{sanitize_name_component(condition_b)}"
                    ),
                    contrast_label=f"{group}:{condition_a}_vs_{condition_b}",
                    study_contrast="condition_within_group",
                    sample_ids_a=sample_ids_a,
                    sample_ids_b=sample_ids_b,
                    condition_a=condition_a,
                    condition_b=condition_b,
                    group_label=group,
                    output_subdir=f"group={sanitize_name_component(group)}",
                )
            )
        return specs, qc_rows

    raise ValueError(f"Unsupported study_contrast: {cfg.study_contrast}")


def _protein_key(row: PTMRow, protein_accession_column: str | None, gene_id_column: str | None, gene_symbol_column: str | None) -> str | None:
    candidates: list[str | None] = [
        protein_accession_column,
        gene_id_column,
        gene_symbol_column,
        "protein_accession",
        "accession",
        "uniprot_accession",
        "gene_id",
        "gene_symbol",
        "gene_name",
    ]
    for column in candidates:
        if not column:
            continue
        value = _clean(row.values.get(column))
        if value:
            return value
    return None


def _compute_protein_contrast_map(
    *,
    protein_rows: list[PTMRow] | None,
    contrast: ContrastSpec,
    cfg: PTMMatrixWorkflowConfig,
) -> dict[str, dict[str, float]]:
    if not protein_rows:
        return {}
    out: dict[str, dict[str, float]] = {}
    for row in protein_rows:
        key = _protein_key(row, cfg.protein_accession_column, cfg.gene_id_column, cfg.gene_symbol_column)
        if not key:
            continue
        if contrast.study_contrast == "baseline":
            stat = _baseline_statistic(
                row,
                sample_ids=contrast.sample_ids_a,
                missing_value_policy=cfg.missing_value_policy,
                min_present_per_condition=cfg.min_present_per_condition,
            )
        else:
            stat = _site_statistic(
                row,
                sample_ids_a=contrast.sample_ids_a,
                sample_ids_b=contrast.sample_ids_b,
                effect_metric=cfg.effect_metric,
                missing_value_policy=cfg.missing_value_policy,
                min_samples_per_condition=cfg.min_samples_per_condition,
                min_present_per_condition=cfg.min_present_per_condition,
            )
        if stat is None:
            continue
        out[key] = {
            "protein_log2fc": float(stat["log2fc"]),
            "protein_stat": float(stat["stat"]),
        }
    return out


def _contrast_rows_for_sites(
    *,
    site_rows: list[PTMRow],
    contrast: ContrastSpec,
    cfg: PTMMatrixWorkflowConfig,
    protein_map: dict[str, dict[str, float]],
) -> tuple[list[PTMRow], dict[str, int]]:
    contrast_rows: list[PTMRow] = []
    summary = {
        "n_input_sites": len(site_rows),
        "n_sites_retained": 0,
        "n_sites_dropped_missing": 0,
        "n_sites_with_protein_contrast": 0,
    }
    for row in site_rows:
        if contrast.study_contrast == "baseline":
            stat = _baseline_statistic(
                row,
                sample_ids=contrast.sample_ids_a,
                missing_value_policy=cfg.missing_value_policy,
                min_present_per_condition=cfg.min_present_per_condition,
            )
        else:
            stat = _site_statistic(
                row,
                sample_ids_a=contrast.sample_ids_a,
                sample_ids_b=contrast.sample_ids_b,
                effect_metric=cfg.effect_metric,
                missing_value_policy=cfg.missing_value_policy,
                min_samples_per_condition=cfg.min_samples_per_condition,
                min_present_per_condition=cfg.min_present_per_condition,
            )
        if stat is None:
            summary["n_sites_dropped_missing"] += 1
            continue
        values = {k: str(v) for k, v in row.values.items()}
        values["log2fc"] = str(stat["log2fc"])
        values["stat"] = str(stat["stat"])
        values["pvalue"] = str(stat["pvalue"])
        if contrast.study_contrast == "baseline":
            values["padj"] = ""
        values["n_group_a"] = str(int(stat["n_a"]))
        values["n_group_b"] = str(int(stat["n_b"]))
        protein_key = _protein_key(row, cfg.protein_accession_column, cfg.gene_id_column, cfg.gene_symbol_column)
        protein_stat = protein_map.get(protein_key or "")
        if protein_stat is not None:
            values["protein_log2fc"] = str(protein_stat["protein_log2fc"])
            values["protein_stat"] = str(protein_stat["protein_stat"])
            summary["n_sites_with_protein_contrast"] += 1
        contrast_rows.append(PTMRow(line_no=row.line_no, values=values))
        summary["n_sites_retained"] += 1
    return contrast_rows, summary


def _make_child_cfg(cfg: PTMMatrixWorkflowConfig, contrast: ContrastSpec, out_dir: Path) -> PTMWorkflowConfig:
    signature_bits = [
        cfg.signature_name,
        f"study_contrast={contrast.study_contrast}",
        f"contrast={contrast.contrast_id}",
    ]
    if contrast.group_label:
        signature_bits.append(f"group={contrast.group_label}")
    signature_name = "__".join(signature_bits)
    dataset_label = f"{cfg.dataset_label}::{contrast.contrast_label}"
    return PTMWorkflowConfig(
        converter_name=cfg.converter_name,
        out_dir=out_dir,
        organism=cfg.organism,
        genome_build=cfg.genome_build,
        dataset_label=dataset_label,
        signature_name=signature_name,
        ptm_type=cfg.ptm_type,
        site_id_column=cfg.site_id_column,
        site_group_column=cfg.site_group_column,
        gene_id_column=cfg.gene_id_column,
        gene_symbol_column=cfg.gene_symbol_column,
        protein_accession_column=cfg.protein_accession_column,
        residue_column=cfg.residue_column,
        position_column=cfg.position_column,
        score_column=cfg.score_column,
        stat_column=cfg.stat_column,
        logfc_column=cfg.logfc_column,
        padj_column=cfg.padj_column,
        pvalue_column=cfg.pvalue_column,
        localization_prob_column=cfg.localization_prob_column,
        peptide_count_column=cfg.peptide_count_column,
        protein_logfc_column=cfg.protein_logfc_column,
        protein_stat_column=cfg.protein_stat_column,
        score_mode=cfg.score_mode,
        score_transform=cfg.score_transform,
        protein_adjustment=cfg.protein_adjustment,
        protein_adjustment_lambda=cfg.protein_adjustment_lambda,
        confidence_weight_mode=cfg.confidence_weight_mode,
        min_localization_prob=cfg.min_localization_prob,
        site_dup_policy=cfg.site_dup_policy,
        gene_aggregation=cfg.gene_aggregation,
        gene_topk_sites=cfg.gene_topk_sites,
        ambiguous_gene_policy=cfg.ambiguous_gene_policy,
        select=cfg.select,
        top_k=cfg.top_k,
        quantile=cfg.quantile,
        min_score=cfg.min_score,
        normalize=cfg.normalize,
        emit_full=cfg.emit_full,
        emit_gmt=cfg.emit_gmt,
        gmt_out=cfg.gmt_out,
        gmt_prefer_symbol=cfg.gmt_prefer_symbol,
        gmt_require_symbol=cfg.gmt_require_symbol,
        gmt_biotype_allowlist=cfg.gmt_biotype_allowlist,
        gmt_min_genes=cfg.gmt_min_genes,
        gmt_max_genes=cfg.gmt_max_genes,
        gmt_topk_list=cfg.gmt_topk_list,
        gmt_mass_list=cfg.gmt_mass_list,
        gmt_split_signed=cfg.gmt_split_signed,
        emit_small_gene_sets=cfg.emit_small_gene_sets,
        neglog10p_cap=cfg.neglog10p_cap,
        neglog10p_eps=cfg.neglog10p_eps,
    )


def run_ptm_site_matrix_workflow(
    *,
    cfg: PTMMatrixWorkflowConfig,
    alias_map: dict[str, object] | None,
    ubiquity_map: dict[str, object] | None,
    input_files: list[dict[str, str]],
    resources_info: dict[str, object] | None,
) -> dict[str, object]:
    if cfg.matrix_format != "wide_sites_by_sample":
        raise ValueError(f"Unsupported matrix_format: {cfg.matrix_format}")

    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    matrix_fieldnames, matrix_rows = read_tsv_rows(cfg.ptm_matrix_tsv)
    _meta_fieldnames, meta_rows = read_sample_metadata_tsv(cfg.sample_metadata_tsv)
    sample_rows = _sample_rows_by_id(meta_rows, cfg.sample_id_column)
    available_sample_ids = [sid for sid in sample_rows if sid in matrix_fieldnames]
    if not available_sample_ids:
        raise ValueError("No sample IDs from sample_metadata_tsv were found as columns in ptm_matrix_tsv.")
    sample_rows = {sid: sample_rows[sid] for sid in available_sample_ids}

    protein_rows: list[PTMRow] | None = None
    if cfg.protein_matrix_tsv:
        _protein_fieldnames, protein_rows = read_tsv_rows(cfg.protein_matrix_tsv)

    contrasts, _qc_seed = _build_contrast_specs(cfg, sample_rows)
    if not contrasts:
        raise ValueError("No PTM matrix contrasts could be constructed from the supplied metadata.")

    files = list(input_files)
    manifest_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    combined_gmt_sets: list[tuple[str, list[str]]] = []
    multiple = len(contrasts) > 1
    emitted = 0

    for contrast in contrasts:
        child_out_dir = out_dir if not multiple else out_dir / str(contrast.output_subdir or f"contrast={sanitize_name_component(contrast.contrast_id)}")
        protein_map = _compute_protein_contrast_map(protein_rows=protein_rows, contrast=contrast, cfg=cfg)
        contrast_rows, contrast_summary = _contrast_rows_for_sites(
            site_rows=matrix_rows,
            contrast=contrast,
            cfg=cfg,
            protein_map=protein_map,
        )
        qc_row = {
            "contrast_id": contrast.contrast_id,
            "contrast_label": contrast.contrast_label,
            "study_contrast": contrast.study_contrast,
            "group_label": contrast.group_label or "",
            "condition_a": contrast.condition_a or "",
            "condition_b": contrast.condition_b or "",
            "n_samples_a": len(contrast.sample_ids_a),
            "n_samples_b": len(contrast.sample_ids_b),
            "n_input_sites": contrast_summary["n_input_sites"],
            "n_sites_retained": contrast_summary["n_sites_retained"],
            "n_sites_dropped_missing": contrast_summary["n_sites_dropped_missing"],
            "n_sites_with_protein_contrast": contrast_summary["n_sites_with_protein_contrast"],
            "status": "pending",
            "reason_if_skipped": "",
            "path": "",
        }
        if not contrast_rows:
            qc_row["status"] = "skipped"
            qc_row["reason_if_skipped"] = (
                "No sites passed contrast QC; likely too many missing values or too few samples per condition."
            )
            qc_rows.append(qc_row)
            print(
                f"warning: skipped PTM contrast {contrast.contrast_label}; no sites passed matrix contrast QC. "
                f"Try lowering --min_samples_per_condition or --min_present_per_condition.",
                file=sys.stderr,
            )
            continue

        child_cfg = _make_child_cfg(cfg, contrast, child_out_dir)
        result = run_ptm_site_diff_workflow(
            cfg=child_cfg,
            fieldnames=list(matrix_fieldnames) + ["log2fc", "stat", "pvalue", "padj", "n_group_a", "n_group_b", "protein_log2fc", "protein_stat"],
            rows=contrast_rows,
            alias_map=alias_map,
            ubiquity_map=ubiquity_map,
            input_files=files,
            resources_info=resources_info,
        )
        emitted += 1
        qc_row["status"] = "emitted"
        qc_row["path"] = "." if not multiple else str(child_out_dir.relative_to(out_dir))
        qc_rows.append(qc_row)
        if multiple:
            manifest_rows.append(
                {
                    "contrast_id": contrast.contrast_id,
                    "study_contrast": contrast.study_contrast,
                    "group_label": contrast.group_label or "",
                    "condition_a": contrast.condition_a or "",
                    "condition_b": contrast.condition_b or "",
                    "path": str(child_out_dir.relative_to(out_dir)),
                }
            )
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
            writer = csv.DictWriter(
                fh,
                delimiter="\t",
                fieldnames=["contrast_id", "study_contrast", "group_label", "condition_a", "condition_b", "path"],
            )
            writer.writeheader()
            for row in manifest_rows:
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
                "n_input_sites",
                "n_sites_retained",
                "n_sites_dropped_missing",
                "n_sites_with_protein_contrast",
                "status",
                "reason_if_skipped",
                "path",
            ],
        )
        writer.writeheader()
        for row in qc_rows:
            writer.writerow(row)

    if multiple and cfg.emit_gmt and combined_gmt_sets:
        write_gmt(combined_gmt_sets, out_dir / "genesets.gmt")

    summary_payload = {
        "converter": cfg.converter_name,
        "dataset_label": cfg.dataset_label,
        "signature_name": cfg.signature_name,
        "study_contrast": cfg.study_contrast,
        "matrix_format": cfg.matrix_format,
        "effect_metric": cfg.effect_metric,
        "missing_value_policy": cfg.missing_value_policy,
        "min_samples_per_condition": cfg.min_samples_per_condition,
        "min_present_per_condition": cfg.min_present_per_condition,
        "n_samples_total": len(sample_rows),
        "n_matrix_rows": len(matrix_rows),
        "n_contrasts_total": len(contrasts),
        "n_contrasts_emitted": emitted,
        "n_contrasts_skipped": len(contrasts) - emitted,
        "resources": resources_info,
        "contrast_qc_path": "contrast_qc.tsv",
    }
    write_run_summary_files(out_dir, summary_payload)

    if emitted == 0:
        raise ValueError("No PTM matrix contrasts produced emit-ready outputs.")

    return {
        "out_dir": str(out_dir),
        "n_contrasts_total": len(contrasts),
        "n_contrasts_emitted": emitted,
        "grouped_output": multiple,
    }
