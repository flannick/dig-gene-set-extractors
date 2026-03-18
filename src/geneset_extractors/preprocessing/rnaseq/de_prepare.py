from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
from typing import Any

import numpy as np

from geneset_extractors.converters import rna_deg_multi
from geneset_extractors.preprocessing.rnaseq.de_backends.lightweight import LightweightContrastResult, run_lightweight_contrast
from geneset_extractors.preprocessing.rnaseq.de_backends import r_dream, r_limma_voom
from geneset_extractors.preprocessing.rnaseq.de_design import ComparisonSpec, build_cli_comparisons, comparison_manifest_rows, parse_csv_columns, read_comparisons_tsv
from geneset_extractors.preprocessing.rnaseq.de_io import MatrixData, join_metadata, read_matrix_tsv, read_tsv_rows, subset_matrix, write_matrix_tsv, write_tsv
from geneset_extractors.preprocessing.rnaseq.de_pseudobulk import build_pseudobulk


STANDARD_DE_FIELDS = [
    "comparison_id",
    "gene_id",
    "gene_symbol",
    "logFC",
    "stat",
    "pvalue",
    "padj",
    "group_a",
    "group_b",
    "stratum",
    "backend",
    "n_group_a",
    "n_group_b",
    "mean_expr",
    "model_formula",
]


@dataclass(frozen=True)
class ComparisonSampleSelection:
    spec: ComparisonSpec
    group_b_label: str
    selected_rows: list[dict[str, str]]
    filter_scope_rows: list[dict[str, str]]
    audit_row: dict[str, Any]
    selected_sample_rows: list[dict[str, Any]]



def _clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()



def _unique_preserving(values: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for value in values:
        if value in seen or not value:
            continue
        out.append(value)
        seen.add(value)
    return out


def _metadata_columns_present(rows: list[dict[str, str]], columns: list[str]) -> tuple[list[str], list[str]]:
    if not columns:
        return [], []
    present: set[str] = set()
    populated: set[str] = set()
    for row in rows:
        for key, value in row.items():
            if key:
                key_str = str(key)
                present.add(key_str)
                if _clean(value):
                    populated.add(key_str)
    missing = [col for col in columns if col not in present]
    blank_only = [col for col in columns if col in present and col not in populated]
    return missing, blank_only


def _stable_seed(*parts: object) -> int:
    digest = hashlib.sha256("|".join(str(part) for part in parts).encode("utf-8")).hexdigest()
    return int(digest[:8], 16)



def _is_integer_like_matrix(counts: np.ndarray) -> bool:
    if counts.size == 0:
        return True
    rounded = np.rint(counts)
    return bool(np.allclose(counts, rounded))



def _validate_count_matrix(counts: np.ndarray, *, allow_non_count_input: bool) -> dict[str, Any]:
    if np.any(counts < 0):
        raise ValueError("Count matrix contains negative values, which are not valid for the RNA DE workflow")
    is_integer_like = _is_integer_like_matrix(counts)
    if not allow_non_count_input and not is_integer_like:
        raise ValueError(
            "Matrix contains non-integer values. Pass --allow_non_count_input true to run an approximate analysis on normalized values."
        )
    return {
        "allow_non_count_input": bool(allow_non_count_input),
        "non_count_input_detected": bool(not is_integer_like),
    }



def _standardize_metadata_rows(rows: list[dict[str, str]], *, source_id_column: str) -> list[dict[str, str]]:
    standardized: list[dict[str, str]] = []
    for row in rows:
        sample_id = _clean(row.get(source_id_column))
        if not sample_id:
            continue
        out = {str(k): _clean(v) for k, v in row.items()}
        out["sample_id"] = sample_id
        standardized.append(out)
    if not standardized:
        raise ValueError(f"No metadata rows contained non-empty IDs in column '{source_id_column}'")
    return standardized



def _align_matrix_and_metadata(matrix: MatrixData, metadata_rows: list[dict[str, str]]) -> tuple[MatrixData, list[dict[str, str]], dict[str, int]]:
    meta_by_sample = {str(row.get("sample_id", "")).strip(): row for row in metadata_rows if str(row.get("sample_id", "")).strip()}
    keep_samples = [sample_id for sample_id in matrix.sample_ids if sample_id in meta_by_sample]
    dropped_matrix = len(matrix.sample_ids) - len(keep_samples)
    dropped_meta = sum(1 for row in metadata_rows if str(row.get("sample_id", "")).strip() not in set(keep_samples))
    if not keep_samples:
        raise ValueError("No overlapping sample IDs were found between the count matrix and metadata table")
    aligned_matrix = subset_matrix(matrix, keep_samples)
    aligned_metadata = [meta_by_sample[sample_id] for sample_id in keep_samples]
    return aligned_matrix, aligned_metadata, {
        "n_samples_matrix": len(matrix.sample_ids),
        "n_samples_metadata": len(metadata_rows),
        "n_samples_overlap": len(keep_samples),
        "n_samples_dropped_from_matrix": dropped_matrix,
        "n_rows_dropped_from_metadata": dropped_meta,
    }



def _materialize_comparison_samples(metadata_rows: list[dict[str, str]], spec: ComparisonSpec) -> list[dict[str, str]]:
    rows = metadata_rows
    for col, value in spec.strata.items():
        rows = [row for row in rows if _clean(row.get(col)) == value]
    if spec.comparison_kind == "group_vs_rest" or spec.group_b is None:
        rows = [row for row in rows if _clean(row.get(spec.group_column))]
    else:
        rows = [row for row in rows if _clean(row.get(spec.group_column)) in {spec.group_a, spec.group_b}]
    return rows


def _comparison_group_label(row: dict[str, str], spec: ComparisonSpec) -> str | None:
    group_value = _clean(row.get(spec.group_column))
    if not group_value:
        return None
    if group_value == spec.group_a:
        return spec.group_a
    if spec.comparison_kind == "group_vs_rest" or spec.group_b is None:
        return "rest"
    if group_value == spec.group_b:
        return spec.group_b
    return None


def _select_balanced_rows(
    rows: list[dict[str, str]],
    *,
    spec: ComparisonSpec,
    balance_groups: bool,
    balance_seed: int,
    gene_filter_scope: str,
    de_mode: str,
) -> ComparisonSampleSelection:
    group_b_label = spec.group_b if spec.group_b else "rest"
    grouped_rows: dict[str, list[dict[str, str]]] = {
        spec.group_a: [],
        group_b_label: [],
    }
    for row in rows:
        label = _comparison_group_label(row, spec)
        if label in grouped_rows:
            grouped_rows[label].append(row)
    pre_a = len(grouped_rows[spec.group_a])
    pre_b = len(grouped_rows[group_b_label])
    target_size = min(pre_a, pre_b) if balance_groups else 0
    balance_applied = False
    selected_groups: dict[str, list[dict[str, str]]] = {}
    effective_seeds: dict[str, int] = {}
    balance_sampler = "none"
    for label, group_rows in grouped_rows.items():
        original_rows = list(group_rows)
        ordered_rows = sorted(group_rows, key=lambda row: _clean(row.get("sample_id")))
        if not balance_groups or len(ordered_rows) <= target_size:
            selected_groups[label] = list(ordered_rows)
            continue
        balance_applied = True
        if de_mode == "harmonizome":
            balance_sampler = "harmonizome_random_state"
            seed_value = int(balance_seed)
            effective_seeds[label] = seed_value
            rng = np.random.RandomState(seed_value)
            picked_idx = list(rng.choice(len(original_rows), size=target_size, replace=False))
            selected_groups[label] = [original_rows[idx] for idx in picked_idx]
        else:
            balance_sampler = "hashed_default_rng"
            seed_value = _stable_seed(balance_seed, spec.comparison_id, label)
            effective_seeds[label] = seed_value
            rng = np.random.default_rng(seed_value)
            picked_idx = sorted(int(idx) for idx in rng.choice(len(ordered_rows), size=target_size, replace=False))
            selected_groups[label] = [ordered_rows[idx] for idx in picked_idx]
    if de_mode == "harmonizome":
        selected_rows = selected_groups[spec.group_a] + selected_groups[group_b_label]
    else:
        selected_rows = sorted(
            selected_groups[spec.group_a] + selected_groups[group_b_label],
            key=lambda row: (_comparison_group_label(row, spec) or "", _clean(row.get("sample_id"))),
        )
    selected_sample_rows: list[dict[str, Any]] = []
    for rank, row in enumerate(selected_rows, start=1):
        label = _comparison_group_label(row, spec) or ""
        sample_row: dict[str, Any] = {
            "comparison_id": spec.comparison_id,
            "sample_id": _clean(row.get("sample_id")),
            "group_label": label,
            "group_value": _clean(row.get(spec.group_column)),
            "selected_rank": rank,
            "balance_requested": bool(balance_groups),
            "balance_applied": bool(balance_applied),
            "balance_seed": int(balance_seed),
            "balance_seed_effective": effective_seeds.get(label, ""),
        }
        for key, value in sorted(spec.strata.items()):
            sample_row[key] = value
        selected_sample_rows.append(sample_row)
    filter_scope_rows = list(rows) if gene_filter_scope == "stratum" else list(selected_rows)
    audit_row: dict[str, Any] = {
        "comparison_id": spec.comparison_id,
        "comparison_kind": spec.comparison_kind,
        "group_column": spec.group_column,
        "group_a": spec.group_a,
        "group_b": group_b_label,
        "gene_filter_scope": gene_filter_scope,
        "n_filter_scope_samples": len(filter_scope_rows),
        "balance_requested": bool(balance_groups),
        "balance_applied": bool(balance_applied),
        "balance_sampler": balance_sampler,
        "balance_seed": int(balance_seed),
        "balance_seed_group_a": effective_seeds.get(spec.group_a, ""),
        "balance_seed_group_b": effective_seeds.get(group_b_label, ""),
        "n_group_a_pre_balance": pre_a,
        "n_group_b_pre_balance": pre_b,
        "n_group_a_post_balance": len(selected_groups[spec.group_a]),
        "n_group_b_post_balance": len(selected_groups[group_b_label]),
        "n_samples_selected": len(selected_rows),
    }
    for key, value in sorted(spec.strata.items()):
        audit_row[key] = value
    return ComparisonSampleSelection(
        spec=spec,
        group_b_label=group_b_label,
        selected_rows=selected_rows,
        filter_scope_rows=filter_scope_rows,
        audit_row=audit_row,
        selected_sample_rows=selected_sample_rows,
    )


def _resolve_de_mode_settings(
    *,
    de_mode: str,
    modality: str,
    balance_groups: bool,
    balance_groups_explicit: bool,
    balance_seed: int,
    balance_seed_explicit: bool,
    gene_filter_scope: str,
    gene_filter_scope_explicit: bool,
    covariate_cols: list[str],
    batch_cols: list[str],
    repeated_measures: bool,
) -> tuple[str, bool, int, str]:
    resolved_mode = _clean(de_mode) or "modern"
    if resolved_mode not in {"modern", "harmonizome"}:
        raise ValueError(f"Unsupported de_mode: {resolved_mode}")
    resolved_balance = bool(balance_groups)
    resolved_seed = int(balance_seed)
    resolved_gene_filter_scope = _clean(gene_filter_scope) or "contrast"
    if resolved_gene_filter_scope not in {"contrast", "stratum"}:
        raise ValueError(f"Unsupported gene_filter_scope: {resolved_gene_filter_scope}")
    if resolved_mode != "harmonizome":
        return resolved_mode, resolved_balance, resolved_seed, resolved_gene_filter_scope
    if modality != "bulk":
        raise ValueError("de_mode=harmonizome currently supports bulk RNA-seq only")
    if batch_cols or repeated_measures:
        raise ValueError(
            "de_mode=harmonizome is a conservative bulk preset and does not support batch columns or repeated-measures flags"
        )
    if not balance_groups_explicit:
        resolved_balance = True
    if not balance_seed_explicit:
        resolved_seed = 1
    if not gene_filter_scope_explicit:
        resolved_gene_filter_scope = "stratum"
    return resolved_mode, resolved_balance, resolved_seed, resolved_gene_filter_scope



def _n_unique_units(rows: list[dict[str, str]], unit_column: str | None, group_column: str, level: str | None) -> int:
    if level is None:
        subset = rows
    elif level == "rest":
        subset = [row for row in rows if _clean(row.get(group_column))]
    else:
        subset = [row for row in rows if _clean(row.get(group_column)) == level]
    if not unit_column:
        return len(subset)
    return len({str(row.get(unit_column, "")).strip() for row in subset if str(row.get(unit_column, "")).strip()})



def _write_backend_inputs(
    work_dir: Path,
    matrix: MatrixData,
    metadata_rows: list[dict[str, str]],
    specs: list[ComparisonSpec],
    selected_sample_rows: list[dict[str, Any]],
    backend_sample_ids: list[str] | None = None,
) -> tuple[Path, Path, Path, Path]:
    work_dir.mkdir(parents=True, exist_ok=True)
    selected_sample_ids = backend_sample_ids or _unique_preserving([_clean(row.get("sample_id")) for row in selected_sample_rows])
    if selected_sample_ids:
        selected_sample_id_set = set(selected_sample_ids)
        matrix = subset_matrix(matrix, [sample_id for sample_id in matrix.sample_ids if sample_id in selected_sample_id_set])
        metadata_rows = [row for row in metadata_rows if _clean(row.get("sample_id")) in selected_sample_id_set]
    counts_path = work_dir / "counts_gene_by_sample.tsv"
    counts_rows: list[dict[str, Any]] = []
    for feature_id, gene_symbol, values in zip(matrix.feature_ids, matrix.gene_symbols, matrix.counts.T, strict=False):
        row: dict[str, Any] = {matrix.feature_id_column: feature_id, "gene_symbol": gene_symbol or feature_id}
        for sample_id, value in zip(matrix.sample_ids, values, strict=False):
            row[sample_id] = float(value)
        counts_rows.append(row)
    write_tsv(counts_path, counts_rows, fieldnames=[matrix.feature_id_column, "gene_symbol"] + list(matrix.sample_ids))
    metadata_path = work_dir / "metadata.tsv"
    metadata_fieldnames = _unique_preserving(["sample_id"] + [col for row in metadata_rows for col in row.keys()])
    write_tsv(metadata_path, metadata_rows, fieldnames=metadata_fieldnames)
    comparisons_path = work_dir / "comparisons.tsv"
    write_tsv(comparisons_path, comparison_manifest_rows(specs))
    selected_path = work_dir / "comparison_selected_samples.tsv"
    selected_fieldnames = _unique_preserving(
        [
            "comparison_id",
            "sample_id",
            "group_label",
            "group_value",
            "selected_rank",
            "balance_requested",
            "balance_applied",
            "balance_seed",
            "balance_seed_effective",
        ]
        + [str(key) for row in selected_sample_rows for key in row.keys()]
    )
    write_tsv(selected_path, selected_sample_rows, fieldnames=selected_fieldnames)
    return counts_path, metadata_path, comparisons_path, selected_path



def _run_r_backend(
    *,
    backend_name: str,
    out_dir: Path,
    matrix: MatrixData,
    metadata_rows: list[dict[str, str]],
    specs: list[ComparisonSpec],
    selected_sample_rows: list[dict[str, Any]],
    covariates: list[str],
    batch_columns: list[str],
    random_effect_column: str | None,
    gene_filter_scope: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    work_dir = out_dir / "backend_work"
    counts_path, metadata_path, comparisons_path, selected_samples_path = _write_backend_inputs(
        work_dir,
        matrix,
        metadata_rows,
        specs,
        selected_sample_rows,
        backend_sample_ids=(
            _unique_preserving([_clean(row.get("sample_id")) for row in metadata_rows])
            if gene_filter_scope == "stratum"
            else None
        ),
    )
    output_path = work_dir / "deg_long.tsv"
    if backend_name == "r_limma_voom":
        available, detail = r_limma_voom.backend_available()
        if not available:
            raise ValueError(f"Requested backend r_limma_voom is not available: {detail}")
        script_path = r_limma_voom.write_script(
            script_path=work_dir / "run_limma_voom.R",
            counts_tsv=counts_path,
            metadata_tsv=metadata_path,
            comparisons_tsv=comparisons_path,
            selected_samples_tsv=selected_samples_path,
            output_tsv=output_path,
            covariates_csv=",".join(covariates),
            batch_columns_csv=",".join(batch_columns),
            gene_filter_scope=gene_filter_scope,
        )
        rscript = detail
    elif backend_name == "r_dream":
        available, detail = r_dream.backend_available()
        if not available:
            raise ValueError(f"Requested backend r_dream is not available: {detail}")
        if not random_effect_column:
            raise ValueError("backend=r_dream requires a random_effect_column (subject/donor column)")
        script_path = r_dream.write_script(
            script_path=work_dir / "run_dream.R",
            counts_tsv=counts_path,
            metadata_tsv=metadata_path,
            comparisons_tsv=comparisons_path,
            selected_samples_tsv=selected_samples_path,
            output_tsv=output_path,
            covariates_csv=",".join(covariates),
            batch_columns_csv=",".join(batch_columns),
            random_effect_column=random_effect_column,
            gene_filter_scope=gene_filter_scope,
        )
        rscript = detail
    else:
        raise ValueError(f"Unsupported R backend: {backend_name}")
    result = subprocess.run([rscript, "--vanilla", str(script_path)], capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise ValueError(
            f"{backend_name} failed with exit code {result.returncode}. stderr: {result.stderr.strip() or '[empty]'}"
        )
    table = read_tsv_rows(output_path)
    rows = [{field: row.get(field, "") for field in STANDARD_DE_FIELDS} for row in table.rows]
    return rows, {
        "backend": backend_name,
        "backend_script": str(script_path),
        "backend_stdout": result.stdout,
        "backend_stderr": result.stderr,
        "backend_inputs_dir": str(work_dir),
        "n_rows": len(rows),
    }



def _resolve_backend(requested: str, *, repeated_measures: bool) -> str:
    if requested != "auto":
        return requested
    if repeated_measures:
        available, _detail = r_dream.backend_available()
        if available:
            return "r_dream"
    available, _detail = r_limma_voom.backend_available()
    if available:
        return "r_limma_voom"
    return "lightweight"



def run_de_prepare(
    *,
    modality: str,
    counts_tsv: str,
    out_dir: str,
    organism: str,
    genome_build: str,
    matrix_orientation: str,
    feature_id_column: str,
    matrix_gene_symbol_column: str | None,
    matrix_delim: str,
    metadata_delim: str,
    sample_id_column: str | None,
    sample_metadata_tsv: str | None,
    subject_metadata_tsv: str | None,
    subject_join_sample_column: str | None,
    subject_join_metadata_column: str | None,
    subject_column: str | None,
    cell_id_column: str | None,
    cell_metadata_tsv: str | None,
    donor_column: str | None,
    cell_type_column: str | None,
    pseudobulk_within_cell_type: bool,
    min_cells_per_pseudobulk: int,
    min_donors_per_group: int,
    group_column: str | None,
    comparison_mode: str | None,
    condition_a: str | None,
    condition_b: str | None,
    reference_level: str | None,
    comparisons_tsv: str | None,
    stratify_by: str | None,
    covariates: str | None,
    batch_columns: str | None,
    de_mode: str,
    balance_groups: bool,
    balance_seed: int,
    gene_filter_scope: str,
    balance_groups_explicit: bool,
    balance_seed_explicit: bool,
    gene_filter_scope_explicit: bool,
    repeated_measures: bool,
    allow_approximate_repeated_measures: bool,
    backend: str,
    allow_non_count_input: bool,
    write_pseudobulk_artifacts: bool,
    run_extractor: bool,
    extractor_out_dir: str | None,
    extractor_signature_name: str | None,
    extractor_score_mode: str,
    extractor_select: str,
    extractor_top_k: int,
    extractor_quantile: float,
    extractor_min_score: float,
    extractor_normalize: str,
    extractor_padj_max: float | None,
    extractor_pvalue_max: float | None,
    extractor_min_abs_logfc: float | None,
    extractor_emit_gmt: bool,
    extractor_gmt_split_signed: bool,
    extractor_gmt_topk_list: str,
    extractor_gmt_min_genes: int,
    extractor_gmt_max_genes: int,
) -> dict[str, Any]:
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    workflow_warnings: list[str] = []

    stratify_cols = parse_csv_columns(stratify_by)
    covariate_cols = parse_csv_columns(covariates)
    batch_cols = parse_csv_columns(batch_columns)
    resolved_de_mode, resolved_balance_groups, resolved_balance_seed, resolved_gene_filter_scope = _resolve_de_mode_settings(
        de_mode=de_mode,
        modality=modality,
        balance_groups=balance_groups,
        balance_groups_explicit=balance_groups_explicit,
        balance_seed=balance_seed,
        balance_seed_explicit=balance_seed_explicit,
        gene_filter_scope=gene_filter_scope,
        gene_filter_scope_explicit=gene_filter_scope_explicit,
        covariate_cols=covariate_cols,
        batch_cols=batch_cols,
        repeated_measures=repeated_measures,
    )

    if modality == "bulk":
        if not sample_metadata_tsv or not sample_id_column:
            raise ValueError("bulk mode requires --sample_metadata_tsv and --sample_id_column")
        matrix = read_matrix_tsv(
            counts_tsv,
            orientation=matrix_orientation,
            sample_id_column=sample_id_column,
            feature_id_column=feature_id_column,
            gene_symbol_column=matrix_gene_symbol_column,
            delimiter=matrix_delim,
        )
        count_input_summary = _validate_count_matrix(matrix.counts, allow_non_count_input=allow_non_count_input)
        sample_meta = _standardize_metadata_rows(read_tsv_rows(sample_metadata_tsv, delimiter=metadata_delim).rows, source_id_column=sample_id_column)
        join_summary = {"n_joined": 0, "n_missing_join": 0}
        if subject_metadata_tsv:
            if not subject_join_sample_column or not subject_join_metadata_column:
                raise ValueError(
                    "subject_metadata_tsv requires both --subject_join_sample_column and --subject_join_metadata_column"
                )
            subject_rows = read_tsv_rows(subject_metadata_tsv, delimiter=metadata_delim).rows
            sample_meta, join_summary = join_metadata(
                sample_meta,
                primary_key=subject_join_sample_column,
                secondary_rows=subject_rows,
                secondary_key=subject_join_metadata_column,
                suffix="_subject",
            )
        matrix_aligned, metadata_rows, align_summary = _align_matrix_and_metadata(matrix, sample_meta)
        pseudobulk_summary = None
        unit_column = subject_column
        prepared_matrix = matrix_aligned
    elif modality == "scrna":
        if not cell_metadata_tsv or not cell_id_column or not donor_column:
            raise ValueError("scrna mode requires --cell_metadata_tsv, --cell_id_column, and --donor_column")
        matrix = read_matrix_tsv(
            counts_tsv,
            orientation=matrix_orientation,
            sample_id_column=cell_id_column,
            feature_id_column=feature_id_column,
            gene_symbol_column=matrix_gene_symbol_column,
            delimiter=matrix_delim,
        )
        count_input_summary = _validate_count_matrix(matrix.counts, allow_non_count_input=allow_non_count_input)
        cell_meta = _standardize_metadata_rows(read_tsv_rows(cell_metadata_tsv, delimiter=metadata_delim).rows, source_id_column=cell_id_column)
        cell_matrix, cell_metadata_rows, align_summary = _align_matrix_and_metadata(matrix, cell_meta)
        pb_group_columns = _unique_preserving(
            ([cell_type_column] if pseudobulk_within_cell_type and cell_type_column else [])
            + ([group_column] if group_column else [])
            + stratify_cols
            + covariate_cols
            + batch_cols
        )
        pb_result = build_pseudobulk(
            cell_matrix,
            cell_metadata_rows=cell_metadata_rows,
            cell_id_column="sample_id",
            donor_column=donor_column,
            group_columns=[col for col in pb_group_columns if col and col != donor_column],
            min_cells_per_pseudobulk=min_cells_per_pseudobulk,
        )
        prepared_matrix = pb_result.matrix
        metadata_rows = []
        for row in pb_result.metadata_rows:
            out = dict(row)
            out["sample_id"] = str(row.get("sample_id", "")).strip()
            metadata_rows.append(out)
        join_summary = {"n_joined": 0, "n_missing_join": 0}
        pseudobulk_summary = {
            "n_cells_input": pb_result.n_cells_input,
            "n_cells_used": pb_result.n_cells_used,
            "n_pseudobulks": pb_result.n_pseudobulks,
            "n_dropped_small_groups": pb_result.n_dropped_small_groups,
        }
        unit_column = donor_column
        if write_pseudobulk_artifacts:
            write_matrix_tsv(output_dir / "pseudobulk_matrix.tsv", prepared_matrix)
            write_tsv(output_dir / "pseudobulk_metadata.tsv", metadata_rows)
    else:
        raise ValueError(f"Unsupported modality: {modality}")

    missing_covariates, blank_covariates = _metadata_columns_present(metadata_rows, covariate_cols)
    missing_batch_columns, blank_batch_columns = _metadata_columns_present(metadata_rows, batch_cols)
    if missing_covariates:
        raise ValueError("Requested covariates were not found in the prepared metadata: " + ", ".join(missing_covariates))
    if missing_batch_columns:
        raise ValueError(
            "Requested batch columns were not found in the prepared metadata: " + ", ".join(missing_batch_columns)
        )
    if blank_covariates:
        workflow_warnings.append(
            "Some requested covariates were present but blank for all aligned samples: " + ", ".join(blank_covariates)
        )
    if blank_batch_columns:
        workflow_warnings.append(
            "Some requested batch columns were present but blank for all aligned samples: "
            + ", ".join(blank_batch_columns)
        )
    if resolved_de_mode == "harmonizome" and not covariate_cols:
        workflow_warnings.append(
            "de_mode=harmonizome was run without explicit covariates. For broad-tissue GTEx/Harmonizome-style analyses, "
            "prefer explicit fixed-effect covariates such as SEX,SMTSD when those columns exist."
        )
    for warning in workflow_warnings:
        print(f"warning: {warning}", file=sys.stderr)

    covariates_used = ",".join(covariate_cols)
    batch_columns_used = ",".join(batch_cols)
    harmonizome_covariate_mode = "explicit" if covariate_cols else "none"

    if comparisons_tsv:
        specs = read_comparisons_tsv(comparisons_tsv, default_group_column=group_column, stratify_by=stratify_cols)
    else:
        if not group_column:
            raise ValueError("--group_column is required for comparison generation")
        if not comparison_mode:
            raise ValueError("Provide either --comparisons_tsv or --comparison_mode")
        specs = build_cli_comparisons(
            metadata_rows,
            comparison_mode=comparison_mode,
            group_column=group_column,
            condition_a=condition_a,
            condition_b=condition_b,
            reference_level=reference_level,
            stratify_by=stratify_cols,
        )
    if resolved_de_mode == "harmonizome":
        bad_specs = [
            spec.comparison_id
            for spec in specs
            if spec.comparison_kind == "group_vs_rest" or not spec.group_b
        ]
        if bad_specs:
            raise ValueError(
                "de_mode=harmonizome requires explicit two-group contrasts with group_a and group_b. "
                f"Unsupported comparisons: {', '.join(bad_specs[:5])}"
            )

    resolved_backend = _resolve_backend(backend, repeated_measures=repeated_measures)
    if resolved_de_mode == "harmonizome" and resolved_backend == "r_dream":
        raise ValueError("de_mode=harmonizome requires a simple two-group backend and does not support backend=r_dream")
    approximate_repeated_measures = False
    random_effect_column = donor_column if modality == "scrna" else subject_column
    if repeated_measures and resolved_backend == "lightweight":
        if not allow_approximate_repeated_measures:
            raise ValueError(
                "Repeated-measures design requested but only lightweight backend is available. "
                "Use backend=r_dream or rerun with --allow_approximate_repeated_measures true."
            )
        approximate_repeated_measures = True

    prepared_sample_id_set = set(prepared_matrix.sample_ids)
    comparison_selections: list[ComparisonSampleSelection] = []
    comparison_audit_rows: list[dict[str, Any]] = []
    selected_sample_rows_all: list[dict[str, Any]] = []
    for spec in specs:
        materialized_rows = _materialize_comparison_samples(metadata_rows, spec)
        materialized_rows = [
            row
            for row in materialized_rows
            if _clean(row.get("sample_id")) in prepared_sample_id_set
        ]
        selection = _select_balanced_rows(
            materialized_rows,
            spec=spec,
            balance_groups=resolved_balance_groups,
            balance_seed=resolved_balance_seed,
            gene_filter_scope=resolved_gene_filter_scope,
            de_mode=resolved_de_mode,
        )
        selected_sample_rows_all.extend(selection.selected_sample_rows)
        audit_row = dict(selection.audit_row)
        audit_row["de_mode"] = resolved_de_mode
        audit_row["gene_filter_scope"] = resolved_gene_filter_scope
        audit_row["covariates_used"] = covariates_used
        audit_row["batch_columns_used"] = batch_columns_used
        audit_row["harmonizome_covariate_mode"] = harmonizome_covariate_mode
        audit_row["workflow_warning_count"] = len(workflow_warnings)
        audit_row["unit_column"] = unit_column or "sample_id"
        n_group_a_units = _n_unique_units(selection.selected_rows, unit_column, spec.group_column, spec.group_a)
        if spec.group_b is None or spec.comparison_kind == "group_vs_rest":
            group_b_rows = [row for row in selection.selected_rows if _comparison_group_label(row, spec) == selection.group_b_label]
            n_group_b_units = _n_unique_units(group_b_rows, unit_column, spec.group_column, None)
        else:
            n_group_b_units = _n_unique_units(selection.selected_rows, unit_column, spec.group_column, selection.group_b_label)
        audit_row["n_group_a_units_post"] = n_group_a_units
        audit_row["n_group_b_units_post"] = n_group_b_units
        if n_group_a_units < int(min_donors_per_group) or n_group_b_units < int(min_donors_per_group):
            audit_row["selected_for_fit"] = False
            audit_row["skip_reason"] = "insufficient_units_after_selection"
            comparison_audit_rows.append(audit_row)
            continue
        audit_row["selected_for_fit"] = True
        audit_row["skip_reason"] = ""
        comparison_audit_rows.append(audit_row)
        comparison_selections.append(selection)

    all_rows: list[dict[str, Any]] = []
    backend_summary: dict[str, Any]
    emitted_specs = [selection.spec for selection in comparison_selections]
    backend_selected_rows = [row for selection in comparison_selections for row in selection.selected_sample_rows]
    if resolved_backend in {"r_limma_voom", "r_dream"}:
        all_rows, backend_summary = _run_r_backend(
            backend_name=resolved_backend,
            out_dir=output_dir,
            matrix=prepared_matrix,
            metadata_rows=metadata_rows,
            specs=emitted_specs,
            selected_sample_rows=backend_selected_rows,
            covariates=covariate_cols,
            batch_columns=batch_cols,
            random_effect_column=random_effect_column,
            gene_filter_scope=resolved_gene_filter_scope,
        )
    else:
        backend_summary = {"backend": "lightweight", "n_rows": 0}
        for selection in comparison_selections:
            spec = selection.spec
            sample_ids = [_clean(row.get("sample_id")) for row in selection.selected_rows if _clean(row.get("sample_id")) in prepared_sample_id_set]
            if not sample_ids:
                continue
            matrix_sub = subset_matrix(prepared_matrix, sample_ids)
            result: LightweightContrastResult = run_lightweight_contrast(
                counts=matrix_sub.counts,
                feature_ids=matrix_sub.feature_ids,
                gene_symbols=matrix_sub.gene_symbols,
                sample_ids=matrix_sub.sample_ids,
                metadata_rows=selection.selected_rows,
                spec=spec,
                covariates=covariate_cols,
                batch_columns=batch_cols,
                low_expression_filter=True,
            )
            all_rows.extend(result.rows)
        backend_summary["n_rows"] = len(all_rows)

    if not all_rows:
        raise ValueError("The DE workflow did not emit any contrast rows")

    metrics_by_comparison: dict[str, dict[str, Any]] = {}
    for row in all_rows:
        comparison_id = _clean(row.get("comparison_id"))
        if not comparison_id:
            continue
        entry = metrics_by_comparison.setdefault(
            comparison_id,
            {
                "n_tested_genes": 0,
                "n_significant_genes_padj_0_05": 0,
                "n_positive_significant_genes_padj_0_05": 0,
                "n_negative_significant_genes_padj_0_05": 0,
            },
        )
        entry["n_tested_genes"] += 1
        try:
            padj_value = float(row.get("padj", ""))
        except (TypeError, ValueError):
            padj_value = float("nan")
        try:
            logfc_value = float(row.get("logFC", ""))
        except (TypeError, ValueError):
            logfc_value = float("nan")
        if np.isfinite(padj_value) and padj_value <= 0.05:
            entry["n_significant_genes_padj_0_05"] += 1
            if np.isfinite(logfc_value) and logfc_value > 0:
                entry["n_positive_significant_genes_padj_0_05"] += 1
            elif np.isfinite(logfc_value) and logfc_value < 0:
                entry["n_negative_significant_genes_padj_0_05"] += 1
    for audit_row in comparison_audit_rows:
        if not audit_row.get("selected_for_fit"):
            continue
        audit_row["backend"] = resolved_backend
        audit_row.update(metrics_by_comparison.get(_clean(audit_row.get("comparison_id")), {}))

    deg_long_path = output_dir / "deg_long.tsv"
    write_tsv(deg_long_path, all_rows, fieldnames=STANDARD_DE_FIELDS)
    comparison_rows: list[dict[str, Any]] = []
    for audit_row in comparison_audit_rows:
        if not audit_row.get("selected_for_fit"):
            continue
        comparison_id = _clean(audit_row.get("comparison_id"))
        merged = dict(audit_row)
        merged["backend"] = resolved_backend
        merged.update(metrics_by_comparison.get(comparison_id, {}))
        comparison_rows.append(merged)
    comparison_manifest_path = output_dir / "comparison_manifest.tsv"
    write_tsv(comparison_manifest_path, comparison_rows)
    comparison_audit_path = output_dir / "comparison_audit.tsv"
    audit_fieldnames = _unique_preserving(
        [
            "comparison_id",
            "comparison_kind",
            "group_column",
            "group_a",
            "group_b",
            "de_mode",
            "gene_filter_scope",
            "n_filter_scope_samples",
            "covariates_used",
            "batch_columns_used",
            "harmonizome_covariate_mode",
            "workflow_warning_count",
            "balance_requested",
            "balance_sampler",
            "balance_applied",
            "balance_seed",
            "balance_seed_group_a",
            "balance_seed_group_b",
            "n_group_a_pre_balance",
            "n_group_b_pre_balance",
            "n_group_a_post_balance",
            "n_group_b_post_balance",
            "n_group_a_units_post",
            "n_group_b_units_post",
            "selected_for_fit",
            "skip_reason",
            "n_tested_genes",
            "n_significant_genes_padj_0_05",
            "n_positive_significant_genes_padj_0_05",
            "n_negative_significant_genes_padj_0_05",
            "backend",
        ]
        + [str(key) for row in comparison_audit_rows for key in row.keys()]
    )
    write_tsv(comparison_audit_path, comparison_audit_rows, fieldnames=audit_fieldnames)
    selected_samples_path = output_dir / "comparison_selected_samples.tsv"
    selected_sample_fieldnames = _unique_preserving(
        [
            "comparison_id",
            "sample_id",
            "group_label",
            "group_value",
            "selected_rank",
            "balance_requested",
            "balance_applied",
            "balance_seed",
            "balance_seed_effective",
        ]
        + [str(key) for row in selected_sample_rows_all for key in row.keys()]
    )
    write_tsv(selected_samples_path, selected_sample_rows_all, fieldnames=selected_sample_fieldnames)

    extractor_result: dict[str, Any] | None = None
    if run_extractor:
        extractor_dir = Path(extractor_out_dir) if extractor_out_dir else output_dir / "extractor"
        extractor_args = SimpleNamespace(
            deg_tsv=str(deg_long_path),
            comparison_column="comparison_id",
            out_dir=str(extractor_dir),
            organism=organism,
            genome_build=genome_build,
            signature_name=extractor_signature_name or "contrast",
            gene_id_column="gene_id",
            gene_symbol_column="gene_symbol",
            stat_column="stat",
            logfc_column="logFC",
            padj_column="padj",
            pvalue_column="pvalue",
            score_column=None,
            score_mode=extractor_score_mode,
            postprocess_mode="harmonizome",
            duplicate_gene_policy="max_abs",
            neglog10p_cap=50.0,
            neglog10p_eps=1e-300,
            exclude_gene_regex=None,
            disable_default_excludes=False,
            gtf=None,
            gtf_gene_id_field="gene_id",
            gtf_source=None,
            select=extractor_select,
            top_k=extractor_top_k,
            quantile=extractor_quantile,
            min_score=extractor_min_score,
            normalize=extractor_normalize,
            emit_full=True,
            emit_gmt=extractor_emit_gmt,
            gmt_out=None,
            gmt_prefer_symbol=True,
            gmt_require_symbol=False,
            gmt_biotype_allowlist="protein_coding",
            gmt_min_genes=extractor_gmt_min_genes,
            gmt_max_genes=extractor_gmt_max_genes,
            gmt_topk_list=extractor_gmt_topk_list,
            gmt_mass_list="",
            gmt_split_signed=extractor_gmt_split_signed,
            gmt_emit_abs=False,
            gmt_source="full",
            emit_small_gene_sets=False,
            provenance_overlay_json=None,
            padj_max=extractor_padj_max,
            pvalue_max=extractor_pvalue_max,
            min_abs_logfc=extractor_min_abs_logfc,
        )
        extractor_result = rna_deg_multi.run(extractor_args)

    summary = {
        "workflow": "rna_de_prepare",
        "modality": modality,
        "organism": organism,
        "genome_build": genome_build,
        "counts_tsv": str(counts_tsv),
        "de_mode": resolved_de_mode,
        "gene_filter_scope": resolved_gene_filter_scope,
        "covariates_used": covariate_cols,
        "batch_columns_used": batch_cols,
        "harmonizome_covariate_mode": harmonizome_covariate_mode,
        "balance_groups": bool(resolved_balance_groups),
        "balance_sampler": ("harmonizome_random_state" if resolved_de_mode == "harmonizome" and resolved_balance_groups else ("hashed_default_rng" if resolved_balance_groups else "none")),
        "balance_seed": int(resolved_balance_seed),
        "resolved_backend": resolved_backend,
        "backend_requested": backend,
        "approximate_repeated_measures": bool(approximate_repeated_measures),
        "repeated_measures": bool(repeated_measures),
        "count_input_summary": count_input_summary,
        "alignment": align_summary,
        "metadata_join": join_summary,
        "pseudobulk": pseudobulk_summary,
        "n_comparisons_requested": len(specs),
        "n_comparisons_emitted": len(comparison_rows),
        "n_deg_rows": len(all_rows),
        "comparison_manifest_path": str(comparison_manifest_path),
        "comparison_audit_path": str(comparison_audit_path),
        "comparison_selected_samples_path": str(selected_samples_path),
        "warnings": workflow_warnings,
        "deg_long_path": str(deg_long_path),
        "backend_summary": backend_summary,
        "extractor": extractor_result,
    }
    summary_path = output_dir / "prepare_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "out_dir": str(output_dir),
        "deg_long_tsv": str(deg_long_path),
        "comparison_manifest_tsv": str(comparison_manifest_path),
        "comparison_audit_tsv": str(comparison_audit_path),
        "comparison_selected_samples_tsv": str(selected_samples_path),
        "n_comparisons": len(comparison_rows),
        "n_deg_rows": len(all_rows),
        "resolved_backend": resolved_backend,
        "extractor_out_dir": str(Path(extractor_result["out_dir"])) if extractor_result else None,
    }
