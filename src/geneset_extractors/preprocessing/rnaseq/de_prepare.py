from __future__ import annotations

import json
from pathlib import Path
import subprocess
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



def _write_backend_inputs(work_dir: Path, matrix: MatrixData, metadata_rows: list[dict[str, str]], specs: list[ComparisonSpec]) -> tuple[Path, Path, Path]:
    work_dir.mkdir(parents=True, exist_ok=True)
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
    return counts_path, metadata_path, comparisons_path



def _run_r_backend(
    *,
    backend_name: str,
    out_dir: Path,
    matrix: MatrixData,
    metadata_rows: list[dict[str, str]],
    specs: list[ComparisonSpec],
    covariates: list[str],
    batch_columns: list[str],
    random_effect_column: str | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    work_dir = out_dir / "backend_work"
    counts_path, metadata_path, comparisons_path = _write_backend_inputs(work_dir, matrix, metadata_rows, specs)
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
            output_tsv=output_path,
            covariates_csv=",".join(covariates),
            batch_columns_csv=",".join(batch_columns),
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
            output_tsv=output_path,
            covariates_csv=",".join(covariates),
            batch_columns_csv=",".join(batch_columns),
            random_effect_column=random_effect_column,
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
    manifest_rows = comparison_manifest_rows(specs)
    qc_rows = [{"comparison_id": spec.comparison_id, "backend": backend_name} for spec in specs]
    return rows, manifest_rows, {
        "backend": backend_name,
        "backend_script": str(script_path),
        "backend_stdout": result.stdout,
        "backend_stderr": result.stderr,
        "backend_inputs_dir": str(work_dir),
        "n_rows": len(rows),
        "qc_rows": qc_rows,
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

    stratify_cols = parse_csv_columns(stratify_by)
    covariate_cols = parse_csv_columns(covariates)
    batch_cols = parse_csv_columns(batch_columns)

    if modality == "bulk":
        if not sample_metadata_tsv or not sample_id_column:
            raise ValueError("bulk mode requires --sample_metadata_tsv and --sample_id_column")
        matrix = read_matrix_tsv(
            counts_tsv,
            orientation=matrix_orientation,
            sample_id_column=sample_id_column,
            feature_id_column=feature_id_column,
            gene_symbol_column=None,
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
            gene_symbol_column=None,
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

    if not group_column:
        raise ValueError("--group_column is required for comparison generation")

    if comparisons_tsv:
        specs = read_comparisons_tsv(comparisons_tsv, default_group_column=group_column, stratify_by=stratify_cols)
    else:
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

    resolved_backend = _resolve_backend(backend, repeated_measures=repeated_measures)
    approximate_repeated_measures = False
    random_effect_column = donor_column if modality == "scrna" else subject_column
    if repeated_measures and resolved_backend == "lightweight":
        if not allow_approximate_repeated_measures:
            raise ValueError(
                "Repeated-measures design requested but only lightweight backend is available. "
                "Use backend=r_dream or rerun with --allow_approximate_repeated_measures true."
            )
        approximate_repeated_measures = True

    all_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    backend_summary: dict[str, Any]
    if resolved_backend in {"r_limma_voom", "r_dream"}:
        all_rows, manifest_rows, backend_summary = _run_r_backend(
            backend_name=resolved_backend,
            out_dir=output_dir,
            matrix=prepared_matrix,
            metadata_rows=metadata_rows,
            specs=specs,
            covariates=covariate_cols,
            batch_columns=batch_cols,
            random_effect_column=random_effect_column,
        )
        comparison_rows = manifest_rows
    else:
        backend_summary = {"backend": "lightweight", "n_rows": 0}
        for spec in specs:
            materialized_rows = _materialize_comparison_samples(metadata_rows, spec)
            group_b_label = spec.group_b if spec.group_b else "rest"
            n_group_a_units = _n_unique_units(materialized_rows, unit_column, spec.group_column, spec.group_a)
            n_group_b_units = _n_unique_units(
                [row for row in materialized_rows if _clean(row.get(spec.group_column)) != spec.group_a]
                if spec.group_b is None or spec.comparison_kind == "group_vs_rest"
                else materialized_rows,
                unit_column if unit_column else None,
                spec.group_column,
                group_b_label if spec.group_b is not None else None,
            )
            if n_group_a_units < int(min_donors_per_group) or n_group_b_units < int(min_donors_per_group):
                continue
            sample_ids = [str(row.get("sample_id", "")).strip() for row in materialized_rows if str(row.get("sample_id", "")).strip() in set(prepared_matrix.sample_ids)]
            if not sample_ids:
                continue
            matrix_sub = subset_matrix(prepared_matrix, sample_ids)
            result: LightweightContrastResult = run_lightweight_contrast(
                counts=matrix_sub.counts,
                feature_ids=matrix_sub.feature_ids,
                gene_symbols=matrix_sub.gene_symbols,
                sample_ids=matrix_sub.sample_ids,
                metadata_rows=materialized_rows,
                spec=spec,
                covariates=covariate_cols,
                batch_columns=batch_cols,
                low_expression_filter=True,
            )
            all_rows.extend(result.rows)
            comparison_row = {
                "comparison_id": spec.comparison_id,
                "comparison_kind": spec.comparison_kind,
                "group_column": spec.group_column,
                "group_a": spec.group_a,
                "group_b": spec.group_b or "rest",
                "backend": result.summary["backend"],
                "n_group_a": result.summary["n_group_a"],
                "n_group_b": result.summary["n_group_b"],
                "n_features_input": result.summary["n_features_input"],
                "n_features_retained": result.summary["n_features_retained"],
            }
            for key, value in sorted(spec.strata.items()):
                comparison_row[key] = value
            comparison_rows.append(comparison_row)
        backend_summary["n_rows"] = len(all_rows)

    if not all_rows:
        raise ValueError("The DE workflow did not emit any contrast rows")

    deg_long_path = output_dir / "deg_long.tsv"
    write_tsv(deg_long_path, all_rows, fieldnames=STANDARD_DE_FIELDS)
    comparison_manifest_path = output_dir / "comparison_manifest.tsv"
    write_tsv(comparison_manifest_path, comparison_rows)

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
        "n_comparisons": len(comparison_rows),
        "n_deg_rows": len(all_rows),
        "resolved_backend": resolved_backend,
        "extractor_out_dir": str(Path(extractor_result["out_dir"])) if extractor_result else None,
    }
