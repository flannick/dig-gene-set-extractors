from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from geneset_extractors.preprocessing.rnaseq.de_prepare import run_de_prepare
from geneset_extractors.workflows.gtex_runtime_common import (
    build_age_binned_comparisons,
    build_sample_metadata_rows,
    open_maybe_gzip,
    parse_gct_header,
    read_tsv,
    write_filtered_counts,
    write_tsv,
    write_workflow_provenance_graph,
)


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _augment_deg_long_with_labels(*, workflow_out: Path, comparisons_path: Path) -> Path:
    deg_path = workflow_out / "deg_long.tsv"
    deg_rows = _read_tsv_rows(deg_path)
    if not deg_rows:
        return deg_path
    fieldnames = list(deg_rows[0].keys())
    if "gmt_comparison_label" in fieldnames:
        return deg_path
    comparison_rows = _read_tsv_rows(comparisons_path)
    label_by_comparison = {
        str(row.get("comparison_id", "")).strip(): str(row.get("gmt_comparison_label", "")).strip()
        for row in comparison_rows
        if str(row.get("comparison_id", "")).strip()
    }
    augmented_rows: list[dict[str, str]] = []
    for row in deg_rows:
        updated = dict(row)
        updated["gmt_comparison_label"] = label_by_comparison.get(str(row.get("comparison_id", "")).strip(), "")
        augmented_rows.append(updated)
    _write_tsv_rows(deg_path, augmented_rows, [*fieldnames, "gmt_comparison_label"])
    return deg_path


def run(args) -> dict[str, object]:
    expression_gct = Path(args.expression_gct).resolve()
    sample_attributes_tsv = Path(args.sample_attributes_tsv).resolve()
    subject_phenotypes_tsv = Path(args.subject_phenotypes_tsv).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_rows = read_tsv(sample_attributes_tsv)
    subject_rows = read_tsv(subject_phenotypes_tsv)
    _header, gct_sample_ids = parse_gct_header(expression_gct)
    age_bins = [token.strip() for token in str(args.age_bins).split(",") if token.strip()]
    tissue_label = str(args.tissue_label).strip()
    tissue_id = str(args.tissue_id).strip()
    tissue_column = str(getattr(args, "tissue_column", "") or "").strip() or None
    tissue_value = str(getattr(args, "tissue_value", "") or "").strip() or None

    prepared_meta = build_sample_metadata_rows(
        sample_rows=sample_rows,
        subject_rows=subject_rows,
        gct_sample_ids=gct_sample_ids,
        tissue_label=tissue_label,
        tissue_id=tissue_id,
        tissue_column=tissue_column,
        tissue_value=tissue_value,
        age_bins=age_bins,
    )
    retained_sample_ids = [row["sample_id"] for row in prepared_meta]
    retained_sample_id_set = set(retained_sample_ids)
    sample_index = [idx for idx, sample_id in enumerate(gct_sample_ids) if sample_id in retained_sample_id_set]
    sample_columns = [gct_sample_ids[idx] for idx in sample_index]

    counts_path = out_dir / "tissue_counts.tsv"
    sample_metadata_path = out_dir / "sample_metadata.tsv"
    comparisons_path = out_dir / "comparisons.tsv"

    n_genes_retained = write_filtered_counts(
        counts_gct=expression_gct,
        sample_columns=sample_columns,
        sample_index=sample_index,
        out_path=counts_path,
    )
    write_tsv(
        sample_metadata_path,
        prepared_meta,
        ["sample_id", "subject_id", "age_bin", "SEX", "primary_tissue", "detailed_tissue", "tissue_id", "tissue_label"],
    )
    comparisons, age_counts = build_age_binned_comparisons(
        prepared_meta=prepared_meta,
        reference_age_bin=str(args.reference_age_bin),
        age_bins=age_bins,
        min_samples_per_group=int(args.min_samples_per_group),
    )
    write_tsv(
        comparisons_path,
        comparisons,
        ["comparison_id", "gmt_comparison_label", "comparison_kind", "group_column", "group_a", "group_b"],
    )

    result = run_de_prepare(
        modality="bulk",
        counts_tsv=str(counts_path),
        out_dir=str(out_dir),
        organism=str(args.organism),
        genome_build=str(args.genome_build),
        matrix_orientation="gene_by_sample",
        feature_id_column="gene_id",
        matrix_gene_symbol_column="gene_symbol",
        matrix_delim="\t",
        metadata_delim="\t",
        sample_id_column="sample_id",
        sample_metadata_tsv=str(sample_metadata_path),
        subject_metadata_tsv=None,
        subject_join_sample_column=None,
        subject_join_metadata_column=None,
        subject_column=None,
        cell_id_column=None,
        cell_metadata_tsv=None,
        donor_column=None,
        cell_type_column=None,
        pseudobulk_within_cell_type=True,
        min_cells_per_pseudobulk=10,
        min_donors_per_group=2,
        group_column="age_bin",
        comparison_mode=None,
        condition_a=None,
        condition_b=None,
        reference_level=None,
        comparisons_tsv=str(comparisons_path),
        stratify_by=None,
        covariates=str(getattr(args, "covariates", "") or ""),
        batch_columns=None,
        de_mode=str(args.de_mode),
        balance_groups=bool(args.balance_groups),
        balance_seed=int(args.balance_seed),
        gene_filter_scope=str(args.gene_filter_scope),
        feature_mapping_tsv=None,
        feature_mapping_from_column=None,
        feature_mapping_to_column=None,
        feature_mapping_strip_version=False,
        drop_unmapped_features=False,
        balance_groups_explicit=True,
        balance_seed_explicit=True,
        gene_filter_scope_explicit=True,
        repeated_measures=False,
        allow_approximate_repeated_measures=False,
        backend=str(args.backend),
        allow_non_count_input=False,
        write_pseudobulk_artifacts=True,
        run_extractor=False,
        extractor_out_dir=None,
        extractor_signature_name=None,
        extractor_score_mode="auto",
        extractor_select="top_k",
        extractor_top_k=200,
        extractor_quantile=0.01,
        extractor_min_score=0.0,
        extractor_normalize="within_set_l1",
        extractor_padj_max=None,
        extractor_pvalue_max=None,
        extractor_min_abs_logfc=None,
        extractor_emit_gmt=True,
        extractor_gmt_split_signed=True,
        extractor_gmt_topk_list="200",
        extractor_gmt_min_genes=100,
        extractor_gmt_max_genes=500,
    )
    deg_long_path = _augment_deg_long_with_labels(workflow_out=out_dir, comparisons_path=comparisons_path)
    comparison_manifest_path = out_dir / "comparison_manifest.tsv"
    comparison_audit_path = out_dir / "comparison_audit.tsv"
    selected_samples_path = out_dir / "comparison_selected_samples.tsv"
    write_workflow_provenance_graph(
        workflow_name="gtex_age_binned",
        module_name="geneset_extractors.workflows.gtex_age_binned",
        output_dir=out_dir,
        focus_output_path=deg_long_path,
        output_paths=[
            (deg_long_path, "deg_tsv"),
            (counts_path, "counts_tsv"),
            (sample_metadata_path, "sample_metadata_tsv"),
            (comparisons_path, "comparisons_tsv"),
            (comparison_manifest_path, "comparison_manifest"),
            (comparison_audit_path, "comparison_audit"),
            (selected_samples_path, "comparison_selected_samples"),
        ],
        input_paths=[
            (expression_gct, "expression_gct"),
            (sample_attributes_tsv, "sample_attributes_tsv"),
            (subject_phenotypes_tsv, "subject_phenotypes_tsv"),
        ],
        parameters={
            "tissue_id": tissue_id,
            "tissue_label": tissue_label,
            "tissue_column": tissue_column,
            "tissue_value": tissue_value,
            "reference_age_bin": str(args.reference_age_bin),
            "de_mode": str(args.de_mode),
            "backend": str(args.backend),
            "covariates": str(getattr(args, "covariates", "") or ""),
            "balance_groups": bool(args.balance_groups),
            "balance_seed": int(args.balance_seed),
            "gene_filter_scope": str(args.gene_filter_scope),
            "min_samples_per_group": int(args.min_samples_per_group),
            "n_samples_retained": len(prepared_meta),
            "n_genes_retained": n_genes_retained,
            "n_comparisons": len(comparisons),
            "age_bin_counts": {key: age_counts.get(key, 0) for key in age_bins},
        },
    )
    result["deg_long_path"] = str(deg_long_path)
    result["n_samples_retained"] = len(prepared_meta)
    result["n_genes_retained"] = n_genes_retained
    return result
