from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

from geneset_extractors.preprocessing.rnaseq.de_prepare import run_de_prepare
from geneset_extractors.workflows.gtex_runtime_common import write_tsv, write_workflow_provenance_graph


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _slugify(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "-", str(value or "").strip().lower()).strip("-")
    return text or "tissue"


def _augment_metadata_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    augmented: list[dict[str, str]] = []
    for row in rows:
        updated = {str(key): str(value) for key, value in row.items()}
        updated["tissue_slug"] = _slugify(updated.get("tissue", ""))
        updated["tissue_code_no"] = str(updated.get("tissue_code_no", "")).strip().lower()
        augmented.append(updated)
    return augmented


def _build_comparison_rows(
    *,
    metadata_rows: list[dict[str, str]],
    min_samples_per_group: int,
) -> tuple[list[dict[str, str]], list[dict[str, Any]]]:
    counts_by_stratum: dict[tuple[str, str, str, str], dict[str, int]] = {}
    for row in metadata_rows:
        sex_label = str(row.get("sex_label", "")).strip().lower()
        timepoint_label = str(row.get("timepoint_label", "")).strip()
        tissue_code_no = str(row.get("tissue_code_no", "")).strip().lower()
        tissue_slug = str(row.get("tissue_slug", "")).strip()
        intervention = str(row.get("intervention", "")).strip().lower()
        if not sex_label or not timepoint_label or not tissue_code_no or not tissue_slug or intervention not in {"control", "training"}:
            continue
        key = (sex_label, timepoint_label, tissue_code_no, tissue_slug)
        counts_by_stratum.setdefault(key, {"control": 0, "training": 0})
        counts_by_stratum[key][intervention] += 1

    comparisons: list[dict[str, str]] = []
    summary_rows: list[dict[str, Any]] = []
    for sex_label, timepoint_label, tissue_code_no, tissue_slug in sorted(counts_by_stratum):
        counts = counts_by_stratum[(sex_label, timepoint_label, tissue_code_no, tissue_slug)]
        comparison_id = f"{tissue_code_no}-{tissue_slug}_{sex_label}_{timepoint_label}"
        summary_rows.append(
            {
                "comparison_id": comparison_id,
                "sex_label": sex_label,
                "timepoint_label": timepoint_label,
                "n_control": counts["control"],
                "n_training": counts["training"],
            }
        )
        if counts["control"] < int(min_samples_per_group) or counts["training"] < int(min_samples_per_group):
            continue
        comparisons.append(
            {
                "comparison_id": comparison_id,
                "comparison_kind": "condition_a_vs_b",
                "group_column": "intervention",
                "group_a": "training",
                "group_b": "control",
                "sex_label": sex_label,
                "timepoint_label": timepoint_label,
                "tissue_code_no": tissue_code_no,
                "tissue_slug": tissue_slug,
            }
        )
    return comparisons, summary_rows


def run(args) -> dict[str, object]:
    counts_tsv = Path(args.counts_tsv).resolve()
    sample_metadata_tsv = Path(args.sample_metadata_tsv).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata_rows = _augment_metadata_rows(_read_tsv_rows(sample_metadata_tsv))
    augmented_metadata_path = out_dir / "sample_metadata.tsv"
    comparisons_path = out_dir / "comparisons.tsv"
    comparison_summary_path = out_dir / "comparison_summary.tsv"
    write_tsv(
        augmented_metadata_path,
        metadata_rows,
        [
            "sample_id",
            "pid",
            "bid",
            "sex",
            "sex_label",
            "intervention",
            "tissue",
            "transcript_tissue",
            "tissue_code_no",
            "timepoint_label",
            "tissue_slug",
        ],
    )
    comparisons, comparison_summary = _build_comparison_rows(
        metadata_rows=metadata_rows,
        min_samples_per_group=int(args.min_samples_per_group),
    )
    if not comparisons:
        raise ValueError("No runnable MoTrPAC timewise comparisons were identified from the prepared sample metadata.")
    write_tsv(
        comparisons_path,
        comparisons,
        [
            "comparison_id",
            "comparison_kind",
            "group_column",
            "group_a",
            "group_b",
            "sex_label",
            "timepoint_label",
            "tissue_code_no",
            "tissue_slug",
        ],
    )
    write_tsv(
        comparison_summary_path,
        comparison_summary,
        ["comparison_id", "sex_label", "timepoint_label", "n_control", "n_training"],
    )

    result = run_de_prepare(
        modality="bulk",
        counts_tsv=str(counts_tsv),
        out_dir=str(out_dir),
        organism=str(args.organism),
        genome_build=str(args.genome_build),
        matrix_orientation="gene_by_sample",
        feature_id_column="gene_id",
        matrix_gene_symbol_column="gene_symbol",
        matrix_delim="\t",
        metadata_delim="\t",
        sample_id_column="sample_id",
        sample_metadata_tsv=str(augmented_metadata_path),
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
        group_column="intervention",
        comparison_mode=None,
        condition_a=None,
        condition_b=None,
        reference_level=None,
        comparisons_tsv=str(comparisons_path),
        stratify_by="sex_label,timepoint_label,tissue_code_no,tissue_slug",
        covariates="",
        batch_columns=None,
        de_mode="modern",
        balance_groups=False,
        balance_seed=0,
        gene_filter_scope="stratum",
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
        backend="auto",
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
    deg_long_path = out_dir / "deg_long.tsv"
    comparison_manifest_path = out_dir / "comparison_manifest.tsv"
    comparison_audit_path = out_dir / "comparison_audit.tsv"
    selected_samples_path = out_dir / "comparison_selected_samples.tsv"
    write_workflow_provenance_graph(
        workflow_name="motrpac_timewise",
        module_name="geneset_extractors.workflows.motrpac_timewise",
        output_dir=out_dir,
        focus_output_path=deg_long_path,
        output_paths=[
            (deg_long_path, "deg_tsv"),
            (augmented_metadata_path, "sample_metadata_tsv"),
            (comparisons_path, "comparisons_tsv"),
            (comparison_summary_path, "comparison_summary"),
            (comparison_manifest_path, "comparison_manifest"),
            (comparison_audit_path, "comparison_audit"),
            (selected_samples_path, "comparison_selected_samples"),
        ],
        input_paths=[
            (counts_tsv, "counts_tsv"),
            (sample_metadata_tsv, "sample_metadata_tsv"),
        ],
        parameters={
            "group_column": "intervention",
            "group_a": "training",
            "group_b": "control",
            "stratify_by": ["sex_label", "timepoint_label", "tissue_code_no", "tissue_slug"],
            "de_mode": "modern",
            "gene_filter_scope": "stratum",
            "backend": "auto",
            "min_samples_per_group": int(args.min_samples_per_group),
            "n_samples": len(metadata_rows),
            "n_comparisons": len(comparisons),
        },
    )
    result["deg_long_path"] = str(deg_long_path)
    result["n_comparisons"] = len(comparisons)
    return result
