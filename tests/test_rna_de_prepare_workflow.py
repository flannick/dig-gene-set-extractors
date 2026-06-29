from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from geneset_extractors.cli import main
from geneset_extractors.preprocessing.rnaseq.de_backends import r_dream, r_limma_voom
from geneset_extractors.preprocessing.rnaseq.de_design import build_cli_comparisons
from geneset_extractors.preprocessing.rnaseq.de_pseudobulk import build_pseudobulk
from geneset_extractors.preprocessing.rnaseq.de_io import MatrixData
from geneset_extractors.workflows.rna_de_prepare import run as run_rna_de_prepare


class BulkArgs(SimpleNamespace):
    modality = "bulk"
    counts_tsv = ""
    out_dir = ""
    organism = "human"
    genome_build = "hg38"
    matrix_orientation = "sample_by_gene"
    feature_id_column = "gene_id"
    matrix_gene_symbol_column = None
    matrix_delim = "\t"
    metadata_delim = "\t"
    sample_id_column = "sample_id"
    sample_metadata_tsv = ""
    subject_metadata_tsv = None
    subject_join_sample_column = None
    subject_join_metadata_column = None
    subject_column = None
    cell_id_column = None
    cell_metadata_tsv = None
    donor_column = None
    cell_type_column = None
    pseudobulk_within_cell_type = True
    min_cells_per_pseudobulk = 2
    min_donors_per_group = 2
    group_column = "condition"
    comparison_mode = "condition_a_vs_b"
    condition_a = "treated"
    condition_b = "control"
    reference_level = None
    comparisons_tsv = None
    stratify_by = None
    covariates = None
    batch_columns = None
    de_mode = "modern"
    balance_groups = False
    balance_seed = 0
    gene_filter_scope = "contrast"
    feature_mapping_tsv = None
    feature_mapping_from_column = None
    feature_mapping_to_column = None
    feature_mapping_strip_version = False
    drop_unmapped_features = False
    balance_groups_explicit = False
    balance_seed_explicit = False
    gene_filter_scope_explicit = False
    repeated_measures = False
    allow_approximate_repeated_measures = False
    backend = "lightweight"
    allow_non_count_input = False
    write_pseudobulk_artifacts = True
    run_extractor = False
    extractor_out_dir = None
    extractor_signature_name = "contrast"
    extractor_score_mode = "auto"
    extractor_select = "top_k"
    extractor_top_k = 50
    extractor_quantile = 0.01
    extractor_min_score = 0.0
    extractor_normalize = "within_set_l1"
    extractor_padj_max = None
    extractor_pvalue_max = None
    extractor_min_abs_logfc = None
    extractor_emit_gmt = True
    extractor_gmt_split_signed = True
    extractor_gmt_topk_list = "50"
    extractor_gmt_min_genes = 1
    extractor_gmt_max_genes = 200


class ScrnaArgs(BulkArgs):
    modality = "scrna"
    sample_id_column = None
    sample_metadata_tsv = None
    cell_id_column = "cell_id"
    cell_metadata_tsv = ""
    donor_column = "donor_id"
    cell_type_column = "cell_type"
    min_cells_per_pseudobulk = 2



def _write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)



def _make_bulk_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    counts_path = tmp_path / "bulk_counts.tsv"
    meta_path = tmp_path / "bulk_meta.tsv"
    subject_path = tmp_path / "subject_meta.tsv"
    _write_tsv(
        counts_path,
        [
            {"sample_id": "S1", "G_UP": 120, "G_DOWN": 10, "G_FLAT": 50},
            {"sample_id": "S2", "G_UP": 110, "G_DOWN": 8, "G_FLAT": 48},
            {"sample_id": "S3", "G_UP": 130, "G_DOWN": 12, "G_FLAT": 51},
            {"sample_id": "S4", "G_UP": 10, "G_DOWN": 100, "G_FLAT": 52},
            {"sample_id": "S5", "G_UP": 12, "G_DOWN": 105, "G_FLAT": 47},
            {"sample_id": "S6", "G_UP": 9, "G_DOWN": 95, "G_FLAT": 49},
        ],
        ["sample_id", "G_UP", "G_DOWN", "G_FLAT"],
    )
    _write_tsv(
        meta_path,
        [
            {"sample_id": "S1", "condition": "treated", "tissue": "lung", "subject_id": "D1"},
            {"sample_id": "S2", "condition": "treated", "tissue": "lung", "subject_id": "D2"},
            {"sample_id": "S3", "condition": "treated", "tissue": "lung", "subject_id": "D3"},
            {"sample_id": "S4", "condition": "control", "tissue": "lung", "subject_id": "D4"},
            {"sample_id": "S5", "condition": "control", "tissue": "lung", "subject_id": "D5"},
            {"sample_id": "S6", "condition": "control", "tissue": "lung", "subject_id": "D6"},
        ],
        ["sample_id", "condition", "tissue", "subject_id"],
    )
    _write_tsv(
        subject_path,
        [
            {"subject_id": "D1", "sex": "F"},
            {"subject_id": "D2", "sex": "M"},
            {"subject_id": "D3", "sex": "F"},
            {"subject_id": "D4", "sex": "F"},
            {"subject_id": "D5", "sex": "M"},
            {"subject_id": "D6", "sex": "M"},
        ],
        ["subject_id", "sex"],
    )
    return counts_path, meta_path, subject_path


def _make_bulk_inputs_unbalanced(tmp_path: Path) -> tuple[Path, Path]:
    counts_path = tmp_path / "bulk_counts_unbalanced.tsv"
    meta_path = tmp_path / "bulk_meta_unbalanced.tsv"
    _write_tsv(
        counts_path,
        [
            {"sample_id": "S1", "G_UP": 120, "G_DOWN": 10, "G_FLAT": 50},
            {"sample_id": "S2", "G_UP": 118, "G_DOWN": 12, "G_FLAT": 48},
            {"sample_id": "S3", "G_UP": 125, "G_DOWN": 11, "G_FLAT": 49},
            {"sample_id": "S4", "G_UP": 122, "G_DOWN": 9, "G_FLAT": 51},
            {"sample_id": "S5", "G_UP": 12, "G_DOWN": 110, "G_FLAT": 50},
            {"sample_id": "S6", "G_UP": 10, "G_DOWN": 105, "G_FLAT": 47},
        ],
        ["sample_id", "G_UP", "G_DOWN", "G_FLAT"],
    )
    _write_tsv(
        meta_path,
        [
            {"sample_id": "S1", "condition": "treated", "tissue": "lung"},
            {"sample_id": "S2", "condition": "treated", "tissue": "lung"},
            {"sample_id": "S3", "condition": "treated", "tissue": "lung"},
            {"sample_id": "S4", "condition": "treated", "tissue": "lung"},
            {"sample_id": "S5", "condition": "control", "tissue": "lung"},
            {"sample_id": "S6", "condition": "control", "tissue": "lung"},
        ],
        ["sample_id", "condition", "tissue"],
    )
    return counts_path, meta_path


def _make_bulk_gene_by_sample_with_mapping_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    counts_path = tmp_path / "bulk_gene_by_sample.tsv"
    meta_path = tmp_path / "bulk_gene_by_sample_meta.tsv"
    mapping_path = tmp_path / "feature_map.tsv"
    _write_tsv(
        counts_path,
        [
            {"gene_id": "ENSG1.1", "gene_symbol": "OLD1", "S1": 50, "S2": 48, "S3": 5, "S4": 4},
            {"gene_id": "ENSG2.2", "gene_symbol": "OLD2", "S1": 4, "S2": 5, "S3": 55, "S4": 60},
            {"gene_id": "ENSG3.3", "gene_symbol": "OLD3", "S1": 20, "S2": 20, "S3": 20, "S4": 20},
        ],
        ["gene_id", "gene_symbol", "S1", "S2", "S3", "S4"],
    )
    _write_tsv(
        meta_path,
        [
            {"sample_id": "S1", "condition": "treated", "tissue": "lung"},
            {"sample_id": "S2", "condition": "treated", "tissue": "lung"},
            {"sample_id": "S3", "condition": "control", "tissue": "lung"},
            {"sample_id": "S4", "condition": "control", "tissue": "lung"},
        ],
        ["sample_id", "condition", "tissue"],
    )
    _write_tsv(
        mapping_path,
        [
            {"ensembl": "ENSG1", "symbol": "NEW1"},
            {"ensembl": "ENSG2", "symbol": "NEW2"},
        ],
        ["ensembl", "symbol"],
    )
    return counts_path, meta_path, mapping_path


def _make_scrna_inputs(tmp_path: Path) -> tuple[Path, Path]:
    counts_path = tmp_path / "scrna_counts.tsv"
    meta_path = tmp_path / "scrna_meta.tsv"
    _write_tsv(
        counts_path,
        [
            {"cell_id": "C1", "G_UP": 20, "G_DOWN": 1, "G_FLAT": 5},
            {"cell_id": "C2", "G_UP": 18, "G_DOWN": 2, "G_FLAT": 4},
            {"cell_id": "C3", "G_UP": 16, "G_DOWN": 1, "G_FLAT": 6},
            {"cell_id": "C4", "G_UP": 17, "G_DOWN": 2, "G_FLAT": 5},
            {"cell_id": "C5", "G_UP": 1, "G_DOWN": 20, "G_FLAT": 5},
            {"cell_id": "C6", "G_UP": 2, "G_DOWN": 18, "G_FLAT": 4},
            {"cell_id": "C7", "G_UP": 1, "G_DOWN": 16, "G_FLAT": 6},
            {"cell_id": "C8", "G_UP": 2, "G_DOWN": 17, "G_FLAT": 5},
        ],
        ["cell_id", "G_UP", "G_DOWN", "G_FLAT"],
    )
    _write_tsv(
        meta_path,
        [
            {"cell_id": "C1", "donor_id": "D1", "condition": "treated", "cell_type": "Tcell"},
            {"cell_id": "C2", "donor_id": "D1", "condition": "treated", "cell_type": "Tcell"},
            {"cell_id": "C3", "donor_id": "D2", "condition": "treated", "cell_type": "Tcell"},
            {"cell_id": "C4", "donor_id": "D2", "condition": "treated", "cell_type": "Tcell"},
            {"cell_id": "C5", "donor_id": "D3", "condition": "control", "cell_type": "Tcell"},
            {"cell_id": "C6", "donor_id": "D3", "condition": "control", "cell_type": "Tcell"},
            {"cell_id": "C7", "donor_id": "D4", "condition": "control", "cell_type": "Tcell"},
            {"cell_id": "C8", "donor_id": "D4", "condition": "control", "cell_type": "Tcell"},
        ],
        ["cell_id", "donor_id", "condition", "cell_type"],
    )
    return counts_path, meta_path



def test_build_cli_comparisons_reference_level_with_stratum():
    rows = [
        {"sample_id": "S1", "age_bin": "20-29", "tissue": "lung"},
        {"sample_id": "S2", "age_bin": "30-39", "tissue": "lung"},
        {"sample_id": "S3", "age_bin": "40-49", "tissue": "lung"},
        {"sample_id": "S4", "age_bin": "20-29", "tissue": "heart"},
        {"sample_id": "S5", "age_bin": "30-39", "tissue": "heart"},
    ]
    specs = build_cli_comparisons(
        rows,
        comparison_mode="reference_level",
        group_column="age_bin",
        condition_a=None,
        condition_b=None,
        reference_level="20-29",
        stratify_by=["tissue"],
    )
    assert {spec.comparison_id for spec in specs} == {
        "age_bin=30-39_vs_20-29__tissue=heart",
        "age_bin=30-39_vs_20-29__tissue=lung",
        "age_bin=40-49_vs_20-29__tissue=lung",
    }



def test_build_pseudobulk_aggregates_cells_by_donor_and_cell_type():
    matrix = MatrixData(
        sample_ids=["C1", "C2", "C3", "C4"],
        feature_ids=["G1", "G2"],
        counts=np.asarray([[1, 2], [3, 4], [10, 20], [30, 40]], dtype=float),
        gene_symbols=["G1", "G2"],
        matrix_orientation="sample_by_gene",
        sample_id_column="cell_id",
        feature_id_column="gene_id",
        source_path="inline",
    )
    rows = [
        {"sample_id": "C1", "donor_id": "D1", "cell_type": "A", "condition": "treated"},
        {"sample_id": "C2", "donor_id": "D1", "cell_type": "A", "condition": "treated"},
        {"sample_id": "C3", "donor_id": "D2", "cell_type": "A", "condition": "control"},
        {"sample_id": "C4", "donor_id": "D2", "cell_type": "A", "condition": "control"},
    ]
    result = build_pseudobulk(
        matrix,
        cell_metadata_rows=rows,
        cell_id_column="sample_id",
        donor_column="donor_id",
        group_columns=["cell_type", "condition"],
        min_cells_per_pseudobulk=2,
    )
    assert result.n_pseudobulks == 2
    assert result.matrix.sample_ids == [
        "donor_id=D1__cell_type=A__condition=treated",
        "donor_id=D2__cell_type=A__condition=control",
    ]
    assert result.matrix.counts.tolist() == [[4.0, 6.0], [40.0, 60.0]]


def test_r_backend_scripts_drop_gene_symbol_metadata_column(tmp_path: Path):
    limma_script = r_limma_voom.write_script(
        script_path=tmp_path / "run_limma_voom.R",
        counts_tsv=tmp_path / "counts.tsv",
        metadata_tsv=tmp_path / "meta.tsv",
        comparisons_tsv=tmp_path / "comparisons.tsv",
        selected_samples_tsv=tmp_path / "selected.tsv",
        output_tsv=tmp_path / "deg.tsv",
        covariates_csv="",
        batch_columns_csv="",
        gene_filter_scope="contrast",
    )
    dream_script = r_dream.write_script(
        script_path=tmp_path / "run_dream.R",
        counts_tsv=tmp_path / "counts.tsv",
        metadata_tsv=tmp_path / "meta.tsv",
        comparisons_tsv=tmp_path / "comparisons.tsv",
        selected_samples_tsv=tmp_path / "selected.tsv",
        output_tsv=tmp_path / "deg.tsv",
        covariates_csv="",
        batch_columns_csv="",
        random_effect_column="subject_id",
        gene_filter_scope="contrast",
    )
    for script_path in (limma_script, dream_script):
        text = script_path.read_text(encoding="utf-8")
        assert 'gene_symbols <- if ("gene_symbol" %in% colnames(counts))' in text
        assert 'count_cols <- setdiff(colnames(counts), c(colnames(counts)[1], "gene_symbol"))' in text
        assert 'storage.mode(count_mat) <- "numeric"' in text



def test_rna_de_prepare_applies_feature_mapping_before_de(tmp_path: Path):
    counts_path, meta_path, mapping_path = _make_bulk_gene_by_sample_with_mapping_inputs(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_feature_mapping"),
        matrix_orientation="gene_by_sample",
        feature_id_column="gene_id",
        matrix_gene_symbol_column="gene_symbol",
        feature_mapping_tsv=str(mapping_path),
        feature_mapping_from_column="ensembl",
        feature_mapping_to_column="symbol",
        feature_mapping_strip_version=True,
        drop_unmapped_features=True,
        backend="lightweight",
    )
    result = run_rna_de_prepare(args)
    assert result["n_comparisons"] == 1
    rows = list(csv.DictReader((Path(args.out_dir) / "deg_long.tsv").open("r", encoding="utf-8"), delimiter="	"))
    genes = {row["gene_symbol"] for row in rows}
    assert genes == {"NEW1", "NEW2"}
    summary = json.loads((Path(args.out_dir) / "prepare_summary.json").read_text(encoding="utf-8"))
    feature_summary = summary["feature_preprocessing"]
    assert feature_summary["mapping_applied"] is True
    assert feature_summary["n_features_input"] == 3
    assert feature_summary["n_features_output"] == 2
    assert feature_summary["n_features_mapped"] == 2
    assert feature_summary["n_features_unmapped"] == 1


def test_rna_de_prepare_bulk_writes_deg_long_and_summary(tmp_path: Path):
    counts_path, meta_path, subject_path = _make_bulk_inputs(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        subject_metadata_tsv=str(subject_path),
        subject_join_sample_column="subject_id",
        subject_join_metadata_column="subject_id",
        out_dir=str(tmp_path / "bulk_out"),
        covariates="sex",
    )
    result = run_rna_de_prepare(args)
    assert result["n_comparisons"] == 1
    deg_long = Path(args.out_dir) / "deg_long.tsv"
    rows = list(csv.DictReader(deg_long.open("r", encoding="utf-8"), delimiter="\t"))
    assert rows
    top = max(rows, key=lambda row: abs(float(row["stat"])))
    assert top["gene_id"] in {"G_UP", "G_DOWN"}
    summary = json.loads((Path(args.out_dir) / "prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["workflow"] == "rna_de_prepare"
    assert summary["metadata_join"]["n_joined"] == 6



def test_rna_de_prepare_harmonizome_mode_defaults_gene_filter_scope_to_stratum(tmp_path: Path):
    counts_path, meta_path = _make_bulk_inputs_unbalanced(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_harmonizome_filter_scope"),
        de_mode="harmonizome",
        backend="lightweight",
    )
    run_rna_de_prepare(args)
    audit_rows = list(
        csv.DictReader((Path(args.out_dir) / "comparison_audit.tsv").open("r", encoding="utf-8"), delimiter="	")
    )
    audit = audit_rows[0]
    assert audit["gene_filter_scope"] == "stratum"
    assert audit["n_filter_scope_samples"] == "6"
    summary = json.loads((Path(args.out_dir) / "prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["gene_filter_scope"] == "stratum"


def test_rna_de_prepare_harmonizome_mode_balances_bulk_contrast_and_writes_audit(tmp_path: Path):
    counts_path, meta_path = _make_bulk_inputs_unbalanced(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_harmonizome"),
        de_mode="harmonizome",
        backend="lightweight",
    )
    result = run_rna_de_prepare(args)
    assert result["n_comparisons"] == 1
    audit_rows = list(
        csv.DictReader((Path(args.out_dir) / "comparison_audit.tsv").open("r", encoding="utf-8"), delimiter="\t")
    )
    assert len(audit_rows) == 1
    audit = audit_rows[0]
    assert audit["de_mode"] == "harmonizome"
    assert audit["balance_requested"] == "True"
    assert audit["balance_sampler"] == "harmonizome_random_state"
    assert audit["n_group_a_pre_balance"] == "4"
    assert audit["n_group_b_pre_balance"] == "2"
    assert audit["n_group_a_post_balance"] == "2"
    assert audit["n_group_b_post_balance"] == "2"
    selected_rows = [
        row
        for row in csv.DictReader((Path(args.out_dir) / "comparison_selected_samples.tsv").open("r", encoding="utf-8"), delimiter="\t")
        if row["comparison_id"] == "condition=treated_vs_control"
    ]
    assert len(selected_rows) == 4
    treated_ids = [row["sample_id"] for row in selected_rows if row["group_label"] == "treated"]
    control_ids = [row["sample_id"] for row in selected_rows if row["group_label"] == "control"]
    assert treated_ids == ["S4", "S3"]
    assert control_ids == ["S5", "S6"]
    summary = json.loads((Path(args.out_dir) / "prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["de_mode"] == "harmonizome"
    assert summary["balance_groups"] is True
    assert summary["balance_sampler"] == "harmonizome_random_state"
    assert summary["balance_seed"] == 1


def test_rna_de_prepare_modern_mode_respects_explicit_gene_filter_scope(tmp_path: Path):
    counts_path, meta_path = _make_bulk_inputs_unbalanced(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_modern_filter_scope"),
        backend="lightweight",
        gene_filter_scope="stratum",
        gene_filter_scope_explicit=True,
    )
    run_rna_de_prepare(args)
    summary = json.loads((Path(args.out_dir) / "prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["gene_filter_scope"] == "stratum"


def test_rna_de_prepare_modern_mode_keeps_all_samples_on_unbalanced_bulk(tmp_path: Path):
    counts_path, meta_path = _make_bulk_inputs_unbalanced(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_modern"),
        backend="lightweight",
    )
    run_rna_de_prepare(args)
    audit_rows = list(
        csv.DictReader((Path(args.out_dir) / "comparison_audit.tsv").open("r", encoding="utf-8"), delimiter="\t")
    )
    audit = audit_rows[0]
    assert audit["de_mode"] == "modern"
    assert audit["balance_requested"] == "False"
    assert audit["balance_applied"] == "False"
    assert audit["balance_sampler"] == "none"
    assert audit["n_group_a_pre_balance"] == "4"
    assert audit["n_group_b_pre_balance"] == "2"
    assert audit["n_group_a_post_balance"] == "4"
    assert audit["n_group_b_post_balance"] == "2"


def test_rna_de_prepare_harmonizome_mode_accepts_and_records_covariates(tmp_path: Path):
    counts_path, meta_path, _subject_path = _make_bulk_inputs(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_harmonizome_covariates"),
        de_mode="harmonizome",
        covariates="tissue",
        backend="lightweight",
    )
    result = run_rna_de_prepare(args)
    assert result["n_comparisons"] == 1
    summary = json.loads((Path(args.out_dir) / "prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["de_mode"] == "harmonizome"
    assert summary["covariates_used"] == ["tissue"]
    assert summary["harmonizome_covariate_mode"] == "explicit"
    audit_rows = list(
        csv.DictReader((Path(args.out_dir) / "comparison_audit.tsv").open("r", encoding="utf-8"), delimiter="\t")
    )
    assert audit_rows[0]["covariates_used"] == "tissue"
    assert audit_rows[0]["harmonizome_covariate_mode"] == "explicit"


def test_rna_de_prepare_harmonizome_mode_warns_when_covariates_missing(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    counts_path, meta_path = _make_bulk_inputs_unbalanced(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_harmonizome_no_covariates"),
        de_mode="harmonizome",
        backend="lightweight",
    )
    result = run_rna_de_prepare(args)
    assert result["n_comparisons"] == 1
    captured = capsys.readouterr()
    assert "without explicit covariates" in captured.err
    summary = json.loads((Path(args.out_dir) / "prepare_summary.json").read_text(encoding="utf-8"))
    assert any("without explicit covariates" in warning for warning in summary["warnings"])
    assert summary["harmonizome_covariate_mode"] == "none"


def test_rna_de_prepare_modern_mode_warns_on_strong_imbalance(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    counts_path, meta_path = _make_bulk_inputs_unbalanced(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_modern_warn_imbalance"),
        backend="lightweight",
    )
    run_rna_de_prepare(args)
    captured = capsys.readouterr()
    assert "strongly imbalanced before fit" in captured.err
    summary = json.loads((Path(args.out_dir) / "prepare_summary.json").read_text(encoding="utf-8"))
    assert any("strongly imbalanced before fit" in warning for warning in summary["warnings"])
    audit_rows = list(csv.DictReader((Path(args.out_dir) / "comparison_audit.tsv").open("r", encoding="utf-8"), delimiter="	"))
    assert audit_rows[0]["imbalance_ratio_pre_balance"] == "2.000"
    assert audit_rows[0]["imbalance_ratio_post_balance"] == "2.000"


def test_rna_de_prepare_modern_bulk_stratified_warns_without_covariates(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    counts_path, meta_path = _make_bulk_inputs_unbalanced(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_modern_stratified_no_covariates"),
        backend="lightweight",
        stratify_by="tissue",
    )
    run_rna_de_prepare(args)
    captured = capsys.readouterr()
    assert "without explicit covariates while stratification is enabled" in captured.err


def test_rna_de_prepare_fails_when_requested_covariate_is_missing(tmp_path: Path):
    counts_path, meta_path = _make_bulk_inputs_unbalanced(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_missing_covariate"),
        covariates="SEX",
        backend="lightweight",
    )
    with pytest.raises(ValueError, match="Requested covariates were not found"):
        run_rna_de_prepare(args)


def test_rna_de_prepare_uses_group_column_from_comparisons_tsv(tmp_path: Path):
    counts_path, meta_path = _make_bulk_inputs_unbalanced(tmp_path)
    comparisons_path = tmp_path / "comparisons.tsv"
    _write_tsv(
        comparisons_path,
        [
            {
                "comparison_id": "condition=treated_vs_control",
                "group_column": "condition",
                "group_a": "treated",
                "group_b": "control",
                "tissue": "lung",
            }
        ],
        ["comparison_id", "group_column", "group_a", "group_b", "tissue"],
    )
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_from_comparisons"),
        group_column=None,
        comparison_mode=None,
        comparisons_tsv=str(comparisons_path),
        stratify_by="tissue",
        backend="lightweight",
    )
    result = run_rna_de_prepare(args)
    assert result["n_comparisons"] == 1
    audit_rows = list(
        csv.DictReader((Path(args.out_dir) / "comparison_audit.tsv").open("r", encoding="utf-8"), delimiter="\t")
    )
    assert audit_rows[0]["comparison_id"] == "condition=treated_vs_control"


def test_rna_de_prepare_scrna_writes_pseudobulk_and_deg_long(tmp_path: Path):
    counts_path, meta_path = _make_scrna_inputs(tmp_path)
    args = ScrnaArgs(
        counts_tsv=str(counts_path),
        cell_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "scrna_out"),
    )
    result = run_rna_de_prepare(args)
    assert result["n_comparisons"] == 1
    assert (Path(args.out_dir) / "pseudobulk_matrix.tsv").exists()
    assert (Path(args.out_dir) / "pseudobulk_metadata.tsv").exists()
    rows = list(csv.DictReader((Path(args.out_dir) / "deg_long.tsv").open("r", encoding="utf-8"), delimiter="\t"))
    assert rows
    assert {row["comparison_id"] for row in rows} == {"condition=treated_vs_control"}



def test_rna_de_prepare_repeated_measures_requires_explicit_approx_fallback(tmp_path: Path):
    counts_path, meta_path, _subject_path = _make_bulk_inputs(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_repeated"),
        repeated_measures=True,
        backend="lightweight",
    )
    with pytest.raises(ValueError, match="Repeated-measures design requested"):
        run_rna_de_prepare(args)



def test_rna_de_prepare_can_call_internal_extractor(tmp_path: Path):
    counts_path, meta_path, _subject_path = _make_bulk_inputs(tmp_path)
    args = BulkArgs(
        counts_tsv=str(counts_path),
        sample_metadata_tsv=str(meta_path),
        out_dir=str(tmp_path / "bulk_with_extractor"),
        run_extractor=True,
        extractor_padj_max=0.2,
        extractor_min_abs_logfc=0.5,
        extractor_gmt_min_genes=1,
        extractor_gmt_max_genes=20,
        extractor_gmt_topk_list="2",
    )
    result = run_rna_de_prepare(args)
    extractor_out = Path(result["extractor_out_dir"])
    assert (extractor_out / "manifest.tsv").exists()
    assert (extractor_out / "genesets.gmt").exists()
    workflow_graph = Path(args.out_dir) / "deg_long.provenance_graph.json"
    assert workflow_graph.exists()

    with (extractor_out / "manifest.tsv").open("r", encoding="utf-8") as fh:
        first_row = next(csv.DictReader(fh, delimiter="\t"))
    provenance = json.loads((extractor_out / str(first_row["path"]) / "geneset.provenance.json").read_text(encoding="utf-8"))
    graph = next(iter(provenance.values()))
    analysis_nodes = [node for node in graph["nodes"] if node.get("type") == "AnalysisType"]
    assert {node["id"].split(":")[1] for node in analysis_nodes} >= {"rna_de_prepare", "rna_deg_multi"}
    deg_long_nodes = [node for node in graph["nodes"] if node.get("type") == "File" and node.get("name") == "deg_long.tsv"]
    assert deg_long_nodes
    deg_long_id = deg_long_nodes[0]["id"]
    assert any(edge["source"].startswith("analysis:rna_de_prepare:") and edge["target"] == deg_long_id for edge in graph["edges"])
    assert any(edge["source"] == deg_long_id and edge["target"].startswith("analysis:rna_deg_multi:") for edge in graph["edges"])



def test_rna_de_prepare_cli_entrypoint_bulk(tmp_path: Path):
    counts_path, meta_path, _subject_path = _make_bulk_inputs(tmp_path)
    out_dir = tmp_path / "bulk_cli"
    code = main(
        [
            "workflows",
            "rna_de_prepare",
            "--modality",
            "bulk",
            "--counts_tsv",
            str(counts_path),
            "--sample_metadata_tsv",
            str(meta_path),
            "--sample_id_column",
            "sample_id",
            "--feature_id_column",
            "gene_id",
            "--group_column",
            "condition",
            "--comparison_mode",
            "condition_a_vs_b",
            "--condition_a",
            "treated",
            "--condition_b",
            "control",
            "--backend",
            "lightweight",
            "--out_dir",
            str(out_dir),
            "--organism",
            "human",
            "--genome_build",
            "hg38",
        ]
    )
    assert code == 0
    assert (out_dir / "deg_long.tsv").exists()
