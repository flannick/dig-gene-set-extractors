from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from types import SimpleNamespace

from geneset_extractors.cli import main
from geneset_extractors.converters import rna_deg_multi
from geneset_extractors.workflows.rna_prepare_bulk_tissue import run as run_rna_prepare_bulk_tissue
from geneset_extractors.workflows.rna_de_prepare import run as run_rna_de_prepare


class Args(SimpleNamespace):
    counts_gct = ""
    sample_metadata_tsv = ""
    subject_metadata_tsv = ""
    tissue_label = "Adipose - Subcutaneous"
    out_dir = ""
    sample_id_column = "SAMPID"
    subject_id_column_sample = "SUBJID"
    subject_id_column_subject = "SUBJID"
    age_column = "AGE"
    sex_column = "SEX"
    primary_tissue_column = "SMTS"
    detailed_tissue_column = "SMTSD"
    reference_age_bin = "20-29"
    age_bins = "20-29,30-39,40-49,50-59,60-69,70-79"
    min_samples_per_group = 2


class DeArgs(SimpleNamespace):
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
    group_column = None
    comparison_mode = None
    condition_a = None
    condition_b = None
    reference_level = None
    comparisons_tsv = ""
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
    run_extractor = True
    extractor_out_dir = None
    extractor_signature_name = "AB1"
    extractor_score_mode = "auto"
    extractor_select = "top_k"
    extractor_top_k = 50
    extractor_quantile = 0.01
    extractor_min_score = 0.0
    extractor_normalize = "within_set_l1"
    extractor_padj_max = 0.2
    extractor_pvalue_max = None
    extractor_min_abs_logfc = 0.5
    extractor_emit_gmt = True
    extractor_gmt_split_signed = True
    extractor_gmt_topk_list = "2"
    extractor_gmt_min_genes = 1
    extractor_gmt_max_genes = 20


class DegMultiArgs(SimpleNamespace):
    deg_tsv = ""
    comparison_column = "comparison_id"
    out_dir = ""
    organism = "human"
    genome_build = "hg38"
    signature_name = "AB1"
    gene_id_column = "gene_id"
    gene_symbol_column = "gene_symbol"
    stat_column = "stat"
    logfc_column = "logFC"
    padj_column = "padj"
    pvalue_column = "pvalue"
    score_column = None
    score_mode = "auto"
    padj_max = None
    pvalue_max = None
    min_abs_logfc = None
    postprocess_mode = "legacy"
    duplicate_gene_policy = "max_abs"
    neglog10p_cap = 50.0
    neglog10p_eps = 1e-300
    exclude_gene_regex = None
    disable_default_excludes = False
    gtf = None
    gtf_gene_id_field = "gene_id"
    gtf_source = None
    select = "top_k"
    top_k = 50
    quantile = 0.01
    min_score = 0.0
    normalize = "within_set_l1"
    emit_full = True
    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = False
    gmt_biotype_allowlist = ""
    gmt_min_genes = 1
    gmt_max_genes = 20
    gmt_topk_list = "2"
    gmt_mass_list = ""
    gmt_split_signed = True
    gmt_emit_abs = False
    gmt_source = "full"
    emit_small_gene_sets = True


def _write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_gct(path: Path) -> None:
    rows = [
        ["GENE1", "GeneOne", "500", "520", "510", "5", "6", "4"],
        ["GENE2", "GeneTwo", "7", "6", "5", "450", "470", "460"],
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("#1.2\n")
        handle.write("2\t6\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Name", "Description", "GTEX-01-A", "GTEX-02-A", "GTEX-03-A", "GTEX-04-A", "GTEX-05-A", "GTEX-06-A"])
        writer.writerows(rows)


def _make_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    counts_gct = tmp_path / "gene_reads.gct.gz"
    sample_metadata_tsv = tmp_path / "sample_attributes.tsv"
    subject_metadata_tsv = tmp_path / "subject_phenotypes.tsv"
    _write_gct(counts_gct)
    _write_tsv(
        sample_metadata_tsv,
        [
            {"SAMPID": "GTEX-01-A", "SMTS": "Adipose Tissue", "SMTSD": "Adipose - Subcutaneous"},
            {"SAMPID": "GTEX-02-A", "SMTS": "Adipose Tissue", "SMTSD": "Adipose - Subcutaneous"},
            {"SAMPID": "GTEX-03-A", "SMTS": "Adipose Tissue", "SMTSD": "Adipose - Subcutaneous"},
            {"SAMPID": "GTEX-04-A", "SMTS": "Adipose Tissue", "SMTSD": "Adipose - Subcutaneous"},
            {"SAMPID": "GTEX-05-A", "SMTS": "Adipose Tissue", "SMTSD": "Adipose - Subcutaneous"},
            {"SAMPID": "GTEX-06-A", "SMTS": "Adipose Tissue", "SMTSD": "Adipose - Subcutaneous"},
            {"SAMPID": "GTEX-99-A", "SMTS": "Adipose Tissue", "SMTSD": "Adipose - Subcutaneous"},
        ],
        ["SAMPID", "SMTS", "SMTSD"],
    )
    _write_tsv(
        subject_metadata_tsv,
        [
            {"SUBJID": "GTEX-01", "AGE": "1", "SEX": "1"},
            {"SUBJID": "GTEX-02", "AGE": "1", "SEX": "2"},
            {"SUBJID": "GTEX-03", "AGE": "2", "SEX": "1"},
            {"SUBJID": "GTEX-04", "AGE": "2", "SEX": "2"},
            {"SUBJID": "GTEX-05", "AGE": "1", "SEX": "1"},
            {"SUBJID": "GTEX-06", "AGE": "2", "SEX": "2"},
        ],
        ["SUBJID", "AGE", "SEX"],
    )
    return counts_gct, sample_metadata_tsv, subject_metadata_tsv


def test_rna_prepare_bulk_tissue_writes_prepared_bundle(tmp_path: Path):
    counts_gct, sample_metadata_tsv, subject_metadata_tsv = _make_inputs(tmp_path)
    out_dir = tmp_path / "prepared"
    args = Args(
        counts_gct=str(counts_gct),
        sample_metadata_tsv=str(sample_metadata_tsv),
        subject_metadata_tsv=str(subject_metadata_tsv),
        out_dir=str(out_dir),
    )

    result = run_rna_prepare_bulk_tissue(args)

    assert result["workflow"] == "rna_prepare_bulk_tissue"
    assert result["n_samples_retained"] == 6
    assert result["n_comparisons"] == 1
    assert result["comparison_ids"] == ["age30_20"]
    assert Path(result["counts_provenance_graph_path"]).exists()

    counts_rows = list(csv.DictReader((out_dir / "tissue_counts.tsv").open("r", encoding="utf-8"), delimiter="\t"))
    assert counts_rows[0]["gene_id"] == "GENE1"
    assert counts_rows[0]["gene_symbol"] == "GeneOne"
    assert list(counts_rows[0].keys()) == ["gene_id", "gene_symbol", "GTEX-01-A", "GTEX-02-A", "GTEX-03-A", "GTEX-04-A", "GTEX-05-A", "GTEX-06-A"]

    sample_rows = list(csv.DictReader((out_dir / "sample_metadata.tsv").open("r", encoding="utf-8"), delimiter="\t"))
    assert sample_rows[0]["subject_id"] == "GTEX-01"
    assert sample_rows[0]["age_bin"] == "20-29"
    assert sample_rows[0]["SEX"] == "M"
    assert sample_rows[1]["SEX"] == "F"

    comparisons = list(csv.DictReader((out_dir / "comparisons.tsv").open("r", encoding="utf-8"), delimiter="\t"))
    assert comparisons == [
        {
            "comparison_id": "age30_20",
            "comparison_kind": "condition_a_vs_b",
            "group_column": "age_bin",
            "group_a": "30-39",
            "group_b": "20-29",
        }
    ]

    summary = json.loads((out_dir / "prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["age_bin_counts"] == {
        "20-29": 3,
        "30-39": 3,
        "40-49": 0,
        "50-59": 0,
        "60-69": 0,
        "70-79": 0,
    }
    assert summary["counts_provenance_graph_path"].endswith("tissue_counts.provenance_graph.json")
    assert "age30_20" in (out_dir / "naming_reference.md").read_text(encoding="utf-8")


def test_rna_prepare_bulk_tissue_cli_entrypoint(tmp_path: Path):
    counts_gct, sample_metadata_tsv, subject_metadata_tsv = _make_inputs(tmp_path)
    out_dir = tmp_path / "cli_prepared"

    rc = main(
        [
            "workflows",
            "rna_prepare_bulk_tissue",
            "--counts_gct",
            str(counts_gct),
            "--sample_metadata_tsv",
            str(sample_metadata_tsv),
            "--subject_metadata_tsv",
            str(subject_metadata_tsv),
            "--tissue_label",
            "Adipose - Subcutaneous",
            "--out_dir",
            str(out_dir),
        ]
    )

    assert rc == 0
    assert (out_dir / "tissue_counts.tsv").exists()
    assert (out_dir / "sample_metadata.tsv").exists()
    assert (out_dir / "comparisons.tsv").exists()


def test_rna_prepare_bulk_tissue_threads_into_final_provenance(tmp_path: Path):
    counts_gct, sample_metadata_tsv, subject_metadata_tsv = _make_inputs(tmp_path)
    prepared_dir = tmp_path / "prepared"
    prep_args = Args(
        counts_gct=str(counts_gct),
        sample_metadata_tsv=str(sample_metadata_tsv),
        subject_metadata_tsv=str(subject_metadata_tsv),
        out_dir=str(prepared_dir),
    )
    prep_result = run_rna_prepare_bulk_tissue(prep_args)
    assert Path(prep_result["counts_provenance_graph_path"]).exists()

    de_args = DeArgs(
        counts_tsv=str(prepared_dir / "tissue_counts.tsv"),
        sample_metadata_tsv=str(prepared_dir / "sample_metadata.tsv"),
        comparisons_tsv=str(prepared_dir / "comparisons.tsv"),
        out_dir=str(tmp_path / "de_workflow"),
        matrix_orientation="gene_by_sample",
        matrix_gene_symbol_column="gene_symbol",
        run_extractor=False,
    )
    result = run_rna_de_prepare(de_args)
    extractor_out = tmp_path / "extractor"
    deg_long_path = Path(de_args.out_dir) / "deg_long.tsv"
    assert deg_long_path.exists()
    deg_args = DegMultiArgs(
        deg_tsv=str(deg_long_path),
        out_dir=str(extractor_out),
    )
    deg_result = rna_deg_multi.run(deg_args)
    assert deg_result["n_groups"] == 1
    manifest_rows = list(csv.DictReader((extractor_out / "manifest.tsv").open("r", encoding="utf-8"), delimiter="\t"))
    provenance_path = extractor_out / str(manifest_rows[0]["path"]) / "geneset.provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    graph = provenance["c2m2"]
    analysis_nodes = [node for node in graph["nodes"] if node.get("type") == "AnalysisType"]
    assert {node["id"].split(":")[1] for node in analysis_nodes} >= {
        "rna_prepare_bulk_tissue",
        "rna_de_prepare",
        "rna_deg_multi",
    }
    tissue_counts_id = next(
        node["id"]
        for node in graph["nodes"]
        if node.get("type") == "File" and node.get("name") == "tissue_counts.tsv"
    )
    deg_long_id = next(
        node["id"]
        for node in graph["nodes"]
        if node.get("type") == "File" and node.get("name") == "deg_long.tsv"
    )
    assert any(edge["source"].startswith("analysis:rna_prepare_bulk_tissue:") and edge["target"] == tissue_counts_id for edge in graph["edges"])
    assert any(edge["source"] == tissue_counts_id and edge["target"].startswith("analysis:rna_de_prepare:") for edge in graph["edges"])
    assert any(edge["source"].startswith("analysis:rna_de_prepare:") and edge["target"] == deg_long_id for edge in graph["edges"])
    assert any(edge["source"] == deg_long_id and edge["target"].startswith("analysis:rna_deg_multi:") for edge in graph["edges"])
