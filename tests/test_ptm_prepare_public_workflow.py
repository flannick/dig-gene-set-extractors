import csv
import json
import shutil
from pathlib import Path

from geneset_extractors.converters import ptm_site_matrix
from geneset_extractors.extractors.proteomics.public_prepare import run_public_prepare
from geneset_extractors.workflows.ptm_prepare_reference_bundle import run as run_ptm_prepare_reference_bundle


class MatrixArgs:
    ptm_matrix_tsv = ""
    sample_metadata_tsv = ""
    protein_matrix_tsv = ""
    out_dir = ""
    organism = "human"
    genome_build = "human"
    signature_name = "prepared_public_ptm"
    dataset_label = "prepared_public_ptm_dataset"
    ptm_type = "phospho"
    matrix_format = "wide_sites_by_sample"
    sample_id_column = "sample_id"
    group_column = "group"
    condition_column = "condition"
    study_contrast = "condition_a_vs_b"
    condition_a = "case"
    condition_b = "control"
    min_samples_per_condition = 2
    effect_metric = "welch_t"
    missing_value_policy = "min_present"
    min_present_per_condition = 2
    site_id_column = None
    site_group_column = None
    gene_id_column = None
    gene_symbol_column = None
    protein_accession_column = None
    residue_column = None
    position_column = None
    score_column = None
    stat_column = "stat"
    logfc_column = "log2fc"
    padj_column = "padj"
    pvalue_column = "pvalue"
    localization_prob_column = None
    peptide_count_column = None
    protein_logfc_column = "protein_log2fc"
    protein_stat_column = "protein_stat"
    score_mode = "auto"
    score_transform = "signed"
    protein_adjustment = "subtract"
    protein_adjustment_run_mode = "compare_if_protein"
    protein_adjustment_lambda = 1.0
    confidence_weight_mode = "combined"
    min_localization_prob = 0.75
    site_dup_policy = "highest_confidence"
    gene_aggregation = "signed_topk_mean"
    gene_topk_sites = 3
    emit_gene_topk_site_comparison = False
    gene_topk_site_compare_to = 1
    ambiguous_gene_policy = "drop"
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    use_reference_bundle = False
    site_alias_resource_id = None
    site_ubiquity_resource_id = None
    select = "top_k"
    top_k = 200
    quantile = 0.01
    min_score = 0.0
    normalize = "within_set_l1"
    emit_full = True
    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = True
    gmt_biotype_allowlist = ""
    gmt_min_genes = 1
    gmt_max_genes = 20
    gmt_topk_list = "3"
    gmt_mass_list = ""
    gmt_split_signed = True
    gmt_format = "dig2col"
    emit_small_gene_sets = True
    neglog10p_cap = 50.0
    neglog10p_eps = 1e-300


class BundleArgs:
    sources_tsv = ""
    out_dir = ""
    organism = "human"
    ptm_type = "phospho"
    bundle_id = "toy_public_ptm_bundle_v1"


def _copy_resources(resources_dir: Path) -> None:
    resources_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile("tests/data/phosphosite_aliases_human_v1.tsv.gz", resources_dir / "phosphosite_aliases_human_v1.tsv.gz")
    shutil.copyfile("tests/data/phosphosite_ubiquity_human_v1.tsv.gz", resources_dir / "phosphosite_ubiquity_human_v1.tsv.gz")


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_ptm_prepare_public_cdap_files_end_to_end(tmp_path: Path):
    out_dir = tmp_path / "prepared_public"
    summary = run_public_prepare(
        input_mode="cdap_files",
        ptm_report_tsv="tests/data/toy_cdap_phosphosite.tmt11.tsv",
        protein_report_tsv="tests/data/toy_cdap_proteome.tmt11.tsv",
        sample_design_tsv="tests/data/toy_cdap.sample.txt",
        sample_annotations_tsv="tests/data/toy_ptm_public_sample_annotations.tsv",
        pdc_manifest_tsv=None,
        source_dir=None,
        out_dir=out_dir,
        organism="human",
        ptm_type="phospho",
        study_id="PTM_STUDY_1",
        study_label="Toy PTM Public Study",
    )

    assert summary["workflow"] == "ptm_prepare_public"
    assert summary["ptm_rows_emitted"] == 3
    assert summary["single_resolved_sites"] == 1
    assert summary["grouped_sites"] == 1
    assert summary["unparsed_fallback_sites"] == 1
    assert summary["n_samples"] == 4
    assert summary["n_samples_with_condition_labels"] == 4

    ptm_rows = _read_rows(out_dir / "ptm_matrix.tsv")
    assert ptm_rows[0]["site_id"] == "NP_000001|S|123|phospho"
    assert ptm_rows[1]["site_id"] == ""
    assert ptm_rows[1]["site_group_id"] == "NP_000002|S50;T52|phospho"
    assert ptm_rows[2]["site_parser_status"] == "unparsed_fallback"
    assert {"Tumor-A", "Tumor-B", "Normal-A", "Normal-B"}.issubset(ptm_rows[0].keys())
    assert "POOL-REF" not in ptm_rows[0]

    protein_rows = _read_rows(out_dir / "protein_matrix.tsv")
    assert len(protein_rows) == 3
    assert protein_rows[0]["Tumor-A"] == "0.8"
    assert "Tumor-A Unshared Log Ratio" not in protein_rows[0]

    sample_rows = _read_rows(out_dir / "sample_metadata.tsv")
    assert [row["sample_id"] for row in sample_rows] == ["Tumor-A", "Tumor-B", "Normal-A", "Normal-B"]
    assert sample_rows[0]["reporter_ion"] == "126"
    assert sample_rows[1]["condition"] == "case"
    assert sample_rows[3]["condition"] == "control"
    assert sample_rows[0]["study_label"] == "Toy PTM Public Study"
    assert sample_rows[0]["cohort"] == "cohortA"

    sample_id_map = _read_rows(out_dir / "sample_id_map.tsv")
    pool_row = next(row for row in sample_id_map if row["sample_id_raw"] == "POOL-REF")
    assert pool_row["drop_status"] == "dropped_pool_reference"

    site_map = _read_rows(out_dir / "site_id_map.tsv")
    assert {row["site_parser_status"] for row in site_map} == {"single_site", "site_group", "unparsed_fallback"}

    prepared = json.loads((out_dir / "prepare_summary.json").read_text(encoding="utf-8"))
    assert prepared["outputs"]["ptm_matrix_tsv"].endswith("ptm_matrix.tsv")
    assert prepared["source_files"]["ptm_report_tsv"]["path"].endswith("toy_cdap_phosphosite.tmt11.tsv")


def test_ptm_prepare_public_pdc_manifest_and_downstream_handoffs(tmp_path: Path):
    out_dir = tmp_path / "prepared_from_manifest"
    summary = run_public_prepare(
        input_mode="pdc_manifest",
        ptm_report_tsv=None,
        protein_report_tsv=None,
        sample_design_tsv=None,
        sample_annotations_tsv="tests/data/toy_ptm_public_sample_annotations.tsv",
        pdc_manifest_tsv="tests/data/toy_pdc_manifest.tsv",
        source_dir="tests/data",
        out_dir=out_dir,
        organism="human",
        ptm_type="phospho",
        study_id="PTM_STUDY_1",
        study_label="Toy PTM Public Study",
    )
    assert summary["n_samples"] == 4

    bundle_args = BundleArgs()
    bundle_args.sources_tsv = str(out_dir / "bundle_source_row.tsv")
    bundle_args.out_dir = str(tmp_path / "ptm_bundle_from_public")
    bundle_result = run_ptm_prepare_reference_bundle(bundle_args)
    assert bundle_result["n_canonical_sites"] == 3
    assert (Path(bundle_args.out_dir) / "phosphosite_aliases_human_v1.tsv.gz").exists()
    assert (Path(bundle_args.out_dir) / "phosphosite_ubiquity_human_v1.tsv.gz").exists()

    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)
    matrix_args = MatrixArgs()
    matrix_args.ptm_matrix_tsv = str(out_dir / "ptm_matrix.tsv")
    matrix_args.sample_metadata_tsv = str(out_dir / "sample_metadata.tsv")
    matrix_args.protein_matrix_tsv = str(out_dir / "protein_matrix.tsv")
    matrix_args.out_dir = str(tmp_path / "ptm_from_public_matrix")
    matrix_args.resources_dir = str(resources_dir)
    matrix_args.use_reference_bundle = True
    result = ptm_site_matrix.run(matrix_args)
    assert result["n_contrasts_emitted"] == 2
    assert result["grouped_output"]
    manifest_rows = _read_rows(Path(matrix_args.out_dir) / "manifest.tsv")
    assert {row["protein_adjustment"] for row in manifest_rows} == {"none", "subtract"}
    for row in manifest_rows:
        child_dir = Path(matrix_args.out_dir) / row["path"]
        assert (child_dir / "geneset.tsv").exists()
    assert (Path(matrix_args.out_dir) / "genesets.gmt").exists()


def test_ptm_prepare_public_assay_type_qc_warns_for_lysine_dominant(tmp_path: Path, capsys):
    out_dir = tmp_path / "prepared_lysine_warn"
    summary = run_public_prepare(
        input_mode="cdap_files",
        ptm_report_tsv="tests/data/toy_cdap_lysine_ptm.tmt11.tsv",
        protein_report_tsv=None,
        sample_design_tsv="tests/data/toy_cdap.sample.txt",
        sample_annotations_tsv="tests/data/toy_ptm_public_sample_annotations.tsv",
        pdc_manifest_tsv=None,
        source_dir=None,
        out_dir=out_dir,
        organism="human",
        ptm_type="phospho",
        study_id="K_STUDY",
        study_label="Lysine Dominant Study",
        assay_type_policy="warn",
        min_phospho_like_fraction=0.6,
        max_k_fraction=0.25,
    )
    captured = capsys.readouterr()
    assert "weakly phospho-like" in captured.err
    assert summary["assay_type_qc"]["status"] == "warn"
    assert summary["assay_type_qc"]["dominant_residue_family"] == "lysine_dominant"


def test_ptm_prepare_public_assay_type_qc_can_fail(tmp_path: Path):
    out_dir = tmp_path / "prepared_lysine_fail"
    try:
        run_public_prepare(
            input_mode="cdap_files",
            ptm_report_tsv="tests/data/toy_cdap_lysine_ptm.tmt11.tsv",
            protein_report_tsv=None,
            sample_design_tsv="tests/data/toy_cdap.sample.txt",
            sample_annotations_tsv="tests/data/toy_ptm_public_sample_annotations.tsv",
            pdc_manifest_tsv=None,
            source_dir=None,
            out_dir=out_dir,
            organism="human",
            ptm_type="phospho",
            study_id="K_STUDY",
            study_label="Lysine Dominant Study",
            assay_type_policy="fail",
            min_phospho_like_fraction=0.6,
            max_k_fraction=0.25,
        )
    except ValueError as exc:
        assert "assay-type QC failed" in str(exc)
    else:
        raise AssertionError("Expected assay-type QC fail policy to raise")
