import csv
import json
from pathlib import Path

from geneset_extractors.converters import splice_event_matrix
from geneset_extractors.extractors.splicing.public_prepare import run_public_prepare
from geneset_extractors.workflows.splice_prepare_reference_bundle import run as run_splice_prepare_reference_bundle


class MatrixArgs:
    psi_matrix_tsv = ""
    sample_metadata_tsv = ""
    event_metadata_tsv = ""
    coverage_matrix_tsv = None
    out_dir = ""
    organism = "human"
    genome_build = "human"
    signature_name = "prepared_public_splice"
    dataset_label = "prepared_public_splice_dataset"
    matrix_format = "wide_events_by_sample"
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
    tool_family = "tcga_spliceseq"
    event_id_column = None
    event_group_column = None
    event_type_column = None
    gene_id_column = None
    gene_symbol_column = None
    chrom_column = None
    start_column = None
    end_column = None
    strand_column = None
    score_column = None
    stat_column = "stat"
    delta_psi_column = "delta_psi"
    psi_column = None
    padj_column = "padj"
    pvalue_column = "pvalue"
    probability_column = None
    read_support_column = "read_support"
    novel_flag_column = None
    annotation_status_column = "annotation_status"
    score_mode = "auto"
    score_transform = "signed"
    confidence_weight_mode = "combined"
    min_probability = 0.8
    min_read_support = 0.0
    neglog10p_cap = 50.0
    neglog10p_eps = 1e-300
    event_dup_policy = "highest_confidence"
    gene_aggregation = "signed_topk_mean"
    gene_topk_events = 3
    ambiguous_gene_policy = "drop"
    impact_mode = "conservative"
    impact_min = 0.75
    impact_max = 1.35
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    use_reference_bundle = False
    event_alias_resource_id = None
    event_ubiquity_resource_id = None
    event_impact_resource_id = None
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
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 1
    gmt_max_genes = 20
    gmt_topk_list = "3"
    gmt_mass_list = ""
    gmt_split_signed = True
    gmt_format = "dig2col"
    emit_small_gene_sets = True
    cluster_stats_tsv = None


class BundleArgs:
    sources_tsv = ""
    out_dir = ""
    organism = "human"
    bundle_id = "toy_splice_public_bundle_v1"
    min_ref_read_support = 0.0



def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))



def test_splice_prepare_public_tcga_spliceseq_end_to_end(tmp_path: Path):
    out_dir = tmp_path / "prepared_public"
    summary = run_public_prepare(
        input_mode="tcga_spliceseq",
        psi_tsv="tests/data/toy_spliceseq_public.tsv",
        sample_annotations_tsv="tests/data/toy_spliceseq_sample_annotations.tsv",
        out_dir=out_dir,
        organism="human",
        genome_build="hg38",
        study_id="TCGA_TOY",
        study_label="Toy TCGA SpliceSeq",
    )

    assert summary["workflow"] == "splice_prepare_public"
    assert summary["n_samples"] == 4
    assert summary["n_events"] == 3
    assert (out_dir / "bundle_source_row.tsv").exists()

    event_rows = _read_rows(out_dir / "event_metadata.tsv")
    assert event_rows[0]["canonical_event_key"]
    bundle_rows = _read_rows(out_dir / "bundle_source_row.tsv")
    assert len(bundle_rows) == 12
    prepared = json.loads((out_dir / "prepare_summary.json").read_text(encoding="utf-8"))
    assert prepared["outputs"]["psi_matrix_tsv"].endswith("psi_matrix.tsv")



def test_splice_prepare_public_duplicate_sample_ids_fail(tmp_path: Path):
    try:
        run_public_prepare(
            input_mode="tcga_spliceseq",
            psi_tsv="tests/data/toy_spliceseq_public.tsv",
            sample_annotations_tsv="tests/data/toy_spliceseq_duplicate_samples.tsv",
            out_dir=tmp_path / "dup_prepare",
            organism="human",
            genome_build="hg38",
            study_id="TCGA_TOY",
            study_label="Toy TCGA SpliceSeq",
        )
    except ValueError as exc:
        assert "Duplicate sample_id" in str(exc)
    else:
        raise AssertionError("Expected duplicate sample ids to raise")



def test_splice_prepare_public_to_matrix_and_bundle_handoff(tmp_path: Path):
    prepared_dir = tmp_path / "prepared_public"
    run_public_prepare(
        input_mode="tcga_spliceseq",
        psi_tsv="tests/data/toy_spliceseq_public.tsv",
        sample_annotations_tsv="tests/data/toy_spliceseq_sample_annotations.tsv",
        out_dir=prepared_dir,
        organism="human",
        genome_build="hg38",
        study_id="TCGA_TOY",
        study_label="Toy TCGA SpliceSeq",
    )

    sources_manifest = tmp_path / "sources.tsv"
    sources_manifest.write_text(
        "path\tsource_dataset\n" + f"{prepared_dir / 'bundle_source_row.tsv'}\tTCGA_TOY\n",
        encoding="utf-8",
    )
    bundle_args = BundleArgs()
    bundle_args.sources_tsv = str(sources_manifest)
    bundle_args.out_dir = str(tmp_path / "splice_bundle")
    bundle_result = run_splice_prepare_reference_bundle(bundle_args)
    assert bundle_result["n_canonical_events"] == 3

    matrix_args = MatrixArgs()
    matrix_args.psi_matrix_tsv = str(prepared_dir / "psi_matrix.tsv")
    matrix_args.sample_metadata_tsv = str(prepared_dir / "sample_metadata.tsv")
    matrix_args.event_metadata_tsv = str(prepared_dir / "event_metadata.tsv")
    matrix_args.out_dir = str(tmp_path / "splice_from_public")
    matrix_args.resources_dir = bundle_args.out_dir
    matrix_args.use_reference_bundle = True
    result = splice_event_matrix.run(matrix_args)
    assert result["n_contrasts_emitted"] == 1
