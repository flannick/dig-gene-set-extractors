import csv
import shutil
from pathlib import Path

from geneset_extractors.converters import splice_event_matrix
from geneset_extractors.core.validate import validate_output_dir


class Args:
    psi_matrix_tsv = "tests/data/toy_splice_matrix.tsv"
    sample_metadata_tsv = "tests/data/toy_splice_sample_metadata.tsv"
    event_metadata_tsv = None
    coverage_matrix_tsv = "tests/data/toy_splice_coverage_matrix.tsv"
    out_dir = "tests/tmp/splice_matrix"
    organism = "human"
    genome_build = "human"
    signature_name = "toy_splice_matrix"
    dataset_label = "toy_splice_matrix_dataset"
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
    tool_family = "generic"
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
    min_read_support = 5.0
    neglog10p_cap = 50.0
    neglog10p_eps = 1e-300
    delta_psi_soft_floor = 0.05
    delta_psi_soft_floor_mode = "auto"
    event_dup_policy = "highest_confidence"
    gene_aggregation = "signed_topk_mean"
    gene_topk_events = 3
    gene_burden_penalty_mode = "auto"
    min_gene_burden_penalty = 0.35
    ambiguous_gene_policy = "drop"
    impact_mode = "conservative"
    impact_min = 0.75
    impact_max = 1.35
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    use_reference_bundle = True
    event_alias_resource_id = None
    event_ubiquity_resource_id = None
    event_impact_resource_id = None
    gene_burden_resource_id = None
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


def _copy_resources(resources_dir: Path) -> None:
    resources_dir.mkdir(parents=True, exist_ok=True)
    for filename in (
        "splice_event_aliases_human_v1.tsv.gz",
        "splice_event_ubiquity_human_v1.tsv.gz",
        "splice_event_impact_human_v1.tsv.gz",
        "splice_gene_event_burden_human_v1.tsv.gz",
    ):
        shutil.copyfile(Path("tests/data") / filename, resources_dir / filename)


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _score_map(path: Path) -> dict[str, float]:
    return {row["gene_id"]: float(row["score"]) for row in _read_rows(path)}


def test_splice_event_matrix_single_contrast_end_to_end(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args = Args()
    args.out_dir = str(tmp_path / "splice_matrix_single")
    args.resources_dir = str(resources_dir)
    result = splice_event_matrix.run(args)
    assert result["n_contrasts_emitted"] == 1
    assert not result["grouped_output"]

    out_dir = Path(args.out_dir)
    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(out_dir, schema)
    assert (out_dir / "contrast_qc.tsv").exists()

    rows = _read_rows(out_dir / "geneset.tsv")
    assert rows
    assert rows[0]["gene_symbol"] in {"KCNN4", "MAPK1", "GENEA"}


def test_splice_event_matrix_grouped_contrast_mode(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args = Args()
    args.out_dir = str(tmp_path / "splice_matrix_grouped")
    args.resources_dir = str(resources_dir)
    args.study_contrast = "condition_within_group"
    args.min_samples_per_condition = 1
    args.min_present_per_condition = 1
    result = splice_event_matrix.run(args)
    assert result["grouped_output"]
    assert result["n_contrasts_emitted"] == 2

    out_dir = Path(args.out_dir)
    assert (out_dir / "manifest.tsv").exists()
    assert (out_dir / "contrast_qc.tsv").exists()
    assert (out_dir / "genesets.gmt").exists()

    rows = _read_rows(out_dir / "manifest.tsv")
    assert len(rows) == 2
    for row in rows:
        child_dir = out_dir / row["path"]
        assert (child_dir / "geneset.tsv").exists()
        assert (child_dir / "geneset.meta.json").exists()


def test_splice_event_matrix_effect_metric_change(tmp_path: Path):
    args_t = Args()
    args_t.out_dir = str(tmp_path / "splice_matrix_t")
    args_t.use_reference_bundle = False
    args_t.effect_metric = "welch_t"
    splice_event_matrix.run(args_t)

    args_m = Args()
    args_m.out_dir = str(tmp_path / "splice_matrix_m")
    args_m.use_reference_bundle = False
    args_m.effect_metric = "mean_diff"
    splice_event_matrix.run(args_m)

    t_scores = _score_map(Path(args_t.out_dir) / "geneset.full.tsv")
    m_scores = _score_map(Path(args_m.out_dir) / "geneset.full.tsv")
    assert t_scores["G_KCNN4"] != m_scores["G_KCNN4"]


def test_splice_event_matrix_min_present_missingness(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "splice_matrix_missing")
    args.use_reference_bundle = False
    args.min_present_per_condition = 2
    splice_event_matrix.run(args)
    rows = _read_rows(Path(args.out_dir) / "contrast_qc.tsv")
    assert rows[0]["n_events_dropped_missing"] != "0"


def test_splice_event_matrix_bundle_matched_fixture_changes_scores(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args_bundle = Args()
    args_bundle.psi_matrix_tsv = "tests/data/toy_splice_matrix_matched.tsv"
    args_bundle.sample_metadata_tsv = "tests/data/toy_splice_sample_metadata_matched.tsv"
    args_bundle.coverage_matrix_tsv = None
    args_bundle.out_dir = str(tmp_path / "matched_bundle")
    args_bundle.resources_dir = str(resources_dir)
    splice_event_matrix.run(args_bundle)

    args_none = Args()
    args_none.psi_matrix_tsv = "tests/data/toy_splice_matrix_matched.tsv"
    args_none.sample_metadata_tsv = "tests/data/toy_splice_sample_metadata_matched.tsv"
    args_none.coverage_matrix_tsv = None
    args_none.out_dir = str(tmp_path / "matched_none")
    args_none.use_reference_bundle = False
    splice_event_matrix.run(args_none)

    bundle_scores = _score_map(Path(args_bundle.out_dir) / "geneset.full.tsv")
    none_scores = _score_map(Path(args_none.out_dir) / "geneset.full.tsv")
    assert bundle_scores != none_scores
    assert bundle_scores["G_MAPK1"] != none_scores["G_MAPK1"]
