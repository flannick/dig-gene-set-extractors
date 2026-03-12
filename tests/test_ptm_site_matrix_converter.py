import csv
import shutil
from pathlib import Path

from geneset_extractors.converters import ptm_site_matrix
from geneset_extractors.core.validate import validate_output_dir


class Args:
    ptm_matrix_tsv = "tests/data/toy_ptm_matrix.tsv"
    sample_metadata_tsv = "tests/data/toy_ptm_sample_metadata.tsv"
    protein_matrix_tsv = "tests/data/toy_protein_matrix.tsv"
    out_dir = "tests/tmp/ptm_site_matrix"
    organism = "human"
    genome_build = "human"
    signature_name = "toy_ptm_matrix"
    dataset_label = "toy_ptm_matrix_dataset"
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
    protein_adjustment_run_mode = "single"
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
    use_reference_bundle = True
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
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 1
    gmt_max_genes = 20
    gmt_topk_list = "3"
    gmt_mass_list = ""
    gmt_split_signed = True
    gmt_format = "dig2col"
    emit_small_gene_sets = True
    neglog10p_cap = 50.0
    neglog10p_eps = 1e-300


def _copy_resources(resources_dir: Path) -> None:
    resources_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile("tests/data/phosphosite_aliases_human_v1.tsv.gz", resources_dir / "phosphosite_aliases_human_v1.tsv.gz")
    shutil.copyfile("tests/data/phosphosite_ubiquity_human_v1.tsv.gz", resources_dir / "phosphosite_ubiquity_human_v1.tsv.gz")


def test_ptm_site_matrix_single_contrast_end_to_end(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args = Args()
    args.out_dir = str(tmp_path / "ptm_matrix_single")
    args.resources_dir = str(resources_dir)
    args.study_contrast = "condition_a_vs_b"
    result = ptm_site_matrix.run(args)
    assert result["n_contrasts_emitted"] == 1
    assert not result["grouped_output"]

    out_dir = Path(args.out_dir)
    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(out_dir, schema)
    assert (out_dir / "contrast_qc.tsv").exists()

    with (out_dir / "geneset.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    assert rows[0]["gene_symbol"] in {"KCNN4", "AKT1"}
    assert abs(sum(float(r["weight"]) for r in rows) - 1.0) < 1e-9


def test_ptm_site_matrix_condition_within_group_grouped_output(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args = Args()
    args.out_dir = str(tmp_path / "ptm_matrix_grouped")
    args.resources_dir = str(resources_dir)
    args.study_contrast = "condition_within_group"
    result = ptm_site_matrix.run(args)
    assert result["grouped_output"]
    assert result["n_contrasts_emitted"] == 2

    out_dir = Path(args.out_dir)
    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(out_dir, schema)
    assert (out_dir / "manifest.tsv").exists()
    assert (out_dir / "contrast_qc.tsv").exists()
    assert (out_dir / "genesets.gmt").exists()

    with (out_dir / "manifest.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 2
    for row in rows:
        group_dir = out_dir / str(row["path"])
        assert (group_dir / "geneset.tsv").exists()
        assert (group_dir / "geneset.meta.json").exists()


def test_ptm_site_matrix_protein_adjustment_changes_scores(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args_sub = Args()
    args_sub.out_dir = str(tmp_path / "ptm_matrix_subtract")
    args_sub.resources_dir = str(resources_dir)
    args_sub.protein_adjustment = "subtract"
    ptm_site_matrix.run(args_sub)

    args_none = Args()
    args_none.out_dir = str(tmp_path / "ptm_matrix_none")
    args_none.resources_dir = str(resources_dir)
    args_none.protein_adjustment = "none"
    ptm_site_matrix.run(args_none)

    def _score_map(path: Path) -> dict[str, float]:
        with path.open("r", encoding="utf-8") as fh:
            return {row["gene_id"]: float(row["score"]) for row in csv.DictReader(fh, delimiter="\t")}

    sub_scores = _score_map(Path(args_sub.out_dir) / "geneset.full.tsv")
    none_scores = _score_map(Path(args_none.out_dir) / "geneset.full.tsv")
    assert sub_scores["G_KCNN4"] != none_scores["G_KCNN4"]


def test_ptm_site_matrix_default_compare_if_protein_emits_grouped_variants(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args = Args()
    args.out_dir = str(tmp_path / "ptm_matrix_compare")
    args.resources_dir = str(resources_dir)
    args.protein_adjustment_run_mode = "compare_if_protein"
    result = ptm_site_matrix.run(args)
    assert result["grouped_output"]
    assert result["n_contrasts_emitted"] == 2

    out_dir = Path(args.out_dir)
    with (out_dir / "manifest.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {row["protein_adjustment"] for row in rows} == {"none", "subtract"}


def test_ptm_site_matrix_site_cap_comparison_emits_variant(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args = Args()
    args.out_dir = str(tmp_path / "ptm_matrix_sitecap")
    args.resources_dir = str(resources_dir)
    args.protein_matrix_tsv = None
    args.emit_gene_topk_site_comparison = True
    args.gene_topk_site_compare_to = 1
    result = ptm_site_matrix.run(args)
    assert result["grouped_output"]
    assert result["n_contrasts_emitted"] == 2

    out_dir = Path(args.out_dir)
    with (out_dir / "manifest.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {row["gene_topk_sites"] for row in rows} == {"1", "3"}
