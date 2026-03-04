import csv
from pathlib import Path

import pytest

from geneset_extractors.converters import drug_response_screen
from geneset_extractors.core.validate import validate_output_dir


class Args:
    out_dir = "tests/tmp/drug_response"
    organism = "human"
    genome_build = "hg38"
    dataset_label = "toy_drug_response"

    response_tsv = "tests/data/toy_drug_response.tsv"
    response_delimiter = "\t"
    sample_id_column = "sample_id"
    drug_id_column = "drug_id"
    response_column = "response"

    drug_targets_tsv = "tests/data/toy_drug_targets.tsv"
    targets_delimiter = "\t"
    targets_drug_id_column = "drug_id"
    targets_gene_symbol_column = "gene_symbol"
    targets_weight_column = "weight"
    targets_source_column = "source"

    sample_metadata_tsv = None
    sample_metadata_delimiter = "\t"
    sample_metadata_sample_id_column = "sample_id"
    group_column = "group"

    groups_tsv = "tests/data/toy_drug_groups.tsv"
    groups_delimiter = "\t"
    groups_sample_id_column = "sample_id"
    groups_group_column = "group"

    case_control_tsv = None
    case_control_delimiter = "\t"
    case_control_sample_id_column = "sample_id"
    case_control_is_case_column = "is_case"
    case_control_within_group = False

    prism_matrix_csv = None
    prism_treatment_info_csv = None
    prism_cell_line_info_csv = None
    prism_matrix_sample_id_column = "row_name"
    prism_cell_line_sample_id_column = "row_name"
    prism_column_name_column = "column_name"
    prism_broad_id_column = "broad_id"
    prism_target_column = "target"

    response_metric = "logfold_change"
    response_direction = None
    response_transform = "robust_z_mad"
    contrast_method = "group_vs_rest"
    scoring_model = "target_weighted_sum"
    sparse_alpha = 0.01
    ubiquity_penalty = "fraction_active"
    ubiquity_tau = 1.0
    ubiquity_epsilon = 0.05
    polypharm_downweight = True
    polypharm_t0 = 5
    max_programs = 50

    select = "top_k"
    top_k = 2
    quantile = 0.01
    min_score = 0.0
    normalize = "within_set_l1"
    emit_full = True

    target_aliases_tsv = None
    alias_delimiter = "\t"
    alias_column = "alias"
    alias_gene_symbol_column = "gene_symbol"
    drug_blacklist_tsv = None
    drug_blacklist_delimiter = "\t"
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    target_aliases_resource_id = None
    drug_blacklist_resource_id = None

    gtf = None
    gtf_source = None
    gtf_gene_id_field = "gene_id"

    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = False
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 1
    gmt_max_genes = 10
    gmt_topk_list = "2"
    gmt_mass_list = ""
    gmt_split_signed = None
    emit_small_gene_sets = True


def _read_manifest_rows(out_dir: Path) -> list[dict[str, str]]:
    with (out_dir / "manifest.tsv").open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read_geneset_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_drug_response_generic_group_vs_rest_signed_gmt(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "drug_response_generic")
    result = drug_response_screen.run(args)
    assert result["n_groups"] == 2

    out_dir = Path(args.out_dir)
    rows = _read_manifest_rows(out_dir)
    assert len(rows) == 2
    for row in rows:
        p = out_dir / str(row["path"])
        assert (p / "geneset.tsv").exists()
        gmt_path = p / "genesets.gmt"
        assert gmt_path.exists()
        text = gmt_path.read_text(encoding="utf-8")
        assert "__pos__" in text
        assert "__neg__" in text

    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(out_dir, schema)


def test_drug_response_prism_convenience_mode(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "drug_response_prism")
    args.response_tsv = None
    args.drug_targets_tsv = None
    args.groups_tsv = None
    args.contrast_method = None
    args.prism_matrix_csv = "tests/data/toy_prism_matrix.csv"
    args.prism_treatment_info_csv = "tests/data/toy_prism_treatment_info.csv"
    args.prism_cell_line_info_csv = "tests/data/toy_prism_cell_line_info.csv"
    result = drug_response_screen.run(args)
    assert result["n_groups"] == 2
    out_dir = Path(args.out_dir)
    assert (out_dir / "manifest.tsv").exists()
    assert (out_dir / "genesets.gmt").exists()


def test_drug_response_ubiquity_penalty_changes_top_gene(tmp_path: Path):
    args_none = Args()
    args_none.out_dir = str(tmp_path / "drug_response_no_penalty")
    args_none.response_tsv = "tests/data/toy_drug_response_ubiquity.tsv"
    args_none.drug_targets_tsv = "tests/data/toy_drug_targets_ubiquity.tsv"
    args_none.groups_tsv = None
    args_none.contrast_method = "none"
    args_none.response_metric = "other"
    args_none.response_direction = "higher_is_more_sensitive"
    args_none.response_transform = "none"
    args_none.ubiquity_penalty = "none"
    args_none.max_programs = 1
    args_none.top_k = 1
    args_none.gmt_topk_list = "1"
    drug_response_screen.run(args_none)

    args_pen = Args()
    args_pen.out_dir = str(tmp_path / "drug_response_with_penalty")
    args_pen.response_tsv = "tests/data/toy_drug_response_ubiquity.tsv"
    args_pen.drug_targets_tsv = "tests/data/toy_drug_targets_ubiquity.tsv"
    args_pen.groups_tsv = None
    args_pen.contrast_method = "none"
    args_pen.response_metric = "other"
    args_pen.response_direction = "higher_is_more_sensitive"
    args_pen.response_transform = "none"
    args_pen.ubiquity_penalty = "fraction_active"
    args_pen.max_programs = 1
    args_pen.top_k = 1
    args_pen.gmt_topk_list = "1"
    drug_response_screen.run(args_pen)

    out_none = Path(args_none.out_dir)
    out_pen = Path(args_pen.out_dir)
    row_none = _read_manifest_rows(out_none)[0]
    row_pen = _read_manifest_rows(out_pen)[0]
    gene_none = _read_geneset_rows(out_none / row_none["path"] / "geneset.tsv")[0]["gene_id"]
    gene_pen = _read_geneset_rows(out_pen / row_pen["path"] / "geneset.tsv")[0]["gene_id"]
    assert gene_none == "GENEUB"
    assert gene_pen == "GENESEL"


def test_drug_response_low_target_coverage_warning(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    args = Args()
    args.out_dir = str(tmp_path / "drug_response_low_coverage")
    args.response_tsv = "tests/data/toy_drug_response.tsv"
    args.drug_targets_tsv = "tests/data/toy_drug_targets_low_coverage.tsv"
    drug_response_screen.run(args)
    captured = capsys.readouterr()
    assert "low drug-to-target coverage" in captured.err
