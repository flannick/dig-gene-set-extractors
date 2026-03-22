from pathlib import Path
import json
import shutil
import csv

import pytest

from geneset_extractors.converters import atac_sc_10x
from geneset_extractors.core.validate import validate_output_dir
from tests.provenance_helpers import (
    assert_manifest_has_enriched_columns,
    assert_node_has_structured_resource_metadata,
    file_node_for_role,
    load_provenance,
)


class Args:
    matrix_dir = "tests/data/toy_10x_mtx"
    gtf = "tests/data/toy.gtf"
    out_dir = "tests/tmp/sc"
    organism = "human"
    genome_build = "hg38"
    dataset_label = None
    groups_tsv = None
    cell_metadata_tsv = None
    condition_column = None
    donor_column = "donor"
    min_donors_per_condition = 1
    min_total_donors_per_condition = 1
    cell_imbalance_warn_ratio = 5.0
    baseline_carry_through_corr_warn = 0.9
    condition_a = None
    condition_b = None
    min_cells_per_group = 1
    min_cells_per_condition = 50
    peak_summary = "sum_counts"
    peak_weight_transform = "positive"
    link_method = "all"
    promoter_upstream_bp = 2000
    promoter_downstream_bp = 500
    max_distance_bp = None
    decay_length_bp = 50000
    max_genes_per_peak = 5
    normalize = "within_set_l1"
    program_preset = "connectable"
    program_methods = None
    calibration_methods = "none"
    resources_manifest = None
    resources_dir = None
    use_reference_bundle = True
    resource_policy = "skip"
    ref_ubiquity_resource_id = None
    atlas_resource_id = None
    atlas_metric = "zscore"
    atlas_eps = 1e-6
    atlas_min_raw_quantile = 0.95
    atlas_use_log1p = True
    select = "top_k"
    top_k = 200
    quantile = 0.01
    min_score = 0.0
    emit_full = True
    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = True
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 1
    gmt_max_genes = 10
    gmt_topk_list = "3"
    gmt_mass_list = ""
    gmt_split_signed = False
    emit_small_gene_sets = False
    qc_marker_genes_tsv = None
    region_gene_links_tsv = None
    study_contrast = None
    contrast_metric = "log2fc"
    contrast_pseudocount = None
    gtf_source = "toy"
    provenance_overlay_json = None


def test_sc_converter_without_groups_default_connectable(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc")
    args.groups_tsv = None
    args.resources_manifest = "tests/data/toy_resources_manifest.json"
    args.resources_dir = "tests/data"
    args.calibration_methods = "auto_prefer_ref_ubiquity_else_none"
    atac_sc_10x.run(args)
    with (Path(args.out_dir) / "geneset.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    assert abs(sum(float(r["weight"]) for r in rows) - 1.0) < 1e-9
    gmt_text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__program=linked_activity__calibration_method=ref_ubiquity_penalty__link_method=nearest_tss__topk=3" in gmt_text
    assert "__program=distal_activity__calibration_method=ref_ubiquity_penalty__link_method=distance_decay__topk=3" in gmt_text
    assert "__program=promoter_activity__" not in gmt_text
    assert "__program=enhancer_bias__" not in gmt_text
    assert "__link_method=promoter_overlap__" not in gmt_text
    assert (Path(args.out_dir) / "run_summary.json").exists()
    assert (Path(args.out_dir) / "run_summary.txt").exists()
    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)
    provenance = load_provenance(args.out_dir)
    node = file_node_for_role(provenance, "resource:ccre_ubiquity_hg38")
    assert_node_has_structured_resource_metadata(node)


def test_sc_converter_with_groups_validates_root_and_groups(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_groups")
    args.groups_tsv = "tests/data/barcode_groups.tsv"
    args.top_k = 1
    atac_sc_10x.run(args)
    manifest = Path(args.out_dir) / "manifest.tsv"
    assert manifest.exists()
    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)
    validate_output_dir(Path(args.out_dir) / "group=g1", schema)
    validate_output_dir(Path(args.out_dir) / "group=g2", schema)
    assert (Path(args.out_dir) / "group=g1" / "run_summary.json").exists()
    assert (Path(args.out_dir) / "group=g2" / "run_summary.json").exists()
    combined_text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "group=g1" in combined_text
    assert "group=g2" in combined_text
    with manifest.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert_manifest_has_enriched_columns(rows)


def test_sc_condition_within_group_emits_open_close_and_connectable_sets(tmp_path: Path):
    args = Args()
    args.matrix_dir = "tests/data/toy_10x_mtx_conditions"
    args.groups_tsv = "tests/data/barcode_groups_conditions.tsv"
    args.cell_metadata_tsv = "tests/data/cell_metadata_conditions.tsv"
    args.condition_column = "condition"
    args.condition_a = "treated"
    args.condition_b = "control"
    args.min_cells_per_condition = 1
    args.study_contrast = "condition_within_group"
    args.out_dir = str(tmp_path / "sc_condition_groups")
    args.link_method = "all"
    atac_sc_10x.run(args)

    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)
    combined_gmt = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__study_contrast=condition_within_group" in combined_gmt
    assert "__direction=OPEN" in combined_gmt
    assert "__direction=CLOSE" in combined_gmt
    assert "__program=linked_activity__calibration_method=none__link_method=nearest_tss__topk=3" in combined_gmt

    meta = json.loads((Path(args.out_dir) / "group=g1" / "geneset.meta.json").read_text(encoding="utf-8"))
    contrast = meta["program_extraction"]["study_contrast"]
    assert contrast["mode"] == "condition_within_group"
    assert contrast["condition_column"] == "condition"
    assert contrast["condition_a"] == "treated"
    assert contrast["condition_b"] == "control"
    assert contrast["donor_column"] == "donor"
    assert contrast["n_donors_a"] == 1
    assert contrast["n_donors_b"] == 1
    assert contrast["directions_emitted"] == ["OPEN", "CLOSE"]


def test_sc_condition_within_group_defaults_to_none_calibration_method(tmp_path: Path):
    args = Args()
    args.matrix_dir = "tests/data/toy_10x_mtx_conditions"
    args.groups_tsv = "tests/data/barcode_groups_conditions.tsv"
    args.cell_metadata_tsv = "tests/data/cell_metadata_conditions.tsv"
    args.condition_column = "condition"
    args.condition_a = "treated"
    args.condition_b = "control"
    args.study_contrast = "condition_within_group"
    args.calibration_methods = "auto_prefer_ref_ubiquity_else_none"
    args.min_cells_per_condition = 1
    args.min_cells_per_group = 1
    args.out_dir = str(tmp_path / "sc_condition_none_default")
    atac_sc_10x.run(args)
    text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__calibration_method=none__" in text
    assert "__calibration_method=ref_ubiquity_penalty__" not in text


def test_sc_auto_contrast_warning_when_bundle_disabled(tmp_path: Path, capsys):
    args = Args()
    args.out_dir = str(tmp_path / "sc_no_reference_bundle")
    args.groups_tsv = "tests/data/barcode_groups.tsv"
    args.use_reference_bundle = False
    args.calibration_methods = "auto_prefer_ref_ubiquity_else_none"
    atac_sc_10x.run(args)
    captured = capsys.readouterr()
    assert "auto_prefer_ref_ubiquity_else_none" in captured.err


def test_sc_all_preset_allows_explicit_cross_product(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_all")
    args.groups_tsv = "tests/data/barcode_groups.tsv"
    args.link_method = "all"
    args.program_preset = "all"
    args.program_methods = "linked_activity,promoter_activity,distal_activity,enhancer_bias,tfidf_distal"
    args.calibration_methods = "none"
    atac_sc_10x.run(args)

    text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__program=promoter_activity__calibration_method=none__link_method=promoter_overlap__topk=3" in text
    assert "__program=tfidf_distal__calibration_method=none__link_method=distance_decay__topk=3" in text


def test_sc_converter_supports_features_tsv_coords(tmp_path: Path):
    matrix_dir = tmp_path / "mtx_features"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy("tests/data/toy_10x_mtx/matrix.mtx", matrix_dir / "matrix.mtx")
    shutil.copy("tests/data/toy_10x_mtx/barcodes.tsv", matrix_dir / "barcodes.tsv")
    (matrix_dir / "features.tsv").write_text(
        "chr1:900-1100\tpeak1\tPeaks\nchr1:1800-1900\tpeak2\tPeaks\nchr1:6000-6100\tpeak3\tPeaks\n",
        encoding="utf-8",
    )

    args = Args()
    args.matrix_dir = str(matrix_dir)
    args.out_dir = str(tmp_path / "sc_features")
    args.groups_tsv = None
    atac_sc_10x.run(args)
    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)


def test_sc_external_requires_region_gene_links(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_external_missing")
    args.link_method = "external"
    args.groups_tsv = None
    args.region_gene_links_tsv = None
    with pytest.raises(ValueError, match="region_gene_links_tsv"):
        atac_sc_10x.run(args)


def test_sc_groups_below_min_cells_per_group_are_skipped(tmp_path: Path, capsys):
    groups_path = tmp_path / "groups.tsv"
    groups_path.write_text(
        "barcode\tgroup\n"
        "c1\tg1\n"
        "c2\tg1\n"
        "c3\tg1\n"
        "c4\tg2\n",
        encoding="utf-8",
    )
    args = Args()
    args.matrix_dir = "tests/data/toy_10x_mtx_conditions"
    args.out_dir = str(tmp_path / "sc_group_size_guard")
    args.groups_tsv = str(groups_path)
    args.cell_metadata_tsv = "tests/data/cell_metadata_conditions.tsv"
    args.condition_column = "condition"
    args.condition_a = "treated"
    args.condition_b = "control"
    args.study_contrast = "condition_within_group"
    args.min_donors_per_condition = 1
    args.min_total_donors_per_condition = 1
    args.min_cells_per_group = 2
    args.min_cells_per_condition = 1
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    atac_sc_10x.run(args)
    captured = capsys.readouterr()
    assert "skipping group due to insufficient cells for stable output" in captured.err
    assert "group=g2" in captured.err
    assert (Path(args.out_dir) / "group=g1" / "geneset.tsv").exists()
    assert not (Path(args.out_dir) / "group=g2" / "geneset.tsv").exists()


def test_sc_condition_group_skipped_when_donor_support_low(tmp_path: Path, capsys):
    args = Args()
    args.matrix_dir = "tests/data/toy_10x_mtx_conditions"
    args.out_dir = str(tmp_path / "sc_donor_gate")
    args.groups_tsv = "tests/data/barcode_groups_conditions.tsv"
    args.cell_metadata_tsv = "tests/data/cell_metadata_conditions.tsv"
    args.condition_column = "condition"
    args.condition_a = "treated"
    args.condition_b = "control"
    args.study_contrast = "condition_within_group"
    args.min_cells_per_group = 1
    args.min_cells_per_condition = 1
    args.min_donors_per_condition = 2
    args.min_total_donors_per_condition = 1
    with pytest.raises(ValueError, match="No groups passed condition_within_group QC gates"):
        atac_sc_10x.run(args)
    captured = capsys.readouterr()
    assert "donor support is too low for condition contrast" in captured.err


def test_sc_condition_dataset_donor_gate_blocks_explicit_contrast(tmp_path: Path):
    args = Args()
    args.matrix_dir = "tests/data/toy_10x_mtx_conditions"
    args.out_dir = str(tmp_path / "sc_donor_total_gate")
    args.groups_tsv = "tests/data/barcode_groups_conditions.tsv"
    args.cell_metadata_tsv = "tests/data/cell_metadata_conditions.tsv"
    args.condition_column = "condition"
    args.condition_a = "treated"
    args.condition_b = "control"
    args.study_contrast = "condition_within_group"
    args.min_cells_per_group = 1
    args.min_cells_per_condition = 1
    args.min_donors_per_condition = 1
    args.min_total_donors_per_condition = 5
    with pytest.raises(ValueError, match="donor gate failed at dataset level"):
        atac_sc_10x.run(args)


def test_sc_group_vs_rest_not_filtered_by_condition_donor_gates(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_group_vs_rest_unchanged")
    args.groups_tsv = "tests/data/barcode_groups.tsv"
    args.study_contrast = "group_vs_rest"
    args.min_cells_per_group = 100
    args.min_cells_per_condition = 100
    args.min_donors_per_condition = 100
    args.min_total_donors_per_condition = 100
    atac_sc_10x.run(args)
    assert (Path(args.out_dir) / "group=g1" / "geneset.tsv").exists()
    assert (Path(args.out_dir) / "group=g2" / "geneset.tsv").exists()
