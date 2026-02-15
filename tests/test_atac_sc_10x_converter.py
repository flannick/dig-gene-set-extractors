from pathlib import Path
import json
import shutil
import csv

import pytest

from omics2geneset.converters import atac_sc_10x
from omics2geneset.core.validate import validate_output_dir


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
    condition_a = None
    condition_b = None
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
    program_preset = "none"
    program_methods = None
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    ref_ubiquity_resource_id = None
    atlas_resource_id = None
    atlas_metric = "logratio"
    atlas_eps = 1e-6
    select = "top_k"
    top_k = 200
    quantile = 0.01
    min_score = 0.0
    emit_full = True
    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_min_genes = 100
    gmt_max_genes = 500
    gmt_topk_list = "100,200,500"
    gmt_mass_list = "0.5,0.8,0.9"
    gmt_split_signed = False
    region_gene_links_tsv = None
    contrast = None
    contrast_metric = "log2fc"
    contrast_pseudocount = None
    gtf_source = "toy"


def test_sc_converter_without_groups(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc")
    args.groups_tsv = None
    atac_sc_10x.run(args)
    with (Path(args.out_dir) / "geneset.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    total = sum(float(r["weight"]) for r in rows)
    assert abs(total - 1.0) < 1e-9
    assert (Path(args.out_dir) / "genesets.gmt").exists()
    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)


def test_sc_converter_with_groups(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_groups")
    args.groups_tsv = "tests/data/barcode_groups.tsv"
    args.top_k = 1
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    atac_sc_10x.run(args)
    manifest = Path(args.out_dir) / "manifest.tsv"
    assert manifest.exists()
    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)
    validate_output_dir(Path(args.out_dir) / "group=g1", schema)
    validate_output_dir(Path(args.out_dir) / "group=g2", schema)
    for group in ("g1", "g2"):
        with (Path(args.out_dir) / f"group={group}" / "geneset.tsv").open("r", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert len(rows) == 1
        assert abs(sum(float(r["weight"]) for r in rows) - 1.0) < 1e-9
        assert (Path(args.out_dir) / f"group={group}" / "genesets.gmt").exists()
    combined_gmt = Path(args.out_dir) / "genesets.gmt"
    assert combined_gmt.exists()
    combined_text = combined_gmt.read_text(encoding="utf-8")
    assert "group=g1" in combined_text
    assert "group=g2" in combined_text
    meta = json.loads((Path(args.out_dir) / "group=g1" / "geneset.meta.json").read_text(encoding="utf-8"))
    roles = {f["role"] for f in meta["input"]["files"]}
    assert {"matrix", "barcodes", "peaks_or_features", "gtf", "groups_tsv"}.issubset(roles)
    params = meta["converter"]["parameters"]
    assert "matrix_dir" not in params
    assert "out_dir" not in params
    assert meta["gmt"]["written"] is True
    assert meta["gmt"]["plans"]


def test_sc_converter_group_vs_rest_contrast_marks_expected_genes(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_contrast")
    args.groups_tsv = "tests/data/barcode_groups.tsv"
    args.peak_summary = "sum_counts"
    args.contrast = "group_vs_rest"
    args.contrast_metric = "log2fc"
    args.peak_weight_transform = "positive"
    args.select = "top_k"
    args.top_k = 1
    args.emit_full = False
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    atac_sc_10x.run(args)

    with (Path(args.out_dir) / "group=g1" / "geneset.tsv").open("r", encoding="utf-8") as fh:
        g1_rows = list(csv.DictReader(fh, delimiter="\t"))
    with (Path(args.out_dir) / "group=g2" / "geneset.tsv").open("r", encoding="utf-8") as fh:
        g2_rows = list(csv.DictReader(fh, delimiter="\t"))

    assert g1_rows and g1_rows[0]["gene_id"] == "GENE1"
    assert g2_rows and g2_rows[0]["gene_id"] == "GENE2"

    meta = json.loads((Path(args.out_dir) / "group=g1" / "geneset.meta.json").read_text(encoding="utf-8"))
    contrast = meta["program_extraction"]["contrast"]
    assert contrast["mode"] == "group_vs_rest"
    assert contrast["metric"] == "log2fc"
    assert float(contrast["pseudocount"]) == 1.0


def test_sc_condition_within_group_emits_open_close_gmt_and_metadata(tmp_path: Path):
    args = Args()
    args.matrix_dir = "tests/data/toy_10x_mtx_conditions"
    args.groups_tsv = "tests/data/barcode_groups_conditions.tsv"
    args.cell_metadata_tsv = "tests/data/cell_metadata_conditions.tsv"
    args.condition_column = "condition"
    args.condition_a = "treated"
    args.condition_b = "control"
    args.min_cells_per_condition = 1
    args.contrast = "condition_within_group"
    args.out_dir = str(tmp_path / "sc_condition_groups")
    args.link_method = "promoter_overlap"
    args.top_k = 1
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.gmt_mass_list = ""
    atac_sc_10x.run(args)

    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)

    combined_gmt = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__contrast=condition_within_group" in combined_gmt
    assert "__direction=OPEN" in combined_gmt
    assert "__direction=CLOSE" in combined_gmt

    meta = json.loads((Path(args.out_dir) / "group=g1" / "geneset.meta.json").read_text(encoding="utf-8"))
    contrast = meta["program_extraction"]["contrast"]
    assert contrast["mode"] == "condition_within_group"
    assert contrast["condition_column"] == "condition"
    assert contrast["condition_a"] == "treated"
    assert contrast["condition_b"] == "control"
    assert contrast["n_cells_a"] == 1
    assert contrast["n_cells_b"] == 1
    assert contrast["directions_emitted"] == ["OPEN", "CLOSE"]


def test_sc_connectable_preset_triggers_condition_within_group(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_connectable_condition")
    args.groups_tsv = None
    args.cell_metadata_tsv = "tests/data/barcode_conditions.tsv"
    args.condition_column = "condition"
    args.condition_a = "treated"
    args.condition_b = "control"
    args.min_cells_per_condition = 1
    args.program_preset = "connectable"
    args.contrast = None
    args.link_method = "promoter_overlap"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.gmt_mass_list = ""
    atac_sc_10x.run(args)

    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["program_extraction"]["contrast"]["mode"] == "condition_within_group"
    gmt_text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__direction=OPEN" in gmt_text
    assert "__direction=CLOSE" in gmt_text


def test_sc_connectable_preset_falls_back_to_group_vs_rest_without_metadata(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_connectable_fallback")
    args.groups_tsv = "tests/data/barcode_groups.tsv"
    args.program_preset = "connectable"
    args.contrast = None
    args.link_method = "promoter_overlap"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.gmt_mass_list = ""
    atac_sc_10x.run(args)

    meta = json.loads((Path(args.out_dir) / "group=g1" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["program_extraction"]["contrast"]["mode"] == "group_vs_rest"


def test_sc_default_program_preset_emits_tfidf_distal_sets(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_program_families")
    args.groups_tsv = "tests/data/barcode_groups.tsv"
    args.program_preset = "default"
    args.link_method = "promoter_overlap"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.gmt_mass_list = ""
    atac_sc_10x.run(args)

    combined_gmt = Path(args.out_dir) / "genesets.gmt"
    text = combined_gmt.read_text(encoding="utf-8")
    assert "__program=tfidf_distal__topk=3" in text
    meta = json.loads((Path(args.out_dir) / "group=g1" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert "program_methods" in meta["program_extraction"]
    assert "tfidf_distal" in meta["program_extraction"]["program_methods"]


def test_sc_resource_backed_program_methods(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_resource_methods")
    args.groups_tsv = "tests/data/barcode_groups.tsv"
    args.link_method = "promoter_overlap"
    args.program_preset = "none"
    args.program_methods = "ref_ubiquity_penalty,atlas_residual"
    args.resources_manifest = "tests/data/toy_resources_manifest.json"
    args.resources_dir = "tests/data"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.gmt_mass_list = ""
    atac_sc_10x.run(args)

    combined_gmt = Path(args.out_dir) / "genesets.gmt"
    text = combined_gmt.read_text(encoding="utf-8")
    assert "__program=ref_ubiquity_penalty__topk=3" in text
    assert "__program=atlas_residual__topk=3" in text
    meta = json.loads((Path(args.out_dir) / "group=g1" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["resources"]["used"]
    assert meta["resources"]["used"][0]["stable_id"]
    assert meta["resources"]["used"][0]["version"]


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
    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)


def test_sc_external_requires_region_gene_links(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_external_missing")
    args.link_method = "external"
    args.groups_tsv = None
    args.region_gene_links_tsv = None
    with pytest.raises(ValueError, match="region_gene_links_tsv"):
        atac_sc_10x.run(args)


def test_sc_all_plus_external_writes_external_gmt(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_external_ok")
    args.groups_tsv = "tests/data/barcode_groups.tsv"
    args.link_method = "all,external"
    args.region_gene_links_tsv = "tests/data/toy_region_gene_links.tsv"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.gmt_mass_list = ""
    atac_sc_10x.run(args)

    combined_gmt = Path(args.out_dir) / "genesets.gmt"
    assert combined_gmt.exists()
    text = combined_gmt.read_text(encoding="utf-8")
    assert "__link_method=external__" in text
