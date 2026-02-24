import csv
import json
from pathlib import Path

from omics2geneset.converters import atac_bulk_matrix
from omics2geneset.core.validate import validate_output_dir


class Args:
    peak_matrix_tsv = "tests/data/toy_bulk_peak_matrix.tsv"
    peak_bed = "tests/data/toy_peaks.bed"
    sample_metadata_tsv = "tests/data/toy_bulk_sample_metadata.tsv"
    gtf = "tests/data/toy.gtf"
    out_dir = "tests/tmp/bulk_matrix"
    organism = "human"
    genome_build = "hg38"
    dataset_label = "toy_bulk"
    sample_id_column = "sample_id"
    condition_column = "condition"
    condition_a = "case"
    condition_b = "control"
    contrast_metric = "log2fc"
    contrast_pseudocount = 1.0
    link_method = "all"
    region_gene_links_tsv = None
    promoter_upstream_bp = 2000
    promoter_downstream_bp = 500
    max_distance_bp = None
    decay_length_bp = 50000
    max_genes_per_peak = 5
    peak_weight_transform = "positive"
    normalize = "within_set_l1"
    program_preset = "connectable"
    program_methods = None
    contrast_methods = "none"
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
    top_k = 1
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
    qc_marker_genes_tsv = None
    gtf_source = "toy"


def test_bulk_matrix_default_connectable_emits_four_directional_sets(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_matrix")
    atac_bulk_matrix.run(args)

    geneset = Path(args.out_dir) / "geneset.tsv"
    assert geneset.exists()
    assert (Path(args.out_dir) / "run_summary.json").exists()
    assert (Path(args.out_dir) / "run_summary.txt").exists()
    with geneset.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    assert abs(sum(float(r["weight"]) for r in rows) - 1.0) < 1e-9

    gmt_text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    lines = [x for x in gmt_text.splitlines() if x.strip()]
    assert 2 <= len(lines) <= 4
    assert "__program=linked_activity__contrast_method=none__link_method=nearest_tss__topk=3" in gmt_text
    assert "__direction=OPEN" in gmt_text
    assert "__direction=CLOSE" in gmt_text
    assert "__program=promoter_activity__" not in gmt_text
    assert "__program=enhancer_bias__" not in gmt_text
    assert "__contrast_method=atlas_residual__" not in gmt_text
    assert "__link_method=promoter_overlap__" not in gmt_text

    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    contrast = meta["program_extraction"]["contrast"]
    assert contrast["mode"] == "condition_between_samples"
    assert contrast["condition_a"] == "case"
    assert contrast["condition_b"] == "control"
    assert contrast["n_samples_a"] == 2
    assert contrast["n_samples_b"] == 2

    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)


def test_bulk_matrix_auto_contrast_warning_when_bundle_disabled(tmp_path: Path, capsys):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_matrix_no_bundle")
    args.use_reference_bundle = False
    args.contrast_methods = "auto_prefer_ref_ubiquity_else_none"
    atac_bulk_matrix.run(args)
    captured = capsys.readouterr()
    assert "auto_prefer_ref_ubiquity_else_none" in captured.err


def test_bulk_matrix_all_preset_allows_explicit_cross_product(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_matrix_all")
    args.program_preset = "all"
    args.link_method = "all"
    args.program_methods = "linked_activity,promoter_activity,distal_activity,enhancer_bias"
    args.contrast_methods = "none"
    atac_bulk_matrix.run(args)

    text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__program=promoter_activity__contrast_method=none__link_method=promoter_overlap__topk=3" in text
    assert "__program=linked_activity__contrast_method=none__link_method=nearest_tss__topk=3" in text
