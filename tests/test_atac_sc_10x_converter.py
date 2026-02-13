from pathlib import Path
import json
import shutil
import csv

from omics2geneset.converters import atac_sc_10x
from omics2geneset.core.validate import validate_output_dir


class Args:
    matrix_dir = "tests/data/toy_10x_mtx"
    gtf = "tests/data/toy.gtf"
    out_dir = "tests/tmp/sc"
    organism = "human"
    genome_build = "hg38"
    groups_tsv = None
    peak_summary = "sum_counts"
    peak_weight_transform = "positive"
    link_method = "all"
    promoter_upstream_bp = 2000
    promoter_downstream_bp = 500
    max_distance_bp = None
    decay_length_bp = 50000
    max_genes_per_peak = 5
    normalize = "within_set_l1"
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
