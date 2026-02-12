from pathlib import Path
import json
import shutil

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
    link_method = "promoter_overlap"
    promoter_upstream_bp = 2000
    promoter_downstream_bp = 500
    max_distance_bp = None
    decay_length_bp = 50000
    max_genes_per_peak = 5
    normalize = "l1"
    gtf_source = "toy"


def test_sc_converter_without_groups(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc")
    args.groups_tsv = None
    atac_sc_10x.run(args)
    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)


def test_sc_converter_with_groups(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_groups")
    args.groups_tsv = "tests/data/barcode_groups.tsv"
    atac_sc_10x.run(args)
    manifest = Path(args.out_dir) / "manifest.tsv"
    assert manifest.exists()
    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)
    validate_output_dir(Path(args.out_dir) / "group=g1", schema)
    validate_output_dir(Path(args.out_dir) / "group=g2", schema)
    meta = json.loads((Path(args.out_dir) / "group=g1" / "geneset.meta.json").read_text(encoding="utf-8"))
    roles = {f["role"] for f in meta["input"]["files"]}
    assert {"matrix", "barcodes", "peaks_or_features", "gtf", "groups_tsv"}.issubset(roles)
    params = meta["converter"]["parameters"]
    assert "matrix_dir" not in params
    assert "out_dir" not in params


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
