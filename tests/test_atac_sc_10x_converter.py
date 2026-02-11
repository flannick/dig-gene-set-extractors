from pathlib import Path

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
    max_distance_bp = 100000
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
    validate_output_dir(Path(args.out_dir) / "group=g1", schema)
    validate_output_dir(Path(args.out_dir) / "group=g2", schema)
