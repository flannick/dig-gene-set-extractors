from pathlib import Path

from omics2geneset.converters import sc_rna_marker
from omics2geneset.core.validate import validate_output_dir


class Args:
    counts_tsv = "tests/data/toy_sc_rna_counts.tsv"
    out_dir = "tests/tmp/sc_rna"
    organism = "human"
    genome_build = "hg38"
    groups_tsv = None
    gene_id_column = "gene_id"
    barcode_column = "barcode"
    value_column = "count"
    group_barcode_column = "barcode"
    group_column = "group"
    peak_summary = "sum_counts"
    normalize = "l1"


def test_sc_rna_marker_without_groups(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_rna")
    args.groups_tsv = None
    sc_rna_marker.run(args)
    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)


def test_sc_rna_marker_with_groups(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_rna_groups")
    args.groups_tsv = "tests/data/toy_sc_rna_groups.tsv"
    sc_rna_marker.run(args)
    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir) / "group=g1", schema)
    validate_output_dir(Path(args.out_dir) / "group=g2", schema)
    assert (Path(args.out_dir) / "manifest.tsv").exists()
