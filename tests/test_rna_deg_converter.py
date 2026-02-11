from pathlib import Path

from omics2geneset.converters import rna_deg
from omics2geneset.core.validate import validate_output_dir


class Args:
    deg_tsv = "tests/data/toy_rna_deg.tsv"
    out_dir = "tests/tmp/rna_deg"
    organism = "human"
    genome_build = "hg38"
    gene_id_column = "gene_id"
    score_column = "log2fc"
    score_transform = "signed"
    normalize = "l1"


def test_rna_deg_converter(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "rna_deg")
    rna_deg.run(args)
    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)
