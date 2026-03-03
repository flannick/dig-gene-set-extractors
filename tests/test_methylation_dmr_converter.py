from pathlib import Path

from geneset_extractors.converters import methylation_dmr
from geneset_extractors.core.validate import validate_output_dir


class Args:
    dmr_tsv = "tests/data/toy_methylation_dmr.tsv"
    out_dir = "tests/tmp/meth"
    organism = "human"
    genome_build = "hg38"
    gene_id_column = "gene_id"
    score_column = "delta_methylation"
    score_transform = "abs"
    normalize = "l1"


def test_methylation_dmr_converter(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "meth")
    methylation_dmr.run(args)
    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)
