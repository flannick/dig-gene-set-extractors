import json
from pathlib import Path

from geneset_extractors.converters import proteomics_diff
from geneset_extractors.core.validate import validate_output_dir


class Args:
    proteomics_tsv = "tests/data/toy_proteomics.tsv"
    out_dir = "tests/tmp/prot"
    organism = "human"
    genome_build = "hg38"
    gene_id_column = "gene_id"
    score_column = "log2fc"
    score_transform = "abs"
    normalize = "l1"


def test_proteomics_diff_converter(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "prot")
    proteomics_diff.run(args)
    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)
    assert (Path(args.out_dir) / "geneset.provenance.json").exists()
    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["provenance"]["focus_node_id"]
