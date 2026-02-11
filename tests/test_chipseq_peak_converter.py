from pathlib import Path

from omics2geneset.converters import chipseq_peak
from omics2geneset.core.validate import validate_output_dir


class Args:
    peaks = "tests/data/toy_chipseq_peaks.bed"
    gtf = "tests/data/toy.gtf"
    out_dir = "tests/tmp/chipseq"
    organism = "human"
    genome_build = "hg38"
    weight_column = 5
    peak_weight_transform = "positive"
    normalize = "l1"
    gtf_source = "toy"
    link_method = "promoter_overlap"
    promoter_upstream_bp = 2000
    promoter_downstream_bp = 500
    max_distance_bp = None
    decay_length_bp = 50000
    max_genes_per_peak = 5


def test_chipseq_peak_converter(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "chipseq")
    chipseq_peak.run(args)
    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)
