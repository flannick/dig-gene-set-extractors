import csv
import gzip
import json
from pathlib import Path

from omics2geneset.converters import atac_bulk
from omics2geneset.core.validate import validate_output_dir


class Args:
    peaks = "tests/data/toy_peaks.bed"
    gtf = "tests/data/toy.gtf"
    out_dir = "tests/tmp/bulk"
    organism = "human"
    genome_build = "hg38"
    peaks_weight_column = 5
    peak_weights_tsv = "tests/data/toy_peak_weights.tsv"
    link_method = "promoter_overlap"
    promoter_upstream_bp = 2000
    promoter_downstream_bp = 500
    max_distance_bp = None
    decay_length_bp = 50000
    max_genes_per_peak = 5
    peak_weight_transform = "abs"
    normalize = "l1"
    gtf_source = "toy"


def test_bulk_converter_end_to_end(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk")
    atac_bulk.run(args)

    geneset = Path(args.out_dir) / "geneset.tsv"
    meta = Path(args.out_dir) / "geneset.meta.json"
    assert geneset.exists()
    assert meta.exists()

    total = 0.0
    with geneset.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    assert "gene_id" in reader.fieldnames
    assert "weight" in reader.fieldnames
    for r in rows:
        w = float(r["weight"])
        assert w >= 0
        total += w
    assert abs(total - 1.0) < 1e-9

    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)

    payload = json.loads(meta.read_text(encoding="utf-8"))
    assert payload["summary"]["fraction_features_assigned"] == 2.0 / 3.0
    params = payload["converter"]["parameters"]
    assert "peaks" not in params
    assert "out_dir" not in params


def test_bulk_nearest_tss_uses_method_default(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_nearest")
    args.link_method = "nearest_tss"
    args.max_distance_bp = None
    atac_bulk.run(args)
    payload = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    assert payload["converter"]["parameters"]["max_distance_bp"] == 100000


def test_bulk_distance_decay_uses_method_default(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_decay")
    args.link_method = "distance_decay"
    args.max_distance_bp = None
    atac_bulk.run(args)
    payload = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    assert payload["converter"]["parameters"]["max_distance_bp"] == 500000


def test_bulk_converter_accepts_gzipped_peaks(tmp_path: Path):
    gz_peaks = tmp_path / "toy_peaks.bed.gz"
    with Path("tests/data/toy_peaks.bed").open("rb") as src, gzip.open(gz_peaks, "wb") as dst:
        dst.write(src.read())

    args = Args()
    args.peaks = str(gz_peaks)
    args.out_dir = str(tmp_path / "bulk_gz")
    args.peak_weights_tsv = "tests/data/toy_peak_weights.tsv"
    atac_bulk.run(args)

    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)
