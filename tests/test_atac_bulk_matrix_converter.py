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
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    ref_ubiquity_resource_id = None
    atlas_resource_id = None
    atlas_metric = "logratio"
    atlas_eps = 1e-6
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
    gtf_source = "toy"


def test_bulk_matrix_converter_emits_open_close(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_matrix")
    atac_bulk_matrix.run(args)

    geneset = Path(args.out_dir) / "geneset.tsv"
    assert geneset.exists()
    with geneset.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    assert abs(sum(float(r["weight"]) for r in rows) - 1.0) < 1e-9

    gmt_path = Path(args.out_dir) / "genesets.gmt"
    text = gmt_path.read_text(encoding="utf-8")
    assert "__direction=OPEN" in text
    assert "__direction=CLOSE" in text
    assert "__contrast=condition_between_samples" in text

    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    contrast = meta["program_extraction"]["contrast"]
    assert contrast["mode"] == "condition_between_samples"
    assert contrast["condition_a"] == "case"
    assert contrast["condition_b"] == "control"
    assert contrast["n_samples_a"] == 2
    assert contrast["n_samples_b"] == 2
    assert contrast["directions_emitted"] == ["OPEN", "CLOSE"]

    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)
