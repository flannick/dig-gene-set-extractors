import csv
import gzip
import json
from pathlib import Path

import pytest

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
    link_method = "all"
    promoter_upstream_bp = 2000
    promoter_downstream_bp = 500
    max_distance_bp = None
    decay_length_bp = 50000
    max_genes_per_peak = 5
    peak_weight_transform = "abs"
    normalize = "within_set_l1"
    program_preset = "none"
    program_methods = None
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    ref_ubiquity_resource_id = None
    atlas_resource_id = None
    atlas_metric = "logratio"
    atlas_eps = 1e-6
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
    region_gene_links_tsv = None
    gtf_source = "toy"


def test_bulk_converter_end_to_end(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk")
    args.top_k = 1
    atac_bulk.run(args)

    geneset = Path(args.out_dir) / "geneset.tsv"
    geneset_full = Path(args.out_dir) / "geneset.full.tsv"
    meta = Path(args.out_dir) / "geneset.meta.json"
    assert geneset.exists()
    assert geneset_full.exists()
    assert meta.exists()

    total = 0.0
    with geneset.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    assert len(rows) == 1
    assert "gene_id" in reader.fieldnames
    assert "score" in reader.fieldnames
    assert "weight" in reader.fieldnames
    assert "rank" in reader.fieldnames
    for r in rows:
        s = float(r["score"])
        assert s >= 0
        w = float(r["weight"])
        assert w >= 0
        total += w
    assert abs(total - 1.0) < 1e-9

    with geneset_full.open("r", encoding="utf-8") as fh:
        full_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(full_rows) >= len(rows)

    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)

    payload = json.loads(meta.read_text(encoding="utf-8"))
    assert payload["summary"]["fraction_features_assigned"] == 2.0 / 3.0
    params = payload["converter"]["parameters"]
    assert "peaks" not in params
    assert "out_dir" not in params
    assert payload["program_extraction"]["n_selected_genes"] == 1
    output_roles = {f["role"] for f in payload["output"]["files"]}
    assert output_roles == {"selected_program", "full_scores", "gmt"}
    assert payload["gmt"]["written"] is True
    assert payload["gmt"]["plans"]


def test_bulk_converter_writes_gmt(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_gmt")
    args.emit_gmt = True
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.gmt_mass_list = ""
    atac_bulk.run(args)

    gmt_path = Path(args.out_dir) / "genesets.gmt"
    assert gmt_path.exists()
    lines = gmt_path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 3
    for line in lines:
        assert line.count("\t") == 1
        name, genes = line.split("\t")
        assert name.endswith("__topk=3")
        assert "__link_method=" in name
        assert set(genes.split(" ")) == {"G1", "G2"}


def test_bulk_default_program_preset_emits_program_family_sets(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_program_families")
    args.program_preset = "default"
    args.link_method = "promoter_overlap"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.gmt_mass_list = ""
    atac_bulk.run(args)

    gmt_text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__program=promoter_activity__topk=3" in gmt_text
    assert "__program=distal_activity__topk=3" in gmt_text
    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    assert "program_methods" in meta["program_extraction"]
    assert "enhancer_bias" in meta["program_extraction"]["program_methods"]


def test_bulk_resource_backed_program_methods(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_resource_methods")
    args.link_method = "promoter_overlap"
    args.program_preset = "none"
    args.program_methods = "ref_ubiquity_penalty,atlas_residual"
    args.resources_manifest = "tests/data/toy_resources_manifest.json"
    args.resources_dir = "tests/data"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.gmt_mass_list = ""
    atac_bulk.run(args)

    gmt_text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__program=ref_ubiquity_penalty__topk=3" in gmt_text
    assert "__program=atlas_residual__topk=3" in gmt_text
    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["resources"]["used"]
    assert not meta["program_extraction"]["program_methods_skipped"]


def test_bulk_resource_policy_skip_skips_missing_method(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_resource_skip")
    args.program_preset = "none"
    args.program_methods = "ref_ubiquity_penalty"
    args.resources_manifest = "tests/data/toy_resources_manifest.json"
    args.resources_dir = str(tmp_path / "missing_resources")
    args.resource_policy = "skip"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.gmt_mass_list = ""
    atac_bulk.run(args)
    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    assert "ref_ubiquity_penalty" in meta["program_extraction"]["program_methods_skipped"]


def test_bulk_resource_policy_fail_raises_on_missing_resource(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_resource_fail")
    args.program_preset = "none"
    args.program_methods = "atlas_residual"
    args.resources_manifest = "tests/data/toy_resources_manifest.json"
    args.resources_dir = str(tmp_path / "missing_resources")
    args.resource_policy = "fail"
    with pytest.raises(FileNotFoundError, match="atlas"):
        atac_bulk.run(args)


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


def test_bulk_external_requires_region_gene_links(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_external_missing")
    args.link_method = "external"
    args.region_gene_links_tsv = None
    with pytest.raises(ValueError, match="region_gene_links_tsv"):
        atac_bulk.run(args)


def test_bulk_all_plus_external_writes_external_gmt(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_external_ok")
    args.link_method = "all,external"
    args.region_gene_links_tsv = "tests/data/toy_region_gene_links.tsv"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.gmt_mass_list = ""
    atac_bulk.run(args)

    gmt_path = Path(args.out_dir) / "genesets.gmt"
    assert gmt_path.exists()
    lines = gmt_path.read_text(encoding="utf-8").splitlines()
    assert any("__link_method=external__" in line for line in lines)


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
