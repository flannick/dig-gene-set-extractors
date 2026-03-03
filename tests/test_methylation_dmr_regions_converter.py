import csv
from pathlib import Path

import pytest

from geneset_extractors.converters import methylation_dmr_regions
from geneset_extractors.core.validate import validate_output_dir


class Args:
    dmr_tsv = "tests/data/toy_methylation_dmr_regions.tsv"
    out_dir = "tests/tmp/methylation_dmr_regions"
    organism = "human"
    genome_build = "hg38"
    dataset_label = "toy_dmr"
    gtf = "tests/data/toy.gtf"
    gtf_source = None
    chrom_column = "chrom"
    start_column = "start"
    end_column = "end"
    delta_column = "delta_methylation"
    padj_column = None
    pvalue_column = "pvalue"
    drop_sex_chrom = True
    score_mode = "delta_times_neglog10p"
    delta_orientation = "activity_oriented"
    neglog10p_cap = 50.0
    neglog10p_eps = 1e-300
    link_method = "promoter_overlap"
    region_gene_links_tsv = None
    promoter_upstream_bp = 2000
    promoter_downstream_bp = 500
    max_distance_bp = None
    decay_length_bp = 50000
    max_genes_per_peak = 5
    program_preset = "connectable"
    program_methods = "promoter_activity"
    select = "top_k"
    top_k = 200
    quantile = 0.01
    min_score = 0.0
    normalize = "within_set_l1"
    aggregation = "weighted_mean"
    score_transform = "signed"
    emit_full = True
    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = True
    gmt_biotype_allowlist = ""
    gmt_min_genes = 1
    gmt_max_genes = 10
    gmt_topk_list = "2"
    gmt_mass_list = ""
    gmt_split_signed = True
    emit_small_gene_sets = True
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    enhancer_resource_id = None
    enhancer_bed = None
    distal_mode = "nonpromoter"
    exclude_gene_symbol_regex = None
    exclude_gene_symbols_tsv = None


def test_methylation_dmr_regions_basic(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "meth_dmr_regions")
    result = methylation_dmr_regions.run(args)
    assert result["n_peaks"] > 0
    assert result["n_genes"] > 0

    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)

    rows = list(csv.DictReader((Path(args.out_dir) / "geneset.tsv").open("r", encoding="utf-8"), delimiter="\t"))
    score_by_gene = {str(row["gene_id"]): float(row["score"]) for row in rows}
    assert "GENE1" in score_by_gene
    assert score_by_gene["GENE1"] < 0.0


def test_methylation_dmr_regions_delta_orientation_raw(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "meth_dmr_regions_raw")
    args.delta_orientation = "raw"
    methylation_dmr_regions.run(args)
    rows = list(csv.DictReader((Path(args.out_dir) / "geneset.tsv").open("r", encoding="utf-8"), delimiter="\t"))
    score_by_gene = {str(row["gene_id"]): float(row["score"]) for row in rows}
    assert "GENE1" in score_by_gene
    assert score_by_gene["GENE1"] > 0.0


def test_methylation_dmr_regions_small_set_warning_has_actionable_hint(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
):
    args = Args()
    args.out_dir = str(tmp_path / "meth_dmr_regions_small_hint")
    args.gmt_min_genes = 100
    args.gmt_max_genes = 500
    args.gmt_topk_list = "200"
    args.emit_small_gene_sets = False
    methylation_dmr_regions.run(args)
    captured = capsys.readouterr()
    assert "--emit_small_gene_sets true" in captured.err
    assert "--gmt_min_genes" in captured.err
