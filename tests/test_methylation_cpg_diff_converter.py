from pathlib import Path

import pytest

from omics2geneset.converters import methylation_cpg_diff
from omics2geneset.core.validate import validate_output_dir


class Args:
    cpg_tsv = "tests/data/toy_methylation_cpg.tsv"
    out_dir = "tests/tmp/methylation_cpg"
    organism = "human"
    genome_build = "hg38"
    dataset_label = "toy_cpg"
    gtf = "tests/data/toy.gtf"
    gtf_source = None
    array_type = "450k"
    probe_id_column = "probe_id"
    chrom_column = "chrom"
    pos_column = "pos"
    start_column = None
    end_column = None
    input_pos_is_0based = False
    delta_column = "delta_beta"
    padj_column = None
    pvalue_column = "pvalue"
    probe_manifest_tsv = "tests/data/toy_probe_manifest.tsv"
    probe_manifest_resource_id = None
    manifest_probe_id_column = "probe_id"
    manifest_chrom_column = "chrom"
    manifest_pos_column = "pos"
    manifest_start_column = None
    manifest_end_column = None
    manifest_pos_is_0based = False
    probe_blacklist_tsv = None
    drop_sex_chrom = True
    score_mode = "delta_times_neglog10p"
    neglog10p_cap = 50.0
    neglog10p_eps = 1e-300
    link_method = "all"
    region_gene_links_tsv = None
    promoter_upstream_bp = 2000
    promoter_downstream_bp = 500
    max_distance_bp = None
    decay_length_bp = 50000
    max_genes_per_peak = 5
    program_preset = "connectable"
    program_methods = None
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


def test_methylation_cpg_diff_runs_and_emits_gmt(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "meth_cpg")
    result = methylation_cpg_diff.run(args)
    assert result["n_peaks"] > 0
    assert result["n_genes"] > 0

    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)

    gmt_text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__pos__" in gmt_text
    assert "__neg__" in gmt_text


def test_methylation_cpg_diff_manifest_missing_warns_or_fails(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
):
    args_skip = Args()
    args_skip.cpg_tsv = "tests/data/toy_methylation_cpg_mixed.tsv"
    args_skip.probe_manifest_tsv = None
    args_skip.out_dir = str(tmp_path / "meth_cpg_skip")
    args_skip.resource_policy = "skip"
    methylation_cpg_diff.run(args_skip)
    captured = capsys.readouterr()
    assert "could not be resolved to coordinates" in captured.err

    args_fail = Args()
    args_fail.cpg_tsv = "tests/data/toy_methylation_cpg_mixed.tsv"
    args_fail.probe_manifest_tsv = None
    args_fail.out_dir = str(tmp_path / "meth_cpg_fail")
    args_fail.resource_policy = "fail"
    with pytest.raises(ValueError):
        methylation_cpg_diff.run(args_fail)
