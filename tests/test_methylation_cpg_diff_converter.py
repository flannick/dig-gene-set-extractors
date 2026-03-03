import csv
import gzip
import json
from pathlib import Path

import pytest

from geneset_extractors.converters import methylation_cpg_diff
from geneset_extractors.core.validate import validate_output_dir


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
    delta_orientation = "activity_oriented"
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
    exclude_gene_symbol_regex = None
    exclude_gene_symbols_tsv = None


def test_methylation_cpg_diff_runs_and_emits_gmt(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "meth_cpg")
    result = methylation_cpg_diff.run(args)
    assert result["n_peaks"] > 0
    assert result["n_genes"] > 0

    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
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


def test_methylation_cpg_diff_auto_manifest_resolution_from_resources_dir(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    resources_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = resources_dir / "methylation_probe_manifest_450k_hg38.tsv.gz"
    with gzip.open(manifest_path, "wt", encoding="utf-8") as fh:
        fh.write("probe_id\tchrom\tpos\n")
        fh.write("cg_auto\tchr1\t1001\n")

    cpg_path = tmp_path / "cpg.tsv"
    cpg_path.write_text("probe_id\tdelta_beta\tpvalue\ncg_auto\t-0.25\t1e-6\n", encoding="utf-8")

    args = Args()
    args.cpg_tsv = str(cpg_path)
    args.probe_manifest_tsv = None
    args.probe_manifest_resource_id = None
    args.resources_dir = str(resources_dir)
    args.out_dir = str(tmp_path / "meth_cpg_auto_resource")
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "1"
    args.emit_small_gene_sets = True
    methylation_cpg_diff.run(args)

    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    used = meta.get("resources", {}).get("used", [])
    assert any(str(entry.get("id")) == "methylation_probe_manifest_450k_hg38" for entry in used)
    parse_summary = meta.get("summary", {}).get("parse_summary", {})
    assert str(parse_summary.get("probe_manifest_source", "")).startswith("resource:")
    assert str(parse_summary.get("probe_manifest_path", "")).endswith(
        "methylation_probe_manifest_450k_hg38.tsv.gz"
    )


def test_methylation_cpg_diff_delta_orientation_controls_score_sign(tmp_path: Path):
    manifest_path = tmp_path / "manifest.tsv"
    manifest_path.write_text("probe_id\tchrom\tpos\ncg_sign\tchr1\t1001\n", encoding="utf-8")
    cpg_path = tmp_path / "cpg.tsv"
    cpg_path.write_text("probe_id\tdelta_beta\tpvalue\ncg_sign\t0.40\t1e-6\n", encoding="utf-8")

    args_default = Args()
    args_default.cpg_tsv = str(cpg_path)
    args_default.probe_manifest_tsv = str(manifest_path)
    args_default.program_methods = "promoter_activity"
    args_default.link_method = "promoter_overlap"
    args_default.top_k = 1
    args_default.gmt_min_genes = 1
    args_default.gmt_max_genes = 10
    args_default.gmt_topk_list = "1"
    args_default.emit_small_gene_sets = True
    args_default.out_dir = str(tmp_path / "meth_cpg_default_orientation")
    methylation_cpg_diff.run(args_default)

    args_raw = Args()
    args_raw.cpg_tsv = str(cpg_path)
    args_raw.probe_manifest_tsv = str(manifest_path)
    args_raw.program_methods = "promoter_activity"
    args_raw.link_method = "promoter_overlap"
    args_raw.top_k = 1
    args_raw.delta_orientation = "raw"
    args_raw.gmt_min_genes = 1
    args_raw.gmt_max_genes = 10
    args_raw.gmt_topk_list = "1"
    args_raw.emit_small_gene_sets = True
    args_raw.out_dir = str(tmp_path / "meth_cpg_raw_orientation")
    methylation_cpg_diff.run(args_raw)

    rows_default = list(
        csv.DictReader((Path(args_default.out_dir) / "geneset.tsv").open("r", encoding="utf-8"), delimiter="\t")
    )
    rows_raw = list(csv.DictReader((Path(args_raw.out_dir) / "geneset.tsv").open("r", encoding="utf-8"), delimiter="\t"))
    score_default = float(rows_default[0]["score"])
    score_raw = float(rows_raw[0]["score"])
    assert score_default < 0.0
    assert score_raw > 0.0
    assert abs(abs(score_default) - abs(score_raw)) < 1e-9


def test_methylation_cpg_diff_warnings_and_resources_missing_recorded(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
):
    gtf_path = tmp_path / "warn.gtf"
    gtf_path.write_text(
        'chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "GENE1"; gene_name "G1";\n'
        'chrX\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "GENEX"; gene_name "GX";\n',
        encoding="utf-8",
    )
    cpg_path = tmp_path / "warn.tsv"
    cpg_path.write_text(
        "chrom\tpos\tdelta_beta\tpvalue\n"
        "chr1\t1001\t-0.2\t1e-6\n"
        "chrX\t1001\t-0.2\t1e-6\n"
        "chrM\t1001\t-0.2\t1e-6\n",
        encoding="utf-8",
    )
    args = Args()
    args.gtf = str(gtf_path)
    args.cpg_tsv = str(cpg_path)
    args.probe_manifest_tsv = None
    args.enhancer_resource_id = "encode_ccre_hg38"
    args.resources_dir = str(tmp_path / "empty_resources")
    Path(args.resources_dir).mkdir(parents=True, exist_ok=True)
    args.out_dir = str(tmp_path / "meth_cpg_warn")
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "1"
    args.emit_small_gene_sets = True
    methylation_cpg_diff.run(args)
    captured = capsys.readouterr()
    assert "drop_sex_chrom=true removed" in captured.err
    assert "many rows could not be mapped to GTF chromosome names/genome build" in captured.err
    assert "Missing enhancer_bed resource file" in captured.err

    run_summary_txt = (Path(args.out_dir) / "run_summary.txt").read_text(encoding="utf-8")
    assert "n_dropped_sex_chrom: 1" in run_summary_txt
    assert "n_rows_skipped_chrom_mismatch: 1" in run_summary_txt

    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    missing = meta.get("resources", {}).get("missing", [])
    assert missing
    assert any(str(entry.get("resource_id")) == "encode_ccre_hg38" for entry in missing)


def test_methylation_cpg_diff_artifact_warning_and_exclude_regex(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
):
    gtf_path = tmp_path / "artifact.gtf"
    gtf_path.write_text(
        'chr1\ttest\tgene\t1001\t1400\t.\t+\t.\tgene_id "OR1A1"; gene_name "OR1A1";\n'
        'chr1\ttest\tgene\t3001\t3400\t.\t+\t.\tgene_id "OR2A1"; gene_name "OR2A1";\n'
        'chr1\ttest\tgene\t5001\t5400\t.\t+\t.\tgene_id "GENEA"; gene_name "GENEA";\n',
        encoding="utf-8",
    )
    cpg_path = tmp_path / "artifact.tsv"
    cpg_path.write_text(
        "chrom\tpos\tdelta_beta\tpvalue\n"
        "chr1\t1001\t-0.6\t1e-6\n"
        "chr1\t3001\t-0.5\t1e-6\n"
        "chr1\t5001\t-0.1\t1e-6\n",
        encoding="utf-8",
    )

    args = Args()
    args.gtf = str(gtf_path)
    args.cpg_tsv = str(cpg_path)
    args.probe_manifest_tsv = None
    args.program_methods = "promoter_activity"
    args.link_method = "promoter_overlap"
    args.out_dir = str(tmp_path / "artifact_warn")
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.emit_small_gene_sets = True
    methylation_cpg_diff.run(args)
    captured = capsys.readouterr()
    assert "potential gene-family dominance" in captured.err
    assert "pattern=^OR[0-9A-Z]+" in captured.err

    args_filtered = Args()
    args_filtered.gtf = str(gtf_path)
    args_filtered.cpg_tsv = str(cpg_path)
    args_filtered.probe_manifest_tsv = None
    args_filtered.program_methods = "promoter_activity"
    args_filtered.link_method = "promoter_overlap"
    args_filtered.out_dir = str(tmp_path / "artifact_filtered")
    args_filtered.exclude_gene_symbol_regex = ["^OR"]
    args_filtered.gmt_min_genes = 1
    args_filtered.gmt_max_genes = 10
    args_filtered.gmt_topk_list = "3"
    args_filtered.emit_small_gene_sets = True
    methylation_cpg_diff.run(args_filtered)

    gmt_text = (Path(args_filtered.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "OR1A1" not in gmt_text
    assert "OR2A1" not in gmt_text
    assert "GENEA" in gmt_text
