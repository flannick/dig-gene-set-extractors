import csv
import gzip
import json
from pathlib import Path

from geneset_extractors.converters import splice_event_diff
from geneset_extractors.workflows.splice_prepare_reference_bundle import run as run_splice_prepare_reference_bundle


class DiffArgs:
    splice_tsv = "tests/data/toy_splice_event_diff.tsv"
    cluster_stats_tsv = None
    out_dir = ""
    organism = "human"
    genome_build = "human"
    signature_name = "bundle_runtime"
    dataset_label = "bundle_runtime"
    tool_family = "generic"
    event_id_column = None
    event_group_column = None
    event_type_column = None
    gene_id_column = None
    gene_symbol_column = None
    chrom_column = None
    start_column = None
    end_column = None
    strand_column = None
    score_column = None
    stat_column = None
    delta_psi_column = None
    psi_column = None
    padj_column = None
    pvalue_column = None
    probability_column = None
    read_support_column = None
    novel_flag_column = None
    annotation_status_column = None
    score_mode = "auto"
    score_transform = "signed"
    confidence_weight_mode = "combined"
    min_probability = 0.8
    min_read_support = 10.0
    neglog10p_cap = 50.0
    neglog10p_eps = 1e-300
    event_dup_policy = "highest_confidence"
    gene_aggregation = "signed_topk_mean"
    gene_topk_events = 3
    ambiguous_gene_policy = "drop"
    impact_mode = "conservative"
    impact_min = 0.75
    impact_max = 1.35
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    use_reference_bundle = True
    event_alias_resource_id = None
    event_ubiquity_resource_id = None
    event_impact_resource_id = None
    select = "top_k"
    top_k = 200
    quantile = 0.01
    min_score = 0.0
    normalize = "within_set_l1"
    emit_full = True
    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = True
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 1
    gmt_max_genes = 10
    gmt_topk_list = "2"
    gmt_mass_list = ""
    gmt_split_signed = True
    gmt_format = "dig2col"
    emit_small_gene_sets = True
    delta_psi_soft_floor = 0.05
    delta_psi_soft_floor_mode = "auto"
    gene_burden_penalty_mode = "auto"
    min_gene_burden_penalty = 0.35
    gene_burden_resource_id = None


class BundleArgs:
    sources_tsv = ""
    out_dir = ""
    organism = "human"
    bundle_id = "toy_splice_bundle_v1"
    min_ref_read_support = 0.0
    exclude_source_datasets = None



def test_splice_prepare_reference_bundle_and_runtime_auto_resolution(tmp_path: Path):
    source_rows = tmp_path / "bundle_source_row.tsv"
    source_rows.write_text(
        "\n".join(
            [
                "source_dataset\tsample_id\tinput_event_key\tcanonical_event_key\tcanonicalization_status\tcanonicalization_confidence\tgene_id\tgene_symbol\tevent_type\tpsi\tread_support\tannotation_status\tchrom\tstart\tend\tstrand\ttool_family",
                "toy\tS1\tE1\tchr19:100-200:+:exon_skip:CL1\tcoordinate_canonical\thigh\tG_KCNN4\tKCNN4\texon_skip\t0.8\t30\tannotated_coding\tchr19\t100\t200\t+\tgeneric",
                "toy\tS2\tE2\tE2\traw_id_fallback\tlow\tG_MAPK1\tMAPK1\tretained_intron\t0.6\t15\tnovel\tchr22\t300\t420\t-\tgeneric",
                "toy\tS3\tE4\tchr19:240-310:+:alt_donor:CL4\tcoordinate_canonical\thigh\tG_KCNN4\tKCNN4\talt_donor\t0.7\t18\tannotated_coding\tchr19\t240\t310\t+\tgeneric",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sources_manifest = tmp_path / "sources.tsv"
    sources_manifest.write_text(f"path\tsource_dataset\n{source_rows}\ttoy\n", encoding="utf-8")

    args = BundleArgs()
    args.sources_tsv = str(sources_manifest)
    args.out_dir = str(tmp_path / "bundle")
    result = run_splice_prepare_reference_bundle(args)
    assert result["n_canonical_events"] == 3

    bundle_dir = Path(args.out_dir)
    assert (bundle_dir / "splice_event_aliases_human_v1.tsv.gz").exists()
    assert (bundle_dir / "splice_event_ubiquity_human_v1.tsv.gz").exists()
    assert (bundle_dir / "splice_event_impact_human_v1.tsv.gz").exists()
    assert (bundle_dir / "splice_gene_event_burden_human_v1.tsv.gz").exists()
    assert (bundle_dir / "bundle_provenance.json").exists()
    assert (bundle_dir / "local_resources_manifest.json").exists()

    with gzip.open(bundle_dir / "splice_event_ubiquity_human_v1.tsv.gz", "rt", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 3
    assert rows[0]["canonical_event_key"]
    assert {row["canonicalization_confidence"] for row in rows} == {"high", "low"}
    assert "fraction_datasets_ref" in rows[0]

    with gzip.open(bundle_dir / "splice_gene_event_burden_human_v1.tsv.gz", "rt", encoding="utf-8") as fh:
        burden_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert "n_unique_event_groups_ref" in burden_rows[0]
    assert "n_studies_ref" in burden_rows[0]
    assert "median_unique_groups_per_study" in burden_rows[0]

    diff_args = DiffArgs()
    diff_args.out_dir = str(tmp_path / "runtime")
    diff_args.resources_dir = str(bundle_dir)
    diff_args.resources_manifest = None
    result = splice_event_diff.run(diff_args)
    assert result["n_genes"] >= 1
    run_summary = json.loads((Path(diff_args.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    resources = run_summary["resources"]
    assert resources is not None
    assert resources["missing"] == []
    used_ids = {str(item["id"]) for item in resources["used"]}
    assert used_ids == {
        "splice_event_aliases_human_v1",
        "splice_event_ubiquity_human_v1",
        "splice_event_impact_human_v1",
        "splice_gene_event_burden_human_v1",
    }


def test_splice_prepare_reference_bundle_exclude_source_datasets(tmp_path: Path):
    source_a = tmp_path / "source_a.tsv"
    source_b = tmp_path / "source_b.tsv"
    header = "source_dataset\tsample_id\tinput_event_key\tcanonical_event_key\tcanonicalization_status\tcanonicalization_confidence\tgene_id\tgene_symbol\tevent_group\tevent_type\tpsi\tread_support\tannotation_status\tchrom\tstart\tend\tstrand\ttool_family\n"
    source_a.write_text(
        header
        + "A\tA1\tE1\tchr1:100-150:+:exon_skip:G1\tcoordinate_canonical\thigh\tG1\tGENE1\tG1\texon_skip\t0.8\t20\tannotated_coding\tchr1\t100\t150\t+\tgeneric\n",
        encoding="utf-8",
    )
    source_b.write_text(
        header
        + "B\tB1\tE2\tchr2:200-250:+:exon_skip:G2\tcoordinate_canonical\thigh\tG2\tGENE2\tG2\texon_skip\t0.7\t20\tannotated_coding\tchr2\t200\t250\t+\tgeneric\n",
        encoding="utf-8",
    )
    sources_manifest = tmp_path / "sources.tsv"
    sources_manifest.write_text(
        f"path\tsource_dataset\n{source_a}\tA\n{source_b}\tB\n",
        encoding="utf-8",
    )
    args = BundleArgs()
    args.sources_tsv = str(sources_manifest)
    args.out_dir = str(tmp_path / "bundle_excluded")
    args.exclude_source_datasets = "B"
    result = run_splice_prepare_reference_bundle(args)
    assert result["n_source_datasets"] == 1
    provenance = json.loads((Path(args.out_dir) / "bundle_provenance.json").read_text(encoding="utf-8"))
    assert provenance["excluded_source_datasets"] == ["B"]
