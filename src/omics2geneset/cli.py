from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from omics2geneset.core.validate import validate_output_dir
from omics2geneset.resource_manager import (
    describe_resource,
    fetch_resources,
    load_manifest,
    resource_status_rows,
    resolve_requested_resource_ids,
)
from omics2geneset.registry import get_converter, get_converter_spec, list_converters


def _parse_bool(value: str) -> bool:
    v = value.strip().lower()
    if v in {"true", "1", "yes", "y"}:
        return True
    if v in {"false", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError("expected true or false")


def _add_linking_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--link_method",
        default="all",
        help=(
            "Linking mode expression. Tokens: all,promoter_overlap,nearest_tss,distance_decay,external. "
            "Use comma-separated tokens to combine methods, e.g. all,external."
        ),
    )
    parser.add_argument("--region_gene_links_tsv", help="Standardized external linkage TSV with columns chrom,start,end,gene_id,link_weight")
    parser.add_argument("--promoter_upstream_bp", type=int, default=2000)
    parser.add_argument("--promoter_downstream_bp", type=int, default=500)
    parser.add_argument("--max_distance_bp", type=int)
    parser.add_argument("--decay_length_bp", type=int, default=50000)
    parser.add_argument("--max_genes_per_peak", type=int, default=5)


def _add_program_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--program_preset",
        choices=["none", "default", "connectable", "qc", "experimental", "all"],
        default="connectable",
        help=(
            "Program family preset for GMT emission and metadata. "
            "connectable/default emit recommended high-value sets; qc/experimental/all are opt-in."
        ),
    )
    parser.add_argument(
        "--program_methods",
        help=(
            "Optional comma-separated program methods overriding preset. "
            "Bulk: linked_activity,promoter_activity,distal_activity,enhancer_bias. "
            "scATAC adds tfidf_distal. Program methods are evaluated per selected link_method. "
            "Legacy reference methods in this list are still accepted."
        ),
    )
    parser.add_argument(
        "--calibration_methods",
        dest="calibration_methods",
        default="auto_prefer_ref_ubiquity_else_none",
        help=(
            "Comma-separated calibration methods applied independently of link_method. "
            "Tokens: all,auto_prefer_ref_ubiquity_else_none,none,ref_ubiquity_penalty,atlas_residual. "
            "Default auto_prefer_ref_ubiquity_else_none."
        ),
    )
    parser.add_argument(
        "--contrast_methods",
        dest="calibration_methods",
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument("--emit_full", type=_parse_bool, default=True)
    parser.add_argument(
        "--qc_marker_genes_tsv",
        help=(
            "Optional marker list for QC hit counting (one gene symbol/id per line or first column TSV). "
            "When provided, marker hit counts are written to run_summary and metadata."
        ),
    )


def _add_resource_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--resources_manifest", help="Optional resources manifest path (defaults to bundled manifest)")
    parser.add_argument("--resources_dir", help="Optional resources cache directory (defaults to OMICS2GENESET_RESOURCES_DIR or ~/.cache)")
    parser.add_argument(
        "--use_reference_bundle",
        type=_parse_bool,
        default=True,
        help="If true (default), allow reference-backed calibration methods; if false, calibration behavior is forced to none.",
    )
    parser.add_argument("--resource_policy", choices=["skip", "fail"], default="skip")
    parser.add_argument("--ref_ubiquity_resource_id", help="Resource id used by calibration_method=ref_ubiquity_penalty")
    parser.add_argument("--atlas_resource_id", help="Resource id used by calibration_method=atlas_residual")
    parser.add_argument("--atlas_metric", choices=["logratio", "zscore"], default="zscore")
    parser.add_argument("--atlas_eps", type=float, default=1e-6)
    parser.add_argument(
        "--atlas_min_raw_quantile",
        type=float,
        default=0.95,
        help="When atlas_residual is active, drop genes below this raw-score quantile before residual ranking.",
    )
    parser.add_argument(
        "--atlas_use_log1p",
        type=_parse_bool,
        default=True,
        help="When atlas_residual is active, apply log1p stabilization before computing residuals.",
    )


def _add_gmt_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--emit_gmt", type=_parse_bool, default=True)
    parser.add_argument("--gmt_out")
    parser.add_argument("--gmt_prefer_symbol", type=_parse_bool, default=True)
    parser.add_argument("--gmt_require_symbol", type=_parse_bool, default=True)
    parser.add_argument("--gmt_biotype_allowlist", default="protein_coding")
    parser.add_argument("--gmt_min_genes", type=int, default=100)
    parser.add_argument("--gmt_max_genes", type=int, default=500)
    parser.add_argument("--gmt_topk_list", default="200")
    parser.add_argument("--gmt_mass_list", default="")
    parser.add_argument("--gmt_split_signed", type=_parse_bool, default=False)
    parser.add_argument(
        "--emit_small_gene_sets",
        type=_parse_bool,
        default=False,
        help="If true, emit GMT sets smaller than gmt_min_genes (still warns).",
    )


def _add_transform_flags(parser: argparse.ArgumentParser, default: str) -> None:
    parser.add_argument("--score_transform", choices=["signed", "abs", "positive", "negative"], default=default)
    parser.add_argument("--normalize", choices=["l1", "none"], default="l1")


def _add_rna_deg_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--signature_name",
        default="contrast",
        help="Signature label for metadata/GMT naming. If left as 'contrast', uses deg_tsv filename stem.",
    )
    parser.add_argument("--gene_id_column", default="gene_id")
    parser.add_argument("--gene_symbol_column")
    parser.add_argument("--stat_column")
    parser.add_argument("--logfc_column")
    parser.add_argument("--padj_column")
    parser.add_argument("--pvalue_column")
    parser.add_argument("--score_column")
    parser.add_argument(
        "--score_mode",
        choices=["auto", "stat", "logfc_times_neglog10p", "custom_column"],
        default="auto",
    )
    parser.add_argument(
        "--duplicate_gene_policy",
        choices=["sum", "max_abs", "mean", "last"],
        default="max_abs",
        help="How to aggregate duplicate gene_id rows in DEG tables.",
    )
    parser.add_argument("--neglog10p_cap", type=float, default=50.0)
    parser.add_argument("--neglog10p_eps", type=float, default=1e-300)
    parser.add_argument("--exclude_gene_regex", action="append")
    parser.add_argument("--disable_default_excludes", action="store_true")
    parser.add_argument("--gtf")
    parser.add_argument("--gtf_gene_id_field", default="gene_id")
    parser.add_argument("--gtf_source")

    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    parser.add_argument("--emit_full", type=_parse_bool, default=True)

    _add_gmt_flags(parser)
    parser.add_argument(
        "--gmt_emit_abs",
        type=_parse_bool,
        default=False,
        help="If true, emit an additional abs(score)-ranked GMT set family.",
    )
    parser.add_argument(
        "--gmt_source",
        choices=["full", "selected"],
        default="full",
        help="Source table for GMT export: full ranked scores or selected geneset.tsv rows.",
    )
    parser.set_defaults(
        emit_gmt=True,
        gmt_split_signed=True,
        gmt_topk_list="200",
        gmt_min_genes=100,
        gmt_max_genes=500,
    )


def _add_rna_sc_program_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--program_loadings_tsv")
    parser.add_argument(
        "--loadings_format",
        choices=[
            "auto",
            "wide_genes_by_program",
            "wide_programs_by_gene",
            "long_tidy",
            "cnmf_gene_spectra_tpm",
            "cnmf_gene_spectra_score",
            "schpf_gene_scores",
        ],
        default="auto",
    )
    parser.add_argument("--gene_id_column", default="gene_id")
    parser.add_argument("--program_id_column", default="program_id")
    parser.add_argument("--loading_column", default="loading")
    parser.add_argument("--transpose", type=_parse_bool, default=False)

    parser.add_argument("--cnmf_gene_spectra_tsv")
    parser.add_argument("--cnmf_kind", choices=["tpm", "score"])
    parser.add_argument("--schpf_gene_scores_tsv")

    parser.add_argument("--dataset_label")
    parser.add_argument("--signature_name", default="programs")
    parser.add_argument("--program_id_prefix")

    parser.add_argument("--exclude_gene_regex", action="append")
    parser.add_argument("--disable_default_excludes", action="store_true")
    parser.add_argument("--gtf")
    parser.add_argument("--gtf_gene_id_field", default="gene_id")
    parser.add_argument("--gtf_source")

    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument(
        "--score_transform",
        choices=["signed", "abs", "positive", "negative"],
        default="positive",
    )
    parser.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    parser.add_argument("--emit_full", type=_parse_bool, default=True)

    _add_gmt_flags(parser)
    parser.set_defaults(
        emit_gmt=True,
        gmt_split_signed=False,
        gmt_topk_list="200",
        gmt_min_genes=100,
        gmt_max_genes=500,
        gmt_require_symbol=False,
    )


def _add_scrna_cnmf_prepare_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--matrix_tsv", required=True, help="Cell x gene dense TSV (header required).")
    parser.add_argument(
        "--matrix_orientation",
        choices=["auto", "cell_by_gene", "gene_by_cell"],
        default="auto",
        help="Input matrix orientation (default auto).",
    )
    parser.add_argument(
        "--matrix_cell_id_column",
        help="Column name containing cell IDs in matrix_tsv (default: first column).",
    )
    parser.add_argument(
        "--matrix_gene_id_column",
        help="Gene ID column for gene_by_cell orientation (default: first column).",
    )
    parser.add_argument("--matrix_delim", default="\t", help="Matrix delimiter (default: tab).")
    parser.add_argument("--meta_tsv", required=True, help="Cell metadata TSV.")
    parser.add_argument("--meta_cell_id_column", default="cell_id")
    parser.add_argument(
        "--donor_column",
        default="donor_id",
        help="Donor column for donor-aware balancing (default: donor_id if present).",
    )
    parser.add_argument("--phenotype_column", help="Optional phenotype column in metadata.")
    parser.add_argument("--cell_type_column", help="Optional cell-type column in metadata.")
    parser.add_argument(
        "--split_by_cell_type",
        type=_parse_bool,
        default=None,
        help="If true, emit one subset per cell_type. Default: true when --cell_type_column is provided, else false.",
    )
    parser.add_argument("--cell_type_allowlist", help="Optional comma-separated cell types to keep.")
    parser.add_argument("--min_cells_per_cell_type", type=int, default=500)
    parser.add_argument(
        "--bucket_columns",
        help=(
            "Optional comma-separated bucket columns for downsampling. "
            "Default: (cell_type,donor) when split, donor when not split, else global."
        ),
    )
    parser.add_argument("--max_cells_per_bucket", type=int, default=200)
    parser.add_argument("--max_cells_total", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--matrix_value_type", choices=["counts", "logcounts"], default="logcounts")
    parser.add_argument("--min_total_per_cell", type=float)
    parser.add_argument("--min_total_per_gene", type=float)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--keep_tmp", type=_parse_bool, default=False)

    parser.add_argument("--cnmf_name", help="Optional base cNMF run name per subset.")
    parser.add_argument(
        "--cnmf_k_list",
        default="auto",
        help=(
            "K grid for cnmf prepare. Use 'auto' (default) for subset-size-dependent K grid, "
            "or provide explicit values, e.g. '10 15 20'."
        ),
    )
    parser.add_argument(
        "--cnmf_k",
        default="auto",
        help="K for consensus script. Use 'auto' (default) to call workflows cnmf_select_k, or an integer.",
    )
    parser.add_argument(
        "--cnmf_select_strategy",
        choices=["largest_stable", "max_stability", "manual"],
        default="largest_stable",
        help="Auto-K strategy used by generated consensus/conversion scripts when --cnmf_k=auto.",
    )
    parser.add_argument("--cnmf_select_stability_frac_of_max", type=float, default=0.95)
    parser.add_argument("--cnmf_select_min_stability_abs", type=float, default=0.0)
    parser.add_argument("--cnmf_select_require_local_max", type=_parse_bool, default=False)
    parser.add_argument("--cnmf_n_iter", type=int, default=100)
    parser.add_argument("--cnmf_numgenes", type=int, default=2000)
    parser.add_argument("--cnmf_total_workers", type=int, default=1)
    parser.add_argument(
        "--cnmf_densify",
        type=_parse_bool,
        default=True,
        help="Whether generated cnmf prepare command includes --densify.",
    )
    parser.add_argument("--write_postprocess_template", type=_parse_bool, default=True)
    parser.add_argument("--execute", type=_parse_bool, default=False)
    parser.add_argument("--cnmf_local_density_threshold", type=float, default=0.5)
    parser.add_argument("--cnmf_local_neighborhood_size", type=float, default=0.3)
    parser.add_argument("--cnmf_show_clustering", type=_parse_bool, default=False)
    parser.add_argument("--cnmf_export_kind", choices=["tpm", "score"], default="tpm")
    parser.add_argument("--organism", choices=["human", "mouse"], default="human")
    parser.add_argument("--genome_build", default="hg38")


def _add_cnmf_select_k_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--cnmf_output_dir", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--stats_tsv", help="Optional TSV fallback with columns k,stability,prediction_error.")
    parser.add_argument("--strategy", choices=["largest_stable", "max_stability", "manual"], default="largest_stable")
    parser.add_argument("--stability_frac_of_max", type=float, default=0.95)
    parser.add_argument("--min_stability_abs", type=float, default=0.0)
    parser.add_argument("--require_local_max", type=_parse_bool, default=False)
    parser.add_argument("--fixed_k", type=int)


def _add_methylation_program_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--program_preset",
        choices=["none", "default", "connectable", "qc", "experimental", "all"],
        default="connectable",
        help="Methylation program preset (connectable/default recommended).",
    )
    parser.add_argument(
        "--program_methods",
        help=(
            "Optional comma-separated methylation program methods overriding preset. "
            "Methods: promoter_activity,distal_activity,linked_activity."
        ),
    )
    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    parser.add_argument(
        "--aggregation",
        choices=["weighted_mean", "sum", "mean", "max_abs"],
        default="weighted_mean",
        help="Gene-score aggregation across linked CpGs/DMRs.",
    )
    parser.add_argument("--score_transform", choices=["signed", "abs", "positive", "negative"], default="signed")


def _add_methylation_resource_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--resources_manifest", help="Optional resources manifest path")
    parser.add_argument("--resources_dir", help="Optional resources directory")
    parser.add_argument("--resource_policy", choices=["skip", "fail"], default="skip")
    parser.add_argument("--enhancer_resource_id", help="Optional resource id for enhancer BED")


def _add_methylation_common_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--dataset_label")
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--gtf_source")
    parser.add_argument("--distal_mode", choices=["nonpromoter", "enhancer_only"], default="nonpromoter")
    parser.add_argument("--enhancer_bed")
    parser.add_argument(
        "--score_mode",
        choices=["delta_times_neglog10p", "delta_only"],
        default="delta_times_neglog10p",
        help="Methylation feature scoring mode before program extraction.",
    )
    parser.add_argument(
        "--delta_orientation",
        choices=["activity_oriented", "raw"],
        default="activity_oriented",
        help=(
            "How delta is interpreted before scoring. "
            "activity_oriented uses -delta (default); raw uses +delta."
        ),
    )
    parser.add_argument("--neglog10p_cap", type=float, default=50.0)
    parser.add_argument("--neglog10p_eps", type=float, default=1e-300)
    parser.add_argument("--drop_sex_chrom", type=_parse_bool, default=True)
    parser.add_argument(
        "--exclude_gene_symbol_regex",
        action="append",
        help="Repeatable regex filter for gene_symbol/gene_id exclusion before selection and GMT.",
    )
    parser.add_argument(
        "--exclude_gene_symbols_tsv",
        help="Optional one-symbol-per-line exclusion list applied before selection and GMT.",
    )
    _add_linking_flags(parser)
    _add_methylation_program_flags(parser)
    _add_methylation_resource_flags(parser)
    _add_gmt_flags(parser)
    parser.add_argument("--emit_full", type=_parse_bool, default=True)
    parser.set_defaults(
        emit_gmt=True,
        gmt_split_signed=True,
        gmt_topk_list="200",
        gmt_min_genes=100,
        gmt_max_genes=500,
        emit_small_gene_sets=False,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="omics2geneset")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("list")

    p_desc = sub.add_parser("describe")
    p_desc.add_argument("converter")

    p_val = sub.add_parser("validate")
    p_val.add_argument("output_dir")

    p_convert = sub.add_parser("convert")
    conv = p_convert.add_subparsers(dest="converter", required=True)

    p_workflows = sub.add_parser("workflows")
    wf_sub = p_workflows.add_subparsers(dest="workflow_command", required=True)
    p_scrna_cnmf_prepare = wf_sub.add_parser("scrna_cnmf_prepare")
    _add_scrna_cnmf_prepare_flags(p_scrna_cnmf_prepare)
    p_cnmf_select_k = wf_sub.add_parser("cnmf_select_k")
    _add_cnmf_select_k_flags(p_cnmf_select_k)

    p_resources = sub.add_parser("resources")
    res_sub = p_resources.add_subparsers(dest="resources_command", required=True)
    p_res_list = res_sub.add_parser("list")
    p_res_list.add_argument("--manifest")
    p_res_list.add_argument(
        "--manifest_mode",
        choices=["overlay", "replace"],
        default="overlay",
        help="When --manifest is provided, overlay merges with bundled manifest (default) or replace uses only provided manifest.",
    )

    p_res_describe = res_sub.add_parser("describe")
    p_res_describe.add_argument("resource_id")
    p_res_describe.add_argument("--manifest")
    p_res_describe.add_argument(
        "--manifest_mode",
        choices=["overlay", "replace"],
        default="overlay",
        help="When --manifest is provided, overlay merges with bundled manifest (default) or replace uses only provided manifest.",
    )
    p_res_describe.add_argument("--resources_dir")
    p_res_describe.add_argument("--verify", action="store_true", help="Compute checksum for local file when describing status")

    p_res_status = res_sub.add_parser("status")
    p_res_status.add_argument("--manifest")
    p_res_status.add_argument(
        "--manifest_mode",
        choices=["overlay", "replace"],
        default="overlay",
        help="When --manifest is provided, overlay merges with bundled manifest (default) or replace uses only provided manifest.",
    )
    p_res_status.add_argument("--resources_dir")
    p_res_status.add_argument("--verify", action="store_true", help="Compute checksums for local files")
    p_res_status.add_argument("--fast", action="store_true", help="Alias for default fast mode (existence + size; checksums optional)")
    p_res_status.add_argument(
        "--check_schema",
        action="store_true",
        help="Run schema checks for known resource ids (for example ccre_ubiquity_* and atac_reference_profiles_*)",
    )

    p_res_fetch = res_sub.add_parser("fetch")
    p_res_fetch.add_argument("resource_ids", nargs="*")
    p_res_fetch.add_argument("--preset")
    p_res_fetch.add_argument("--manifest")
    p_res_fetch.add_argument(
        "--manifest_mode",
        choices=["overlay", "replace"],
        default="overlay",
        help="When --manifest is provided, overlay merges with bundled manifest (default) or replace uses only provided manifest.",
    )
    p_res_fetch.add_argument("--resources_dir")
    p_res_fetch.add_argument("--overwrite", action="store_true")
    p_res_fetch.add_argument(
        "--skip_missing_urls",
        type=_parse_bool,
        default=True,
        help="Skip resources without URL as manual entries (default true).",
    )

    p_res_validate = res_sub.add_parser("manifest-validate")
    p_res_validate.add_argument("--manifest")
    p_res_validate.add_argument(
        "--manifest_mode",
        choices=["overlay", "replace"],
        default="overlay",
        help="When --manifest is provided, overlay merges with bundled manifest (default) or replace uses only provided manifest.",
    )
    p_res_validate.add_argument(
        "--strict",
        action="store_true",
        help="Treat manifest warnings (e.g. download URL without sha256) as errors.",
    )

    p_bulk = conv.add_parser("atac_bulk")
    p_bulk.add_argument("--peaks", required=True)
    p_bulk.add_argument("--gtf", required=True)
    p_bulk.add_argument("--out_dir", required=True)
    p_bulk.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_bulk.add_argument("--genome_build", required=True)
    p_bulk.add_argument("--dataset_label")
    p_bulk.add_argument("--peaks_weight_column", type=int, default=5)
    p_bulk.add_argument("--peak_weights_tsv")
    p_bulk.add_argument("--peak_weight_transform", choices=["signed", "abs", "positive", "negative"], default="abs")
    p_bulk.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    p_bulk.add_argument("--gtf_source")
    _add_linking_flags(p_bulk)
    _add_program_flags(p_bulk)
    _add_resource_flags(p_bulk)
    _add_gmt_flags(p_bulk)

    p_bulk_matrix = conv.add_parser("atac_bulk_matrix")
    p_bulk_matrix.add_argument("--peak_matrix_tsv", required=True)
    p_bulk_matrix.add_argument("--peak_bed", required=True)
    p_bulk_matrix.add_argument("--sample_metadata_tsv", required=True)
    p_bulk_matrix.add_argument("--gtf", required=True)
    p_bulk_matrix.add_argument("--out_dir", required=True)
    p_bulk_matrix.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_bulk_matrix.add_argument("--genome_build", required=True)
    p_bulk_matrix.add_argument("--dataset_label")
    p_bulk_matrix.add_argument("--sample_id_column", default="sample_id")
    p_bulk_matrix.add_argument("--condition_column", default="condition")
    p_bulk_matrix.add_argument("--condition_a")
    p_bulk_matrix.add_argument("--condition_b")
    p_bulk_matrix.add_argument("--contrast_metric", choices=["log2fc", "diff"], default="log2fc")
    p_bulk_matrix.add_argument("--contrast_pseudocount", type=float, default=1.0)
    p_bulk_matrix.add_argument("--peak_weight_transform", choices=["signed", "abs", "positive", "negative"], default="positive")
    p_bulk_matrix.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    p_bulk_matrix.add_argument("--gtf_source")
    _add_linking_flags(p_bulk_matrix)
    _add_program_flags(p_bulk_matrix)
    _add_resource_flags(p_bulk_matrix)
    _add_gmt_flags(p_bulk_matrix)

    p_sc = conv.add_parser("atac_sc_10x")
    p_sc.add_argument("--matrix_dir", required=True)
    p_sc.add_argument("--gtf", required=True)
    p_sc.add_argument("--out_dir", required=True)
    p_sc.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_sc.add_argument("--genome_build", required=True)
    p_sc.add_argument("--dataset_label")
    p_sc.add_argument("--groups_tsv")
    p_sc.add_argument("--cell_metadata_tsv")
    p_sc.add_argument("--condition_column")
    p_sc.add_argument("--donor_column", default="donor")
    p_sc.add_argument("--min_donors_per_condition", type=int, default=3)
    p_sc.add_argument("--min_total_donors_per_condition", type=int, default=3)
    p_sc.add_argument("--cell_imbalance_warn_ratio", type=float, default=5.0)
    p_sc.add_argument("--baseline_carry_through_corr_warn", type=float, default=0.9)
    p_sc.add_argument("--condition_a")
    p_sc.add_argument("--condition_b")
    p_sc.add_argument("--min_cells_per_group", type=int, default=100)
    p_sc.add_argument("--min_cells_per_condition", type=int, default=50)
    p_sc.add_argument("--peak_summary", choices=["sum_counts", "mean_counts", "frac_cells_nonzero"], default="sum_counts")
    p_sc.add_argument("--peak_weight_transform", choices=["signed", "abs", "positive", "negative"], default="positive")
    p_sc.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    p_sc.add_argument(
        "--study_contrast",
        dest="study_contrast",
        choices=["none", "group_vs_rest", "condition_within_group"],
        help="Study-design contrast used to define peak statistics before calibration.",
    )
    p_sc.add_argument("--contrast", dest="study_contrast", default=argparse.SUPPRESS, help=argparse.SUPPRESS)
    p_sc.add_argument("--contrast_metric", choices=["log2fc", "diff"], default="log2fc")
    p_sc.add_argument("--contrast_pseudocount", type=float)
    p_sc.add_argument("--gtf_source")
    _add_linking_flags(p_sc)
    _add_program_flags(p_sc)
    _add_resource_flags(p_sc)
    _add_gmt_flags(p_sc)

    p_rna = conv.add_parser("rna_deg")
    p_rna.add_argument("--deg_tsv", required=True)
    p_rna.add_argument("--out_dir", required=True)
    p_rna.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_rna.add_argument("--genome_build", required=True)
    _add_rna_deg_flags(p_rna)

    p_rna_multi = conv.add_parser("rna_deg_multi")
    p_rna_multi.add_argument("--deg_tsv", required=True)
    p_rna_multi.add_argument("--comparison_column", required=True)
    p_rna_multi.add_argument("--out_dir", required=True)
    p_rna_multi.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_rna_multi.add_argument("--genome_build", required=True)
    _add_rna_deg_flags(p_rna_multi)

    p_rna_sc_programs = conv.add_parser("rna_sc_programs")
    p_rna_sc_programs.add_argument("--out_dir", required=True)
    p_rna_sc_programs.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_rna_sc_programs.add_argument("--genome_build", required=True)
    _add_rna_sc_program_flags(p_rna_sc_programs)

    p_cnv = conv.add_parser("cnv_gene_extractor")
    p_cnv.add_argument("--segments_tsv", required=True)
    p_cnv.add_argument("--out_dir", required=True)
    p_cnv.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_cnv.add_argument("--genome_build", required=True)
    p_cnv.add_argument("--dataset_label")
    p_cnv.add_argument("--gtf", required=True)
    p_cnv.add_argument("--gtf_source")
    p_cnv.add_argument("--gtf_gene_id_field", default="gene_id")
    p_cnv.add_argument("--chrom_column")
    p_cnv.add_argument("--start_column")
    p_cnv.add_argument("--end_column")
    p_cnv.add_argument("--amplitude_column")
    p_cnv.add_argument("--sample_id_column")
    p_cnv.add_argument(
        "--coord_system",
        choices=["one_based_closed", "zero_based_half_open"],
        default="one_based_closed",
    )
    p_cnv.add_argument(
        "--chrom_prefix_mode",
        choices=["auto", "add_chr", "drop_chr", "none"],
        default="auto",
    )
    p_cnv.add_argument("--purity_tsv")
    p_cnv.add_argument("--purity_sample_id_column", default="sample_id")
    p_cnv.add_argument("--purity_value_column", default="purity")
    p_cnv.add_argument(
        "--use_purity_correction",
        nargs="?",
        const=True,
        type=_parse_bool,
        default=False,
        help="Enable purity correction (supports --use_purity_correction or --use_purity_correction true/false).",
    )
    p_cnv.add_argument(
        "--use-purity-correction",
        dest="use_purity_correction",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    p_cnv.add_argument(
        "--no-use-purity-correction",
        dest="use_purity_correction",
        action="store_false",
        help="Disable purity correction.",
    )
    p_cnv.add_argument("--purity_floor", type=float, default=0.1)
    p_cnv.add_argument("--max_abs_amplitude", type=float, default=3.0)
    p_cnv.add_argument("--min_abs_amplitude", type=float, default=0.10)
    p_cnv.add_argument("--focal_length_scale_bp", type=float, default=10_000_000)
    p_cnv.add_argument("--focal_length_alpha", type=float, default=1.0)
    p_cnv.add_argument("--gene_count_penalty", choices=["none", "inv_sqrt", "inv"], default="inv_sqrt")
    p_cnv.add_argument("--aggregation", choices=["weighted_mean", "sum", "mean", "max_abs"], default="weighted_mean")
    p_cnv.add_argument(
        "--program_preset",
        choices=["default", "connectable", "broad", "qc", "all", "none"],
        default="default",
    )
    p_cnv.add_argument("--program_methods")
    p_cnv.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    p_cnv.add_argument("--top_k", type=int, default=200)
    p_cnv.add_argument("--quantile", type=float, default=0.01)
    p_cnv.add_argument("--min_score", type=float, default=0.0)
    p_cnv.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    p_cnv.add_argument("--emit_full", type=_parse_bool, default=True)
    p_cnv.add_argument(
        "--emit_cohort_sets",
        nargs="?",
        const=True,
        type=_parse_bool,
        default=False,
        help="Enable cohort recurrence gene sets (supports --emit_cohort_sets or --emit_cohort_sets true/false).",
    )
    p_cnv.add_argument(
        "--emit-cohort-sets",
        dest="emit_cohort_sets",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    p_cnv.add_argument(
        "--no-emit-cohort-sets",
        dest="emit_cohort_sets",
        action="store_false",
        help="Disable cohort recurrence gene sets.",
    )
    p_cnv.add_argument("--cohort_score_threshold", type=float, default=0.15)
    p_cnv.add_argument("--cohort_min_fraction", type=float, default=0.05)
    p_cnv.add_argument("--cohort_min_samples", type=int, default=5)
    _add_gmt_flags(p_cnv)
    p_cnv.set_defaults(
        emit_gmt=True,
        gmt_split_signed=False,
        gmt_topk_list="200",
        gmt_min_genes=100,
        gmt_max_genes=500,
        emit_small_gene_sets=False,
    )

    p_chip = conv.add_parser("chipseq_peak")
    p_chip.add_argument("--peaks", required=True)
    p_chip.add_argument("--gtf", required=True)
    p_chip.add_argument("--out_dir", required=True)
    p_chip.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_chip.add_argument("--genome_build", required=True)
    p_chip.add_argument("--weight_column", type=int, default=5)
    p_chip.add_argument("--peak_weight_transform", choices=["signed", "abs", "positive", "negative"], default="positive")
    p_chip.add_argument("--normalize", choices=["l1", "none"], default="l1")
    p_chip.add_argument("--gtf_source")
    _add_linking_flags(p_chip)
    p_chip.set_defaults(link_method="promoter_overlap")

    p_meth_cpg = conv.add_parser("methylation_cpg_diff")
    p_meth_cpg.add_argument("--cpg_tsv", required=True)
    p_meth_cpg.add_argument("--out_dir", required=True)
    p_meth_cpg.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_meth_cpg.add_argument("--genome_build", required=True)
    p_meth_cpg.add_argument("--array_type", choices=["450k", "epic"], default="450k")
    p_meth_cpg.add_argument("--probe_id_column", default="probe_id")
    p_meth_cpg.add_argument("--chrom_column", default="chrom")
    p_meth_cpg.add_argument("--pos_column", default="pos")
    p_meth_cpg.add_argument("--start_column")
    p_meth_cpg.add_argument("--end_column")
    p_meth_cpg.add_argument("--input_pos_is_0based", type=_parse_bool, default=False)
    p_meth_cpg.add_argument("--delta_column", default="delta_beta")
    p_meth_cpg.add_argument("--padj_column")
    p_meth_cpg.add_argument("--pvalue_column")
    p_meth_cpg.add_argument("--probe_manifest_tsv")
    p_meth_cpg.add_argument("--probe_manifest_resource_id")
    p_meth_cpg.add_argument("--manifest_probe_id_column", default="probe_id")
    p_meth_cpg.add_argument("--manifest_chrom_column", default="chrom")
    p_meth_cpg.add_argument("--manifest_pos_column", default="pos")
    p_meth_cpg.add_argument("--manifest_start_column")
    p_meth_cpg.add_argument("--manifest_end_column")
    p_meth_cpg.add_argument("--manifest_pos_is_0based", type=_parse_bool, default=False)
    p_meth_cpg.add_argument("--probe_blacklist_tsv")
    _add_methylation_common_flags(p_meth_cpg)

    p_meth_dmr_regions = conv.add_parser("methylation_dmr_regions")
    p_meth_dmr_regions.add_argument("--dmr_tsv", required=True)
    p_meth_dmr_regions.add_argument("--out_dir", required=True)
    p_meth_dmr_regions.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_meth_dmr_regions.add_argument("--genome_build", required=True)
    p_meth_dmr_regions.add_argument("--chrom_column", default="chrom")
    p_meth_dmr_regions.add_argument("--start_column", default="start")
    p_meth_dmr_regions.add_argument("--end_column", default="end")
    p_meth_dmr_regions.add_argument("--delta_column", default="delta_methylation")
    p_meth_dmr_regions.add_argument("--padj_column")
    p_meth_dmr_regions.add_argument("--pvalue_column")
    _add_methylation_common_flags(p_meth_dmr_regions)

    p_meth = conv.add_parser("methylation_dmr")
    p_meth.add_argument("--dmr_tsv", required=True)
    p_meth.add_argument("--out_dir", required=True)
    p_meth.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_meth.add_argument("--genome_build", required=True)
    p_meth.add_argument("--gene_id_column", default="gene_id")
    p_meth.add_argument("--score_column", default="delta_methylation")
    _add_transform_flags(p_meth, default="abs")

    p_prot = conv.add_parser("proteomics_diff")
    p_prot.add_argument("--proteomics_tsv", required=True)
    p_prot.add_argument("--out_dir", required=True)
    p_prot.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_prot.add_argument("--genome_build", required=True)
    p_prot.add_argument("--gene_id_column", default="gene_id")
    p_prot.add_argument("--score_column", default="log2fc")
    _add_transform_flags(p_prot, default="abs")

    p_scrna = conv.add_parser("sc_rna_marker")
    p_scrna.add_argument("--counts_tsv", required=True)
    p_scrna.add_argument("--out_dir", required=True)
    p_scrna.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_scrna.add_argument("--genome_build", required=True)
    p_scrna.add_argument("--groups_tsv")
    p_scrna.add_argument("--gene_id_column", default="gene_id")
    p_scrna.add_argument("--barcode_column", default="barcode")
    p_scrna.add_argument("--value_column", default="count")
    p_scrna.add_argument("--group_barcode_column", default="barcode")
    p_scrna.add_argument("--group_column", default="group")
    p_scrna.add_argument("--peak_summary", choices=["sum_counts", "mean_counts", "frac_cells_nonzero"], default="sum_counts")
    p_scrna.add_argument("--normalize", choices=["l1", "none"], default="l1")

    return parser


def _flag_used(argv: list[str], flag: str) -> bool:
    for token in argv:
        if token == flag or token.startswith(flag + "="):
            return True
    return False


def _emit_cli_alias_warnings(argv: list[str]) -> None:
    if _flag_used(argv, "--contrast_methods"):
        print(
            "warning: --contrast_methods is deprecated; use --calibration_methods.",
            file=sys.stderr,
        )
    if _flag_used(argv, "--contrast"):
        print(
            "warning: --contrast is deprecated; use --study_contrast.",
            file=sys.stderr,
        )


def main(argv: list[str] | None = None) -> int:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(raw_argv)
    _emit_cli_alias_warnings(raw_argv)
    try:
        if args.command == "list":
            for name, desc in list_converters():
                print(f"{name}\t{desc}")
            return 0

        if args.command == "describe":
            spec = get_converter_spec(args.converter)
            print(json.dumps(spec, indent=2, sort_keys=True))
            return 0

        if args.command == "validate":
            schema = Path(__file__).parent / "schemas" / "geneset_metadata.schema.json"
            result = validate_output_dir(args.output_dir, schema)
            if result.get("mode") == "grouped":
                print(f"ok: validated grouped output with n_groups={result.get('n_groups')}")
            else:
                print("ok")
            return 0

        if args.command == "resources":
            merge_with_bundled = bool(getattr(args, "manifest_mode", "overlay") == "overlay")
            manifest_path, resources, presets, manifest_warnings = load_manifest(
                args.manifest,
                merge_with_bundled=merge_with_bundled,
            )
            for warning in manifest_warnings:
                print(f"warning: {warning}", file=sys.stderr)

            if args.resources_command == "list":
                for rid in sorted(resources):
                    entry = resources[rid]
                    print(
                        "\t".join(
                            [
                                rid,
                                str(entry.get("genome_build", "")),
                                str(entry.get("provider", "")),
                                str(entry.get("filename", "")),
                                ("download" if str(entry.get("url", "")).strip() else "manual"),
                            ]
                        )
                    )
                print(f"# manifest={manifest_path}", file=sys.stderr)
                return 0

            if args.resources_command == "describe":
                payload = describe_resource(
                    resources=resources,
                    resource_id=args.resource_id,
                    resources_dir=args.resources_dir,
                    verify=bool(args.verify),
                )
                print(json.dumps(payload, indent=2, sort_keys=True))
                print(f"# manifest={manifest_path}", file=sys.stderr)
                return 0

            if args.resources_command == "status":
                if args.fast and args.verify:
                    raise ValueError("status accepts either --fast or --verify, not both")
                verify = bool(args.verify)
                rows = resource_status_rows(
                    resources,
                    args.resources_dir,
                    verify=verify,
                    check_schema=bool(args.check_schema),
                )
                for row in rows:
                    schema_error = str(row.get("schema_error", "") or "").replace("\t", " ").replace("\n", " ")
                    print(
                        "\t".join(
                            [
                                str(row["id"]),
                                str(row["status"]),
                                str(row["path"]),
                                str(row.get("availability", "")),
                                str(row.get("size_bytes", "")),
                                str(row.get("schema_status", "")),
                                str(row.get("schema_kind", "")),
                                str(row.get("schema_rows", "")),
                                schema_error,
                            ]
                        )
                    )
                print(f"# manifest={manifest_path}", file=sys.stderr)
                return 0

            if args.resources_command == "fetch":
                selected = resolve_requested_resource_ids(resources, presets, list(args.resource_ids), args.preset)
                results = fetch_resources(
                    resources=resources,
                    resource_ids=selected,
                    resources_dir=args.resources_dir,
                    overwrite=bool(args.overwrite),
                    skip_missing_urls=bool(args.skip_missing_urls),
                )
                for row in results:
                    print(
                        "\t".join(
                            [
                                str(row["id"]),
                                str(row["status"]),
                                str(row["path"]),
                                str(row.get("reason", "")),
                            ]
                        )
                    )
                manual = [row for row in results if str(row.get("status")) == "manual"]
                if manual:
                    ids = ",".join(sorted(str(row["id"]) for row in manual))
                    print(
                        "warning: some resources are manual (no url in manifest) and were skipped: "
                        f"{ids}",
                        file=sys.stderr,
                    )
                print(f"# manifest={manifest_path}", file=sys.stderr)
                return 0

            if args.resources_command == "manifest-validate":
                if manifest_warnings:
                    for warning in manifest_warnings:
                        print(f"warning: {warning}", file=sys.stderr)
                if manifest_warnings and bool(args.strict):
                    raise ValueError("manifest validation failed in strict mode due to warnings")
                print("ok")
                print(f"# manifest={manifest_path}", file=sys.stderr)
                return 0

        if args.command == "workflows":
            if args.workflow_command == "scrna_cnmf_prepare":
                from omics2geneset.workflows.scrna_cnmf_prepare import run as run_scrna_cnmf_prepare

                result = run_scrna_cnmf_prepare(args)
                print(
                    "workflow_completed "
                    f"workflow=scrna_cnmf_prepare n_subsets={result.get('n_subsets')} "
                    f"out={result.get('out_dir')}",
                    file=sys.stderr,
                )
                return 0
            if args.workflow_command == "cnmf_select_k":
                from omics2geneset.workflows.cnmf_select_k import run as run_cnmf_select_k

                _ = run_cnmf_select_k(args)
                return 0

        if args.command == "convert":
            converter = get_converter(args.converter)
            result = converter.run(args)
            if int(result.get("n_groups", 1)) > 1:
                print(
                    "converted "
                    f"peaks={result.get('n_peaks')} "
                    f"groups={result.get('n_groups')} "
                    f"genes_per_group={result.get('n_genes_min')}-{result.get('n_genes_max')} "
                    f"unique_genes={result.get('n_genes_unique')} "
                    f"out={result.get('out_dir')}",
                    file=sys.stderr,
                )
            else:
                print(
                    f"converted peaks={result.get('n_peaks')} genes={result.get('n_genes')} out={result.get('out_dir')}",
                    file=sys.stderr,
                )
            active_methods = result.get("program_methods")
            if isinstance(active_methods, list) and active_methods:
                print("program_methods_active=" + ",".join(str(x) for x in active_methods), file=sys.stderr)
                skipped_methods = result.get("program_methods_skipped")
                if isinstance(skipped_methods, dict):
                    if skipped_methods:
                        print("program_methods_skipped=" + ",".join(sorted(str(k) for k in skipped_methods)), file=sys.stderr)
                    else:
                        print("program_methods_skipped=none", file=sys.stderr)
            active_contrasts = result.get("calibration_methods")
            if isinstance(active_contrasts, list) and active_contrasts:
                print("calibration_methods_active=" + ",".join(str(x) for x in active_contrasts), file=sys.stderr)
                skipped_contrasts = result.get("calibration_methods_skipped")
                if isinstance(skipped_contrasts, dict):
                    if skipped_contrasts:
                        print(
                            "calibration_methods_skipped=" + ",".join(sorted(str(k) for k in skipped_contrasts)),
                            file=sys.stderr,
                        )
                    else:
                        print("calibration_methods_skipped=none", file=sys.stderr)
            return 0
    except Exception as exc:  # pragma: no cover
        print(f"error: {exc}", file=sys.stderr)
        return 1

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
