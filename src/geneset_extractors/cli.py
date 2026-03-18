from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from geneset_extractors.core.validate import validate_output_dir
from geneset_extractors.resource_manager import (
    describe_resource,
    fetch_resources,
    load_manifest,
    resource_status_rows,
    resolve_requested_resource_ids,
)
from geneset_extractors.registry import get_converter, get_converter_spec, list_converters


def _parse_bool(value: str) -> bool:
    v = value.strip().lower()
    if v in {"true", "1", "yes", "y"}:
        return True
    if v in {"false", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError("expected true or false")


class _StoreWithExplicitFlag(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):  # type: ignore[override]
        setattr(namespace, self.dest, values)
        setattr(namespace, f"{self.dest}_explicit", True)


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
    parser.add_argument(
        "--resources_dir",
        help=(
            "Optional resources cache directory (defaults to "
            "GENESET_EXTRACTORS_RESOURCES_DIR, then OMICS2GENESET_RESOURCES_DIR, then ~/.cache)"
        ),
    )
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
    parser.add_argument("--gmt_format", choices=["dig2col", "classic"], default="dig2col")
    parser.add_argument(
        "--emit_small_gene_sets",
        type=_parse_bool,
        default=False,
        help="If true, emit GMT sets smaller than gmt_min_genes (still warns).",
    )


def _add_provenance_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--provenance_overlay_json",
        help="Optional JSON overlay that adds canonical URIs, public URLs, and operation replay metadata to emitted provenance.",
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
        choices=["auto", "stat", "logfc_times_neglog10p", "signed_neglog10padj", "custom_column"],
        default="auto",
    )
    parser.add_argument(
        "--postprocess_mode",
        choices=["harmonizome", "legacy"],
        default="harmonizome",
        help="RNA DEG post-processing preset. harmonizome is the default; legacy preserves prior behavior.",
    )
    parser.add_argument("--padj_max", type=float, help="Optional adjusted p-value filter applied before gene aggregation.")
    parser.add_argument("--pvalue_max", type=float, help="Optional p-value filter applied before gene aggregation.")
    parser.add_argument("--min_abs_logfc", type=float, help="Optional abs(logFC) filter applied before gene aggregation.")
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


def _add_ptm_site_diff_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--signature_name")
    parser.add_argument("--dataset_label")
    parser.add_argument("--ptm_type", choices=["phospho", "acetyl", "ub", "generic"], default="phospho")
    parser.add_argument("--site_id_column")
    parser.add_argument("--site_group_column")
    parser.add_argument("--gene_id_column")
    parser.add_argument("--gene_symbol_column")
    parser.add_argument("--protein_accession_column")
    parser.add_argument("--residue_column")
    parser.add_argument("--position_column")
    parser.add_argument("--score_column")
    parser.add_argument("--stat_column")
    parser.add_argument("--logfc_column")
    parser.add_argument("--padj_column")
    parser.add_argument("--pvalue_column")
    parser.add_argument("--localization_prob_column")
    parser.add_argument("--peptide_count_column")
    parser.add_argument("--protein_logfc_column")
    parser.add_argument("--protein_stat_column")
    parser.add_argument(
        "--score_mode",
        choices=["auto", "stat", "logfc_times_neglog10p", "custom_column"],
        default="auto",
    )
    parser.add_argument("--score_transform", choices=["signed", "abs", "positive", "negative"], default="signed")
    parser.add_argument("--protein_adjustment", choices=["none", "subtract", "residual"], default="subtract")
    parser.add_argument("--protein_adjustment_lambda", type=float, default=1.0)
    parser.add_argument(
        "--confidence_weight_mode",
        choices=["none", "pvalue", "localization", "combined"],
        default="combined",
    )
    parser.add_argument("--min_localization_prob", type=float, default=0.75)
    parser.add_argument(
        "--site_dup_policy",
        choices=["highest_confidence", "max_abs", "mean", "sum"],
        default="highest_confidence",
    )
    parser.add_argument(
        "--gene_aggregation",
        choices=["signed_topk_mean", "max_abs", "sum", "mean"],
        default="signed_topk_mean",
    )
    parser.add_argument("--gene_topk_sites", type=int, default=3)
    parser.add_argument("--ambiguous_gene_policy", choices=["drop", "split_equal", "first"], default="drop")
    parser.add_argument("--neglog10p_cap", type=float, default=50.0)
    parser.add_argument("--neglog10p_eps", type=float, default=1e-300)
    parser.add_argument("--resources_manifest")
    parser.add_argument("--resources_dir")
    parser.add_argument("--resource_policy", choices=["skip", "fail"], default="skip")
    parser.add_argument("--use_reference_bundle", type=_parse_bool, default=True)
    parser.add_argument("--site_alias_resource_id")
    parser.add_argument("--site_ubiquity_resource_id")
    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    parser.add_argument("--emit_full", type=_parse_bool, default=True)
    _add_gmt_flags(parser)
    parser.set_defaults(
        emit_gmt=True,
        gmt_split_signed=True,
        gmt_topk_list="200",
        gmt_min_genes=100,
        gmt_max_genes=500,
    )


def _add_ptm_site_matrix_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--signature_name")
    parser.add_argument("--dataset_label")
    parser.add_argument("--ptm_type", choices=["phospho", "acetyl", "ub", "generic"], default="phospho")
    parser.add_argument("--matrix_format", choices=["wide_sites_by_sample"], default="wide_sites_by_sample")
    parser.add_argument("--sample_metadata_tsv", required=True)
    parser.add_argument("--protein_matrix_tsv")
    parser.add_argument("--sample_id_column", default="sample_id")
    parser.add_argument("--group_column", default="group")
    parser.add_argument("--condition_column", default="condition")
    parser.add_argument(
        "--study_contrast",
        choices=["condition_a_vs_b", "group_vs_rest", "condition_within_group", "baseline"],
        default="condition_a_vs_b",
    )
    parser.add_argument("--condition_a")
    parser.add_argument("--condition_b")
    parser.add_argument("--min_samples_per_condition", type=int, default=2)
    parser.add_argument("--effect_metric", choices=["welch_t", "mean_diff"], default="welch_t")
    parser.add_argument("--missing_value_policy", choices=["drop", "min_present"], default="min_present")
    parser.add_argument("--min_present_per_condition", type=int, default=2)
    parser.add_argument("--site_id_column")
    parser.add_argument("--site_group_column")
    parser.add_argument("--gene_id_column")
    parser.add_argument("--gene_symbol_column")
    parser.add_argument("--protein_accession_column")
    parser.add_argument("--residue_column")
    parser.add_argument("--position_column")
    parser.add_argument("--score_column")
    parser.add_argument("--stat_column", default="stat")
    parser.add_argument("--logfc_column", default="log2fc")
    parser.add_argument("--padj_column", default="padj")
    parser.add_argument("--pvalue_column", default="pvalue")
    parser.add_argument("--localization_prob_column")
    parser.add_argument("--peptide_count_column")
    parser.add_argument("--protein_logfc_column", default="protein_log2fc")
    parser.add_argument("--protein_stat_column", default="protein_stat")
    parser.add_argument(
        "--score_mode",
        choices=["auto", "stat", "logfc_times_neglog10p", "custom_column"],
        default="auto",
    )
    parser.add_argument("--score_transform", choices=["signed", "abs", "positive", "negative"], default="signed")
    parser.add_argument("--protein_adjustment", choices=["none", "subtract", "residual"], default="subtract")
    parser.add_argument(
        "--protein_adjustment_run_mode",
        choices=["single", "compare_if_protein", "compare"],
        default="compare_if_protein",
    )
    parser.add_argument("--protein_adjustment_lambda", type=float, default=1.0)
    parser.add_argument(
        "--confidence_weight_mode",
        choices=["none", "pvalue", "localization", "combined"],
        default="combined",
    )
    parser.add_argument("--min_localization_prob", type=float, default=0.75)
    parser.add_argument(
        "--site_dup_policy",
        choices=["highest_confidence", "max_abs", "mean", "sum"],
        default="highest_confidence",
    )
    parser.add_argument(
        "--gene_aggregation",
        choices=["signed_topk_mean", "max_abs", "sum", "mean"],
        default="signed_topk_mean",
    )
    parser.add_argument("--gene_topk_sites", type=int, default=3)
    parser.add_argument("--emit_gene_topk_site_comparison", type=_parse_bool, default=False)
    parser.add_argument("--gene_topk_site_compare_to", type=int, default=1)
    parser.add_argument("--ambiguous_gene_policy", choices=["drop", "split_equal", "first"], default="drop")
    parser.add_argument("--neglog10p_cap", type=float, default=50.0)
    parser.add_argument("--neglog10p_eps", type=float, default=1e-300)
    parser.add_argument("--resources_manifest")
    parser.add_argument("--resources_dir")
    parser.add_argument("--resource_policy", choices=["skip", "fail"], default="skip")
    parser.add_argument("--use_reference_bundle", type=_parse_bool, default=True)
    parser.add_argument("--site_alias_resource_id")
    parser.add_argument("--site_ubiquity_resource_id")
    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    parser.add_argument("--emit_full", type=_parse_bool, default=True)
    _add_gmt_flags(parser)
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


def _add_rna_de_prepare_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--modality", choices=["bulk", "scrna"], required=True)
    parser.add_argument("--counts_tsv", required=True, help="Dense count matrix TSV.")
    parser.add_argument("--matrix_orientation", choices=["sample_by_gene", "gene_by_sample"], default="sample_by_gene")
    parser.add_argument("--feature_id_column", required=True, help="Feature/gene identifier column for gene_by_sample inputs.")
    parser.add_argument("--matrix_gene_symbol_column", help="Optional gene symbol column for gene_by_sample matrices.")
    parser.add_argument("--matrix_delim", default="\t")
    parser.add_argument("--metadata_delim", default="\t")
    parser.add_argument("--sample_id_column", help="Sample ID column for bulk count matrices and sample metadata.")
    parser.add_argument("--sample_metadata_tsv", help="Bulk sample metadata TSV.")
    parser.add_argument("--subject_metadata_tsv", help="Optional subject/donor metadata TSV for bulk mode.")
    parser.add_argument("--subject_join_sample_column", help="Column in sample metadata used to join subject metadata.")
    parser.add_argument("--subject_join_metadata_column", help="Column in subject metadata used to join onto samples.")
    parser.add_argument("--subject_column", help="Subject/random-effect column after metadata joins in bulk mode.")
    parser.add_argument("--cell_id_column", help="Cell ID column for scRNA count matrices and cell metadata.")
    parser.add_argument("--cell_metadata_tsv", help="scRNA cell metadata TSV.")
    parser.add_argument("--donor_column", help="Donor column for scRNA pseudobulk and repeated-measures designs.")
    parser.add_argument("--cell_type_column", help="Optional cell-type column for scRNA pseudobulk grouping.")
    parser.add_argument(
        "--pseudobulk_within_cell_type",
        type=_parse_bool,
        default=True,
        help="If true, scRNA pseudobulks split by donor and cell_type when cell_type_column is provided.",
    )
    parser.add_argument("--min_cells_per_pseudobulk", type=int, default=10)
    parser.add_argument("--min_donors_per_group", type=int, default=2)
    parser.add_argument("--group_column", help="Primary condition/group column used to define contrasts.")
    parser.add_argument(
        "--comparison_mode",
        choices=["condition_a_vs_b", "group_vs_rest", "reference_level"],
        help="Simple CLI-driven contrast generation mode.",
    )
    parser.add_argument("--condition_a", help="Group A / focal level for CLI-generated comparisons.")
    parser.add_argument("--condition_b", help="Group B / reference level for condition_a_vs_b.")
    parser.add_argument("--reference_level", help="Reference level for reference_level comparisons.")
    parser.add_argument("--comparisons_tsv", help="Explicit comparisons TSV for advanced contrast definitions.")
    parser.add_argument("--stratify_by", help="Optional comma-separated stratification columns.")
    parser.add_argument("--covariates", help="Optional comma-separated fixed-effect covariates.")
    parser.add_argument("--batch_columns", help="Optional comma-separated batch columns.")
    parser.add_argument(
        "--de_mode",
        choices=["modern", "harmonizome"],
        default="modern",
        help=(
            "DE workflow preset. modern keeps the general-purpose use-all-samples behavior. "
            "harmonizome enables deterministic per-contrast group balancing plus explicit fixed-effect covariates "
            "for conservative GTEx/Harmonizome-style bulk runs."
        ),
    )
    parser.add_argument(
        "--balance_groups",
        type=_parse_bool,
        default=False,
        action=_StoreWithExplicitFlag,
        help="If true, downsample each emitted contrast to equal group sizes using min(n_group_a, n_group_b).",
    )
    parser.add_argument(
        "--balance_seed",
        type=int,
        default=0,
        action=_StoreWithExplicitFlag,
        help="Base seed for deterministic per-comparison balancing. harmonizome mode defaults this to 1 unless explicitly overridden.",
    )
    parser.add_argument("--repeated_measures", type=_parse_bool, default=False)
    parser.add_argument("--allow_approximate_repeated_measures", type=_parse_bool, default=False)
    parser.add_argument(
        "--backend",
        choices=["auto", "lightweight", "r_limma_voom", "r_dream"],
        default="auto",
        help="DE backend. auto prefers implemented/available R backends, else lightweight.",
    )
    parser.add_argument(
        "--allow_non_count_input",
        type=_parse_bool,
        default=False,
        help="Allow approximate analysis on non-integer normalized input and record that the run is approximate.",
    )
    parser.add_argument("--write_pseudobulk_artifacts", type=_parse_bool, default=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--organism", choices=["human", "mouse"], default="human")
    parser.add_argument("--genome_build", default="hg38")
    parser.add_argument("--run_extractor", type=_parse_bool, default=False)
    parser.add_argument("--extractor_out_dir")
    parser.add_argument("--extractor_signature_name", default="contrast")
    parser.add_argument(
        "--extractor_score_mode",
        choices=["auto", "stat", "logfc_times_neglog10p", "signed_neglog10padj", "custom_column"],
        default="auto",
    )
    parser.add_argument("--extractor_select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--extractor_top_k", type=int, default=200)
    parser.add_argument("--extractor_quantile", type=float, default=0.01)
    parser.add_argument("--extractor_min_score", type=float, default=0.0)
    parser.add_argument("--extractor_normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    parser.add_argument("--extractor_padj_max", type=float)
    parser.add_argument("--extractor_pvalue_max", type=float)
    parser.add_argument("--extractor_min_abs_logfc", type=float)
    parser.add_argument("--extractor_emit_gmt", type=_parse_bool, default=True)
    parser.add_argument("--extractor_gmt_split_signed", type=_parse_bool, default=True)
    parser.add_argument("--extractor_gmt_topk_list", default="200")
    parser.add_argument("--extractor_gmt_min_genes", type=int, default=100)
    parser.add_argument("--extractor_gmt_max_genes", type=int, default=500)


def _add_prism_prepare_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--release", default="24Q2")
    parser.add_argument("--matrix_url")
    parser.add_argument("--treatment_info_url")
    parser.add_argument("--cell_line_info_url")
    parser.add_argument("--matrix_file_id", default="20237709")
    parser.add_argument("--treatment_info_file_id", default="20237715")
    parser.add_argument("--cell_line_info_file_id", default="20237718")
    parser.add_argument("--depmap_file_id_url_template", default="https://ndownloader.figshare.com/files/{file_id}")
    parser.add_argument("--max_retries", type=int, default=5)
    parser.add_argument("--retry_backoff_sec", type=float, default=1.0)
    parser.add_argument("--user_agent", default="geneset-extractors/1.0 prism_prepare")
    parser.add_argument("--subset_seed", type=int, default=0)
    parser.add_argument("--max_cell_lines_total", type=int)
    parser.add_argument("--max_compounds_total", type=int)
    parser.add_argument("--max_cell_lines_per_group", type=int)
    parser.add_argument("--balance_by")
    parser.add_argument("--min_per_balance_bin", type=int, default=5)
    parser.add_argument("--keep_raw_downloads", type=_parse_bool, default=False)


def _add_ptm_prepare_reference_bundle_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--sources_tsv", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--organism", choices=["human", "mouse"], required=True)
    parser.add_argument("--ptm_type", choices=["phospho", "acetyl", "ub", "generic"], default="phospho")
    parser.add_argument("--bundle_id", required=True)


def _add_ptm_prepare_public_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input_mode", choices=["cdap_files", "pdc_manifest"], default="cdap_files")
    parser.add_argument("--ptm_report_tsv")
    parser.add_argument("--protein_report_tsv")
    parser.add_argument("--sample_design_tsv")
    parser.add_argument("--sample_annotations_tsv")
    parser.add_argument("--pdc_manifest_tsv")
    parser.add_argument("--source_dir")
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--organism", choices=["human", "mouse"], required=True)
    parser.add_argument("--ptm_type", choices=["phospho", "acetyl", "ub", "generic"], default="phospho")
    parser.add_argument("--study_id")
    parser.add_argument("--study_label")
    parser.add_argument("--assay_type_policy", choices=["off", "warn", "fail"], default="warn")
    parser.add_argument("--min_phospho_like_fraction", type=float, default=0.6)
    parser.add_argument("--max_k_fraction", type=float, default=0.25)


def _add_splice_event_diff_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--signature_name")
    parser.add_argument("--dataset_label")
    parser.add_argument(
        "--tool_family",
        choices=["auto", "generic", "leafcutter", "majiq", "whippet", "tcga_spliceseq"],
        default="auto",
    )
    parser.add_argument("--cluster_stats_tsv")
    parser.add_argument("--event_id_column")
    parser.add_argument("--event_group_column")
    parser.add_argument("--event_type_column")
    parser.add_argument("--gene_id_column")
    parser.add_argument("--gene_symbol_column")
    parser.add_argument("--chrom_column")
    parser.add_argument("--start_column")
    parser.add_argument("--end_column")
    parser.add_argument("--strand_column")
    parser.add_argument("--score_column")
    parser.add_argument("--stat_column")
    parser.add_argument("--delta_psi_column")
    parser.add_argument("--psi_column")
    parser.add_argument("--padj_column")
    parser.add_argument("--pvalue_column")
    parser.add_argument("--probability_column")
    parser.add_argument("--read_support_column")
    parser.add_argument("--novel_flag_column")
    parser.add_argument("--annotation_status_column")
    parser.add_argument(
        "--score_mode",
        choices=["auto", "stat", "delta_psi_times_neglog10p", "delta_psi_times_confidence", "custom_column"],
        default="auto",
    )
    parser.add_argument("--score_transform", choices=["signed", "abs", "positive", "negative"], default="signed")
    parser.add_argument(
        "--confidence_weight_mode",
        choices=["none", "pvalue", "probability", "read_support", "combined"],
        default="combined",
    )
    parser.add_argument("--min_probability", type=float, default=0.8)
    parser.add_argument("--min_read_support", type=float, default=10.0)
    parser.add_argument("--delta_psi_soft_floor", type=float, default=0.05)
    parser.add_argument("--delta_psi_soft_floor_mode", choices=["auto", "off"], default="auto")
    parser.add_argument(
        "--event_dup_policy",
        choices=["highest_confidence", "max_abs", "mean", "sum"],
        default="highest_confidence",
    )
    parser.add_argument(
        "--gene_aggregation",
        choices=["signed_topk_mean", "max_abs", "sum", "mean"],
        default="signed_topk_mean",
    )
    parser.add_argument("--gene_topk_events", type=int, default=3)
    parser.add_argument("--gene_burden_penalty_mode", choices=["none", "current_input", "reference_bundle", "auto"], default="auto")
    parser.add_argument("--min_gene_burden_penalty", type=float, default=0.35)
    parser.add_argument("--gene_support_penalty_mode", choices=["none", "independent_groups", "auto"], default="auto")
    parser.add_argument(
        "--locus_density_penalty_mode",
        choices=["none", "window_diversity", "chromosome_diversity", "auto"],
        default="none",
        action=_StoreWithExplicitFlag,
    )
    parser.add_argument("--locus_density_window_bp", type=int, default=20000000)
    parser.add_argument("--locus_density_top_n", type=int, default=20)
    parser.add_argument("--ambiguous_gene_policy", choices=["drop", "split_equal", "first"], default="drop")
    parser.add_argument("--impact_mode", choices=["none", "conservative", "custom_bundle"], default="conservative")
    parser.add_argument("--impact_min", type=float, default=0.75)
    parser.add_argument("--impact_max", type=float, default=1.35)
    parser.add_argument(
        "--bundle_prior_profile",
        choices=["auto", "full", "nuisance_only", "none"],
        default="auto",
        action=_StoreWithExplicitFlag,
    )
    parser.add_argument(
        "--artifact_action",
        choices=["warn", "prune", "suppress_if_persistent", "fail"],
        default="warn",
        action=_StoreWithExplicitFlag,
    )
    parser.add_argument("--neglog10p_cap", type=float, default=50.0)
    parser.add_argument("--neglog10p_eps", type=float, default=1e-300)
    parser.add_argument("--resources_manifest")
    parser.add_argument("--resources_dir")
    parser.add_argument("--resource_policy", choices=["skip", "fail"], default="skip")
    parser.add_argument("--use_reference_bundle", type=_parse_bool, default=True)
    parser.add_argument("--source_dataset")
    parser.add_argument("--bundle_same_dataset_policy", choices=["exclude", "warn", "fail", "ignore"], default="exclude")
    parser.add_argument("--event_alias_resource_id")
    parser.add_argument("--event_ubiquity_resource_id")
    parser.add_argument("--event_ubiquity_by_dataset_resource_id")
    parser.add_argument("--event_impact_resource_id")
    parser.add_argument("--gene_burden_resource_id")
    parser.add_argument("--gene_burden_by_dataset_resource_id")
    parser.add_argument("--gene_locus_resource_id")
    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    parser.add_argument("--emit_full", type=_parse_bool, default=True)
    _add_gmt_flags(parser)
    parser.set_defaults(
        emit_gmt=True,
        gmt_split_signed=True,
        gmt_topk_list="200",
        gmt_min_genes=100,
        gmt_max_genes=500,
    )


def _add_splice_event_matrix_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--signature_name")
    parser.add_argument("--dataset_label")
    parser.add_argument("--sample_metadata_tsv", required=True)
    parser.add_argument("--event_metadata_tsv")
    parser.add_argument("--coverage_matrix_tsv")
    parser.add_argument("--matrix_format", choices=["wide_events_by_sample"], default="wide_events_by_sample")
    parser.add_argument("--sample_id_column", default="sample_id")
    parser.add_argument("--group_column", default="group")
    parser.add_argument("--condition_column", default="condition")
    parser.add_argument(
        "--study_contrast",
        choices=["condition_a_vs_b", "group_vs_rest", "condition_within_group", "baseline"],
        default="condition_a_vs_b",
    )
    parser.add_argument("--condition_a")
    parser.add_argument("--condition_b")
    parser.add_argument("--min_samples_per_condition", type=int, default=3)
    parser.add_argument("--effect_metric", choices=["welch_t", "mean_diff"], default="welch_t")
    parser.add_argument("--missing_value_policy", choices=["drop", "min_present"], default="min_present")
    parser.add_argument("--min_present_per_condition", type=int, default=3)
    parser.add_argument(
        "--tool_family",
        choices=["auto", "generic", "leafcutter", "majiq", "whippet", "tcga_spliceseq"],
        default="generic",
    )
    parser.add_argument("--event_id_column")
    parser.add_argument("--event_group_column")
    parser.add_argument("--event_type_column")
    parser.add_argument("--gene_id_column")
    parser.add_argument("--gene_symbol_column")
    parser.add_argument("--chrom_column")
    parser.add_argument("--start_column")
    parser.add_argument("--end_column")
    parser.add_argument("--strand_column")
    parser.add_argument("--score_column")
    parser.add_argument("--stat_column", default="stat")
    parser.add_argument("--delta_psi_column", default="delta_psi")
    parser.add_argument("--psi_column")
    parser.add_argument("--padj_column", default="padj")
    parser.add_argument("--pvalue_column", default="pvalue")
    parser.add_argument("--probability_column")
    parser.add_argument("--read_support_column", default="read_support")
    parser.add_argument("--novel_flag_column")
    parser.add_argument("--annotation_status_column", default="annotation_status")
    parser.add_argument(
        "--score_mode",
        choices=["auto", "stat", "delta_psi_times_neglog10p", "delta_psi_times_confidence", "custom_column"],
        default="auto",
    )
    parser.add_argument("--score_transform", choices=["signed", "abs", "positive", "negative"], default="signed")
    parser.add_argument(
        "--confidence_weight_mode",
        choices=["none", "pvalue", "probability", "read_support", "combined"],
        default="combined",
    )
    parser.add_argument("--min_probability", type=float, default=0.8)
    parser.add_argument("--min_read_support", type=float, default=10.0)
    parser.add_argument("--delta_psi_soft_floor", type=float, default=0.05)
    parser.add_argument("--delta_psi_soft_floor_mode", choices=["auto", "off"], default="auto")
    parser.add_argument("--impact_mode", choices=["none", "conservative", "custom_bundle"], default="conservative")
    parser.add_argument("--impact_min", type=float, default=0.75)
    parser.add_argument("--impact_max", type=float, default=1.35)
    parser.add_argument(
        "--event_dup_policy",
        choices=["highest_confidence", "max_abs", "mean", "sum"],
        default="highest_confidence",
    )
    parser.add_argument(
        "--gene_aggregation",
        choices=["signed_topk_mean", "max_abs", "sum", "mean"],
        default="signed_topk_mean",
    )
    parser.add_argument("--gene_topk_events", type=int, default=3)
    parser.add_argument("--gene_burden_penalty_mode", choices=["none", "current_input", "reference_bundle", "auto"], default="auto")
    parser.add_argument("--min_gene_burden_penalty", type=float, default=0.35)
    parser.add_argument("--gene_support_penalty_mode", choices=["none", "independent_groups", "auto"], default="auto")
    parser.add_argument(
        "--locus_density_penalty_mode",
        choices=["none", "window_diversity", "chromosome_diversity", "auto"],
        default="none",
        action=_StoreWithExplicitFlag,
    )
    parser.add_argument("--locus_density_window_bp", type=int, default=20000000)
    parser.add_argument("--locus_density_top_n", type=int, default=20)
    parser.add_argument("--ambiguous_gene_policy", choices=["drop", "split_equal", "first"], default="drop")
    parser.add_argument("--neglog10p_cap", type=float, default=50.0)
    parser.add_argument("--neglog10p_eps", type=float, default=1e-300)
    parser.add_argument(
        "--bundle_prior_profile",
        choices=["auto", "full", "nuisance_only", "none"],
        default="auto",
        action=_StoreWithExplicitFlag,
    )
    parser.add_argument(
        "--artifact_action",
        choices=["warn", "prune", "suppress_if_persistent", "fail"],
        default="warn",
        action=_StoreWithExplicitFlag,
    )
    parser.add_argument("--resources_manifest")
    parser.add_argument("--resources_dir")
    parser.add_argument("--resource_policy", choices=["skip", "fail"], default="skip")
    parser.add_argument("--use_reference_bundle", type=_parse_bool, default=True)
    parser.add_argument("--source_dataset")
    parser.add_argument("--bundle_same_dataset_policy", choices=["exclude", "warn", "fail", "ignore"], default="exclude")
    parser.add_argument("--event_alias_resource_id")
    parser.add_argument("--event_ubiquity_resource_id")
    parser.add_argument("--event_ubiquity_by_dataset_resource_id")
    parser.add_argument("--event_impact_resource_id")
    parser.add_argument("--gene_burden_resource_id")
    parser.add_argument("--gene_burden_by_dataset_resource_id")
    parser.add_argument("--gene_locus_resource_id")
    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    parser.add_argument("--emit_full", type=_parse_bool, default=True)
    _add_gmt_flags(parser)
    parser.set_defaults(
        emit_gmt=True,
        gmt_split_signed=True,
        gmt_topk_list="200",
        gmt_min_genes=100,
        gmt_max_genes=500,
    )


def _add_splice_prepare_public_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input_mode", choices=["tcga_spliceseq"], default="tcga_spliceseq")
    parser.add_argument("--psi_tsv")
    parser.add_argument("--sample_annotations_tsv")
    parser.add_argument("--event_metadata_tsv")
    parser.add_argument("--sample_id_map_tsv")
    parser.add_argument("--event_type_allowlist")
    parser.add_argument("--missing_psi_policy", choices=["retain", "drop"], default="retain")
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--organism", choices=["human", "mouse"], required=True)
    parser.add_argument("--genome_build", required=True)
    parser.add_argument("--study_id")
    parser.add_argument("--study_label")


def _add_splice_prepare_reference_bundle_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--sources_tsv", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--organism", choices=["human", "mouse"], required=True)
    parser.add_argument("--bundle_id", required=True)
    parser.add_argument("--min_ref_read_support", type=float, default=0.0)
    parser.add_argument("--exclude_source_datasets")


def _add_calr_common_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--calr_data_csv", required=True)
    parser.add_argument("--session_csv")
    parser.add_argument("--exclusions_tsv")
    parser.add_argument("--dataset_label")
    parser.add_argument("--signature_name")
    parser.add_argument("--analysis_start_hour", type=float)
    parser.add_argument("--analysis_end_hour", type=float)
    parser.add_argument("--photoperiod_lights_on_hour", type=float)
    parser.add_argument("--photoperiod_hours_light", type=float, default=12.0)
    parser.add_argument("--exploratory_without_session", type=_parse_bool, default=True)
    parser.add_argument("--mass_covariate")
    parser.add_argument("--min_group_size", type=int, default=2)
    parser.add_argument("--resources_manifest")
    parser.add_argument("--resources_dir")
    parser.add_argument("--resource_policy", choices=["skip", "fail"], default="skip")
    parser.add_argument("--reference_bundle_id")
    parser.add_argument("--output_gene_species", choices=["human", "source"], default="human")
    parser.add_argument("--ortholog_policy", choices=["unique_only", "expand_all"], default="unique_only")
    parser.add_argument("--mouse_human_orthologs_tsv")
    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    parser.add_argument("--emit_full", type=_parse_bool, default=True)
    _add_gmt_flags(parser)
    parser.set_defaults(
        emit_gmt=True,
        gmt_topk_list="200",
        gmt_min_genes=100,
        gmt_max_genes=500,
        gmt_split_signed=False,
        gmt_require_symbol=False,
        emit_small_gene_sets=False,
    )


def _add_calr_ontology_mapper_flags(parser: argparse.ArgumentParser) -> None:
    _add_calr_common_flags(parser)
    parser.add_argument("--term_templates_tsv")
    parser.add_argument("--phenotype_gene_edges_tsv")
    parser.add_argument("--term_hierarchy_tsv")


def _add_calr_profile_query_flags(parser: argparse.ArgumentParser) -> None:
    _add_calr_common_flags(parser)
    parser.add_argument("--reference_profiles_tsv")
    parser.add_argument("--reference_metadata_tsv")
    parser.add_argument("--feature_schema_tsv")
    parser.add_argument("--feature_stats_tsv")
    parser.add_argument("--similarity_metric", choices=["cosine", "pearson"], default="cosine")
    parser.add_argument("--similarity_floor", type=float, default=0.0)
    parser.add_argument("--similarity_power", type=float, default=1.0)
    parser.add_argument("--hubness_penalty", choices=["none", "inverse_linear", "inverse_rank"], default="inverse_linear")
    parser.add_argument("--provenance_mismatch_penalty", type=float, default=0.15)


def _add_calr_prepare_reference_bundle_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--reference_profiles_tsv", required=True)
    parser.add_argument("--reference_metadata_tsv", required=True)
    parser.add_argument("--feature_schema_tsv")
    parser.add_argument("--feature_stats_tsv")
    parser.add_argument("--term_templates_tsv")
    parser.add_argument("--phenotype_gene_edges_tsv")
    parser.add_argument("--term_hierarchy_tsv")
    parser.add_argument("--include_packaged_term_hierarchy", type=_parse_bool, default=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--organism", choices=["mouse", "human"], required=True)
    parser.add_argument("--bundle_id", required=True)
    parser.add_argument("--mouse_human_orthologs_tsv")
    parser.add_argument("--write_distribution_artifact", type=_parse_bool, default=True)
    parser.add_argument("--distribution_dir")


def _add_calr_prepare_public_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--studies_tsv", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--organism", choices=["mouse", "human"], required=True)
    parser.add_argument("--output_gene_species", choices=["human", "source"], default="human")
    parser.add_argument("--ortholog_policy", choices=["unique_only", "expand_all"], default="unique_only")
    parser.add_argument("--mouse_human_orthologs_tsv")
    parser.add_argument("--bundle_id")
    parser.add_argument("--build_bundle", type=_parse_bool, default=True)
    parser.add_argument("--write_distribution_artifact", type=_parse_bool, default=True)
    parser.add_argument("--term_templates_tsv")
    parser.add_argument("--phenotype_gene_edges_tsv")
    parser.add_argument("--term_hierarchy_tsv")
    parser.add_argument("--include_packaged_term_hierarchy", type=_parse_bool, default=True)
    parser.add_argument("--exploratory_without_session", type=_parse_bool, default=True)
    parser.add_argument("--min_group_size", type=int, default=2)
    parser.add_argument("--mass_covariate")
    parser.add_argument("--analysis_start_hour", type=float)
    parser.add_argument("--analysis_end_hour", type=float)
    parser.add_argument("--photoperiod_lights_on_hour", type=float)
    parser.add_argument("--photoperiod_hours_light", type=float, default=12.0)


def _add_jump_prepare_reference_bundle_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--profile_paths", required=True, help="Comma-separated morphology profile TSV/CSV paths.")
    parser.add_argument("--experimental_metadata_tsv", required=True)
    parser.add_argument("--compound_targets_tsv", required=True)
    parser.add_argument("--target_annotations_tsv")
    parser.add_argument("--use_default_target_annotations", type=_parse_bool, default=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--bundle_id")
    parser.add_argument("--profile_id_column", default="sample_id")
    parser.add_argument("--metadata_join_id_column", default="sample_id")
    parser.add_argument("--perturbation_id_column", default="perturbation_id")
    parser.add_argument("--perturbation_type_column", default="perturbation_type")
    parser.add_argument("--cell_type_column", default="cell_type_or_line")
    parser.add_argument("--timepoint_column", default="timepoint")
    parser.add_argument("--gene_symbol_column", default="gene_symbol")
    parser.add_argument("--compound_id_column", default="compound_id")
    parser.add_argument("--is_control_column", default="is_control")
    parser.add_argument("--control_type_column", default="control_type")
    parser.add_argument("--cell_type_filter")
    parser.add_argument("--timepoint_filter")
    parser.add_argument("--require_same_timepoint_across_modalities", type=_parse_bool, default=True)
    parser.add_argument("--allow_missing_modalities", type=_parse_bool, default=True)
    parser.add_argument("--allow_mixed_timepoints", type=_parse_bool, default=False)
    parser.add_argument("--hubness_k", type=int, default=50)
    parser.add_argument("--profile_kind", default="normalized_feature_select_negcon_plate")
    parser.add_argument("--consensus_aggregate", choices=["median", "mean"], default="median")
    parser.add_argument("--profile_delimiter", default="\t")
    parser.add_argument("--metadata_delimiter", default="\t")
    parser.add_argument("--compound_targets_delimiter", default="\t")
    parser.add_argument("--compound_target_gene_symbol_column", default="gene_symbol")
    parser.add_argument("--compound_target_weight_column", default="weight")
    parser.add_argument("--target_annotations_delimiter", default="\t")
    parser.add_argument("--target_annotation_gene_symbol_column", default="gene_symbol")
    parser.add_argument("--write_distribution_artifact", type=_parse_bool, default=True)
    parser.add_argument("--distribution_dir")


def _add_morphology_profile_query_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--query_profiles_tsv", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--organism", choices=["human", "mouse"], required=True)
    parser.add_argument("--genome_build", required=True)
    parser.add_argument("--dataset_label")
    parser.add_argument("--signature_name")

    parser.add_argument("--query_id_column", default="sample_id")
    parser.add_argument("--query_profiles_delimiter", default="\t")
    parser.add_argument("--query_metadata_tsv")
    parser.add_argument("--query_metadata_id_column", default="sample_id")
    parser.add_argument("--query_metadata_delimiter", default="\t")
    parser.add_argument("--query_modality_column", default="perturbation_type")
    parser.add_argument("--group_query_by")
    parser.add_argument("--query_aggregate", choices=["none", "median", "mean"], default="median")
    parser.add_argument("--exclude_query_ids_from_reference", type=_parse_bool, default=False)

    parser.add_argument("--reference_profiles_tsv")
    parser.add_argument("--reference_profiles_parquet")
    parser.add_argument("--reference_profiles_delimiter", default="\t")
    parser.add_argument("--reference_id_column", default="perturbation_id")
    parser.add_argument("--reference_metadata_tsv")
    parser.add_argument("--reference_metadata_id_column", default="perturbation_id")
    parser.add_argument("--reference_metadata_delimiter", default="\t")
    parser.add_argument("--compound_targets_tsv")
    parser.add_argument("--target_annotations_tsv")
    parser.add_argument("--compound_targets_delimiter", default="\t")
    parser.add_argument("--compound_id_column", default="compound_id")
    parser.add_argument("--compound_target_gene_symbol_column", default="gene_symbol")
    parser.add_argument("--compound_target_weight_column", default="weight")
    parser.add_argument("--feature_stats_tsv")
    parser.add_argument("--feature_schema_tsv")

    parser.add_argument("--resources_manifest")
    parser.add_argument("--resources_dir")
    parser.add_argument("--resource_policy", choices=["skip", "fail"], default="skip")
    parser.add_argument("--reference_bundle_id")

    parser.add_argument("--mode", choices=["direct_target", "mechanism", "hybrid"], default="direct_target")
    parser.add_argument("--similarity_metric", choices=["cosine", "pearson"], default="cosine")
    parser.add_argument("--similarity_power", type=float, default=1.0)
    parser.add_argument("--polarity", choices=["similar", "opposite", "both"], default="similar")
    parser.add_argument("--same_modality_first", type=_parse_bool, default=True)
    parser.add_argument("--cross_modality_penalty", type=float, default=0.35)
    parser.add_argument("--max_reference_neighbors", type=int, default=20)
    parser.add_argument("--adaptive_neighbors", type=_parse_bool, default=True)
    parser.add_argument("--min_effective_neighbors", type=int, default=5)
    parser.add_argument("--neighbor_evidence_drop_ratio", type=float, default=0.25)
    parser.add_argument("--mutual_neighbor_filter", type=_parse_bool, default=True)
    parser.add_argument("--min_similarity", type=float, default=0.0)
    parser.add_argument("--control_calibration", choices=["none", "mean_center", "residualize_controls"], default="mean_center")
    parser.add_argument("--control_residual_components", type=int, default=2)
    parser.add_argument("--control_min_profiles_for_residualization", type=int, default=5)
    parser.add_argument("--hubness_penalty", choices=["none", "inverse_linear", "inverse_rank"], default="inverse_rank")
    parser.add_argument("--gene_recurrence_penalty", choices=["none", "idf"], default="idf")
    parser.add_argument("--min_specificity_confidence_to_emit_opposite", choices=["low", "medium", "high"], default="medium")
    parser.add_argument("--compound_weight", type=float, default=0.5)
    parser.add_argument("--genetic_weight", type=float, default=0.5)

    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    parser.add_argument("--emit_full", type=_parse_bool, default=True)
    _add_gmt_flags(parser)
    parser.set_defaults(
        emit_gmt=True,
        gmt_topk_list="200",
        gmt_min_genes=100,
        gmt_max_genes=500,
        gmt_split_signed=False,
        gmt_require_symbol=False,
        emit_small_gene_sets=False,
    )


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


def _add_drug_response_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--program_preset",
        choices=["default", "connectable", "broad_pharmacology"],
        default="connectable",
        help="Drug-response preset. connectable (default) uses bundle priors/QC when available; broad_pharmacology is lighter-touch.",
    )
    parser.add_argument("--dataset_label")

    parser.add_argument("--response_tsv")
    parser.add_argument("--response_delimiter", default="\t")
    parser.add_argument("--sample_id_column", default="sample_id")
    parser.add_argument("--drug_id_column", default="drug_id")
    parser.add_argument("--response_column", default="response")

    parser.add_argument("--drug_targets_tsv")
    parser.add_argument("--targets_delimiter", default="\t")
    parser.add_argument("--targets_drug_id_column", default="drug_id")
    parser.add_argument("--targets_gene_symbol_column", default="gene_symbol")
    parser.add_argument("--targets_weight_column", default="weight")
    parser.add_argument("--targets_source_column", default="source")

    parser.add_argument("--sample_metadata_tsv")
    parser.add_argument("--sample_metadata_delimiter", default="\t")
    parser.add_argument("--sample_metadata_sample_id_column", default="sample_id")
    parser.add_argument("--group_column", default="group")

    parser.add_argument("--groups_tsv")
    parser.add_argument("--groups_delimiter", default="\t")
    parser.add_argument("--groups_sample_id_column", default="sample_id")
    parser.add_argument("--groups_group_column", default="group")

    parser.add_argument("--case_control_tsv")
    parser.add_argument("--case_control_delimiter", default="\t")
    parser.add_argument("--case_control_sample_id_column", default="sample_id")
    parser.add_argument("--case_control_is_case_column", default="is_case")
    parser.add_argument("--case_control_within_group", type=_parse_bool, default=False)

    parser.add_argument("--prism_matrix_csv")
    parser.add_argument("--prism_treatment_info_csv")
    parser.add_argument("--prism_cell_line_info_csv")
    parser.add_argument("--prism_matrix_sample_id_column", default="row_name")
    parser.add_argument("--prism_cell_line_sample_id_column", default="row_name")
    parser.add_argument("--prism_column_name_column", default="column_name")
    parser.add_argument("--prism_broad_id_column", default="broad_id")
    parser.add_argument("--prism_target_column", default="target")

    parser.add_argument("--response_metric", choices=["logfold_change", "auc", "ic50", "other"], default="logfold_change")
    parser.add_argument(
        "--response_direction",
        choices=["lower_is_more_sensitive", "higher_is_more_sensitive"],
        help="If omitted, defaults by response_metric (logfold_change/auc/ic50 -> lower_is_more_sensitive).",
    )
    parser.add_argument("--response_transform", choices=["robust_z_mad", "rank_normal_score", "none"], default="robust_z_mad")
    parser.add_argument("--contrast_method", choices=["none", "group_mean", "group_vs_rest", "case_control"])
    parser.add_argument("--min_group_size", type=int, default=5)
    parser.add_argument("--scoring_model", choices=["target_weighted_sum", "sparse_deconvolution"], default="target_weighted_sum")
    parser.add_argument("--sparse_alpha", type=float, default=0.01)
    parser.add_argument(
        "--response_ubiquity_penalty",
        choices=["none", "fraction_active"],
        help="Response-side ubiquity penalty (preferred flag; supersedes --ubiquity_penalty when set).",
    )
    parser.add_argument("--ubiquity_penalty", choices=["none", "fraction_active"], default="fraction_active")
    parser.add_argument("--target_ubiquity_penalty", choices=["none", "idf"], default="idf")
    parser.add_argument("--ubiquity_tau", type=float, default=1.0)
    parser.add_argument("--ubiquity_epsilon", type=float, default=0.05)
    parser.add_argument("--polypharm_downweight", type=_parse_bool, default=True)
    parser.add_argument("--polypharm_t0", type=int, default=5)
    parser.add_argument("--max_targets_per_drug", type=int, default=50)
    parser.add_argument("--target_promiscuity_policy", choices=["warn", "drop", "cap"], default="warn")
    parser.add_argument("--max_programs", type=int, default=50)

    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    parser.add_argument("--emit_full", type=_parse_bool, default=True)

    parser.add_argument("--target_aliases_tsv")
    parser.add_argument("--alias_delimiter", default="\t")
    parser.add_argument("--alias_column", default="alias")
    parser.add_argument("--alias_gene_symbol_column", default="gene_symbol")
    parser.add_argument("--drug_blacklist_tsv")
    parser.add_argument("--drug_blacklist_delimiter", default="\t")
    parser.add_argument("--use_compound_qc_bundle", type=_parse_bool)
    parser.add_argument("--include_blacklisted_compounds", type=_parse_bool, default=False)
    parser.add_argument("--resources_manifest")
    parser.add_argument("--resources_dir")
    parser.add_argument("--resource_policy", choices=["skip", "fail"], default="skip")
    parser.add_argument("--target_aliases_resource_id")
    parser.add_argument("--drug_blacklist_resource_id")
    parser.add_argument("--drug_alias_map_resource_id")
    parser.add_argument("--target_edges_resource_id")
    parser.add_argument("--target_ubiquity_resource_id")
    parser.add_argument("--compound_qc_resource_id")

    parser.add_argument("--gtf")
    parser.add_argument("--gtf_source")
    parser.add_argument("--gtf_gene_id_field", default="gene_id")

    _add_gmt_flags(parser)
    parser.set_defaults(
        emit_gmt=True,
        gmt_split_signed=None,
        gmt_format="classic",
        gmt_topk_list="200",
        gmt_min_genes=100,
        gmt_max_genes=500,
        emit_small_gene_sets=False,
        gmt_require_symbol=False,
    )


def build_parser() -> argparse.ArgumentParser:
    prog_name = Path(sys.argv[0]).name if sys.argv else "geneset-extractors"
    parser = argparse.ArgumentParser(prog=prog_name or "geneset-extractors")
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
    p_rna_de_prepare = wf_sub.add_parser("rna_de_prepare")
    _add_rna_de_prepare_flags(p_rna_de_prepare)
    p_prism_prepare = wf_sub.add_parser("prism_prepare")
    _add_prism_prepare_flags(p_prism_prepare)
    p_ptm_public = wf_sub.add_parser("ptm_prepare_public")
    _add_ptm_prepare_public_flags(p_ptm_public)
    p_ptm_prepare = wf_sub.add_parser("ptm_prepare_reference_bundle")
    _add_ptm_prepare_reference_bundle_flags(p_ptm_prepare)
    p_splice_public = wf_sub.add_parser("splice_prepare_public")
    _add_splice_prepare_public_flags(p_splice_public)
    p_splice_prepare = wf_sub.add_parser("splice_prepare_reference_bundle")
    _add_splice_prepare_reference_bundle_flags(p_splice_prepare)
    p_calr_public = wf_sub.add_parser("calr_prepare_public")
    _add_calr_prepare_public_flags(p_calr_public)
    p_calr_prepare = wf_sub.add_parser("calr_prepare_reference_bundle")
    _add_calr_prepare_reference_bundle_flags(p_calr_prepare)
    p_jump_prepare = wf_sub.add_parser("jump_prepare_reference_bundle")
    _add_jump_prepare_reference_bundle_flags(p_jump_prepare)

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
    _add_provenance_flags(p_bulk)

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
    _add_provenance_flags(p_bulk_matrix)

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
    _add_provenance_flags(p_sc)

    p_rna = conv.add_parser("rna_deg")
    p_rna.add_argument("--deg_tsv", required=True)
    p_rna.add_argument("--out_dir", required=True)
    p_rna.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_rna.add_argument("--genome_build", required=True)
    _add_rna_deg_flags(p_rna)
    _add_provenance_flags(p_rna)

    p_rna_multi = conv.add_parser("rna_deg_multi")
    p_rna_multi.add_argument("--deg_tsv", required=True)
    p_rna_multi.add_argument("--comparison_column", required=True)
    p_rna_multi.add_argument("--out_dir", required=True)
    p_rna_multi.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_rna_multi.add_argument("--genome_build", required=True)
    _add_rna_deg_flags(p_rna_multi)
    _add_provenance_flags(p_rna_multi)

    p_rna_sc_programs = conv.add_parser("rna_sc_programs")
    p_rna_sc_programs.add_argument("--out_dir", required=True)
    p_rna_sc_programs.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_rna_sc_programs.add_argument("--genome_build", required=True)
    _add_rna_sc_program_flags(p_rna_sc_programs)
    _add_provenance_flags(p_rna_sc_programs)

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
        "--segments_format",
        choices=["auto", "gdc_seg", "cbio_seg"],
        default="auto",
        help="Segment table format preset for column autodetection.",
    )
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
        help=(
            "Enable purity correction. Examples: --use_purity_correction ; "
            "--use_purity_correction true ; --use_purity_correction false."
        ),
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
    p_cnv.add_argument(
        "--max_segment_length_bp",
        type=int,
        default=0,
        help="Optional hard cap on segment length after coordinate normalization (0 disables).",
    )
    p_cnv.add_argument("--focal_length_scale_bp", type=float, default=10_000_000)
    p_cnv.add_argument("--focal_length_alpha", type=float, default=1.0)
    p_cnv.add_argument("--gene_count_penalty", choices=["none", "inv_sqrt", "inv"], default="inv_sqrt")
    p_cnv.add_argument("--aggregation", choices=["weighted_mean", "sum", "mean", "max_abs"], default="weighted_mean")
    p_cnv.add_argument(
        "--program_preset",
        choices=["default", "connectable", "focal", "broad", "qc", "all", "none"],
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
        help=(
            "Enable cohort recurrence gene sets. Examples: --emit_cohort_sets ; "
            "--emit_cohort_sets true ; --emit_cohort_sets false."
        ),
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
    _add_provenance_flags(p_cnv)
    p_cnv.set_defaults(
        emit_gmt=True,
        gmt_split_signed=False,
        gmt_topk_list="200",
        gmt_min_genes=100,
        gmt_max_genes=500,
        emit_small_gene_sets=False,
    )

    p_drug = conv.add_parser("drug_response_screen")
    p_drug.add_argument("--out_dir", required=True)
    p_drug.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_drug.add_argument("--genome_build", required=True)
    _add_drug_response_flags(p_drug)
    _add_provenance_flags(p_drug)

    p_calr_onto = conv.add_parser("calr_ontology_mapper")
    p_calr_onto.add_argument("--out_dir", required=True)
    p_calr_onto.add_argument("--organism", choices=["mouse", "human"], required=True)
    p_calr_onto.add_argument("--genome_build", required=True)
    _add_calr_ontology_mapper_flags(p_calr_onto)
    _add_provenance_flags(p_calr_onto)

    p_calr_profile = conv.add_parser("calr_profile_query")
    p_calr_profile.add_argument("--out_dir", required=True)
    p_calr_profile.add_argument("--organism", choices=["mouse", "human"], required=True)
    p_calr_profile.add_argument("--genome_build", required=True)
    _add_calr_profile_query_flags(p_calr_profile)
    _add_provenance_flags(p_calr_profile)

    p_morph = conv.add_parser("morphology_profile_query")
    _add_morphology_profile_query_flags(p_morph)
    _add_provenance_flags(p_morph)

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
    _add_provenance_flags(p_chip)

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
    _add_provenance_flags(p_meth_cpg)

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
    _add_provenance_flags(p_meth_dmr_regions)

    p_meth = conv.add_parser("methylation_dmr")
    p_meth.add_argument("--dmr_tsv", required=True)
    p_meth.add_argument("--out_dir", required=True)
    p_meth.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_meth.add_argument("--genome_build", required=True)
    p_meth.add_argument("--gene_id_column", default="gene_id")
    p_meth.add_argument("--score_column", default="delta_methylation")
    _add_transform_flags(p_meth, default="abs")
    _add_provenance_flags(p_meth)

    p_prot = conv.add_parser("proteomics_diff")
    p_prot.add_argument("--proteomics_tsv", required=True)
    p_prot.add_argument("--out_dir", required=True)
    p_prot.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_prot.add_argument("--genome_build", required=True)
    p_prot.add_argument("--gene_id_column", default="gene_id")
    p_prot.add_argument("--score_column", default="log2fc")
    _add_transform_flags(p_prot, default="abs")
    _add_provenance_flags(p_prot)

    p_ptm = conv.add_parser("ptm_site_diff")
    p_ptm.add_argument("--ptm_tsv", required=True)
    p_ptm.add_argument("--out_dir", required=True)
    p_ptm.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_ptm.add_argument("--genome_build", required=True)
    _add_ptm_site_diff_flags(p_ptm)
    _add_provenance_flags(p_ptm)

    p_ptm_matrix = conv.add_parser("ptm_site_matrix")
    p_ptm_matrix.add_argument("--ptm_matrix_tsv", required=True)
    p_ptm_matrix.add_argument("--out_dir", required=True)
    p_ptm_matrix.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_ptm_matrix.add_argument("--genome_build", required=True)
    _add_ptm_site_matrix_flags(p_ptm_matrix)
    _add_provenance_flags(p_ptm_matrix)

    p_splice = conv.add_parser("splice_event_diff")
    p_splice.add_argument("--splice_tsv", required=True)
    p_splice.add_argument("--out_dir", required=True)
    p_splice.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_splice.add_argument("--genome_build", required=True)
    _add_splice_event_diff_flags(p_splice)
    _add_provenance_flags(p_splice)

    p_splice_matrix = conv.add_parser("splice_event_matrix")
    p_splice_matrix.add_argument("--psi_matrix_tsv", required=True)
    p_splice_matrix.add_argument("--out_dir", required=True)
    p_splice_matrix.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_splice_matrix.add_argument("--genome_build", required=True)
    _add_splice_event_matrix_flags(p_splice_matrix)
    _add_provenance_flags(p_splice_matrix)

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
    _add_provenance_flags(p_scrna)

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
                from geneset_extractors.workflows.scrna_cnmf_prepare import run as run_scrna_cnmf_prepare

                result = run_scrna_cnmf_prepare(args)
                print(
                    "workflow_completed "
                    f"workflow=scrna_cnmf_prepare n_subsets={result.get('n_subsets')} "
                    f"out={result.get('out_dir')}",
                    file=sys.stderr,
                )
                return 0
            if args.workflow_command == "cnmf_select_k":
                from geneset_extractors.workflows.cnmf_select_k import run as run_cnmf_select_k

                _ = run_cnmf_select_k(args)
                return 0
            if args.workflow_command == "rna_de_prepare":
                from geneset_extractors.workflows.rna_de_prepare import run as run_rna_de_prepare

                result = run_rna_de_prepare(args)
                print(
                    "workflow_completed "
                    f"workflow=rna_de_prepare n_comparisons={result.get('n_comparisons')} "
                    f"out={result.get('out_dir')}",
                    file=sys.stderr,
                )
                return 0
            if args.workflow_command == "prism_prepare":
                from geneset_extractors.workflows.prism_prepare import run as run_prism_prepare

                result = run_prism_prepare(args)
                print(
                    "workflow_completed "
                    f"workflow=prism_prepare n_rows={result.get('n_response_rows')} "
                    f"out={result.get('out_dir')}",
                    file=sys.stderr,
                )
                return 0
            if args.workflow_command == "ptm_prepare_public":
                from geneset_extractors.workflows.ptm_prepare_public import run as run_ptm_prepare_public

                result = run_ptm_prepare_public(args)
                print(
                    "workflow_completed "
                    f"workflow=ptm_prepare_public n_sites={result.get('ptm_rows_emitted')} "
                    f"n_samples={result.get('n_samples')} "
                    f"out={result.get('outputs', {}).get('ptm_matrix_tsv')}",
                    file=sys.stderr,
                )
                return 0
            if args.workflow_command == "ptm_prepare_reference_bundle":
                from geneset_extractors.workflows.ptm_prepare_reference_bundle import run as run_ptm_prepare_reference_bundle

                result = run_ptm_prepare_reference_bundle(args)
                print(
                    "workflow_completed "
                    f"workflow=ptm_prepare_reference_bundle n_canonical_sites={result.get('n_canonical_sites')} "
                    f"out={result.get('out_dir')}",
                    file=sys.stderr,
                )
                return 0
            if args.workflow_command == "splice_prepare_public":
                from geneset_extractors.workflows.splice_prepare_public import run as run_splice_prepare_public

                result = run_splice_prepare_public(args)
                print(
                    "workflow_completed "
                    f"workflow=splice_prepare_public n_events={result.get('n_events')} "
                    f"n_samples={result.get('n_samples')} "
                    f"out={result.get('outputs', {}).get('psi_matrix_tsv')}",
                    file=sys.stderr,
                )
                return 0
            if args.workflow_command == "splice_prepare_reference_bundle":
                from geneset_extractors.workflows.splice_prepare_reference_bundle import run as run_splice_prepare_reference_bundle

                result = run_splice_prepare_reference_bundle(args)
                print(
                    "workflow_completed "
                    f"workflow=splice_prepare_reference_bundle n_canonical_events={result.get('n_canonical_events')} "
                    f"out={result.get('out_dir')}",
                    file=sys.stderr,
                )
                return 0
            if args.workflow_command == "calr_prepare_reference_bundle":
                from geneset_extractors.workflows.calr_prepare_reference_bundle import run as run_calr_prepare_reference_bundle

                result = run_calr_prepare_reference_bundle(args)
                print(
                    "workflow_completed "
                    f"workflow=calr_prepare_reference_bundle n_profiles={result.get('n_profiles')} "
                    f"out={result.get('bundle_manifest')} "
                    f"tarball={result.get('tarball')}",
                    file=sys.stderr,
                )
                return 0
            if args.workflow_command == "calr_prepare_public":
                from geneset_extractors.workflows.calr_prepare_public import run as run_calr_prepare_public

                result = run_calr_prepare_public(args)
                print(
                    "workflow_completed "
                    f"workflow=calr_prepare_public n_profiles={result.get('n_profiles')} "
                    f"out={result.get('out_dir')} "
                    f"bundle={result.get('bundle_manifest')}",
                    file=sys.stderr,
                )
                return 0
            if args.workflow_command == "jump_prepare_reference_bundle":
                from geneset_extractors.workflows.jump_prepare_reference_bundle import run as run_jump_prepare_reference_bundle

                result = run_jump_prepare_reference_bundle(args)
                print(
                    "workflow_completed "
                    f"workflow=jump_prepare_reference_bundle n_profiles={result.get('n_profiles')} "
                    f"out={result.get('bundle_manifest')} "
                    f"tarball={result.get('tarball')}",
                    file=sys.stderr,
                )
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
