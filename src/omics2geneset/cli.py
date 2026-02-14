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
        choices=["none", "default", "all"],
        default="default",
        help="Program family preset for GMT emission and metadata.",
    )
    parser.add_argument(
        "--program_methods",
        help=(
            "Optional comma-separated program methods overriding preset. "
            "Bulk: linked_activity,promoter_activity,distal_activity,enhancer_bias,ref_ubiquity_penalty,atlas_residual. "
            "scATAC adds tfidf_distal."
        ),
    )
    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument("--emit_full", type=_parse_bool, default=True)


def _add_resource_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--resources_manifest", help="Optional resources manifest path (defaults to bundled manifest)")
    parser.add_argument("--resources_dir", help="Optional resources cache directory (defaults to OMICS2GENESET_RESOURCES_DIR or ~/.cache)")
    parser.add_argument("--resource_policy", choices=["skip", "fail"], default="skip")
    parser.add_argument("--ref_ubiquity_resource_id", help="Resource id for ref_ubiquity_penalty program")
    parser.add_argument("--atlas_resource_id", help="Resource id for atlas_residual program")
    parser.add_argument("--atlas_metric", choices=["logratio", "zscore"], default="logratio")
    parser.add_argument("--atlas_eps", type=float, default=1e-6)


def _add_gmt_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--emit_gmt", type=_parse_bool, default=True)
    parser.add_argument("--gmt_out")
    parser.add_argument("--gmt_prefer_symbol", type=_parse_bool, default=True)
    parser.add_argument("--gmt_require_symbol", type=_parse_bool, default=True)
    parser.add_argument("--gmt_biotype_allowlist", default="protein_coding")
    parser.add_argument("--gmt_min_genes", type=int, default=100)
    parser.add_argument("--gmt_max_genes", type=int, default=500)
    parser.add_argument("--gmt_topk_list", default="100,200,500")
    parser.add_argument("--gmt_mass_list", default="0.5,0.8,0.9")
    parser.add_argument("--gmt_split_signed", type=_parse_bool, default=False)


def _add_transform_flags(parser: argparse.ArgumentParser, default: str) -> None:
    parser.add_argument("--score_transform", choices=["signed", "abs", "positive", "negative"], default=default)
    parser.add_argument("--normalize", choices=["l1", "none"], default="l1")


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
    p_bulk.add_argument("--peaks_weight_column", type=int, default=5)
    p_bulk.add_argument("--peak_weights_tsv")
    p_bulk.add_argument("--peak_weight_transform", choices=["signed", "abs", "positive", "negative"], default="abs")
    p_bulk.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    p_bulk.add_argument("--gtf_source")
    _add_linking_flags(p_bulk)
    _add_program_flags(p_bulk)
    _add_resource_flags(p_bulk)
    _add_gmt_flags(p_bulk)

    p_sc = conv.add_parser("atac_sc_10x")
    p_sc.add_argument("--matrix_dir", required=True)
    p_sc.add_argument("--gtf", required=True)
    p_sc.add_argument("--out_dir", required=True)
    p_sc.add_argument("--organism", choices=["human", "mouse"], required=True)
    p_sc.add_argument("--genome_build", required=True)
    p_sc.add_argument("--groups_tsv")
    p_sc.add_argument("--peak_summary", choices=["sum_counts", "mean_counts", "frac_cells_nonzero"], default="sum_counts")
    p_sc.add_argument("--peak_weight_transform", choices=["signed", "abs", "positive", "negative"], default="positive")
    p_sc.add_argument("--normalize", choices=["none", "l1", "within_set_l1"], default="within_set_l1")
    p_sc.add_argument("--contrast", choices=["none", "group_vs_rest"])
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
    p_rna.add_argument("--gene_id_column", default="gene_id")
    p_rna.add_argument("--score_column", default="log2fc")
    _add_transform_flags(p_rna, default="signed")

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


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
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
                rows = resource_status_rows(resources, args.resources_dir, verify=verify)
                for row in rows:
                    print(
                        "\t".join(
                            [
                                str(row["id"]),
                                str(row["status"]),
                                str(row["path"]),
                                str(row.get("availability", "")),
                                str(row.get("size_bytes", "")),
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
            return 0
    except Exception as exc:  # pragma: no cover
        print(f"error: {exc}", file=sys.stderr)
        return 1

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
