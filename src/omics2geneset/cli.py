from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from omics2geneset.core.validate import validate_output_dir
from omics2geneset.registry import get_converter, get_converter_spec, list_converters


def _parse_bool(value: str) -> bool:
    v = value.strip().lower()
    if v in {"true", "1", "yes", "y"}:
        return True
    if v in {"false", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError("expected true or false")


def _add_linking_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--link_method", choices=["promoter_overlap", "nearest_tss", "distance_decay"], default="promoter_overlap")
    parser.add_argument("--promoter_upstream_bp", type=int, default=2000)
    parser.add_argument("--promoter_downstream_bp", type=int, default=500)
    parser.add_argument("--max_distance_bp", type=int)
    parser.add_argument("--decay_length_bp", type=int, default=50000)
    parser.add_argument("--max_genes_per_peak", type=int, default=5)


def _add_program_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--select", choices=["none", "top_k", "quantile", "threshold"], default="top_k")
    parser.add_argument("--top_k", type=int, default=200)
    parser.add_argument("--quantile", type=float, default=0.01)
    parser.add_argument("--min_score", type=float, default=0.0)
    parser.add_argument("--emit_full", type=_parse_bool, default=True)


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
