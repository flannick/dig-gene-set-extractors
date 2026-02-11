from __future__ import annotations

import csv
from pathlib import Path

from omics2geneset.core.metadata import input_file_record, make_metadata, write_metadata
from omics2geneset.core.models import GeneWeights
from omics2geneset.core.normalization import normalize
from omics2geneset.core.peak_to_gene import link_distance_decay, link_nearest_tss, link_promoter_overlap
from omics2geneset.core.scoring import score_genes
from omics2geneset.io.gtf import read_genes_from_gtf
from omics2geneset.io.mtx_10x import make_group_indices, read_10x_matrix_dir, read_groups_tsv, summarize_peaks


def _resolve_max_distance(args) -> int:
    if args.max_distance_bp is not None:
        return int(args.max_distance_bp)
    if args.link_method == "nearest_tss":
        return 100000
    if args.link_method == "distance_decay":
        return 500000
    return 100000


def _resolved_parameters(args, group_name: str) -> dict[str, object]:
    params = vars(args).copy()
    params["max_distance_bp"] = _resolve_max_distance(args)
    params["group"] = group_name
    return params


def _link(peaks: list[dict[str, object]], genes, args):
    max_distance_bp = _resolve_max_distance(args)
    if args.link_method == "promoter_overlap":
        return link_promoter_overlap(peaks, genes, args.promoter_upstream_bp, args.promoter_downstream_bp)
    if args.link_method == "nearest_tss":
        return link_nearest_tss(peaks, genes, max_distance_bp)
    if args.link_method == "distance_decay":
        return link_distance_decay(peaks, genes, max_distance_bp, args.decay_length_bp, args.max_genes_per_peak)
    raise ValueError(f"Unsupported link_method: {args.link_method}")


def _safe_group_name(name: str) -> str:
    return name.replace("/", "_").replace(" ", "_")


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    peaks, barcodes, entries, matrix_files = read_10x_matrix_dir(args.matrix_dir)
    genes = read_genes_from_gtf(args.gtf)
    groups = read_groups_tsv(args.groups_tsv) if args.groups_tsv else None
    group_indices = make_group_indices(barcodes, groups)
    links = _link(peaks, genes, args)

    files = [
        input_file_record(matrix_files["matrix"], "matrix"),
        input_file_record(matrix_files["barcodes"], "barcodes"),
        input_file_record(matrix_files["peaks_or_features"], "peaks_or_features"),
        input_file_record(args.gtf, "gtf"),
    ]
    if args.groups_tsv:
        files.append(input_file_record(args.groups_tsv, "groups_tsv"))
    manifest_rows: list[tuple[str, str]] = []

    gene_symbol_by_id = {g.gene_id: g.gene_symbol for g in genes}

    for group_name, cell_indices in group_indices.items():
        peak_weights = summarize_peaks(entries, len(peaks), len(barcodes), cell_indices, args.peak_summary)
        raw_gene_weights = score_genes(peak_weights, links, "positive")
        final_gene_weights = normalize(raw_gene_weights, args.normalize)

        rows = [
            {"gene_id": gid, "weight": w, "gene_symbol": gene_symbol_by_id.get(gid)}
            for gid, w in final_gene_weights.items()
        ]
        gw = GeneWeights(rows).sort_desc()

        if args.groups_tsv:
            group_dir = out_dir / f"group={_safe_group_name(group_name)}"
        else:
            group_dir = out_dir
        group_dir.mkdir(parents=True, exist_ok=True)

        gw.to_tsv(group_dir / "geneset.tsv")
        assigned_peaks = len({int(link["peak_index"]) for link in links})

        meta = make_metadata(
            converter_name="atac_sc_10x",
            parameters=_resolved_parameters(args, group_name),
            data_type="atac_seq",
            assay="single_cell",
            organism=args.organism,
            genome_build=args.genome_build,
            files=files,
            gene_annotation={
                "gtf_path": str(args.gtf),
                "source": args.gtf_source or "user",
                "gene_id_field": "gene_id",
            },
            weights={
                "weight_type": "nonnegative",
                "normalization": {"method": args.normalize, "target_sum": 1.0 if args.normalize == "l1" else None},
                "aggregation": "sum",
            },
            summary={
                "n_input_features": len(peaks),
                "n_genes": len(rows),
                "n_features_assigned": assigned_peaks,
                "fraction_features_assigned": assigned_peaks / len(peaks) if peaks else 0.0,
            },
        )
        write_metadata(group_dir / "geneset.meta.json", meta)
        manifest_rows.append((group_name, str(group_dir)))

    if args.groups_tsv:
        with (out_dir / "manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["group", "path"])
            writer.writerows(manifest_rows)

    return {"n_peaks": len(peaks), "n_genes": len(genes), "out_dir": str(out_dir), "n_groups": len(group_indices)}
