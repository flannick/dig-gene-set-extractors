from __future__ import annotations

import csv
from pathlib import Path

from omics2geneset.core.metadata import input_file_record, make_metadata, write_metadata
from omics2geneset.core.peak_to_gene import link_distance_decay, link_nearest_tss, link_promoter_overlap
from omics2geneset.core.scoring import score_genes
from omics2geneset.core.selection import (
    global_l1_weights,
    ranked_gene_ids,
    select_quantile,
    select_threshold,
    select_top_k,
    within_set_l1_weights,
)
from omics2geneset.io.bed import read_bed, read_peak_weights_tsv
from omics2geneset.io.gtf import read_genes_from_gtf


def _arg(args, name: str, default):
    return getattr(args, name, default)


def _extract_peak_weights(peaks: list[dict[str, object]], peaks_weight_column: int | None, peak_weights_tsv: str | None) -> list[float]:
    if peak_weights_tsv:
        m = read_peak_weights_tsv(peak_weights_tsv)
        out: list[float] = []
        for p in peaks:
            key = (str(p["chrom"]), int(p["start"]), int(p["end"]))
            if key not in m:
                raise ValueError(f"Peak missing in peak_weights_tsv: {key}")
            out.append(float(m[key]))
        return out
    if peaks_weight_column is None:
        peaks_weight_column = 5
    idx = peaks_weight_column - 1
    out = []
    for p in peaks:
        cols = p.get("columns", [])
        if idx >= len(cols):
            raise ValueError("peaks_weight_column out of range for peaks file")
        out.append(float(cols[idx]))
    return out


def _resolve_max_distance(args) -> int:
    max_distance_bp = _arg(args, "max_distance_bp", None)
    if max_distance_bp is not None:
        return int(max_distance_bp)
    link_method = _arg(args, "link_method", "promoter_overlap")
    if link_method == "nearest_tss":
        return 100000
    if link_method == "distance_decay":
        return 500000
    return 100000


def _resolved_parameters(args) -> dict[str, object]:
    return {
        "link_method": _arg(args, "link_method", "promoter_overlap"),
        "promoter_upstream_bp": _arg(args, "promoter_upstream_bp", 2000),
        "promoter_downstream_bp": _arg(args, "promoter_downstream_bp", 500),
        "max_distance_bp": _resolve_max_distance(args),
        "decay_length_bp": _arg(args, "decay_length_bp", 50000),
        "max_genes_per_peak": _arg(args, "max_genes_per_peak", 5),
        "peak_weight_transform": _arg(args, "peak_weight_transform", "abs"),
        "normalize": _arg(args, "normalize", "within_set_l1"),
        "selection_method": _arg(args, "select", "top_k"),
        "top_k": _arg(args, "top_k", 200),
        "quantile": _arg(args, "quantile", 0.01),
        "min_score": _arg(args, "min_score", 0.0),
        "emit_full": bool(_arg(args, "emit_full", True)),
        "aggregation": "sum",
    }


def _link(peaks: list[dict[str, object]], genes, args):
    max_distance_bp = _resolve_max_distance(args)
    link_method = _arg(args, "link_method", "promoter_overlap")
    if link_method == "promoter_overlap":
        return link_promoter_overlap(
            peaks,
            genes,
            _arg(args, "promoter_upstream_bp", 2000),
            _arg(args, "promoter_downstream_bp", 500),
        )
    if link_method == "nearest_tss":
        return link_nearest_tss(peaks, genes, max_distance_bp)
    if link_method == "distance_decay":
        return link_distance_decay(
            peaks,
            genes,
            max_distance_bp,
            _arg(args, "decay_length_bp", 50000),
            _arg(args, "max_genes_per_peak", 5),
        )
    raise ValueError(f"Unsupported link_method: {link_method}")


def _select_gene_ids(scores: dict[str, float], args) -> list[str]:
    select = _arg(args, "select", "top_k")
    if select == "none":
        return ranked_gene_ids(scores)
    if select == "top_k":
        return select_top_k(scores, int(_arg(args, "top_k", 200)))
    if select == "quantile":
        return select_quantile(scores, float(_arg(args, "quantile", 0.01)))
    if select == "threshold":
        return select_threshold(scores, float(_arg(args, "min_score", 0.0)))
    raise ValueError(f"Unsupported selection method: {select}")


def _selected_weights(full_scores: dict[str, float], selected_gene_ids: list[str], normalize: str) -> dict[str, float]:
    if normalize == "none":
        return {g: float(full_scores[g]) for g in selected_gene_ids}
    if normalize == "l1":
        global_weights = global_l1_weights(full_scores)
        return {g: float(global_weights.get(g, 0.0)) for g in selected_gene_ids}
    if normalize == "within_set_l1":
        return within_set_l1_weights(full_scores, selected_gene_ids)
    raise ValueError(f"Unsupported normalization method: {normalize}")


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["gene_id", "score", "rank"]
    if any("weight" in row for row in rows):
        fieldnames.append("weight")
    if any("gene_symbol" in row and row["gene_symbol"] is not None for row in rows):
        fieldnames.append("gene_symbol")
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    peaks = read_bed(args.peaks)
    genes = read_genes_from_gtf(args.gtf)
    peak_weights = _extract_peak_weights(peaks, args.peaks_weight_column, args.peak_weights_tsv)
    links = _link(peaks, genes, args)

    peak_weight_transform = _arg(args, "peak_weight_transform", "abs")
    normalize = _arg(args, "normalize", "within_set_l1")
    emit_full = bool(_arg(args, "emit_full", True))

    raw_scores = score_genes(peak_weights, links, peak_weight_transform)
    full_scores = {g: float(s) for g, s in raw_scores.items() if float(s) != 0.0}
    selected_gene_ids = _select_gene_ids(full_scores, args)
    selected_weights = _selected_weights(full_scores, selected_gene_ids, normalize)

    gene_symbol_by_id = {g.gene_id: g.gene_symbol for g in genes}

    selected_rows: list[dict[str, object]] = []
    for rank, gene_id in enumerate(selected_gene_ids, start=1):
        selected_rows.append(
            {
                "gene_id": gene_id,
                "score": float(full_scores[gene_id]),
                "weight": float(selected_weights.get(gene_id, 0.0)),
                "rank": rank,
                "gene_symbol": gene_symbol_by_id.get(gene_id),
            }
        )
    _write_rows(out_dir / "geneset.tsv", selected_rows)

    full_gene_ids = ranked_gene_ids(full_scores)
    full_rows: list[dict[str, object]] = []
    for rank, gene_id in enumerate(full_gene_ids, start=1):
        full_rows.append(
            {
                "gene_id": gene_id,
                "score": float(full_scores[gene_id]),
                "rank": rank,
                "gene_symbol": gene_symbol_by_id.get(gene_id),
            }
        )
    if emit_full:
        _write_rows(out_dir / "geneset.full.tsv", full_rows)

    assigned_peaks = len({int(link["peak_index"]) for link in links})
    files = [input_file_record(args.peaks, "peaks"), input_file_record(args.gtf, "gtf")]
    if args.peak_weights_tsv:
        files.append(input_file_record(args.peak_weights_tsv, "peak_weights_tsv"))

    output_files = [{"path": str(out_dir / "geneset.tsv"), "role": "selected_program"}]
    if emit_full:
        output_files.append({"path": str(out_dir / "geneset.full.tsv"), "role": "full_scores"})

    params = _resolved_parameters(args)
    meta = make_metadata(
        converter_name="atac_bulk",
        parameters=params,
        data_type="atac_seq",
        assay="bulk",
        organism=args.organism,
        genome_build=args.genome_build,
        files=files,
        gene_annotation={
            "mode": "gtf",
            "gtf_path": str(args.gtf),
            "source": args.gtf_source or "user",
            "gene_id_field": "gene_id",
        },
        weights={
            "weight_type": "signed" if peak_weight_transform == "signed" else "nonnegative",
            "normalization": {
                "method": normalize,
                "target_sum": 1.0 if normalize == "within_set_l1" else None,
            },
            "aggregation": "sum",
        },
        summary={
            "n_input_features": len(peaks),
            "n_genes": len(selected_rows),
            "n_features_assigned": assigned_peaks,
            "fraction_features_assigned": assigned_peaks / len(peaks) if peaks else 0.0,
        },
        program_extraction={
            "selection_method": _arg(args, "select", "top_k"),
            "selection_params": {
                "k": _arg(args, "top_k", 200),
                "quantile": _arg(args, "quantile", 0.01),
                "min_score": _arg(args, "min_score", 0.0),
            },
            "normalize": normalize,
            "n_selected_genes": len(selected_rows),
            "score_definition": "sum over peaks of transformed peak_weight times link_weight",
        },
        output_files=output_files,
    )
    write_metadata(out_dir / "geneset.meta.json", meta)

    return {
        "n_peaks": len(peaks),
        "n_genes": len(selected_rows),
        "out_dir": str(out_dir),
    }
