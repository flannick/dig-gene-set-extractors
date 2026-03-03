from __future__ import annotations

from pathlib import Path

from omics2geneset.core.metadata import input_file_record, make_metadata, write_metadata
from omics2geneset.core.models import GeneWeights
from omics2geneset.core.normalization import normalize
from omics2geneset.core.peak_to_gene import link_distance_decay, link_nearest_tss, link_promoter_overlap
from omics2geneset.core.scoring import score_genes
from omics2geneset.io.bed import read_bed
from omics2geneset.io.gtf import read_genes_from_gtf


def _resolve_max_distance(args) -> int:
    if args.max_distance_bp is not None:
        return int(args.max_distance_bp)
    if args.link_method == "nearest_tss":
        return 100000
    if args.link_method == "distance_decay":
        return 500000
    return 100000


def _link(peaks: list[dict[str, object]], genes, args):
    max_distance_bp = _resolve_max_distance(args)
    if args.link_method == "promoter_overlap":
        return link_promoter_overlap(peaks, genes, args.promoter_upstream_bp, args.promoter_downstream_bp)
    if args.link_method == "nearest_tss":
        return link_nearest_tss(peaks, genes, max_distance_bp)
    if args.link_method == "distance_decay":
        return link_distance_decay(peaks, genes, max_distance_bp, args.decay_length_bp, args.max_genes_per_peak)
    raise ValueError(f"Unsupported link_method: {args.link_method}")


def _extract_peak_weights(peaks: list[dict[str, object]], weight_column: int) -> list[float]:
    idx = weight_column - 1
    out: list[float] = []
    for p in peaks:
        cols = p.get("columns", [])
        if idx >= len(cols):
            raise ValueError("weight_column out of range for peaks file")
        out.append(float(cols[idx]))
    return out


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    peaks = read_bed(args.peaks)
    genes = read_genes_from_gtf(args.gtf)
    weights = _extract_peak_weights(peaks, args.weight_column)
    links = _link(peaks, genes, args)

    raw_gene_weights = score_genes(weights, links, args.peak_weight_transform)
    final = normalize(raw_gene_weights, args.normalize)

    gene_symbol_by_id = {g.gene_id: g.gene_symbol for g in genes}
    rows = [{"gene_id": gid, "weight": w, "gene_symbol": gene_symbol_by_id.get(gid)} for gid, w in final.items()]
    GeneWeights(rows).sort_desc().to_tsv(out_dir / "geneset.tsv")

    assigned_peaks = len({int(link["peak_index"]) for link in links})
    params = {
        "link_method": args.link_method,
        "promoter_upstream_bp": args.promoter_upstream_bp,
        "promoter_downstream_bp": args.promoter_downstream_bp,
        "max_distance_bp": _resolve_max_distance(args),
        "decay_length_bp": args.decay_length_bp,
        "max_genes_per_peak": args.max_genes_per_peak,
        "weight_column": args.weight_column,
        "peak_weight_transform": args.peak_weight_transform,
        "normalize": args.normalize,
        "aggregation": "sum",
    }

    meta = make_metadata(
        converter_name="chipseq_peak",
        parameters=params,
        data_type="chip_seq",
        assay="bulk",
        organism=args.organism,
        genome_build=args.genome_build,
        files=[input_file_record(args.peaks, "peaks"), input_file_record(args.gtf, "gtf")],
        gene_annotation={"mode": "gtf", "gtf_path": str(args.gtf), "source": args.gtf_source or "user", "gene_id_field": "gene_id"},
        weights={
            "weight_type": "signed" if args.peak_weight_transform == "signed" else "nonnegative",
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
    write_metadata(out_dir / "geneset.meta.json", meta)
    return {"n_peaks": len(peaks), "n_genes": len(rows), "out_dir": str(out_dir)}
