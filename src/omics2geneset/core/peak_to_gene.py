from __future__ import annotations

import math

from omics2geneset.core.models import Gene


def _overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start < b_end and b_start < a_end


def link_promoter_overlap(
    peaks: list[dict[str, object]],
    genes: list[Gene],
    promoter_upstream_bp: int,
    promoter_downstream_bp: int,
) -> list[dict[str, object]]:
    links: list[dict[str, object]] = []
    by_chrom: dict[str, list[Gene]] = {}
    for g in genes:
        by_chrom.setdefault(g.chrom, []).append(g)
    for pi, p in enumerate(peaks):
        chrom = str(p["chrom"])
        p_start = int(p["start"])
        p_end = int(p["end"])
        for g in by_chrom.get(chrom, []):
            if g.strand == "+":
                prom_start = g.tss - promoter_upstream_bp
                prom_end = g.tss + promoter_downstream_bp
            else:
                prom_start = g.tss - promoter_downstream_bp
                prom_end = g.tss + promoter_upstream_bp
            if _overlap(p_start, p_end, max(0, prom_start), prom_end):
                links.append({"peak_index": pi, "gene_id": g.gene_id, "distance": 0, "link_weight": 1.0})
    return links


def link_nearest_tss(
    peaks: list[dict[str, object]], genes: list[Gene], max_distance_bp: int
) -> list[dict[str, object]]:
    links: list[dict[str, object]] = []
    by_chrom: dict[str, list[Gene]] = {}
    for g in genes:
        by_chrom.setdefault(g.chrom, []).append(g)
    for pi, p in enumerate(peaks):
        chrom = str(p["chrom"])
        peak_center = (int(p["start"]) + int(p["end"])) // 2
        best_gene = None
        best_dist = None
        for g in by_chrom.get(chrom, []):
            d = abs(peak_center - g.tss)
            if best_dist is None or d < best_dist:
                best_dist = d
                best_gene = g
        if best_gene is not None and best_dist is not None and best_dist <= max_distance_bp:
            links.append({"peak_index": pi, "gene_id": best_gene.gene_id, "distance": best_dist, "link_weight": 1.0})
    return links


def link_distance_decay(
    peaks: list[dict[str, object]],
    genes: list[Gene],
    max_distance_bp: int,
    decay_length_bp: int,
    max_genes_per_peak: int,
) -> list[dict[str, object]]:
    links: list[dict[str, object]] = []
    by_chrom: dict[str, list[Gene]] = {}
    for g in genes:
        by_chrom.setdefault(g.chrom, []).append(g)

    for pi, p in enumerate(peaks):
        chrom = str(p["chrom"])
        peak_center = (int(p["start"]) + int(p["end"])) // 2
        candidates: list[dict[str, object]] = []
        for g in by_chrom.get(chrom, []):
            d = abs(peak_center - g.tss)
            if d <= max_distance_bp:
                w = math.exp(-float(d) / float(decay_length_bp)) if decay_length_bp > 0 else 0.0
                candidates.append({"peak_index": pi, "gene_id": g.gene_id, "distance": d, "link_weight": w})
        candidates.sort(key=lambda x: float(x["link_weight"]), reverse=True)
        links.extend(candidates[:max_genes_per_peak])
    return links
