from __future__ import annotations

from bisect import bisect_left, bisect_right
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
    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for g in genes:
        if g.strand == "+":
            prom_start = max(0, g.tss - promoter_upstream_bp)
            prom_end = g.tss + promoter_downstream_bp
        else:
            prom_start = max(0, g.tss - promoter_downstream_bp)
            prom_end = g.tss + promoter_upstream_bp
        by_chrom.setdefault(g.chrom, []).append((prom_start, prom_end, g.gene_id))

    peaks_by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    for pi, p in enumerate(peaks):
        peaks_by_chrom.setdefault(str(p["chrom"]), []).append((int(p["start"]), int(p["end"]), pi))

    for chrom, peak_list in peaks_by_chrom.items():
        intervals = sorted(by_chrom.get(chrom, []), key=lambda x: x[0])
        if not intervals:
            continue
        peak_list_sorted = sorted(peak_list, key=lambda x: x[0])
        active: list[tuple[int, int, str]] = []
        interval_idx = 0
        for p_start, p_end, pi in peak_list_sorted:
            while interval_idx < len(intervals) and intervals[interval_idx][0] < p_end:
                active.append(intervals[interval_idx])
                interval_idx += 1
            active = [it for it in active if it[1] > p_start]
            for prom_start, prom_end, gene_id in active:
                if _overlap(p_start, p_end, prom_start, prom_end):
                    links.append({"peak_index": pi, "gene_id": gene_id, "distance": 0, "link_weight": 1.0})
    return links


def _index_genes_by_chrom(genes: list[Gene]) -> dict[str, tuple[list[int], list[Gene]]]:
    by_chrom: dict[str, list[Gene]] = {}
    for g in genes:
        by_chrom.setdefault(g.chrom, []).append(g)
    indexed: dict[str, tuple[list[int], list[Gene]]] = {}
    for chrom, g_list in by_chrom.items():
        g_sorted = sorted(g_list, key=lambda g: g.tss)
        indexed[chrom] = ([g.tss for g in g_sorted], g_sorted)
    return indexed


def link_nearest_tss(
    peaks: list[dict[str, object]], genes: list[Gene], max_distance_bp: int
) -> list[dict[str, object]]:
    links: list[dict[str, object]] = []
    indexed = _index_genes_by_chrom(genes)
    for pi, p in enumerate(peaks):
        chrom = str(p["chrom"])
        if chrom not in indexed:
            continue
        tss_positions, genes_sorted = indexed[chrom]
        peak_center = (int(p["start"]) + int(p["end"])) // 2
        pos = bisect_left(tss_positions, peak_center)
        candidate_indices = []
        if pos < len(genes_sorted):
            candidate_indices.append(pos)
        if pos > 0:
            candidate_indices.append(pos - 1)
        best_gene = None
        best_dist = None
        for idx in candidate_indices:
            d = abs(peak_center - genes_sorted[idx].tss)
            if best_dist is None or d < best_dist:
                best_dist = d
                best_gene = genes_sorted[idx]
        if best_gene is not None and best_dist is not None and int(best_dist) <= max_distance_bp:
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
    indexed = _index_genes_by_chrom(genes)

    for pi, p in enumerate(peaks):
        chrom = str(p["chrom"])
        if chrom not in indexed:
            continue
        tss_positions, genes_sorted = indexed[chrom]
        peak_center = (int(p["start"]) + int(p["end"])) // 2
        left = bisect_left(tss_positions, peak_center - max_distance_bp)
        right = bisect_right(tss_positions, peak_center + max_distance_bp)
        candidates: list[dict[str, object]] = []
        for i in range(left, right):
            g = genes_sorted[i]
            d = abs(peak_center - g.tss)
            if d <= max_distance_bp:
                w = math.exp(-float(d) / float(decay_length_bp)) if decay_length_bp > 0 else 0.0
                candidates.append({"peak_index": pi, "gene_id": g.gene_id, "distance": d, "link_weight": w})
        candidates.sort(key=lambda x: (-float(x["link_weight"]), str(x["gene_id"])))
        links.extend(candidates[:max_genes_per_peak])
    return links
