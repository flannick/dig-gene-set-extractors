from __future__ import annotations


def _overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start < b_end and b_start < a_end


def peak_ref_idf_by_overlap(
    peaks: list[dict[str, object]],
    ref_intervals: list[dict[str, object]],
    default_idf: float = 1.0,
) -> list[float]:
    out = [float(default_idf)] * len(peaks)

    by_chrom: dict[str, list[tuple[int, int, float]]] = {}
    for row in ref_intervals:
        by_chrom.setdefault(str(row["chrom"]), []).append(
            (int(row["start"]), int(row["end"]), float(row["idf_ref"]))
        )

    peaks_by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    for pi, p in enumerate(peaks):
        peaks_by_chrom.setdefault(str(p["chrom"]), []).append((int(p["start"]), int(p["end"]), pi))

    for chrom, peak_list in peaks_by_chrom.items():
        intervals = sorted(by_chrom.get(chrom, []), key=lambda x: x[0])
        if not intervals:
            continue
        peak_list_sorted = sorted(peak_list, key=lambda x: x[0])
        active: list[tuple[int, int, float]] = []
        interval_idx = 0
        for p_start, p_end, pi in peak_list_sorted:
            while interval_idx < len(intervals) and intervals[interval_idx][0] < p_end:
                active.append(intervals[interval_idx])
                interval_idx += 1
            active = [it for it in active if it[1] > p_start]
            best = float(default_idf)
            for r_start, r_end, idf in active:
                if _overlap(p_start, p_end, r_start, r_end):
                    if idf > best:
                        best = float(idf)
            out[pi] = best
    return out


def apply_peak_idf(
    peak_weights: list[float],
    peak_idf: list[float],
) -> list[float]:
    if len(peak_weights) != len(peak_idf):
        raise ValueError("peak_weights and peak_idf length mismatch")
    return [float(w) * float(i) for w, i in zip(peak_weights, peak_idf)]
