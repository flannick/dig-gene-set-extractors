from __future__ import annotations

import gzip
import json
from pathlib import Path
from statistics import median


_HEADER_TOKENS = {
    "gene",
    "genes",
    "gene_symbol",
    "symbol",
    "marker",
    "markers",
    "gene_id",
}


def summarize_numeric(values: list[float]) -> dict[str, float] | None:
    if not values:
        return None
    numeric = [float(v) for v in values]
    return {
        "min": float(min(numeric)),
        "median": float(median(numeric)),
        "max": float(max(numeric)),
    }


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def load_marker_genes(path: str | None) -> set[str]:
    if path is None:
        return set()
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Marker gene file not found: {p}")
    markers: set[str] = set()
    first_data = True
    with _open_text(p) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            token = line
            if "\t" in token:
                token = token.split("\t", 1)[0]
            elif "," in token:
                token = token.split(",", 1)[0]
            else:
                token = token.split()[0]
            value = token.strip()
            if not value:
                continue
            lowered = value.lower()
            if first_data and lowered in _HEADER_TOKENS:
                first_data = False
                continue
            first_data = False
            markers.add(value.upper())
    return markers


def marker_hit_summary(
    rows: list[dict[str, object]],
    marker_genes: set[str],
    cutoffs: tuple[int, ...] = (10, 20, 50),
) -> dict[str, object] | None:
    if not marker_genes:
        return None

    normalized_markers = {m.upper() for m in marker_genes if m}

    def _row_tokens(row: dict[str, object]) -> list[str]:
        out: list[str] = []
        for key in ("gene_symbol", "gene_id"):
            value = row.get(key)
            if value is None:
                continue
            token = str(value).strip()
            if token:
                out.append(token.upper())
        return out

    def _row_label(row: dict[str, object]) -> str:
        symbol = str(row.get("gene_symbol", "")).strip()
        if symbol:
            return symbol
        return str(row.get("gene_id", "")).strip()

    rows_ranked = list(rows)
    if rows_ranked and "rank" in rows_ranked[0]:
        rows_ranked = sorted(rows_ranked, key=lambda r: int(r.get("rank", 0) or 0))

    by_cutoff: dict[str, dict[str, object]] = {}
    for cutoff in cutoffs:
        top_rows = rows_ranked[:cutoff]
        matched_markers: set[str] = set()
        matched_rows: list[str] = []
        for row in top_rows:
            tokens = _row_tokens(row)
            hit = False
            for token in tokens:
                if token in normalized_markers:
                    matched_markers.add(token)
                    hit = True
            if hit:
                label = _row_label(row)
                if label:
                    matched_rows.append(label)
        by_cutoff[f"top_{cutoff}"] = {
            "n_rows_considered": len(top_rows),
            "n_marker_hits": len(matched_rows),
            "matched_rows": matched_rows,
            "matched_markers_unique": sorted(matched_markers),
        }

    selected_markers: set[str] = set()
    for row in rows_ranked:
        for token in _row_tokens(row):
            if token in normalized_markers:
                selected_markers.add(token)

    return {
        "n_markers_provided": len(normalized_markers),
        "matched_markers_in_selected_program": len(selected_markers),
        "matched_markers_in_selected_program_unique": sorted(selected_markers),
        "top_hits": by_cutoff,
    }


def collect_emitted_method_combinations(gmt_plans: list[dict[str, object]]) -> list[dict[str, str]]:
    seen: set[tuple[str, str, str, str]] = set()
    out: list[dict[str, str]] = []
    for plan in gmt_plans:
        params = plan.get("parameters", {})
        if not isinstance(params, dict):
            continue
        program_method = str(params.get("program_method", ""))
        contrast_method = str(params.get("contrast_method", ""))
        link_method = str(params.get("link_method", ""))
        direction = str(params.get("direction", ""))
        key = (program_method, contrast_method, link_method, direction)
        if key in seen:
            continue
        seen.add(key)
        record: dict[str, str] = {
            "program_method": program_method,
            "contrast_method": contrast_method,
            "link_method": link_method,
        }
        if direction:
            record["direction"] = direction
        out.append(record)
    return out


def render_run_summary_text(payload: dict[str, object]) -> str:
    lines: list[str] = []
    lines.append("ATAC Run Summary")
    lines.append("================")
    for key in (
        "converter",
        "dataset_label",
        "group",
        "program_preset",
        "primary_program_method",
        "primary_link_method",
        "primary_contrast_method",
        "selected_direction",
    ):
        if key in payload:
            lines.append(f"{key}: {payload[key]}")

    if "n_input_peaks" in payload:
        lines.append(f"n_input_peaks: {payload['n_input_peaks']}")
    if "link_assignment" in payload:
        lines.append("link_assignment:")
        assignment = payload.get("link_assignment", {})
        if isinstance(assignment, dict):
            for method in sorted(assignment):
                stats = assignment.get(method, {})
                if not isinstance(stats, dict):
                    continue
                lines.append(
                    "  "
                    + f"{method}: n_assigned={stats.get('n_features_assigned')} "
                    + f"fraction={stats.get('fraction_features_assigned')}"
                )
    if "promoter_peak_count" in payload or "distal_peak_count" in payload:
        lines.append(
            "peak_partition: "
            + f"promoter={payload.get('promoter_peak_count', 0)} "
            + f"distal={payload.get('distal_peak_count', 0)}"
        )

    ref = payload.get("ref_ubiquity")
    if isinstance(ref, dict):
        lines.append("ref_ubiquity:")
        for k in ("n_overlapping_peaks", "fraction_overlapping_peaks", "idf_stats", "adjusted_peak_weight_stats"):
            if k in ref:
                lines.append(f"  {k}: {ref[k]}")

    emitted = payload.get("emitted_method_combinations")
    if isinstance(emitted, list):
        lines.append(f"emitted_method_combinations: {len(emitted)}")
        for row in emitted:
            if not isinstance(row, dict):
                continue
            lines.append(
                "  "
                + f"program={row.get('program_method','')}"
                + f" contrast_method={row.get('contrast_method','')}"
                + f" link_method={row.get('link_method','')}"
                + (f" direction={row.get('direction','')}" if row.get("direction") else "")
            )

    for key in ("program_methods_skipped", "contrast_methods_skipped"):
        value = payload.get(key)
        if isinstance(value, dict):
            if value:
                lines.append(f"{key}:")
                for k in sorted(value):
                    lines.append(f"  {k}: {value[k]}")
            else:
                lines.append(f"{key}: none")

    marker = payload.get("marker_qc")
    if isinstance(marker, dict):
        lines.append("marker_qc:")
        lines.append(
            "  "
            + f"n_markers_provided={marker.get('n_markers_provided', 0)} "
            + f"matched_markers_in_selected_program={marker.get('matched_markers_in_selected_program', 0)}"
        )
        top_hits = marker.get("top_hits", {})
        if isinstance(top_hits, dict):
            for cutoff in sorted(top_hits):
                entry = top_hits.get(cutoff, {})
                if not isinstance(entry, dict):
                    continue
                lines.append(
                    "  "
                    + f"{cutoff}: n_marker_hits={entry.get('n_marker_hits', 0)} "
                    + f"matched_rows={entry.get('matched_rows', [])}"
                )

    return "\n".join(lines) + "\n"


def write_run_summary_files(out_dir: Path, payload: dict[str, object]) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "run_summary.json"
    txt_path = out_dir / "run_summary.txt"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    txt_path.write_text(render_run_summary_text(payload), encoding="utf-8")
    return json_path, txt_path
