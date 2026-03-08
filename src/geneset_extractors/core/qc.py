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
        calibration_method = str(params.get("calibration_method", params.get("contrast_method", "")))
        link_method = str(params.get("link_method", ""))
        direction = str(params.get("direction", ""))
        key = (program_method, calibration_method, link_method, direction)
        if key in seen:
            continue
        seen.add(key)
        record: dict[str, str] = {
            "program_method": program_method,
            "calibration_method": calibration_method,
            "link_method": link_method,
        }
        if direction:
            record["direction"] = direction
        out.append(record)
    return out


def render_run_summary_text(payload: dict[str, object]) -> str:
    lines: list[str] = []
    lines.append("Run Summary")
    lines.append("===========")
    for key in (
        "converter",
        "dataset_label",
        "group",
        "study_contrast",
        "program_preset",
        "primary_program_method",
        "primary_link_method",
        "primary_calibration_method",
        "selected_direction",
    ):
        if key in payload:
            lines.append(f"{key}: {payload[key]}")
    if "primary_calibration_method" not in payload and "primary_contrast_method" in payload:
        lines.append(f"primary_calibration_method: {payload['primary_contrast_method']}")

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
                + f" calibration_method={row.get('calibration_method','')}"
                + f" link_method={row.get('link_method','')}"
                + (f" direction={row.get('direction','')}" if row.get("direction") else "")
            )
    for key in (
        "retrieval_confidence",
        "specificity_confidence",
        "top10_gene_mass",
        "neighbor_target_concentration",
        "neighbor_primary_target_agreement",
        "neighbor_top3_target_agreement",
        "high_hub_mass_fraction",
        "n_positive_neighbors",
        "n_negative_neighbors",
        "n_retained_neighbors",
        "neg_similarity_fraction",
    ):
        if key in payload:
            lines.append(f"{key}: {payload[key]}")
    skipped_programs = payload.get("skipped_programs")
    if isinstance(skipped_programs, list) and skipped_programs:
        lines.append(f"skipped_programs: {len(skipped_programs)}")
        for row in skipped_programs:
            if not isinstance(row, dict):
                continue
            parts = []
            for key in ("query_id", "polarity", "reason", "specificity_confidence", "min_required"):
                if key in row:
                    parts.append(f"{key}={row[key]}")
            lines.append("  " + " ".join(parts))

    for key in ("program_methods_skipped", "calibration_methods_skipped", "contrast_methods_skipped"):
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

    parse_summary = payload.get("parse_summary")
    if isinstance(parse_summary, dict):
        lines.append("parse_summary:")
        for key in sorted(parse_summary):
            lines.append(f"  {key}: {parse_summary[key]}")

    sign_semantics = payload.get("sign_semantics")
    if isinstance(sign_semantics, dict):
        lines.append("sign_semantics:")
        for key in ("delta_orientation", "pos_set_meaning", "neg_set_meaning"):
            if key in sign_semantics:
                lines.append(f"  {key}: {sign_semantics[key]}")

    artifact = payload.get("artifact_diagnostics")
    if isinstance(artifact, dict):
        patterns = artifact.get("patterns")
        threshold = artifact.get("pattern_threshold")
        warnings = artifact.get("warnings")
        lines.append("artifact_diagnostics:")
        if isinstance(patterns, list):
            lines.append(f"  patterns: {patterns}")
        if threshold is not None:
            lines.append(f"  pattern_threshold: {threshold}")
        if isinstance(warnings, list):
            lines.append(f"  n_warnings: {len(warnings)}")
            for warning in warnings:
                if not isinstance(warning, dict):
                    continue
                lines.append(
                    "  "
                    + f"set={warning.get('set_name','')} "
                    + f"pattern={warning.get('pattern','')} "
                    + f"fraction={warning.get('fraction','')}"
                )

    symbol_filter = payload.get("symbol_filter")
    if isinstance(symbol_filter, dict):
        lines.append("symbol_filter:")
        for key in ("exclude_gene_symbol_regex", "exclude_gene_symbols_tsv"):
            if key in symbol_filter:
                lines.append(f"  {key}: {symbol_filter[key]}")
        counts_by_output = symbol_filter.get("counts_by_output")
        if isinstance(counts_by_output, dict):
            lines.append("  counts_by_output:")
            for output_key in sorted(counts_by_output):
                lines.append(f"    {output_key}: {counts_by_output[output_key]}")

    resources = payload.get("resources")
    if isinstance(resources, dict):
        lines.append("resources:")
        if "manifest" in resources:
            lines.append(f"  manifest: {resources.get('manifest')}")
        if "resources_dir" in resources:
            lines.append(f"  resources_dir: {resources.get('resources_dir')}")
        used = resources.get("used")
        if isinstance(used, list):
            lines.append(f"  used: {len(used)}")
            for entry in used:
                if not isinstance(entry, dict):
                    continue
                lines.append(
                    "  "
                    + f"resource_id={entry.get('id','')} "
                    + f"method={entry.get('method','')} "
                    + f"path={entry.get('path','')}"
                )
        missing = resources.get("missing")
        if isinstance(missing, list):
            lines.append(f"  missing: {len(missing)}")
            for entry in missing:
                if not isinstance(entry, dict):
                    continue
                lines.append(
                    "  "
                    + f"resource_id={entry.get('resource_id','')} "
                    + f"role={entry.get('role_label','')} "
                    + f"expected_path={entry.get('expected_path','')}"
                )

    for key in ("n_samples", "n_drugs", "n_response_rows", "n_groups_total", "n_groups_emitted", "n_groups_skipped"):
        if key in payload:
            lines.append(f"{key}: {payload[key]}")

    groups = payload.get("groups")
    if isinstance(groups, list):
        lines.append("groups:")
        for row in groups:
            if not isinstance(row, dict):
                continue
            lines.append(
                "  "
                + f"group={row.get('group','')} "
                + f"n_group={row.get('n_group',0)} "
                + f"n_rest={row.get('n_rest',0)} "
                + f"sets_emitted={row.get('sets_emitted',0)} "
                + f"reason_if_skipped={row.get('reason_if_skipped','')}"
            )

    target_table_stats = payload.get("target_table_stats")
    if isinstance(target_table_stats, dict):
        lines.append("target_table_stats:")
        for key in sorted(target_table_stats):
            lines.append(f"  {key}: {target_table_stats[key]}")

    gene_set_emission = payload.get("gene_set_emission")
    if isinstance(gene_set_emission, dict):
        lines.append("gene_set_emission:")
        for key in sorted(gene_set_emission):
            lines.append(f"  {key}: {gene_set_emission[key]}")

    warnings_payload = payload.get("warnings")
    if isinstance(warnings_payload, list):
        lines.append(f"warnings: {len(warnings_payload)}")
        for warning in warnings_payload:
            if isinstance(warning, dict):
                lines.append(f"  {warning.get('code','warning')}: {warning.get('message','')}")
            else:
                lines.append(f"  warning: {warning}")

    return "\n".join(lines) + "\n"


def write_run_summary_files(out_dir: Path, payload: dict[str, object]) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "run_summary.json"
    txt_path = out_dir / "run_summary.txt"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    txt_path.write_text(render_run_summary_text(payload), encoding="utf-8")
    return json_path, txt_path
