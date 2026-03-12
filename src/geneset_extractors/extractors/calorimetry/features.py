from __future__ import annotations

from dataclasses import dataclass
from math import atan2, cos, pi, sin, sqrt
from pathlib import Path
import sys

from .io import (
    AMBIENT_TEMP_CANDIDATES,
    DIET_CANDIDATES,
    FAT_MASS_CANDIDATES,
    LEAN_MASS_CANDIDATES,
    SEX_CANDIDATES,
    STRAIN_CANDIDATES,
    TOTAL_MASS_CANDIDATES,
    CalRRow,
    find_first_column,
    parse_float,
    parse_timestamp,
    resolve_group_column,
    resolve_run_id_column,
    resolve_subject_id_column,
    resolve_time_columns,
    resolve_variable_columns,
)


MASS_DEPENDENT_VARS = {"vo2", "vco2", "ee", "feed", "drink", "feed.acc", "drink.acc"}
MASS_INDEPENDENT_VARS = {"rer", "body.temp", "xytot", "xyamb", "wheel", "pedmeter"}

PROGRAM_VARIABLES = {
    "global": set(MASS_DEPENDENT_VARS | MASS_INDEPENDENT_VARS),
    "thermogenesis": {"vo2", "vco2", "ee", "body.temp"},
    "substrate_use": {"rer", "vo2", "vco2"},
    "intake_balance": {"feed", "feed.acc", "drink", "drink.acc"},
    "activity_circadian": {"xytot", "xyamb", "wheel", "pedmeter", "body.temp"},
}


@dataclass
class CalorimetryFeatureConfig:
    analysis_start_hour: float | None
    analysis_end_hour: float | None
    photoperiod_lights_on_hour: float | None
    photoperiod_hours_light: float
    exclusions_tsv: str | None
    exploratory_without_session: bool


@dataclass
class SubjectSummary:
    subject_id: str
    group: str
    run_id: str
    feature_values: dict[str, float]
    metadata: dict[str, object]


@dataclass
class FeatureExtractionResult:
    subjects: dict[str, SubjectSummary]
    feature_names: list[str]
    metadata_summary: dict[str, object]
    warnings: list[str]


def _mean(values: list[float]) -> float | None:
    if not values:
        return None
    return float(sum(values) / float(len(values)))


def _sd(values: list[float]) -> float | None:
    if len(values) < 2:
        return None
    ctr = _mean(values)
    assert ctr is not None
    return float((sum((v - ctr) ** 2 for v in values) / float(len(values) - 1)) ** 0.5)


def _auc(times: list[float], values: list[float]) -> float | None:
    if len(times) != len(values) or len(times) < 2:
        return None
    total = 0.0
    for left, right, v_left, v_right in zip(times[:-1], times[1:], values[:-1], values[1:]):
        width = max(0.0, float(right) - float(left))
        total += width * (float(v_left) + float(v_right)) / 2.0
    return float(total)


def _circadian(times: list[float], values: list[float]) -> tuple[float | None, float | None, float | None]:
    if len(times) < 3 or len(values) < 3 or len(times) != len(values):
        return None, None, None
    n = float(len(times))
    a = 0.0
    b = 0.0
    for hour_value, measurement in zip(times, values):
        angle = 2.0 * pi * (float(hour_value) % 24.0) / 24.0
        a += float(measurement) * cos(angle)
        b += float(measurement) * sin(angle)
    a *= 2.0 / n
    b *= 2.0 / n
    amp = sqrt(a * a + b * b)
    phase = atan2(b, a)
    return float(amp), float(cos(phase)), float(sin(phase))


def _selected_window_from_session(session_rows: dict[str, dict[str, str]] | None) -> tuple[float | None, float | None]:
    if not session_rows:
        return None, None
    start_candidates = ("analysis_start", "analysis_start_hour", "start_hour", "start")
    end_candidates = ("analysis_end", "analysis_end_hour", "end_hour", "end")
    start = None
    end = None
    for row in session_rows.values():
        if start is None:
            col = find_first_column(list(row.keys()), start_candidates)
            if col:
                start = parse_float(row.get(col))
        if end is None:
            col = find_first_column(list(row.keys()), end_candidates)
            if col:
                end = parse_float(row.get(col))
    return start, end


def _photoperiod_from_session(session_rows: dict[str, dict[str, str]] | None) -> tuple[float | None, str]:
    if not session_rows:
        return None, "default"
    candidates = ("lights_on_hour", "light_start_hour", "light", "lights_on")
    for row in session_rows.values():
        col = find_first_column(list(row.keys()), candidates)
        if col is None:
            continue
        value = parse_float(row.get(col))
        if value is not None:
            return float(value), "session"
    return None, "default"


def _subject_level_metadata(row_values: dict[str, str], session_row: dict[str, str] | None, *, group_col: str | None, run_col: str | None) -> dict[str, object]:
    out: dict[str, object] = {}
    source_keys = list(row_values.keys()) + (list(session_row.keys()) if session_row else [])
    out["group"] = ""
    if group_col:
        if session_row and group_col in session_row and str(session_row.get(group_col, "")).strip():
            out["group"] = str(session_row.get(group_col, "")).strip()
        elif group_col in row_values:
            out["group"] = str(row_values.get(group_col, "")).strip()
    out["run_id"] = str(row_values.get(run_col or "", "")).strip() if run_col else ""
    for label, candidates in (
        ("total_mass", TOTAL_MASS_CANDIDATES),
        ("lean_mass", LEAN_MASS_CANDIDATES),
        ("fat_mass", FAT_MASS_CANDIDATES),
        ("sex", SEX_CANDIDATES),
        ("diet", DIET_CANDIDATES),
        ("strain", STRAIN_CANDIDATES),
        ("ambient_temperature", AMBIENT_TEMP_CANDIDATES),
    ):
        col = find_first_column(source_keys, candidates)
        value = ""
        if session_row and col and col in session_row:
            value = str(session_row.get(col, "")).strip()
        elif col and col in row_values:
            value = str(row_values.get(col, "")).strip()
        if label.endswith("mass") or label == "ambient_temperature":
            out[label] = parse_float(value)
        else:
            out[label] = value
    return out


def _is_light(hour_value: float, lights_on: float, hours_light: float) -> bool:
    phase = (float(hour_value) - float(lights_on)) % 24.0
    return phase < float(hours_light)


def _parse_exclusion_rules(exclusion_rows: list[dict[str, str]]) -> dict[str, list[tuple[float | None, str]]]:
    out: dict[str, list[tuple[float | None, str]]] = {}
    for row in exclusion_rows:
        subject_id = str(row.get("subject_id", row.get("subject.id", ""))).strip()
        if not subject_id:
            continue
        start_hour = parse_float(row.get("start_hour", row.get("hour_start", "")))
        reason = str(row.get("reason", "manual_exclusion")).strip() or "manual_exclusion"
        out.setdefault(subject_id, []).append((start_hour, reason))
    return out


def _resolve_time_value(row: CalRRow, *, hour_col: str | None, day_col: str | None, datetime_col: str | None, fallback_index: int) -> tuple[float, str]:
    hour = parse_float(row.values.get(hour_col, "")) if hour_col else None
    if hour is not None:
        day = parse_float(row.values.get(day_col, "")) if day_col else None
        if day is not None:
            return float((day - 1.0) * 24.0 + hour), "table_hour_day"
        return float(hour), "table_hour"
    if datetime_col:
        dt = parse_timestamp(row.values.get(datetime_col, ""))
        if dt is not None:
            return float(dt.hour) + float(dt.minute) / 60.0 + float(dt.day - 1) * 24.0, "datetime"
    return float(fallback_index), "implicit_index"


def extract_subject_features(
    *,
    data_rows: list[CalRRow],
    data_fieldnames: list[str],
    session_rows: dict[str, dict[str, str]] | None,
    session_fieldnames: list[str] | None,
    exclusion_rows: list[dict[str, str]],
    cfg: CalorimetryFeatureConfig,
) -> FeatureExtractionResult:
    warnings: list[str] = []
    subject_col = resolve_subject_id_column(data_fieldnames)
    hour_col, day_col, datetime_col = resolve_time_columns(data_fieldnames)
    run_col = resolve_run_id_column(data_fieldnames)
    group_col = resolve_group_column(data_fieldnames, session_fieldnames or [])
    variable_cols = resolve_variable_columns(data_fieldnames)
    if "ee" in variable_cols:
        ee_source = "provided"
    else:
        ee_source = "absent"
    if "rer" in variable_cols:
        rer_source = "provided"
    else:
        rer_source = "absent"

    if session_rows:
        session_mode = "explicit"
    else:
        session_mode = "exploratory"
        warning = "session file absent; running in exploratory mode and not fully CalR-equivalent"
        warnings.append(warning)
        print(f"warning: {warning}", file=sys.stderr)
        if not cfg.exploratory_without_session:
            raise ValueError("Session file is required unless exploratory_without_session=true")

    session_start, session_end = _selected_window_from_session(session_rows)
    analysis_start = cfg.analysis_start_hour if cfg.analysis_start_hour is not None else session_start
    analysis_end = cfg.analysis_end_hour if cfg.analysis_end_hour is not None else session_end
    if cfg.analysis_start_hour is not None or cfg.analysis_end_hour is not None:
        window_source = "explicit"
    elif session_start is not None or session_end is not None:
        window_source = "session"
    else:
        analysis_start = None
        analysis_end = None
        window_source = "full_trace"
        warning = "analysis window fell back to the full trace; results may not match a typical CalR selected-window analysis"
        warnings.append(warning)
        print(f"warning: {warning}", file=sys.stderr)

    lights_on, photoperiod_source = _photoperiod_from_session(session_rows)
    if cfg.photoperiod_lights_on_hour is not None:
        lights_on = float(cfg.photoperiod_lights_on_hour)
        photoperiod_source = "explicit"
    if lights_on is None:
        lights_on = 7.0

    exclusion_map = _parse_exclusion_rules(exclusion_rows)
    if exclusion_map:
        warnings.append(f"applied explicit exclusions for {len(exclusion_map)} subject(s)")

    rows_by_subject: dict[str, list[tuple[float, CalRRow]]] = {}
    time_source = "implicit_index"
    exclusion_notices_emitted: set[tuple[str, float | None, str]] = set()
    for idx, row in enumerate(data_rows):
        subject_id = str(row.values.get(subject_col, "")).strip()
        if not subject_id:
            continue
        time_value, source = _resolve_time_value(row, hour_col=hour_col, day_col=day_col, datetime_col=datetime_col, fallback_index=idx)
        time_source = source
        if analysis_start is not None and time_value < analysis_start:
            continue
        if analysis_end is not None and time_value > analysis_end:
            continue
        excluded = False
        for start_hour, reason in exclusion_map.get(subject_id, []):
            if start_hour is None or time_value >= float(start_hour):
                notice_key = (subject_id, start_hour, reason)
                if notice_key not in exclusion_notices_emitted:
                    warnings.append(
                        f"excluded subject={subject_id} from hour={start_hour if start_hour is not None else 'all'} reason={reason}"
                    )
                    exclusion_notices_emitted.add(notice_key)
                excluded = True
                break
        if excluded:
            continue
        rows_by_subject.setdefault(subject_id, []).append((time_value, row))

    subjects: dict[str, SubjectSummary] = {}
    feature_names: set[str] = set()
    run_ids: set[str] = set()
    session_keys = session_fieldnames or []
    for subject_id, entries in sorted(rows_by_subject.items()):
        entries.sort(key=lambda item: item[0])
        first_row = entries[0][1].values
        session_row = session_rows.get(subject_id) if session_rows else None
        meta = _subject_level_metadata(first_row, session_row, group_col=group_col, run_col=run_col)
        group_value = str(meta.get("group", "")).strip() or "all"
        run_id = str(meta.get("run_id", "")).strip()
        if run_id:
            run_ids.add(run_id)
        feature_values: dict[str, float] = {}
        hours = [time_value for time_value, _row in entries]
        for canonical_var, column_name in variable_cols.items():
            values = [parse_float(row.values.get(column_name, "")) for _time, row in entries]
            paired = [(time_value, value) for time_value, value in zip(hours, values) if value is not None]
            if not paired:
                continue
            times_var = [time_value for time_value, _value in paired]
            vals_var = [float(value) for _time_value, value in paired]
            mean_selected = _mean(vals_var)
            if mean_selected is not None:
                feature_values[f"{canonical_var}_mean_selected"] = float(mean_selected)
                feature_names.add(f"{canonical_var}_mean_selected")
            sd_selected = _sd(vals_var)
            if sd_selected is not None:
                feature_values[f"{canonical_var}_sd_selected"] = float(sd_selected)
                feature_names.add(f"{canonical_var}_sd_selected")
            light_vals = [value for time_value, value in paired if _is_light(time_value, lights_on, cfg.photoperiod_hours_light)]
            dark_vals = [value for time_value, value in paired if not _is_light(time_value, lights_on, cfg.photoperiod_hours_light)]
            light_mean = _mean(light_vals)
            dark_mean = _mean(dark_vals)
            if light_mean is not None:
                feature_values[f"{canonical_var}_mean_light"] = float(light_mean)
                feature_names.add(f"{canonical_var}_mean_light")
            if dark_mean is not None:
                feature_values[f"{canonical_var}_mean_dark"] = float(dark_mean)
                feature_names.add(f"{canonical_var}_mean_dark")
            if light_mean is not None and dark_mean is not None:
                feature_values[f"{canonical_var}_dark_minus_light"] = float(dark_mean - light_mean)
                feature_names.add(f"{canonical_var}_dark_minus_light")
            if canonical_var in {"feed", "feed.acc", "drink", "drink.acc"}:
                auc = _auc(times_var, vals_var)
                if auc is not None:
                    feature_values[f"{canonical_var}_auc_selected"] = float(auc)
                    feature_names.add(f"{canonical_var}_auc_selected")
            amp1, phase_cos, phase_sin = _circadian(times_var, vals_var)
            if amp1 is not None:
                feature_values[f"{canonical_var}_amp1"] = float(amp1)
                feature_names.add(f"{canonical_var}_amp1")
            if phase_cos is not None and phase_sin is not None:
                feature_values[f"{canonical_var}_phase1_cos"] = float(phase_cos)
                feature_values[f"{canonical_var}_phase1_sin"] = float(phase_sin)
                feature_names.add(f"{canonical_var}_phase1_cos")
                feature_names.add(f"{canonical_var}_phase1_sin")
        subjects[subject_id] = SubjectSummary(
            subject_id=subject_id,
            group=group_value,
            run_id=run_id,
            feature_values=feature_values,
            metadata=meta,
        )

    if len(run_ids) > 1:
        warnings.append("multiple run_id values detected; runs are preserved in metadata and should not be treated as silently interchangeable")
        print("warning: multiple run_id values detected; combined-run interpretation may be exploratory", file=sys.stderr)

    session_rows_used = int(len(session_rows or {}))
    session_missing = max(0, len(subjects) - session_rows_used) if session_rows else len(subjects)
    summary = {
        "session_mode": session_mode,
        "analysis_start": analysis_start,
        "analysis_end": analysis_end,
        "analysis_window_source": window_source,
        "photoperiod_source": photoperiod_source,
        "lights_on_hour": lights_on,
        "time_source": time_source,
        "ee_source": ee_source,
        "rer_source": rer_source,
        "n_subjects": len(subjects),
        "n_features": len(feature_names),
        "n_rows_input": len(data_rows),
        "n_rows_selected": sum(len(entries) for entries in rows_by_subject.values()),
        "n_subjects_missing_session_metadata": session_missing,
        "run_ids": sorted(run_ids),
        "design_class": "simple_grouped",
        "excluded_subjects": [subject_id for subject_id in sorted(exclusion_map)],
    }
    if not session_rows and window_source == "full_trace":
        summary["acclimation_state"] = "unknown_or_mixed"
    elif window_source == "full_trace":
        summary["acclimation_state"] = "mixed_or_full_run"
    else:
        summary["acclimation_state"] = "post_window_selected"
    if window_source == "full_trace":
        warning = "analysis window includes the apparent full run; results may include acclimation and differ from a typical post-acclimation CalR analysis"
        warnings.append(warning)
        print(f"warning: {warning}", file=sys.stderr)

    return FeatureExtractionResult(
        subjects=subjects,
        feature_names=sorted(feature_names),
        metadata_summary=summary,
        warnings=warnings,
    )


def feature_base_variable(feature_name: str) -> str:
    for base in sorted(PROGRAM_VARIABLES["global"], key=len, reverse=True):
        prefix = f"{base}_"
        if feature_name.startswith(prefix):
            return base
    return feature_name.split("_", 1)[0]
