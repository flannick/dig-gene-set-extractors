from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re
from typing import Iterable


SUBJECT_ID_CANDIDATES = ("subject.id", "subject_id", "id", "animal_id", "mouse_id")
DATE_TIME_CANDIDATES = ("Date.Time", "date_time", "datetime", "timestamp")
HOUR_CANDIDATES = ("exp.hour", "hour", "Hour", "time_hours")
DAY_CANDIDATES = ("exp.day", "day", "Day")
RUN_ID_CANDIDATES = ("run_id", "Run", "run")
GROUP_CANDIDATES = ("group", "group_name", "group1", "group2", "group3", "group4")
DIRECT_GROUP_CANDIDATES = ("group", "group_name")
GROUP_NAME_LAYOUT_CANDIDATES = ("group_names", "group_labels", "group.names")
SESSION_ID_CANDIDATES = ("id", "subject.id", "subject_id")
TOTAL_MASS_CANDIDATES = ("Total.Mass", "total.mass", "total_mass", "subject.mass", "body_mass")
LEAN_MASS_CANDIDATES = ("Lean.Mass", "lean.mass", "lean_mass")
FAT_MASS_CANDIDATES = ("Fat.Mass", "fat.mass", "fat_mass")
SEX_CANDIDATES = ("sex", "Sex")
DIET_CANDIDATES = ("diet", "Diet", "diet_name")
STRAIN_CANDIDATES = ("strain", "Strain", "genotype", "Genotype")
AMBIENT_TEMP_CANDIDATES = ("room.temperature.C", "ambient_temperature", "temp", "temperature")

VARIABLE_ALIASES = {
    "vo2": ("vo2", "VO2"),
    "vco2": ("vco2", "VCO2"),
    "ee": ("ee", "EE"),
    "rer": ("rer", "RER"),
    "feed": ("feed", "food", "food.intake"),
    "feed.acc": ("feed.acc", "food.acc", "food.cumulative", "feed_acc"),
    "drink": ("drink", "water", "water.intake"),
    "drink.acc": ("drink.acc", "water.acc", "water.cumulative", "drink_acc"),
    "xytot": ("xytot", "XYTOT"),
    "xyamb": ("xyamb", "XYAMB"),
    "wheel": ("wheel",),
    "pedmeter": ("pedmeter",),
    "body.temp": ("body.temp", "body_temp", "Body.Temp"),
}

MISSING_TOKENS = {"", "NA", "NaN", "nan", "NULL", "null", "None", "none"}
_GROUP_LAYOUT_RE = re.compile(r"^group[0-9]+$")


@dataclass
class CalRRow:
    line_no: int
    values: dict[str, str]


def _clean(value: object) -> str:
    return str(value).strip()


def _is_missing(value: object) -> bool:
    return _clean(value) in MISSING_TOKENS


def _read_delimited(path: str | Path, delimiter: str) -> tuple[list[str], list[CalRRow]]:
    p = Path(path)
    with p.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Missing header in {p}")
        fieldnames = [str(name) for name in reader.fieldnames]
        rows = [CalRRow(line_no=i + 2, values={str(k): str(v or "") for k, v in row.items()}) for i, row in enumerate(reader)]
    return fieldnames, rows


def read_calr_data_csv(path: str | Path, delimiter: str = ",") -> tuple[list[str], list[CalRRow]]:
    return _read_delimited(path, delimiter)


def read_session_csv(path: str | Path | None, delimiter: str = ",") -> tuple[list[str], dict[str, dict[str, str]] | None]:
    if not path:
        return [], None
    fieldnames, rows = _read_delimited(path, delimiter)
    by_id: dict[str, dict[str, str]] = {}
    id_col = find_first_column(fieldnames, SESSION_ID_CANDIDATES)
    if id_col is None:
        raise ValueError(f"Session CSV missing subject id column. Tried: {', '.join(SESSION_ID_CANDIDATES)}")
    for row in rows:
        subject_id = _clean(row.values.get(id_col, ""))
        if not subject_id or _is_missing(subject_id):
            continue
        by_id[subject_id] = dict(row.values)
    derived_fieldnames = list(fieldnames)
    _augment_session_rows(fieldnames=fieldnames, rows=rows, by_id=by_id, derived_fieldnames=derived_fieldnames)
    return derived_fieldnames, by_id


def read_exclusions_tsv(path: str | Path | None, delimiter: str = "\t") -> list[dict[str, str]]:
    if not path:
        return []
    p = Path(path)
    with p.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Exclusions table missing header: {p}")
        return [{str(k): str(v or "") for k, v in row.items()} for row in reader]


def find_first_column(fieldnames: Iterable[str], candidates: Iterable[str]) -> str | None:
    lookup = {str(name).strip().lower(): str(name) for name in fieldnames}
    for candidate in candidates:
        key = str(candidate).strip().lower()
        if key in lookup:
            return lookup[key]
    return None


def resolve_subject_id_column(fieldnames: list[str]) -> str:
    col = find_first_column(fieldnames, SUBJECT_ID_CANDIDATES)
    if col is None:
        raise ValueError(f"CalR data missing subject id column. Tried: {', '.join(SUBJECT_ID_CANDIDATES)}")
    return col


def resolve_time_columns(fieldnames: list[str]) -> tuple[str | None, str | None, str | None]:
    return (
        find_first_column(fieldnames, HOUR_CANDIDATES),
        find_first_column(fieldnames, DAY_CANDIDATES),
        find_first_column(fieldnames, DATE_TIME_CANDIDATES),
    )


def resolve_run_id_column(fieldnames: list[str]) -> str | None:
    return find_first_column(fieldnames, RUN_ID_CANDIDATES)


def resolve_group_column(fieldnames: list[str], session_fieldnames: list[str] | None = None) -> str | None:
    col = find_first_column(fieldnames, GROUP_CANDIDATES)
    if col is not None:
        return col
    if session_fieldnames:
        return find_first_column(session_fieldnames, GROUP_CANDIDATES)
    return None


def resolve_variable_columns(fieldnames: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for canonical, aliases in VARIABLE_ALIASES.items():
        col = find_first_column(fieldnames, aliases)
        if col is not None:
            out[canonical] = col
    return out


def parse_float(value: object) -> float | None:
    text = _clean(value)
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_timestamp(value: str) -> datetime | None:
    text = _clean(value)
    if not text:
        return None
    for fmt in (
        "%Y-%m-%d %H:%M:%S",
        "%m/%d/%Y %H:%M",
        "%m/%d/%y %H:%M",
        "%Y/%m/%d %H:%M:%S",
    ):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    try:
        return datetime.fromisoformat(text)
    except ValueError:
        return None


def numeric_session_value(session_rows: dict[str, dict[str, str]] | None, candidates: tuple[str, ...]) -> float | None:
    if not session_rows:
        return None
    for row in session_rows.values():
        col = find_first_column(list(row.keys()), candidates)
        if col is None:
            continue
        value = parse_float(row.get(col))
        if value is not None:
            return value
    return None


def unique_session_values(session_rows: dict[str, dict[str, str]] | None, candidates: tuple[str, ...]) -> list[str]:
    if not session_rows:
        return []
    out: list[str] = []
    seen: set[str] = set()
    for row in session_rows.values():
        col = find_first_column(list(row.keys()), candidates)
        if col is None:
            continue
        value = _clean(row.get(col, ""))
        if not value or value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def _ordered_unique_nonmissing(rows: list[CalRRow], column: str | None) -> list[str]:
    if not column:
        return []
    values: list[str] = []
    seen: set[str] = set()
    for row in rows:
        value = _clean(row.values.get(column, ""))
        if _is_missing(value) or value in seen:
            continue
        seen.add(value)
        values.append(value)
    return values


def _ordered_unique_numeric(rows: list[CalRRow], column: str | None) -> list[float]:
    if not column:
        return []
    values: list[float] = []
    seen: set[float] = set()
    for row in rows:
        parsed = parse_float(row.values.get(column, ""))
        if parsed is None:
            continue
        marker = round(float(parsed), 8)
        if marker in seen:
            continue
        seen.add(marker)
        values.append(float(parsed))
    return values


def _session_group_membership_columns(fieldnames: list[str], rows: list[CalRRow], subject_ids: set[str]) -> list[str]:
    out: list[str] = []
    for field in fieldnames:
        field_lower = _clean(field).lower()
        if not _GROUP_LAYOUT_RE.match(field_lower):
            continue
        values = [value for value in (_clean(row.values.get(field, "")) for row in rows) if not _is_missing(value)]
        if not values:
            continue
        subject_matches = sum(1 for value in values if value in subject_ids)
        if subject_matches > 0 and (subject_matches / float(len(values))) >= 0.5:
            out.append(field)
    return out


def _session_group_assignments(fieldnames: list[str], rows: list[CalRRow], subject_ids: set[str]) -> tuple[dict[str, str], str | None]:
    direct_group_col = find_first_column(fieldnames, DIRECT_GROUP_CANDIDATES)
    membership_cols = _session_group_membership_columns(fieldnames, rows, subject_ids)
    if membership_cols:
        group_name_col = find_first_column(fieldnames, GROUP_NAME_LAYOUT_CANDIDATES)
        group_labels = _ordered_unique_nonmissing(rows, group_name_col)
        assignments: dict[str, str] = {}
        for idx, group_col in enumerate(membership_cols):
            label = group_labels[idx] if idx < len(group_labels) else group_col
            for row in rows:
                subject_id = _clean(row.values.get(group_col, ""))
                if _is_missing(subject_id):
                    continue
                assignments.setdefault(subject_id, label)
        return assignments, "wide_membership"
    if direct_group_col:
        assignments = {}
        id_col = find_first_column(fieldnames, SESSION_ID_CANDIDATES)
        if id_col is None:
            return assignments, None
        for row in rows:
            subject_id = _clean(row.values.get(id_col, ""))
            group_label = _clean(row.values.get(direct_group_col, ""))
            if _is_missing(subject_id) or _is_missing(group_label):
                continue
            assignments[subject_id] = group_label
        return assignments, "direct_column"
    fallback_group_col = find_first_column(fieldnames, GROUP_CANDIDATES)
    if fallback_group_col:
        assignments = {}
        id_col = find_first_column(fieldnames, SESSION_ID_CANDIDATES)
        if id_col is None:
            return assignments, None
        for row in rows:
            subject_id = _clean(row.values.get(id_col, ""))
            group_label = _clean(row.values.get(fallback_group_col, ""))
            if _is_missing(subject_id) or _is_missing(group_label):
                continue
            assignments[subject_id] = group_label
        return assignments, "fallback_column"
    return {}, None


def _session_window_bounds(fieldnames: list[str], rows: list[CalRRow]) -> tuple[float | None, float | None, str | None]:
    if find_first_column(fieldnames, ("analysis_start", "analysis_start_hour", "start_hour", "start")) or find_first_column(
        fieldnames, ("analysis_end", "analysis_end_hour", "end_hour", "end")
    ):
        return None, None, None
    xrange_col = find_first_column(fieldnames, ("xrange",))
    exc_col = find_first_column(fieldnames, ("exc",))
    xrange_values = _ordered_unique_numeric(rows, xrange_col)
    exc_values = _ordered_unique_numeric(rows, exc_col)
    if not xrange_values and not exc_values:
        return None, None, None
    start = xrange_values[0] if xrange_values else None
    end = None
    if exc_values:
        end = max(exc_values)
    elif len(xrange_values) >= 2:
        end = xrange_values[1]
    return start, end, "calr_matrix"


def _augment_session_rows(
    *,
    fieldnames: list[str],
    rows: list[CalRRow],
    by_id: dict[str, dict[str, str]],
    derived_fieldnames: list[str],
) -> None:
    if not by_id:
        return
    subject_ids = set(by_id)
    assignments, group_mode = _session_group_assignments(fieldnames, rows, subject_ids)
    if assignments and "group" not in {name.lower() for name in derived_fieldnames}:
        derived_fieldnames.append("group")
    start, end, window_mode = _session_window_bounds(fieldnames, rows)
    if start is not None and "analysis_start_hour" not in {name.lower() for name in derived_fieldnames}:
        derived_fieldnames.append("analysis_start_hour")
    if end is not None and "analysis_end_hour" not in {name.lower() for name in derived_fieldnames}:
        derived_fieldnames.append("analysis_end_hour")
    if end is not None and "session_exclusion_hour" not in {name.lower() for name in derived_fieldnames}:
        derived_fieldnames.append("session_exclusion_hour")
    if group_mode and "_session_group_layout_mode" not in derived_fieldnames:
        derived_fieldnames.append("_session_group_layout_mode")
    if window_mode and "_session_window_layout_mode" not in derived_fieldnames:
        derived_fieldnames.append("_session_window_layout_mode")
    for subject_id, row in by_id.items():
        if subject_id in assignments:
            row["group"] = assignments[subject_id]
        if start is not None:
            row["analysis_start_hour"] = str(start)
        if end is not None:
            row["analysis_end_hour"] = str(end)
            row["session_exclusion_hour"] = str(end)
        if group_mode:
            row["_session_group_layout_mode"] = group_mode
        if window_mode:
            row["_session_window_layout_mode"] = window_mode
