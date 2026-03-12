from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


SUBJECT_ID_CANDIDATES = ("subject.id", "subject_id", "id", "animal_id", "mouse_id")
DATE_TIME_CANDIDATES = ("Date.Time", "date_time", "datetime", "timestamp")
HOUR_CANDIDATES = ("exp.hour", "hour", "Hour", "time_hours")
DAY_CANDIDATES = ("exp.day", "day", "Day")
RUN_ID_CANDIDATES = ("run_id", "Run", "run")
GROUP_CANDIDATES = ("group", "group_name", "group1", "group2", "group3", "group4")
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


@dataclass
class CalRRow:
    line_no: int
    values: dict[str, str]


def _clean(value: object) -> str:
    return str(value).strip()


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
        if not subject_id:
            continue
        by_id[subject_id] = dict(row.values)
    return fieldnames, by_id


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
