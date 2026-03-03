from __future__ import annotations

import csv
from dataclasses import dataclass
import gzip
from pathlib import Path


DEFAULT_CHROM_COLUMNS = ("chromosome", "chrom", "chr")
DEFAULT_START_COLUMNS = ("start", "loc_start")
DEFAULT_END_COLUMNS = ("end", "loc_end")
DEFAULT_AMPLITUDE_COLUMNS = (
    "segment_mean",
    "segmean",
    "log2ratio",
    "log2_copy_ratio",
    "amplitude",
)
DEFAULT_SAMPLE_COLUMNS = ("sample", "sample_id", "tumor_sample_barcode")

CBIO_CHROM_COLUMNS = ("chrom",)
CBIO_START_COLUMNS = ("loc.start",)
CBIO_END_COLUMNS = ("loc.end",)
CBIO_AMPLITUDE_COLUMNS = ("seg.mean",)
CBIO_SAMPLE_COLUMNS = ("ID",)

SUPPORTED_SEGMENT_FORMATS = {"auto", "gdc_seg", "cbio_seg"}


@dataclass
class CNVSegment:
    sample_id: str
    chrom: str
    start: int
    end: int
    amplitude: float
    raw_chrom: str
    raw_start: int
    raw_end: int
    raw_amplitude: float


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _resolve_column(
    fieldnames: list[str],
    explicit: str | None,
    defaults: tuple[str, ...],
    *,
    label: str,
    required: bool,
) -> str | None:
    if explicit is not None and str(explicit).strip():
        chosen = str(explicit).strip()
        if chosen in fieldnames:
            return chosen
        low_map = {x.lower(): x for x in fieldnames}
        mapped = low_map.get(chosen.lower())
        if mapped is not None:
            return mapped
        if required:
            raise ValueError(
                f"Column '{chosen}' for {label} not found. Available columns: {', '.join(fieldnames)}"
            )
        return None
    low_map = {x.lower(): x for x in fieldnames}
    for candidate in defaults:
        if candidate in fieldnames:
            return candidate
        mapped = low_map.get(candidate.lower())
        if mapped is not None:
            return mapped
    if required:
        raise ValueError(
            f"Could not resolve required column for {label}. Tried defaults: {', '.join(defaults)}. "
            f"Available columns: {', '.join(fieldnames)}"
        )
    return None


def _can_resolve(fieldnames: list[str], defaults: tuple[str, ...]) -> bool:
    low_map = {x.lower(): x for x in fieldnames}
    for candidate in defaults:
        if candidate in fieldnames:
            return True
        if low_map.get(candidate.lower()) is not None:
            return True
    return False


def _is_cbio_seg_header(fieldnames: list[str]) -> bool:
    required = {"id", "chrom", "loc.start", "loc.end", "seg.mean"}
    observed = {str(x).strip().lower() for x in fieldnames}
    return required.issubset(observed)


def _detect_segments_format(fieldnames: list[str]) -> str:
    cbio_like = _is_cbio_seg_header(fieldnames)
    gdc_like = (
        _can_resolve(fieldnames, DEFAULT_CHROM_COLUMNS)
        and _can_resolve(fieldnames, DEFAULT_START_COLUMNS)
        and _can_resolve(fieldnames, DEFAULT_END_COLUMNS)
        and _can_resolve(fieldnames, DEFAULT_AMPLITUDE_COLUMNS)
    )
    if cbio_like and gdc_like:
        return "ambiguous"
    if cbio_like:
        return "cbio_seg"
    if gdc_like:
        return "gdc_seg"
    return "unknown"


def _defaults_for_segments_format(segments_format: str) -> tuple[
    tuple[str, ...], tuple[str, ...], tuple[str, ...], tuple[str, ...], tuple[str, ...]
]:
    fmt = str(segments_format).strip().lower()
    if fmt == "cbio_seg":
        return (
            CBIO_CHROM_COLUMNS,
            CBIO_START_COLUMNS,
            CBIO_END_COLUMNS,
            CBIO_AMPLITUDE_COLUMNS,
            CBIO_SAMPLE_COLUMNS,
        )
    if fmt == "gdc_seg":
        return (
            DEFAULT_CHROM_COLUMNS,
            DEFAULT_START_COLUMNS,
            DEFAULT_END_COLUMNS,
            DEFAULT_AMPLITUDE_COLUMNS,
            DEFAULT_SAMPLE_COLUMNS,
        )
    raise ValueError(f"Unsupported segments_format: {segments_format}")


def _parse_int(raw: object) -> int | None:
    if raw is None:
        return None
    text = str(raw).strip()
    if not text:
        return None
    try:
        if "." in text:
            return int(float(text))
        return int(text)
    except ValueError:
        return None


def _parse_float(raw: object) -> float | None:
    if raw is None:
        return None
    text = str(raw).strip()
    if not text:
        return None
    low = text.lower()
    if low in {"na", "nan", "none", "null", "inf", "-inf"}:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if value != value:
        return None
    return float(value)


def _apply_chrom_prefix_mode(chrom: str, mode: str) -> str:
    c = str(chrom).strip()
    if not c:
        return ""
    if mode == "none" or mode == "auto":
        return c
    if mode == "add_chr":
        return c if c.lower().startswith("chr") else f"chr{c}"
    if mode == "drop_chr":
        return c[3:] if c.lower().startswith("chr") else c
    raise ValueError(f"Unsupported chrom_prefix_mode: {mode}")


def _to_zero_based_half_open(start: int, end: int, coord_system: str) -> tuple[int, int]:
    s = int(start)
    e = int(end)
    if coord_system == "one_based_closed":
        return s - 1, e
    if coord_system == "zero_based_half_open":
        return s, e
    raise ValueError(f"Unsupported coord_system: {coord_system}")


def parse_segments_tsv(
    *,
    path: str | Path,
    segments_format: str,
    chrom_column: str | None,
    start_column: str | None,
    end_column: str | None,
    amplitude_column: str | None,
    sample_id_column: str | None,
    coord_system: str,
    chrom_prefix_mode: str,
) -> tuple[list[CNVSegment], dict[str, object]]:
    segments: list[CNVSegment] = []
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Segments table has no header: {path}")
        fieldnames = [str(x) for x in reader.fieldnames]
        requested_format = str(segments_format).strip().lower() or "auto"
        if requested_format not in SUPPORTED_SEGMENT_FORMATS:
            raise ValueError(
                "Unsupported --segments_format "
                f"'{segments_format}'. Expected one of: {', '.join(sorted(SUPPORTED_SEGMENT_FORMATS))}"
            )
        has_explicit_columns = any(
            (x is not None and str(x).strip())
            for x in (chrom_column, start_column, end_column, amplitude_column, sample_id_column)
        )
        detected_format = _detect_segments_format(fieldnames)
        if requested_format == "auto":
            if has_explicit_columns:
                resolved_format = "gdc_seg"
            elif detected_format in {"gdc_seg", "cbio_seg"}:
                resolved_format = detected_format
            elif detected_format == "ambiguous":
                raise ValueError(
                    "Could not auto-detect a unique segment format; header matches multiple known patterns. "
                    "Use --segments_format gdc_seg|cbio_seg or pass explicit --chrom_column/--start_column/"
                    "--end_column/--amplitude_column. "
                    f"Available columns: {', '.join(fieldnames)}"
                )
            else:
                raise ValueError(
                    "Could not auto-detect segment format from header. "
                    "Use --segments_format gdc_seg|cbio_seg or pass explicit --chrom_column/--start_column/"
                    "--end_column/--amplitude_column. "
                    f"Available columns: {', '.join(fieldnames)}"
                )
        else:
            resolved_format = requested_format
        (
            default_chrom_cols,
            default_start_cols,
            default_end_cols,
            default_amp_cols,
            default_sample_cols,
        ) = _defaults_for_segments_format(resolved_format)

        chrom_col = _resolve_column(fieldnames, chrom_column, default_chrom_cols, label="chrom", required=True)
        start_col = _resolve_column(fieldnames, start_column, default_start_cols, label="start", required=True)
        end_col = _resolve_column(fieldnames, end_column, default_end_cols, label="end", required=True)
        amp_col = _resolve_column(fieldnames, amplitude_column, default_amp_cols, label="amplitude", required=True)
        sample_col = _resolve_column(fieldnames, sample_id_column, default_sample_cols, label="sample_id", required=False)

        n_rows = 0
        n_missing_required = 0
        n_invalid_coords = 0
        n_negative_coords = 0
        n_non_numeric_amplitude = 0
        for row in reader:
            n_rows += 1
            raw_chrom = str(row.get(chrom_col or "", "")).strip()
            start_raw = _parse_int(row.get(start_col or ""))
            end_raw = _parse_int(row.get(end_col or ""))
            amp_raw = _parse_float(row.get(amp_col or ""))
            if not raw_chrom or start_raw is None or end_raw is None or amp_raw is None:
                n_missing_required += 1
                if amp_raw is None:
                    n_non_numeric_amplitude += 1
                continue
            start_0, end_0 = _to_zero_based_half_open(start_raw, end_raw, coord_system)
            if start_0 < 0 or end_0 < 0:
                n_negative_coords += 1
                continue
            if end_0 <= start_0:
                n_invalid_coords += 1
                continue
            sample_id = (
                str(row.get(sample_col, "")).strip()
                if sample_col is not None
                else ""
            ) or "sample1"
            segments.append(
                CNVSegment(
                    sample_id=sample_id,
                    chrom=_apply_chrom_prefix_mode(raw_chrom, chrom_prefix_mode),
                    start=int(start_0),
                    end=int(end_0),
                    amplitude=float(amp_raw),
                    raw_chrom=raw_chrom,
                    raw_start=int(start_raw),
                    raw_end=int(end_raw),
                    raw_amplitude=float(amp_raw),
                )
            )

    if not segments:
        raise ValueError("No usable segments parsed from input table.")

    parse_summary: dict[str, object] = {
        "input_path": str(p),
        "n_rows": n_rows,
        "n_segments": len(segments),
        "n_rows_missing_required": n_missing_required,
        "n_rows_invalid_coords": n_invalid_coords,
        "n_rows_negative_coords": n_negative_coords,
        "n_rows_non_numeric_amplitude": n_non_numeric_amplitude,
        "segments_format_requested": requested_format,
        "segments_format_resolved": resolved_format,
        "segments_format_detected": detected_format,
        "coord_system": coord_system,
        "chrom_prefix_mode": chrom_prefix_mode,
        "columns": {
            "chrom": chrom_col,
            "start": start_col,
            "end": end_col,
            "amplitude": amp_col,
            "sample_id": sample_col if sample_col is not None else "sample1 (implicit)",
        },
    }
    return segments, parse_summary


def read_purity_tsv(
    *,
    path: str | Path,
    sample_id_column: str,
    purity_value_column: str,
) -> tuple[dict[str, float], dict[str, object]]:
    purity_by_sample: dict[str, float] = {}
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Purity table has no header: {path}")
        fieldnames = [str(x) for x in reader.fieldnames]
        sample_col = _resolve_column(
            fieldnames,
            sample_id_column,
            ("sample_id", "sample"),
            label="purity sample_id",
            required=True,
        )
        purity_col = _resolve_column(
            fieldnames,
            purity_value_column,
            ("purity",),
            label="purity value",
            required=True,
        )
        n_rows = 0
        n_invalid = 0
        for row in reader:
            n_rows += 1
            sample = str(row.get(sample_col or "", "")).strip()
            purity = _parse_float(row.get(purity_col or ""))
            if not sample or purity is None:
                n_invalid += 1
                continue
            purity_by_sample[sample] = float(purity)
    if not purity_by_sample:
        raise ValueError("No usable purity rows parsed.")
    summary: dict[str, object] = {
        "input_path": str(p),
        "n_rows": n_rows,
        "n_samples_with_purity": len(purity_by_sample),
        "n_rows_invalid": n_invalid,
    }
    return purity_by_sample, summary
