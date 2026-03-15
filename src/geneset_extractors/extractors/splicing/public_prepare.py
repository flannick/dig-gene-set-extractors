from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

from geneset_extractors.extractors.splicing.splice_event_diff_workflow import _normalize_event_type
from geneset_extractors.hashing import sha256_file


RAW_EVENT_ID_COLUMNS = ("event_id", "splice_event_id", "as_id", "Event")
RAW_GENE_SYMBOL_COLUMNS = ("gene_symbol", "symbol", "Gene", "gene_name")
RAW_GENE_ID_COLUMNS = ("gene_id", "GeneID")
RAW_EVENT_GROUP_COLUMNS = ("event_group", "event_group_id", "cluster_id", "cluster", "lsv_id")
RAW_EVENT_TYPE_COLUMNS = ("event_type", "splice_type", "Type")
RAW_CHROM_COLUMNS = ("chrom", "chr", "Chromosome")
RAW_START_COLUMNS = ("start", "Start")
RAW_END_COLUMNS = ("end", "End")
RAW_STRAND_COLUMNS = ("strand", "Strand")
RAW_ANNOTATION_COLUMNS = ("annotation_status", "annotation", "status")
RAW_READ_SUPPORT_COLUMNS = ("read_support", "ReadSupport", "junction_reads")
KNOWN_METADATA_COLUMNS = {
    "event_id",
    "splice_event_id",
    "as_id",
    "Event",
    "event_group",
    "event_group_id",
    "cluster_id",
    "cluster",
    "lsv_id",
    "event_type",
    "splice_type",
    "Type",
    "gene_id",
    "GeneID",
    "gene_symbol",
    "symbol",
    "Gene",
    "gene_name",
    "annotation_status",
    "annotation",
    "status",
    "chrom",
    "chr",
    "Chromosome",
    "start",
    "Start",
    "end",
    "End",
    "strand",
    "Strand",
    "read_support",
    "ReadSupport",
    "junction_reads",
}


def _clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _parse_float(value: object) -> float | None:
    text = _clean(value)
    if not text or text.lower() in {"na", "nan", "none", "null"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _first_value(row: dict[str, str], candidates: tuple[str, ...]) -> str:
    for key in candidates:
        value = _clean(row.get(key))
        if value:
            return value
    return ""


def _read_tsv(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"TSV has no header: {path}")
        return [str(x) for x in reader.fieldnames], [{k: str(v) for k, v in row.items() if k is not None} for row in reader]


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _event_key(row: dict[str, str]) -> str:
    raw = _first_value(row, RAW_EVENT_ID_COLUMNS)
    if raw:
        return raw
    chrom = _first_value(row, RAW_CHROM_COLUMNS)
    start = _first_value(row, RAW_START_COLUMNS)
    end = _first_value(row, RAW_END_COLUMNS)
    strand = _first_value(row, RAW_STRAND_COLUMNS) or "."
    event_type = _normalize_event_type(_first_value(row, RAW_EVENT_TYPE_COLUMNS))
    if chrom and start and end:
        return f"{chrom}:{start}-{end}:{strand}:{event_type}"
    return f"row_{abs(hash(tuple(sorted(row.items()))))}"


def _canonical_event_key(row: dict[str, str], *, input_mode: str) -> tuple[str, str, str, str]:
    chrom = _first_value(row, RAW_CHROM_COLUMNS)
    start = _first_value(row, RAW_START_COLUMNS)
    end = _first_value(row, RAW_END_COLUMNS)
    strand = _first_value(row, RAW_STRAND_COLUMNS) or "."
    event_type = _normalize_event_type(_first_value(row, RAW_EVENT_TYPE_COLUMNS))
    raw_key = _event_key(row)
    group_key = _first_value(row, RAW_EVENT_GROUP_COLUMNS)
    if chrom and start and end:
        suffix = group_key or raw_key
        return (
            f"{chrom}:{start}-{end}:{strand}:{event_type}:{suffix}",
            "coordinate_canonical",
            "high",
            "global_coordinate",
        )
    if input_mode == "tcga_spliceseq":
        stable_event_id = _first_value(row, RAW_EVENT_ID_COLUMNS)
        gene_symbol = _first_value(row, RAW_GENE_SYMBOL_COLUMNS) or _first_value(row, RAW_GENE_ID_COLUMNS)
        if stable_event_id and gene_symbol and event_type and event_type != "unknown":
            return (
                f"tcga_spliceseq_asid::{gene_symbol}::{event_type}::{stable_event_id}",
                "source_family_stable_id",
                "medium",
                "tcga_spliceseq_asid",
            )
    return raw_key, "raw_id_fallback", "low", "raw_id_fallback"


def _load_sample_id_map(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    _fieldnames, rows = _read_tsv(path)
    out: dict[str, str] = {}
    for row in rows:
        raw = _clean(row.get("sample_id_raw")) or _clean(row.get("source_sample_id"))
        resolved = _clean(row.get("sample_id"))
        if raw and resolved:
            out[raw] = resolved
    return out


def _load_sample_annotations(path: str, study_id: str, study_label: str, sample_id_map: dict[str, str]) -> tuple[list[str], list[dict[str, str]], dict[str, object]]:
    fieldnames, rows = _read_tsv(path)
    sample_id_col = "sample_id" if "sample_id" in fieldnames else fieldnames[0]
    out_rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in rows:
        sample_id_raw = _clean(row.get(sample_id_col))
        sample_id = sample_id_map.get(sample_id_raw, sample_id_raw)
        if not sample_id:
            continue
        if sample_id in seen:
            raise ValueError(f"Duplicate sample_id after mapping: {sample_id}")
        seen.add(sample_id)
        record = {k: _clean(v) for k, v in row.items()}
        record["sample_id_raw"] = sample_id_raw
        record["sample_id"] = sample_id
        record.setdefault("condition", "")
        record.setdefault("group", study_id)
        record["study_id"] = study_id
        record["study_label"] = study_label
        record["annotation_source"] = "explicit"
        out_rows.append(record)
    out_fieldnames = [
        "sample_id",
        "sample_id_raw",
        "condition",
        "group",
        "study_id",
        "study_label",
        "annotation_source",
    ] + [f for f in fieldnames if f not in {sample_id_col, "sample_id", "sample_id_raw", "condition", "group", "study_id", "study_label", "annotation_source"}]
    return out_fieldnames, out_rows, {
        "sample_annotation_source": "explicit",
        "sample_inference_rule": "explicit_table",
        "n_tumor_samples": sum(1 for row in out_rows if _clean(row.get("condition")).lower() in {"tumor", "case"}),
        "n_adjacent_normal_samples": sum(1 for row in out_rows if _clean(row.get("condition")).lower() in {"adjacent_normal", "control", "normal"}),
        "barcode_validation": {
            "attempted": False,
            "n_checked": 0,
            "n_agree": 0,
            "n_disagree": 0,
        },
        "warnings": [],
    }


def _infer_tcga_condition_from_barcode(sample_id: str) -> str | None:
    token = _clean(sample_id)
    if not token.startswith("TCGA-"):
        return None
    parts = token.split("-")
    if len(parts) < 4:
        return None
    sample_code = parts[3][:2]
    if not sample_code.isdigit():
        return None
    code = int(sample_code)
    if 1 <= code <= 9:
        return "tumor"
    if 10 <= code <= 19:
        return "adjacent_normal"
    return None


def _infer_sample_annotations_from_header(
    *,
    psi_fieldnames: list[str],
    study_id: str,
    study_label: str,
    sample_id_map: dict[str, str],
) -> tuple[list[str], list[dict[str, str]], dict[str, object]]:
    sample_columns = [col for col in psi_fieldnames if col not in KNOWN_METADATA_COLUMNS]
    if not sample_columns:
        raise ValueError("Could not infer any sample columns from the PSI matrix header.")
    out_rows: list[dict[str, str]] = []
    seen: set[str] = set()
    n_tumor = 0
    n_adjacent = 0
    barcode_checked = 0
    barcode_agree = 0
    barcode_disagree = 0
    warnings: list[str] = []
    for raw_header in sample_columns:
        header = _clean(raw_header)
        is_norm = header.endswith("_Norm")
        base_sample = header[:-5] if is_norm else header
        sample_id = sample_id_map.get(base_sample, base_sample)
        if not sample_id:
            continue
        if sample_id in seen:
            raise ValueError(
                f"Duplicate inferred sample_id after _Norm stripping: {sample_id}. "
                "Provide --sample_annotations_tsv to disambiguate these columns."
            )
        seen.add(sample_id)
        condition = "adjacent_normal" if is_norm else "tumor"
        if condition == "tumor":
            n_tumor += 1
        else:
            n_adjacent += 1
        barcode_condition = _infer_tcga_condition_from_barcode(base_sample)
        if barcode_condition is not None:
            barcode_checked += 1
            if barcode_condition == condition:
                barcode_agree += 1
            else:
                barcode_disagree += 1
                warnings.append(f"sample_inference_barcode_disagree sample_id={sample_id} suffix_condition={condition} barcode_condition={barcode_condition}")
        out_rows.append(
            {
                "sample_id": sample_id,
                "sample_id_raw": header,
                "condition": condition,
                "group": study_id,
                "study_id": study_id,
                "study_label": study_label,
                "annotation_source": "inferred_header_suffix",
                "inference_rule": "tcga_spliceseq_header_norm_suffix",
                "barcode_condition": barcode_condition or "",
            }
        )
    fieldnames = [
        "sample_id",
        "sample_id_raw",
        "condition",
        "group",
        "study_id",
        "study_label",
        "annotation_source",
        "inference_rule",
        "barcode_condition",
    ]
    return fieldnames, out_rows, {
        "sample_annotation_source": "inferred_header",
        "sample_inference_rule": "tcga_spliceseq_norm_suffix",
        "n_tumor_samples": n_tumor,
        "n_adjacent_normal_samples": n_adjacent,
        "barcode_validation": {
            "attempted": True,
            "n_checked": barcode_checked,
            "n_agree": barcode_agree,
            "n_disagree": barcode_disagree,
        },
        "warnings": warnings,
    }


def run_public_prepare(
    *,
    input_mode: str,
    psi_tsv: str | None,
    sample_annotations_tsv: str | None,
    out_dir: str | Path,
    organism: str,
    genome_build: str,
    study_id: str | None,
    study_label: str | None,
    event_metadata_tsv: str | None = None,
    sample_id_map_tsv: str | None = None,
    event_type_allowlist: str | None = None,
    missing_psi_policy: str = "retain",
) -> dict[str, object]:
    if input_mode != "tcga_spliceseq":
        raise ValueError(f"Unsupported input_mode: {input_mode}")
    if not psi_tsv:
        raise ValueError("tcga_spliceseq mode requires --psi_tsv")

    out_root = Path(out_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    study_id_resolved = str(study_id or Path(psi_tsv).stem)
    study_label_resolved = str(study_label or study_id_resolved)
    sample_id_map = _load_sample_id_map(sample_id_map_tsv)

    psi_fieldnames, psi_rows = _read_tsv(psi_tsv)
    if sample_annotations_tsv:
        sample_meta_fields, sample_meta_rows, sample_meta_summary = _load_sample_annotations(
            sample_annotations_tsv,
            study_id_resolved,
            study_label_resolved,
            sample_id_map,
        )
    else:
        sample_meta_fields, sample_meta_rows, sample_meta_summary = _infer_sample_annotations_from_header(
            psi_fieldnames=psi_fieldnames,
            study_id=study_id_resolved,
            study_label=study_label_resolved,
            sample_id_map=sample_id_map,
        )

    sample_ids = [row["sample_id"] for row in sample_meta_rows]
    raw_to_sample = {row["sample_id_raw"]: row["sample_id"] for row in sample_meta_rows}

    metadata_rows_by_event: dict[str, dict[str, str]] = {}
    if event_metadata_tsv:
        meta_fields, meta_rows = _read_tsv(event_metadata_tsv)
        meta_id_col = "event_id" if "event_id" in meta_fields else meta_fields[0]
        metadata_rows_by_event = {_clean(row.get(meta_id_col)): row for row in meta_rows}

    allowlist = {_normalize_event_type(x) for x in str(event_type_allowlist or "").split(",") if _clean(x)}
    unknown_event_types = 0
    parser_status_counts: dict[str, int] = {}
    canonicalization_confidence_counts: dict[str, int] = {"high": 0, "medium": 0, "low": 0}
    sample_id_map_rows: list[dict[str, object]] = []
    for raw_id, sample_id in raw_to_sample.items():
        sample_id_map_rows.append({"sample_id_raw": raw_id, "sample_id": sample_id, "drop_status": "retained"})

    event_metadata_rows: list[dict[str, object]] = []
    event_id_map_rows: list[dict[str, object]] = []
    bundle_source_rows: list[dict[str, object]] = []
    psi_matrix_rows: list[dict[str, object]] = []
    dropped_events = 0
    missing_values = 0
    total_values = 0
    event_type_counts: dict[str, int] = {}

    sample_columns = [col for col in psi_fieldnames if col in raw_to_sample]
    if not sample_columns:
        raise ValueError("No sample columns from explicit or inferred sample annotations were found in psi_tsv")

    for row in psi_rows:
        row_values = dict(row)
        event_id = _first_value(row, RAW_EVENT_ID_COLUMNS)
        if event_id and event_id in metadata_rows_by_event:
            for key, value in metadata_rows_by_event[event_id].items():
                row_values.setdefault(key, value)
        raw_event_key = _event_key(row_values)
        canonical_event_key, canonicalization_status, canonicalization_confidence, event_key_namespace = _canonical_event_key(
            row_values,
            input_mode=input_mode,
        )
        parser_status_counts[canonicalization_status] = parser_status_counts.get(canonicalization_status, 0) + 1
        canonicalization_confidence_counts[canonicalization_confidence] = canonicalization_confidence_counts.get(canonicalization_confidence, 0) + 1
        event_type = _normalize_event_type(_first_value(row_values, RAW_EVENT_TYPE_COLUMNS))
        if not event_type:
            event_type = "unknown"
        event_type_counts[event_type] = event_type_counts.get(event_type, 0) + 1
        if event_type == "unknown":
            unknown_event_types += 1
        if allowlist and event_type not in allowlist:
            dropped_events += 1
            continue
        gene_id = _first_value(row_values, RAW_GENE_ID_COLUMNS) or _first_value(row_values, RAW_GENE_SYMBOL_COLUMNS)
        gene_symbol = _first_value(row_values, RAW_GENE_SYMBOL_COLUMNS) or gene_id
        meta_row = {
            "event_id": raw_event_key,
            "canonical_event_key": canonical_event_key,
            "event_group": _first_value(row_values, RAW_EVENT_GROUP_COLUMNS),
            "event_type": event_type,
            "gene_id": gene_id,
            "gene_symbol": gene_symbol,
            "annotation_status": _first_value(row_values, RAW_ANNOTATION_COLUMNS),
            "chrom": _first_value(row_values, RAW_CHROM_COLUMNS),
            "start": _first_value(row_values, RAW_START_COLUMNS),
            "end": _first_value(row_values, RAW_END_COLUMNS),
            "strand": _first_value(row_values, RAW_STRAND_COLUMNS),
            "canonicalization_status": canonicalization_status,
            "canonicalization_confidence": canonicalization_confidence,
            "event_key_namespace": event_key_namespace,
            "parser_status": canonicalization_status,
            "source_study_id": study_id_resolved,
            "source_file_name": Path(psi_tsv).name,
            "tool_family": "tcga_spliceseq",
        }
        event_metadata_rows.append(meta_row)
        event_id_map_rows.append(
            {
                "input_event_key": raw_event_key,
                "canonical_event_key": canonical_event_key,
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "event_type": event_type,
                "canonicalization_status": canonicalization_status,
                "canonicalization_confidence": canonicalization_confidence,
                "event_key_namespace": event_key_namespace,
                "parser_status": canonicalization_status,
            }
        )
        psi_row = {"event_id": raw_event_key}
        for key in (
            "canonical_event_key",
            "event_group",
            "event_type",
            "gene_id",
            "gene_symbol",
            "annotation_status",
            "chrom",
            "start",
            "end",
            "strand",
            "canonicalization_status",
            "canonicalization_confidence",
            "event_key_namespace",
        ):
            psi_row[key] = meta_row[key]
        for raw_sample in sample_columns:
            sample_id = raw_to_sample[raw_sample]
            value = _clean(row_values.get(raw_sample))
            psi_row[sample_id] = value
            total_values += 1
            parsed_value = _parse_float(value)
            if parsed_value is None:
                missing_values += 1
                if missing_psi_policy == "drop":
                    continue
            bundle_source_rows.append(
                {
                    "source_dataset": study_id_resolved,
                    "sample_id": sample_id,
                    "input_event_key": raw_event_key,
                    "canonical_event_key": canonical_event_key,
                    "canonicalization_status": canonicalization_status,
                    "canonicalization_confidence": canonicalization_confidence,
                    "event_key_namespace": event_key_namespace,
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "event_type": event_type,
                    "psi": value,
                    "read_support": _first_value(row_values, RAW_READ_SUPPORT_COLUMNS),
                    "annotation_status": meta_row["annotation_status"],
                    "chrom": meta_row["chrom"],
                    "start": meta_row["start"],
                    "end": meta_row["end"],
                    "strand": meta_row["strand"],
                    "tool_family": "tcga_spliceseq",
                }
            )
        psi_matrix_rows.append(psi_row)

    outputs = {
        "psi_matrix_tsv": str(out_root / "psi_matrix.tsv"),
        "sample_metadata_tsv": str(out_root / "sample_metadata.tsv"),
        "event_metadata_tsv": str(out_root / "event_metadata.tsv"),
        "sample_id_map_tsv": str(out_root / "sample_id_map.tsv"),
        "event_id_map_tsv": str(out_root / "event_id_map.tsv"),
        "bundle_source_row_tsv": str(out_root / "bundle_source_row.tsv"),
    }
    _write_tsv(
        out_root / "psi_matrix.tsv",
        [
            "event_id",
            "canonical_event_key",
            "event_group",
            "event_type",
            "gene_id",
            "gene_symbol",
            "annotation_status",
            "chrom",
            "start",
            "end",
            "strand",
            "canonicalization_status",
            "canonicalization_confidence",
            "event_key_namespace",
        ] + sample_ids,
        psi_matrix_rows,
    )
    _write_tsv(out_root / "sample_metadata.tsv", sample_meta_fields, sample_meta_rows)
    _write_tsv(
        out_root / "event_metadata.tsv",
        [
            "event_id",
            "canonical_event_key",
            "event_group",
            "event_type",
            "gene_id",
            "gene_symbol",
            "annotation_status",
            "chrom",
            "start",
            "end",
            "strand",
            "canonicalization_status",
            "canonicalization_confidence",
            "event_key_namespace",
            "parser_status",
            "source_study_id",
            "source_file_name",
            "tool_family",
        ],
        event_metadata_rows,
    )
    _write_tsv(out_root / "sample_id_map.tsv", ["sample_id_raw", "sample_id", "drop_status"], sample_id_map_rows)
    _write_tsv(
        out_root / "event_id_map.tsv",
        [
            "input_event_key",
            "canonical_event_key",
            "gene_id",
            "gene_symbol",
            "event_type",
            "canonicalization_status",
            "canonicalization_confidence",
            "event_key_namespace",
            "parser_status",
        ],
        event_id_map_rows,
    )
    _write_tsv(
        out_root / "bundle_source_row.tsv",
        [
            "source_dataset",
            "sample_id",
            "input_event_key",
            "canonical_event_key",
            "canonicalization_status",
            "canonicalization_confidence",
            "event_key_namespace",
            "gene_id",
            "gene_symbol",
            "event_type",
            "psi",
            "read_support",
            "annotation_status",
            "chrom",
            "start",
            "end",
            "strand",
            "tool_family",
        ],
        bundle_source_rows,
    )

    warnings: list[str] = list(sample_meta_summary.get("warnings", []))
    if event_metadata_rows and (unknown_event_types / float(len(event_metadata_rows))) >= 0.8:
        warnings.append("event_types_overwhelmingly_unknown")
    if event_metadata_rows and (len([r for r in event_metadata_rows if r["canonicalization_status"] == "raw_id_fallback"]) / float(len(event_metadata_rows))) >= 0.8:
        warnings.append("most_events_failed_canonicalization")
    if event_metadata_rows and (len([r for r in event_metadata_rows if r["canonicalization_confidence"] == "medium"]) / float(len(event_metadata_rows))) >= 0.5:
        warnings.append("many_events_use_tcga_medium_confidence_ids")
    if total_values and (missing_values / float(total_values)) >= 0.9:
        warnings.append("psi_matrix_mostly_missing")
    if sample_meta_summary.get("sample_annotation_source") == "inferred_header":
        print(
            "warning: sample annotations were inferred from TCGA SpliceSeq header suffixes; explicit sample annotations are preferred when available.",
            file=sys.stderr,
        )
    for warning in warnings:
        if warning.startswith("sample_inference_barcode_disagree"):
            print(f"warning: {warning}; preferring explicit _Norm suffix rule for sample typing.", file=sys.stderr)

    summary = {
        "workflow": "splice_prepare_public",
        "input_mode": input_mode,
        "study_id": study_id_resolved,
        "study_label": study_label_resolved,
        "organism": organism,
        "genome_build": genome_build,
        "n_samples": len(sample_meta_rows),
        "n_events": len(event_metadata_rows),
        "fraction_missing_psi": (float(missing_values) / float(total_values)) if total_values else 0.0,
        "event_type_distribution": event_type_counts,
        "canonicalization_status_counts": parser_status_counts,
        "canonicalization_confidence_counts": canonicalization_confidence_counts,
        "dropped_events": dropped_events,
        "bundle_source_rows": len(bundle_source_rows),
        "bundle_source_eligible_events": len({row["canonical_event_key"] for row in bundle_source_rows}),
        "sample_annotation_source": sample_meta_summary.get("sample_annotation_source"),
        "sample_inference_rule": sample_meta_summary.get("sample_inference_rule"),
        "n_tumor_samples": sample_meta_summary.get("n_tumor_samples", 0),
        "n_adjacent_normal_samples": sample_meta_summary.get("n_adjacent_normal_samples", 0),
        "barcode_validation": sample_meta_summary.get("barcode_validation"),
        "warnings": warnings,
        "outputs": outputs,
        "source_files": {
            "psi_tsv": {"path": str(Path(psi_tsv)), "sha256": sha256_file(Path(psi_tsv))},
        },
    }
    if sample_annotations_tsv:
        summary["source_files"]["sample_annotations_tsv"] = {
            "path": str(Path(sample_annotations_tsv)),
            "sha256": sha256_file(Path(sample_annotations_tsv)),
        }
    if event_metadata_tsv:
        summary["source_files"]["event_metadata_tsv"] = {"path": str(Path(event_metadata_tsv)), "sha256": sha256_file(Path(event_metadata_tsv))}
    if sample_id_map_tsv:
        summary["source_files"]["sample_id_map_tsv"] = {"path": str(Path(sample_id_map_tsv)), "sha256": sha256_file(Path(sample_id_map_tsv))}

    (out_root / "prepare_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary
