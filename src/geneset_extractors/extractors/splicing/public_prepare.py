from __future__ import annotations

import csv
import json
from pathlib import Path

from geneset_extractors.hashing import sha256_file
from geneset_extractors.extractors.splicing.splice_event_diff_workflow import _normalize_event_type



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
    raw = _clean(row.get("event_id")) or _clean(row.get("splice_event_id"))
    if raw:
        return raw
    chrom = _clean(row.get("chrom"))
    start = _clean(row.get("start"))
    end = _clean(row.get("end"))
    strand = _clean(row.get("strand")) or "."
    event_type = _normalize_event_type(_clean(row.get("event_type")))
    if chrom and start and end:
        return f"{chrom}:{start}-{end}:{strand}:{event_type}"
    return f"row_{abs(hash(tuple(sorted(row.items()))))}"



def _canonical_event_key(row: dict[str, str]) -> tuple[str, str]:
    chrom = _clean(row.get("chrom"))
    start = _clean(row.get("start"))
    end = _clean(row.get("end"))
    strand = _clean(row.get("strand")) or "."
    event_type = _normalize_event_type(_clean(row.get("event_type")))
    raw_key = _event_key(row)
    if chrom and start and end:
        return f"{chrom}:{start}-{end}:{strand}:{event_type}:{raw_key}", "coordinate_plus_id"
    return raw_key, "raw_id_fallback"



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



def _load_sample_annotations(path: str, study_id: str, study_label: str, sample_id_map: dict[str, str]) -> tuple[list[str], list[dict[str, str]]]:
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
        out_rows.append(record)
    out_fieldnames = ["sample_id", "sample_id_raw", "condition", "group", "study_id", "study_label"] + [f for f in fieldnames if f not in {sample_id_col, "sample_id", "sample_id_raw", "condition", "group", "study_id", "study_label"}]
    return out_fieldnames, out_rows



def run_public_prepare(*, input_mode: str, psi_tsv: str | None, sample_annotations_tsv: str | None, out_dir: str | Path, organism: str, genome_build: str, study_id: str | None, study_label: str | None, event_metadata_tsv: str | None = None, sample_id_map_tsv: str | None = None, event_type_allowlist: str | None = None, missing_psi_policy: str = "retain") -> dict[str, object]:
    if input_mode != "tcga_spliceseq":
        raise ValueError(f"Unsupported input_mode: {input_mode}")
    if not psi_tsv or not sample_annotations_tsv:
        raise ValueError("tcga_spliceseq mode requires --psi_tsv and --sample_annotations_tsv")

    out_root = Path(out_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    study_id_resolved = str(study_id or Path(psi_tsv).stem)
    study_label_resolved = str(study_label or study_id_resolved)
    sample_id_map = _load_sample_id_map(sample_id_map_tsv)
    sample_meta_fields, sample_meta_rows = _load_sample_annotations(sample_annotations_tsv, study_id_resolved, study_label_resolved, sample_id_map)
    sample_ids = [row["sample_id"] for row in sample_meta_rows]
    raw_to_sample = {row["sample_id_raw"]: row["sample_id"] for row in sample_meta_rows}

    psi_fieldnames, psi_rows = _read_tsv(psi_tsv)
    metadata_rows_by_event: dict[str, dict[str, str]] = {}
    if event_metadata_tsv:
        meta_fields, meta_rows = _read_tsv(event_metadata_tsv)
        meta_id_col = "event_id" if "event_id" in meta_fields else meta_fields[0]
        metadata_rows_by_event = {_clean(row.get(meta_id_col)): row for row in meta_rows}

    allowlist = {_normalize_event_type(x) for x in str(event_type_allowlist or "").split(",") if _clean(x)}
    unknown_event_types = 0
    parser_status_counts: dict[str, int] = {}
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
        raise ValueError("No sample columns from sample_annotations_tsv were found in psi_tsv")

    for row in psi_rows:
        row_values = dict(row)
        event_id = _clean(row.get("event_id")) or _clean(row.get("splice_event_id")) or _clean(row.get("Event"))
        if event_id and event_id in metadata_rows_by_event:
            for key, value in metadata_rows_by_event[event_id].items():
                row_values.setdefault(key, value)
        raw_event_key = _event_key(row_values)
        canonical_event_key, parser_status = _canonical_event_key(row_values)
        parser_status_counts[parser_status] = parser_status_counts.get(parser_status, 0) + 1
        event_type = _normalize_event_type(_clean(row_values.get("event_type")) or _clean(row_values.get("splice_type")))
        if not event_type:
            event_type = "unknown"
        event_type_counts[event_type] = event_type_counts.get(event_type, 0) + 1
        if event_type == "unknown":
            unknown_event_types += 1
        if allowlist and event_type not in allowlist:
            dropped_events += 1
            continue
        gene_id = _clean(row_values.get("gene_id")) or _clean(row_values.get("GeneID")) or _clean(row_values.get("gene_symbol")) or _clean(row_values.get("Gene"))
        gene_symbol = _clean(row_values.get("gene_symbol")) or _clean(row_values.get("Gene")) or gene_id
        meta_row = {
            "event_id": raw_event_key,
            "canonical_event_key": canonical_event_key,
            "event_group": _clean(row_values.get("event_group")) or _clean(row_values.get("event_group_id")),
            "event_type": event_type,
            "gene_id": gene_id,
            "gene_symbol": gene_symbol,
            "annotation_status": _clean(row_values.get("annotation_status")) or _clean(row_values.get("annotation")),
            "chrom": _clean(row_values.get("chrom")) or _clean(row_values.get("chr")),
            "start": _clean(row_values.get("start")),
            "end": _clean(row_values.get("end")),
            "strand": _clean(row_values.get("strand")),
            "parser_status": parser_status,
            "source_study_id": study_id_resolved,
            "source_file_name": Path(psi_tsv).name,
            "tool_family": "tcga_spliceseq",
        }
        event_metadata_rows.append(meta_row)
        event_id_map_rows.append({
            "input_event_key": raw_event_key,
            "canonical_event_key": canonical_event_key,
            "gene_id": gene_id,
            "gene_symbol": gene_symbol,
            "event_type": event_type,
            "parser_status": parser_status,
        })
        psi_row = {"event_id": raw_event_key}
        for key in ("event_group", "event_type", "gene_id", "gene_symbol", "annotation_status", "chrom", "start", "end", "strand"):
            psi_row[key] = meta_row[key]
        for raw_sample in sample_columns:
            sample_id = raw_to_sample[raw_sample]
            value = _clean(row_values.get(raw_sample))
            psi_row[sample_id] = value
            total_values += 1
            if _parse_float(value) is None:
                missing_values += 1
                continue
            bundle_source_rows.append({
                "source_dataset": study_id_resolved,
                "sample_id": sample_id,
                "input_event_key": raw_event_key,
                "canonical_event_key": canonical_event_key,
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "event_type": event_type,
                "psi": value,
                "read_support": _clean(row_values.get("read_support")),
                "annotation_status": meta_row["annotation_status"],
                "chrom": meta_row["chrom"],
                "start": meta_row["start"],
                "end": meta_row["end"],
                "strand": meta_row["strand"],
                "tool_family": "tcga_spliceseq",
            })
        psi_matrix_rows.append(psi_row)

    outputs = {
        "psi_matrix_tsv": str(out_root / "psi_matrix.tsv"),
        "sample_metadata_tsv": str(out_root / "sample_metadata.tsv"),
        "event_metadata_tsv": str(out_root / "event_metadata.tsv"),
        "sample_id_map_tsv": str(out_root / "sample_id_map.tsv"),
        "event_id_map_tsv": str(out_root / "event_id_map.tsv"),
        "bundle_source_row_tsv": str(out_root / "bundle_source_row.tsv"),
    }
    _write_tsv(out_root / "psi_matrix.tsv", ["event_id", "event_group", "event_type", "gene_id", "gene_symbol", "annotation_status", "chrom", "start", "end", "strand"] + sample_ids, psi_matrix_rows)
    _write_tsv(out_root / "sample_metadata.tsv", sample_meta_fields, sample_meta_rows)
    _write_tsv(out_root / "event_metadata.tsv", ["event_id", "canonical_event_key", "event_group", "event_type", "gene_id", "gene_symbol", "annotation_status", "chrom", "start", "end", "strand", "parser_status", "source_study_id", "source_file_name", "tool_family"], event_metadata_rows)
    _write_tsv(out_root / "sample_id_map.tsv", ["sample_id_raw", "sample_id", "drop_status"], sample_id_map_rows)
    _write_tsv(out_root / "event_id_map.tsv", ["input_event_key", "canonical_event_key", "gene_id", "gene_symbol", "event_type", "parser_status"], event_id_map_rows)
    _write_tsv(out_root / "bundle_source_row.tsv", ["source_dataset", "sample_id", "input_event_key", "canonical_event_key", "gene_id", "gene_symbol", "event_type", "psi", "read_support", "annotation_status", "chrom", "start", "end", "strand", "tool_family"], bundle_source_rows)

    warnings: list[str] = []
    if event_metadata_rows and (unknown_event_types / float(len(event_metadata_rows))) >= 0.8:
        warnings.append("event_types_overwhelmingly_unknown")
    if event_metadata_rows and (len([r for r in event_metadata_rows if r["parser_status"] == "raw_id_fallback"]) / float(len(event_metadata_rows))) >= 0.8:
        warnings.append("most_events_failed_canonicalization")
    if total_values and (missing_values / float(total_values)) >= 0.9:
        warnings.append("psi_matrix_mostly_missing")

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
        "parser_status_counts": parser_status_counts,
        "dropped_events": dropped_events,
        "bundle_source_rows": len(bundle_source_rows),
        "bundle_source_eligible_events": len({row["canonical_event_key"] for row in bundle_source_rows}),
        "warnings": warnings,
        "outputs": outputs,
        "source_files": {
            "psi_tsv": {"path": str(Path(psi_tsv)), "sha256": sha256_file(Path(psi_tsv))},
            "sample_annotations_tsv": {"path": str(Path(sample_annotations_tsv)), "sha256": sha256_file(Path(sample_annotations_tsv))},
        },
    }
    if event_metadata_tsv:
        summary["source_files"]["event_metadata_tsv"] = {"path": str(Path(event_metadata_tsv)), "sha256": sha256_file(Path(event_metadata_tsv))}
    if sample_id_map_tsv:
        summary["source_files"]["sample_id_map_tsv"] = {"path": str(Path(sample_id_map_tsv)), "sha256": sha256_file(Path(sample_id_map_tsv))}

    (out_root / "prepare_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary
