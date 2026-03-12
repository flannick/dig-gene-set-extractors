from __future__ import annotations

import csv
import hashlib
import json
from importlib.resources import files
from pathlib import Path
import re
import sys
from typing import Any

from geneset_extractors.extractors.rnaseq.deg_scoring import sanitize_name_component
from geneset_extractors.hashing import sha256_file


SITE_TOKEN_RE = re.compile(r"([A-Z])\s*([0-9]+)")


def _clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _read_tsv(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    p = Path(path)
    with p.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"TSV has no header: {p}")
        fieldnames = [str(x) for x in reader.fieldnames]
        rows = [{k: str(v) for k, v in row.items() if k is not None} for row in reader]
    return fieldnames, rows


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _parse_float(value: object) -> float | None:
    text = _clean(value)
    if not text or text.lower() in {"na", "nan", "null", "none"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _load_profiles() -> dict[str, dict[str, Any]]:
    path = files("geneset_extractors.resources").joinpath("ptm_public_profiles.json")
    return json.loads(path.read_text(encoding="utf-8"))


def _find_first_existing(candidates: list[Path]) -> Path | None:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def discover_cdap_inputs(
    *,
    input_mode: str,
    ptm_report_tsv: str | None,
    protein_report_tsv: str | None,
    sample_design_tsv: str | None,
    pdc_manifest_tsv: str | None,
    source_dir: str | None,
) -> dict[str, Path | None]:
    if input_mode == "cdap_files":
        if not ptm_report_tsv:
            raise ValueError("input_mode=cdap_files requires --ptm_report_tsv")
        return {
            "ptm_report_tsv": Path(ptm_report_tsv),
            "protein_report_tsv": Path(protein_report_tsv) if protein_report_tsv else None,
            "sample_design_tsv": Path(sample_design_tsv) if sample_design_tsv else None,
            "pdc_manifest_tsv": Path(pdc_manifest_tsv) if pdc_manifest_tsv else None,
        }

    if input_mode != "pdc_manifest":
        raise ValueError(f"Unsupported input_mode: {input_mode}")
    if not pdc_manifest_tsv or not source_dir:
        raise ValueError("input_mode=pdc_manifest requires --pdc_manifest_tsv and --source_dir")

    source_root = Path(source_dir)
    manifest_path = Path(pdc_manifest_tsv)
    _fieldnames, manifest_rows = _read_tsv(manifest_path)
    filenames: list[str] = []
    for row in manifest_rows:
        for key in ("file_name", "filename", "name", "path", "file_path"):
            value = _clean(row.get(key))
            if value:
                filenames.append(Path(value).name)
                break

    named_candidates = [source_root / name for name in filenames]
    ptm_path = _find_first_existing([p for p in named_candidates if "phosphosite" in p.name.lower()])
    protein_path = _find_first_existing(
        [p for p in named_candidates if "phosphosite" not in p.name.lower() and ("proteome" in p.name.lower() or "protein" in p.name.lower())]
    )
    sample_design_path = _find_first_existing([p for p in named_candidates if p.name.lower().endswith(".sample.txt") or "sample" in p.name.lower()])

    if ptm_path is None:
        ptm_path = _find_first_existing(sorted(source_root.glob("*phosphosite*.tsv")))
    if protein_path is None:
        protein_path = _find_first_existing(sorted(source_root.glob("*proteome*.tsv")) + sorted(source_root.glob("*protein*.tsv")))
    if sample_design_path is None:
        sample_design_path = _find_first_existing(sorted(source_root.glob("*.sample.txt")) + sorted(source_root.glob("*sample*.txt")))

    if ptm_path is None:
        raise FileNotFoundError("Could not discover a phosphosite report in source_dir using the PDC manifest and filename patterns.")

    return {
        "ptm_report_tsv": ptm_path,
        "protein_report_tsv": protein_path,
        "sample_design_tsv": sample_design_path,
        "pdc_manifest_tsv": manifest_path,
    }


def _detect_value_columns(
    fieldnames: list[str],
    rows: list[dict[str, str]],
    *,
    required_columns: list[str],
    optional_columns: list[str],
    drop_columns_regex: list[str],
) -> list[str]:
    patterns = [re.compile(p) for p in drop_columns_regex]
    excluded = set(required_columns) | set(optional_columns)
    sample_columns: list[str] = []
    preview_rows = rows[: min(25, len(rows))]
    for name in fieldnames:
        if name in excluded:
            continue
        if any(p.search(name) for p in patterns):
            continue
        if any(_parse_float(row.get(name)) is not None for row in preview_rows):
            sample_columns.append(name)
    if not sample_columns:
        raise ValueError("Could not detect any numeric sample columns in the public PTM report.")
    return sample_columns


def _parse_site_label(raw_site_label: str, ptm_type: str) -> dict[str, object]:
    text = _clean(raw_site_label)
    accession = ""
    suffix = text
    match = re.match(r"^([A-Za-z0-9_.-]+)\s*[:|]\s*(.+)$", text)
    if match:
        accession = match.group(1)
        suffix = match.group(2)
    tokens = [(residue, int(position)) for residue, position in SITE_TOKEN_RE.findall(suffix)]
    if tokens:
        tokens = sorted(set(tokens), key=lambda item: (item[1], item[0]))
    if accession and len(tokens) == 1:
        residue, position = tokens[0]
        site_id = f"{accession}|{residue}|{position}|{ptm_type}"
        return {
            "site_id": site_id,
            "site_group_id": site_id,
            "protein_accession": accession,
            "residue": residue,
            "position": str(position),
            "site_parser_status": "single_site",
            "n_sites_in_group": 1,
        }
    if accession and len(tokens) > 1:
        group_label = ";".join(f"{residue}{position}" for residue, position in tokens)
        return {
            "site_id": "",
            "site_group_id": f"{accession}|{group_label}|{ptm_type}",
            "protein_accession": accession,
            "residue": "",
            "position": "",
            "site_parser_status": "site_group",
            "n_sites_in_group": len(tokens),
        }
    fallback_token = sanitize_name_component(text) or "unparsed_site"
    return {
        "site_id": "",
        "site_group_id": f"raw::{fallback_token}|{ptm_type}",
        "protein_accession": accession,
        "residue": "",
        "position": "",
        "site_parser_status": "unparsed_fallback",
        "n_sites_in_group": max(1, len(tokens)),
    }


def _assay_type_qc(
    site_id_map_rows: list[dict[str, object]],
    *,
    ptm_type: str,
    assay_type_policy: str,
    min_phospho_like_fraction: float,
    max_k_fraction: float,
) -> dict[str, object]:
    residue_counts: dict[str, int] = {}
    for row in site_id_map_rows:
        raw_site_label = _clean(row.get("raw_site_label"))
        for residue, _position in SITE_TOKEN_RE.findall(raw_site_label):
            residue_counts[residue] = int(residue_counts.get(residue, 0)) + 1

    total_tokens = int(sum(residue_counts.values()))
    phospho_like_count = int(residue_counts.get("S", 0) + residue_counts.get("T", 0) + residue_counts.get("Y", 0))
    lysine_count = int(residue_counts.get("K", 0))
    phospho_like_fraction = float(phospho_like_count) / float(total_tokens) if total_tokens else 0.0
    lysine_fraction = float(lysine_count) / float(total_tokens) if total_tokens else 0.0

    if total_tokens == 0:
        dominant_residue_family = "unknown"
    elif lysine_count > phospho_like_count:
        dominant_residue_family = "lysine_dominant"
    elif phospho_like_count > 0:
        dominant_residue_family = "phospho_like"
    else:
        dominant_residue_family = "other"

    status = "not_applicable"
    reason = ""
    if ptm_type == "phospho":
        if assay_type_policy == "off":
            status = "unchecked"
            reason = "Assay-type QC disabled."
        elif total_tokens == 0:
            status = "warn"
            reason = "No parsable residue-position tokens were found; phospho assay typing is uncertain."
        elif phospho_like_fraction < float(min_phospho_like_fraction) or lysine_fraction > float(max_k_fraction):
            status = "warn"
            reason = (
                f"Public PTM report looks weakly phospho-like: phospho_like_fraction={phospho_like_fraction:.3f}, "
                f"fraction_k={lysine_fraction:.3f}, dominant_residue_family={dominant_residue_family}."
            )
        else:
            status = "pass"
            reason = "Residue composition is consistent with a phospho-like public report."

    result = {
        "status": status,
        "reason": reason,
        "dominant_residue_family": dominant_residue_family,
        "residue_counts": {key: int(residue_counts[key]) for key in sorted(residue_counts)},
        "n_residue_tokens_parsed": total_tokens,
        "fraction_s": float(residue_counts.get("S", 0)) / float(total_tokens) if total_tokens else 0.0,
        "fraction_t": float(residue_counts.get("T", 0)) / float(total_tokens) if total_tokens else 0.0,
        "fraction_y": float(residue_counts.get("Y", 0)) / float(total_tokens) if total_tokens else 0.0,
        "fraction_k": lysine_fraction,
        "phospho_like_fraction": phospho_like_fraction,
        "thresholds": {
            "min_phospho_like_fraction": float(min_phospho_like_fraction),
            "max_k_fraction": float(max_k_fraction),
        },
    }

    if ptm_type == "phospho" and status == "warn":
        if assay_type_policy == "fail":
            raise ValueError(
                f"PTM public assay-type QC failed: {reason} "
                "Use --assay_type_policy warn to continue for exploratory runs."
            )
        print(f"warning: {reason}", file=sys.stderr)

    return result


def read_cdap_sample_design(path: str | Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    fieldnames, rows = _read_tsv(path)
    if not fieldnames:
        return {}
    analytical_keys = ["AnalyticalSample", "analytical_sample", "Sample", "sample"]
    reporter_keys = ["Reporter", "reporter_ion", "ReporterIon", "Reporter Ion"]
    label_keys = ["LabelReagent", "label_reagent", "Label Reagent"]

    out: dict[str, dict[str, str]] = {}
    for row in rows:
        analytical_sample = ""
        reporter_ion = ""
        label_reagent = ""
        for key in analytical_keys:
            analytical_sample = _clean(row.get(key))
            if analytical_sample:
                break
        for key in reporter_keys:
            reporter_ion = _clean(row.get(key))
            if reporter_ion:
                break
        for key in label_keys:
            label_reagent = _clean(row.get(key))
            if label_reagent:
                break
        if not analytical_sample:
            continue
        out[analytical_sample] = {
            "analytical_sample": analytical_sample,
            "reporter_ion": reporter_ion,
            "label_reagent": label_reagent,
        }
    return out


def _annotations_lookup(path: str | Path | None) -> tuple[dict[str, dict[str, str]], list[str]]:
    if path is None:
        return {}, []
    fieldnames, rows = _read_tsv(path)
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        key = ""
        for candidate in ("sample_id_raw", "sample_id", "analytical_sample", "sample", "Sample"):
            key = _clean(row.get(candidate))
            if key:
                break
        if not key:
            continue
        lookup[key] = row
    extra_columns = [
        name
        for name in fieldnames
        if name
        not in {
            "sample_id_raw",
            "sample_id",
            "analytical_sample",
            "sample",
            "Sample",
        }
    ]
    return lookup, extra_columns


def _derive_condition(sample_type: str) -> str:
    value = _clean(sample_type).lower()
    if value in {"normal", "control", "adjacent_normal"}:
        return "control"
    if value in {"tumor", "case", "treated"}:
        return "case"
    return ""


def normalize_sample_ids(
    *,
    sample_columns: list[str],
    sample_design: dict[str, dict[str, str]],
    sample_annotations: dict[str, dict[str, str]],
    study_id: str,
    study_label: str,
) -> tuple[list[str], dict[str, str], list[dict[str, object]], list[dict[str, object]], list[str]]:
    mapping: dict[str, str] = {}
    metadata_rows: list[dict[str, object]] = []
    id_map_rows: list[dict[str, object]] = []
    kept_sample_ids: list[str] = []
    used_ids: set[str] = set()
    extra_annotation_columns = sorted(
        {
            key
            for row in sample_annotations.values()
            for key in row
            if key
            not in {
                "sample_id_raw",
                "sample_id",
                "analytical_sample",
                "sample",
                "Sample",
                "condition",
                "group",
                "study_id",
                "study_label",
                "case_submitter_id",
                "sample_submitter_id",
                "aliquot_submitter_id",
                "sample_type",
                "tissue_type",
            }
        }
    )

    for raw_name in sample_columns:
        normalized = sanitize_name_component(raw_name) or "sample"
        candidate = normalized
        idx = 2
        while candidate in used_ids:
            candidate = f"{normalized}_{idx}"
            idx += 1
        used_ids.add(candidate)
        mapping[raw_name] = candidate

        design = sample_design.get(raw_name, {})
        anno = sample_annotations.get(raw_name, sample_annotations.get(candidate, {}))
        sample_type = _clean(anno.get("sample_type"))
        condition = _clean(anno.get("condition")) or _derive_condition(sample_type)
        group = _clean(anno.get("group")) or study_id
        drop_status = "dropped_pool_reference" if "pool" in raw_name.lower() or "pool" in _clean(design.get("analytical_sample")).lower() else "kept"

        id_map_rows.append(
            {
                "sample_id_raw": raw_name,
                "sample_id": candidate,
                "source_report": "cdap_report",
                "analytical_sample": _clean(design.get("analytical_sample")) or raw_name,
                "reporter_ion": _clean(design.get("reporter_ion")),
                "label_reagent": _clean(design.get("label_reagent")),
                "drop_status": drop_status,
            }
        )

        if drop_status != "kept":
            continue

        row: dict[str, object] = {
            "sample_id": candidate,
            "sample_id_raw": raw_name,
            "condition": condition,
            "group": group,
            "study_id": _clean(anno.get("study_id")) or study_id,
            "study_label": _clean(anno.get("study_label")) or study_label,
            "analytical_sample": _clean(design.get("analytical_sample")) or raw_name,
            "reporter_ion": _clean(design.get("reporter_ion")),
            "label_reagent": _clean(design.get("label_reagent")),
            "case_submitter_id": _clean(anno.get("case_submitter_id")),
            "sample_submitter_id": _clean(anno.get("sample_submitter_id")),
            "aliquot_submitter_id": _clean(anno.get("aliquot_submitter_id")),
            "sample_type": sample_type,
            "tissue_type": _clean(anno.get("tissue_type")),
        }
        for key in extra_annotation_columns:
            row[key] = _clean(anno.get(key))
        metadata_rows.append(row)
        kept_sample_ids.append(candidate)

    return kept_sample_ids, mapping, metadata_rows, id_map_rows, extra_annotation_columns


def read_cdap_phosphosite_report(
    path: str | Path,
    *,
    profile: dict[str, Any],
    ptm_type: str,
    study_id: str,
    sample_id_map: dict[str, str],
    kept_sample_ids: list[str],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str], dict[str, int]]:
    fieldnames, rows = _read_tsv(path)
    sample_columns_raw = _detect_value_columns(
        fieldnames,
        rows,
        required_columns=list(profile.get("required_columns", [])),
        optional_columns=list(profile.get("optional_columns", [])),
        drop_columns_regex=list(profile.get("drop_columns_regex", [])),
    )
    kept_raw_columns = [raw for raw in sample_columns_raw if raw in sample_id_map and sample_id_map[raw] in kept_sample_ids]

    output_rows: list[dict[str, object]] = []
    site_id_map_rows: list[dict[str, object]] = []
    seen_site_labels: set[str] = set()
    summary = {
        "rows_read": len(rows),
        "rows_emitted": 0,
        "rows_dropped": 0,
        "single_site_rows": 0,
        "site_group_rows": 0,
        "unparsed_rows": 0,
    }
    for row in rows:
        raw_site_label = _clean(row.get("Phosphosite"))
        raw_peptide = _clean(row.get("Peptide"))
        gene = _clean(row.get("Gene"))
        if not raw_site_label:
            summary["rows_dropped"] += 1
            continue
        parsed = _parse_site_label(raw_site_label, ptm_type)
        status = str(parsed["site_parser_status"])
        if status == "single_site":
            summary["single_site_rows"] += 1
        elif status == "site_group":
            summary["site_group_rows"] += 1
        else:
            summary["unparsed_rows"] += 1

        record: dict[str, object] = {
            "site_id": parsed["site_id"],
            "site_group_id": parsed["site_group_id"],
            "gene_id": gene,
            "gene_symbol": gene,
            "protein_accession": parsed["protein_accession"],
            "residue": parsed["residue"],
            "position": parsed["position"],
            "raw_site_label": raw_site_label,
            "raw_peptide": raw_peptide,
            "site_parser_status": status,
            "source_study_id": study_id,
            "source_file_name": Path(path).name,
        }
        for raw_col in kept_raw_columns:
            record[sample_id_map[raw_col]] = _clean(row.get(raw_col))
        output_rows.append(record)
        summary["rows_emitted"] += 1

        if raw_site_label not in seen_site_labels:
            seen_site_labels.add(raw_site_label)
            site_id_map_rows.append(
                {
                    "raw_site_label": raw_site_label,
                    "site_id": parsed["site_id"],
                    "site_group_id": parsed["site_group_id"],
                    "protein_accession": parsed["protein_accession"],
                    "residue": parsed["residue"],
                    "position": parsed["position"],
                    "gene_symbol": gene,
                    "site_parser_status": status,
                    "n_sites_in_group": parsed["n_sites_in_group"],
                }
            )

    sample_columns = [sample_id_map[raw] for raw in kept_raw_columns]
    return output_rows, site_id_map_rows, sample_columns, summary


def read_cdap_proteome_report(
    path: str | Path | None,
    *,
    profile: dict[str, Any],
    study_id: str,
    sample_id_map: dict[str, str],
    kept_sample_ids: list[str],
) -> tuple[list[dict[str, object]], list[str], dict[str, int]]:
    if path is None:
        return [], [], {"rows_read": 0, "rows_emitted": 0, "rows_dropped": 0}
    fieldnames, rows = _read_tsv(path)
    sample_columns_raw = _detect_value_columns(
        fieldnames,
        rows,
        required_columns=list(profile.get("required_columns", [])),
        optional_columns=list(profile.get("optional_columns", [])),
        drop_columns_regex=list(profile.get("drop_columns_regex", [])),
    )
    kept_raw_columns = [raw for raw in sample_columns_raw if raw in sample_id_map and sample_id_map[raw] in kept_sample_ids]
    output_rows: list[dict[str, object]] = []
    summary = {"rows_read": len(rows), "rows_emitted": 0, "rows_dropped": 0}

    for row in rows:
        gene = _clean(row.get("Gene"))
        if gene.lower() in {"mean", "median", "stddev", "standard deviation"}:
            summary["rows_dropped"] += 1
            continue
        values = [_parse_float(row.get(raw_col)) for raw_col in kept_raw_columns]
        if not any(value is not None for value in values):
            summary["rows_dropped"] += 1
            continue
        protein_accession = _clean(row.get("Protein"))
        record: dict[str, object] = {
            "protein_id": gene or protein_accession or "protein",
            "gene_id": gene,
            "gene_symbol": gene,
            "protein_accession": protein_accession,
            "source_study_id": study_id,
            "source_file_name": Path(path).name,
        }
        for raw_col in kept_raw_columns:
            record[sample_id_map[raw_col]] = _clean(row.get(raw_col))
        output_rows.append(record)
        summary["rows_emitted"] += 1
    sample_columns = [sample_id_map[raw] for raw in kept_raw_columns]
    return output_rows, sample_columns, summary


def merge_sample_annotations(
    *,
    kept_sample_ids: list[str],
    sample_metadata_rows: list[dict[str, object]],
    extra_annotation_columns: list[str],
) -> tuple[list[dict[str, object]], int]:
    ordered = {str(row["sample_id"]): row for row in sample_metadata_rows}
    rows = [ordered[sample_id] for sample_id in kept_sample_ids if sample_id in ordered]
    n_with_condition = sum(1 for row in rows if _clean(row.get("condition")))
    return rows, n_with_condition


def _file_md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_qc_outputs(
    *,
    out_dir: Path,
    ptm_rows: list[dict[str, object]],
    protein_rows: list[dict[str, object]],
    sample_metadata_rows: list[dict[str, object]],
    sample_id_map_rows: list[dict[str, object]],
    site_id_map_rows: list[dict[str, object]],
    source_paths: dict[str, Path | None],
    parser_profile: str,
    study_id: str,
    study_label: str,
    ptm_summary: dict[str, int],
    protein_summary: dict[str, int],
    n_samples_with_condition: int,
    assay_qc: dict[str, object],
) -> dict[str, object]:
    ptm_path = out_dir / "ptm_matrix.tsv"
    sample_meta_path = out_dir / "sample_metadata.tsv"
    protein_path = out_dir / "protein_matrix.tsv"
    sample_id_map_path = out_dir / "sample_id_map.tsv"
    site_id_map_path = out_dir / "site_id_map.tsv"
    prepare_summary_path = out_dir / "prepare_summary.json"
    bundle_source_row_path = out_dir / "bundle_source_row.tsv"

    ptm_fieldnames = [
        "site_id",
        "site_group_id",
        "gene_id",
        "gene_symbol",
        "protein_accession",
        "residue",
        "position",
        "raw_site_label",
        "raw_peptide",
        "site_parser_status",
        "source_study_id",
        "source_file_name",
    ] + [
        key for key in ptm_rows[0].keys() if key not in {
            "site_id",
            "site_group_id",
            "gene_id",
            "gene_symbol",
            "protein_accession",
            "residue",
            "position",
            "raw_site_label",
            "raw_peptide",
            "site_parser_status",
            "source_study_id",
            "source_file_name",
        }
    ] if ptm_rows else []
    _write_tsv(ptm_path, ptm_fieldnames, ptm_rows)

    sample_meta_fieldnames = [
        "sample_id",
        "sample_id_raw",
        "condition",
        "group",
        "study_id",
        "study_label",
        "analytical_sample",
        "reporter_ion",
        "label_reagent",
        "case_submitter_id",
        "sample_submitter_id",
        "aliquot_submitter_id",
        "sample_type",
        "tissue_type",
    ] + [
        key
        for key in (sample_metadata_rows[0].keys() if sample_metadata_rows else [])
        if key not in {
            "sample_id",
            "sample_id_raw",
            "condition",
            "group",
            "study_id",
            "study_label",
            "analytical_sample",
            "reporter_ion",
            "label_reagent",
            "case_submitter_id",
            "sample_submitter_id",
            "aliquot_submitter_id",
            "sample_type",
            "tissue_type",
        }
    ]
    _write_tsv(sample_meta_path, sample_meta_fieldnames, sample_metadata_rows)

    if protein_rows:
        protein_fieldnames = [
            "protein_id",
            "gene_id",
            "gene_symbol",
            "protein_accession",
            "source_study_id",
            "source_file_name",
        ] + [
            key for key in protein_rows[0].keys() if key not in {
                "protein_id",
                "gene_id",
                "gene_symbol",
                "protein_accession",
                "source_study_id",
                "source_file_name",
            }
        ]
        _write_tsv(protein_path, protein_fieldnames, protein_rows)

    _write_tsv(
        sample_id_map_path,
        ["sample_id_raw", "sample_id", "source_report", "analytical_sample", "reporter_ion", "label_reagent", "drop_status"],
        sample_id_map_rows,
    )
    _write_tsv(
        site_id_map_path,
        ["raw_site_label", "site_id", "site_group_id", "protein_accession", "residue", "position", "gene_symbol", "site_parser_status", "n_sites_in_group"],
        site_id_map_rows,
    )
    _write_tsv(bundle_source_row_path, ["path", "source_dataset"], [{"path": str(ptm_path), "source_dataset": study_id or study_label or ptm_path.stem}])

    source_hashes = {}
    for label, path in source_paths.items():
        if path is None:
            continue
        source_hashes[label] = {
            "path": str(path),
            "sha256": sha256_file(path),
            "md5": _file_md5(path),
        }

    prepare_summary = {
        "workflow": "ptm_prepare_public",
        "parser_profile": parser_profile,
        "study_id": study_id,
        "study_label": study_label,
        "source_files": source_hashes,
        "ptm_rows_read": ptm_summary["rows_read"],
        "ptm_rows_emitted": ptm_summary["rows_emitted"],
        "ptm_rows_dropped": ptm_summary["rows_dropped"],
        "single_resolved_sites": ptm_summary["single_site_rows"],
        "grouped_sites": ptm_summary["site_group_rows"],
        "unparsed_fallback_sites": ptm_summary["unparsed_rows"],
        "protein_rows_read": protein_summary["rows_read"],
        "protein_rows_emitted": protein_summary["rows_emitted"],
        "protein_rows_dropped": protein_summary["rows_dropped"],
        "n_samples": len(sample_metadata_rows),
        "n_samples_with_condition_labels": n_samples_with_condition,
        "assay_type_qc": assay_qc,
        "outputs": {
            "ptm_matrix_tsv": str(ptm_path),
            "sample_metadata_tsv": str(sample_meta_path),
            "protein_matrix_tsv": str(protein_path) if protein_rows else None,
            "sample_id_map_tsv": str(sample_id_map_path),
            "site_id_map_tsv": str(site_id_map_path),
            "bundle_source_row_tsv": str(bundle_source_row_path),
        },
    }
    prepare_summary_path.write_text(json.dumps(prepare_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return prepare_summary


def run_public_prepare(
    *,
    input_mode: str,
    ptm_report_tsv: str | None,
    protein_report_tsv: str | None,
    sample_design_tsv: str | None,
    sample_annotations_tsv: str | None,
    pdc_manifest_tsv: str | None,
    source_dir: str | None,
    out_dir: str | Path,
    organism: str,
    ptm_type: str,
    study_id: str | None,
    study_label: str | None,
    assay_type_policy: str = "warn",
    min_phospho_like_fraction: float = 0.6,
    max_k_fraction: float = 0.25,
) -> dict[str, object]:
    profiles = _load_profiles()
    phospho_profile = profiles["cdap_phosphosite_v1"]
    proteome_profile = profiles["cdap_proteome_v1"]

    discovered = discover_cdap_inputs(
        input_mode=input_mode,
        ptm_report_tsv=ptm_report_tsv,
        protein_report_tsv=protein_report_tsv,
        sample_design_tsv=sample_design_tsv,
        pdc_manifest_tsv=pdc_manifest_tsv,
        source_dir=source_dir,
    )

    inferred_study_id = _clean(study_id) or Path(discovered["ptm_report_tsv"]).stem
    inferred_study_label = _clean(study_label) or inferred_study_id

    ptm_fieldnames, ptm_rows_raw = _read_tsv(discovered["ptm_report_tsv"])
    ptm_sample_columns_raw = _detect_value_columns(
        ptm_fieldnames,
        ptm_rows_raw,
        required_columns=list(phospho_profile.get("required_columns", [])),
        optional_columns=list(phospho_profile.get("optional_columns", [])),
        drop_columns_regex=list(phospho_profile.get("drop_columns_regex", [])),
    )

    sample_design = read_cdap_sample_design(discovered.get("sample_design_tsv"))
    sample_annotations, _extra_columns = _annotations_lookup(sample_annotations_tsv)
    kept_sample_ids, sample_id_map, sample_metadata_rows, sample_id_map_rows, _extra_annotation_columns = normalize_sample_ids(
        sample_columns=ptm_sample_columns_raw,
        sample_design=sample_design,
        sample_annotations=sample_annotations,
        study_id=inferred_study_id,
        study_label=inferred_study_label,
    )
    sample_metadata_rows, n_samples_with_condition = merge_sample_annotations(
        kept_sample_ids=kept_sample_ids,
        sample_metadata_rows=sample_metadata_rows,
        extra_annotation_columns=[],
    )

    ptm_rows, site_id_map_rows, _ptm_sample_columns, ptm_summary = read_cdap_phosphosite_report(
        discovered["ptm_report_tsv"],
        profile=phospho_profile,
        ptm_type=ptm_type,
        study_id=inferred_study_id,
        sample_id_map=sample_id_map,
        kept_sample_ids=kept_sample_ids,
    )
    protein_rows, _protein_sample_columns, protein_summary = read_cdap_proteome_report(
        discovered.get("protein_report_tsv"),
        profile=proteome_profile,
        study_id=inferred_study_id,
        sample_id_map=sample_id_map,
        kept_sample_ids=kept_sample_ids,
    )
    assay_qc = _assay_type_qc(
        site_id_map_rows,
        ptm_type=ptm_type,
        assay_type_policy=assay_type_policy,
        min_phospho_like_fraction=min_phospho_like_fraction,
        max_k_fraction=max_k_fraction,
    )

    out_root = Path(out_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    prepare_summary = write_qc_outputs(
        out_dir=out_root,
        ptm_rows=ptm_rows,
        protein_rows=protein_rows,
        sample_metadata_rows=sample_metadata_rows,
        sample_id_map_rows=sample_id_map_rows,
        site_id_map_rows=site_id_map_rows,
        source_paths=discovered,
        parser_profile="cdap_phosphosite_v1",
        study_id=inferred_study_id,
        study_label=inferred_study_label,
        ptm_summary=ptm_summary,
        protein_summary=protein_summary,
        n_samples_with_condition=n_samples_with_condition,
        assay_qc=assay_qc,
    )
    prepare_summary["out_dir"] = str(out_root)
    return prepare_summary
