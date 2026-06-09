from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path
from typing import Any


def open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_maybe_gzip(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def normalize_intervention(raw: str) -> str:
    value = str(raw or "").strip().lower()
    if not value:
        return ""
    if value == "control":
        return "control"
    if "training" in value:
        return "training"
    return value


def normalize_sex(raw: str) -> str:
    value = str(raw or "").strip().lower()
    if not value:
        return ""
    if value in {"male", "m"}:
        return "M"
    if value in {"female", "f"}:
        return "F"
    return value


def normalize_sex_label(raw: str) -> str:
    value = str(raw or "").strip().lower()
    if value in {"male", "m"}:
        return "male"
    if value in {"female", "f"}:
        return "female"
    return value


def strip_ensembl_version(value: str) -> str:
    text = str(value or "").strip()
    return text.split(".", 1)[0] if text else ""


def row_variance(values: list[float]) -> float:
    if len(values) <= 1:
        return 0.0
    mean_value = sum(values) / len(values)
    return sum((value - mean_value) ** 2 for value in values) / (len(values) - 1)


def parse_timepoint_label(raw: str) -> str:
    value = str(raw or "").strip().lower()
    match = re.search(r"(\d+)\s*weeks?", value)
    if match:
        return f"{match.group(1)}w"
    return ""


def parse_counts(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with open_maybe_gzip(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = [{str(k): str(v) for k, v in row.items()} for row in reader]
    if len(fieldnames) < 5:
        raise ValueError(f"Unexpected counts schema in {path}")
    return fieldnames[4:], rows


def prepare_tissue_inputs(
    *,
    counts_tsv: Path,
    transcript_metadata_tsv: Path,
    phenotype_metadata_tsv: Path,
    feature_to_gene_tsv: Path,
    rat_to_human_tsv: Path,
    tissue_label: str,
    transcript_tissue_label: str,
) -> dict[str, Any]:
    count_sample_ids, count_rows = parse_counts(counts_tsv)
    count_sample_id_set = set(count_sample_ids)
    transcript_rows = read_tsv(transcript_metadata_tsv)
    phenotype_rows = read_tsv(phenotype_metadata_tsv)
    feature_rows = read_tsv(feature_to_gene_tsv)
    ortholog_rows = read_tsv(rat_to_human_tsv)

    transcript_by_vial: dict[str, dict[str, str]] = {}
    for row in transcript_rows:
        tissue_value = str(row.get("Tissue", "")).strip()
        vial = str(row.get("viallabel", "")).strip() or str(row.get("vial_label", "")).strip()
        if not vial or vial not in count_sample_id_set:
            continue
        if tissue_value != transcript_tissue_label:
            continue
        transcript_by_vial[vial] = row

    phenotype_by_sample_id = {
        str(row.get("sample_id", "")).strip(): row
        for row in phenotype_rows
        if str(row.get("sample_id", "")).strip()
    }
    phenotype_by_pid_bid: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in phenotype_rows:
        key = (str(row.get("pid", "")).strip(), str(row.get("bid", "")).strip())
        if key == ("", ""):
            continue
        phenotype_by_pid_bid.setdefault(key, []).append(row)

    prepared_meta: list[dict[str, str]] = []
    retained_sample_ids: list[str] = []
    skipped_for_join = 0
    for sample_id in count_sample_ids:
        transcript_row = transcript_by_vial.get(sample_id)
        if transcript_row is None:
            skipped_for_join += 1
            continue
        key = (str(transcript_row.get("PID", "")).strip(), str(transcript_row.get("BID", "")).strip())
        phenotype_row = phenotype_by_sample_id.get(sample_id)
        if phenotype_row is None:
            candidate_rows = phenotype_by_pid_bid.get(key, [])
            phenotype_row = next(
                (
                    row
                    for row in candidate_rows
                    if str(row.get("tissue_description", "")).strip().lower() == tissue_label.strip().lower()
                ),
                None,
            )
        if phenotype_row is None:
            skipped_for_join += 1
            continue
        intervention = normalize_intervention(phenotype_row.get("key.intervention", ""))
        sex = normalize_sex(phenotype_row.get("sex", "") or phenotype_row.get("registration.sex", ""))
        sex_label = normalize_sex_label(phenotype_row.get("sex", "") or phenotype_row.get("registration.sex", ""))
        tissue_code_no = str(phenotype_row.get("tissue_code_no", "")).strip()
        timepoint_label = parse_timepoint_label(phenotype_row.get("key.sacrificetime", ""))
        if intervention not in {"control", "training"} or sex not in {"M", "F"} or not tissue_code_no or not timepoint_label:
            skipped_for_join += 1
            continue
        retained_sample_ids.append(sample_id)
        prepared_meta.append(
            {
                "sample_id": sample_id,
                "pid": key[0],
                "bid": key[1],
                "sex": sex,
                "sex_label": sex_label,
                "intervention": intervention,
                "tissue": tissue_label,
                "transcript_tissue": str(transcript_row.get("Tissue", "")).strip(),
                "tissue_code_no": tissue_code_no,
                "timepoint_label": timepoint_label,
            }
        )

    counts_by_group: dict[tuple[str, str], int] = {}
    for row in prepared_meta:
        key = (row["intervention"], row["sex"])
        counts_by_group[key] = counts_by_group.get(key, 0) + 1

    feature_to_gene = {
        strip_ensembl_version(row.get("feature_ID", "")): row
        for row in feature_rows
        if strip_ensembl_version(row.get("feature_ID", ""))
    }
    ortholog_by_rat_ensembl = {
        strip_ensembl_version(row.get("RAT_ENSEMBL_ID", "")): row
        for row in ortholog_rows
        if strip_ensembl_version(row.get("RAT_ENSEMBL_ID", ""))
    }
    ortholog_by_rat_symbol = {
        str(row.get("RAT_SYMBOL", "")).strip(): row
        for row in ortholog_rows
        if str(row.get("RAT_SYMBOL", "")).strip()
    }

    sample_fieldnames = retained_sample_ids
    mapped_rows: dict[str, tuple[float, dict[str, str]]] = {}
    dropped_missing_mapping = 0
    for row in count_rows:
        feature_id = strip_ensembl_version(row.get("feature_ID", ""))
        feature_map = feature_to_gene.get(feature_id)
        if feature_map is None:
            dropped_missing_mapping += 1
            continue
        rat_symbol = str(feature_map.get("gene_symbol", "")).strip()
        rat_ensembl = strip_ensembl_version(feature_map.get("ensembl_gene", "")) or feature_id
        ortholog_row = ortholog_by_rat_ensembl.get(rat_ensembl)
        if ortholog_row is None and rat_symbol:
            ortholog_row = ortholog_by_rat_symbol.get(rat_symbol)
        if ortholog_row is None:
            dropped_missing_mapping += 1
            continue
        human_symbol = str(ortholog_row.get("HUMAN_ORTHOLOG_SYMBOL", "")).strip()
        human_ensembl = strip_ensembl_version(ortholog_row.get("HUMAN_ORTHOLOG_ENSEMBL_ID", ""))
        if not human_symbol:
            dropped_missing_mapping += 1
            continue
        gene_id = human_ensembl or human_symbol
        out_row = {"gene_id": gene_id, "gene_symbol": human_symbol}
        values: list[float] = []
        for sample_id in sample_fieldnames:
            value_text = str(row.get(sample_id, "0")).strip() or "0"
            out_row[sample_id] = value_text
            values.append(float(value_text))
        variance = row_variance(values)
        existing = mapped_rows.get(human_symbol)
        if existing is None or variance >= existing[0]:
            mapped_rows[human_symbol] = (variance, out_row)

    counts_rows_out = [item[1] for _, item in sorted(mapped_rows.items(), key=lambda kv: kv[0])]
    summary = {
        "tissue_label": tissue_label,
        "transcript_tissue_label": transcript_tissue_label,
        "counts_tsv": str(counts_tsv.resolve()),
        "transcript_metadata_tsv": str(transcript_metadata_tsv.resolve()),
        "phenotype_metadata_tsv": str(phenotype_metadata_tsv.resolve()),
        "feature_to_gene_tsv": str(feature_to_gene_tsv.resolve()),
        "rat_to_human_tsv": str(rat_to_human_tsv.resolve()),
        "n_input_samples": len(count_sample_ids),
        "n_retained_samples": len(prepared_meta),
        "n_input_features": len(count_rows),
        "n_genes_retained": len(counts_rows_out),
        "n_features_dropped_missing_mapping": dropped_missing_mapping,
        "n_samples_skipped_for_join": skipped_for_join,
        "group_counts": {
            f"{intervention}_{sex}": count
            for (intervention, sex), count in sorted(counts_by_group.items())
        },
    }
    return {
        "sample_metadata_rows": prepared_meta,
        "counts_rows": counts_rows_out,
        "counts_fieldnames": ["gene_id", "gene_symbol", *sample_fieldnames],
        "sample_metadata_fieldnames": [
            "sample_id",
            "pid",
            "bid",
            "sex",
            "sex_label",
            "intervention",
            "tissue",
            "transcript_tissue",
            "tissue_code_no",
            "timepoint_label",
        ],
        "summary": summary,
    }
