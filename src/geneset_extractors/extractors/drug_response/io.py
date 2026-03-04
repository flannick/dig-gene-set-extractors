from __future__ import annotations

import csv
from dataclasses import dataclass
import gzip
from pathlib import Path


@dataclass(frozen=True)
class ResponseRecord:
    sample_id: str
    drug_id: str
    response: float


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _parse_float(value: object) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    low = text.lower()
    if low in {"na", "nan", "none", "null", "inf", "-inf"}:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    if parsed != parsed:
        return None
    return float(parsed)


def _parse_bool(value: object) -> bool | None:
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y", "case"}:
        return True
    if text in {"0", "false", "f", "no", "n", "control"}:
        return False
    return None


def _read_table(path: str | Path, delimiter: str) -> tuple[list[str], list[dict[str, str]]]:
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Input table has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        rows: list[dict[str, str]] = []
        for row in reader:
            norm: dict[str, str] = {}
            for key, value in row.items():
                if key is None:
                    continue
                norm[str(key)] = "" if value is None else str(value)
            rows.append(norm)
    return fieldnames, rows


def read_response_tsv(
    *,
    path: str | Path,
    sample_id_column: str,
    drug_id_column: str,
    response_column: str,
    delimiter: str = "\t",
) -> tuple[list[ResponseRecord], dict[str, object]]:
    fieldnames, rows = _read_table(path, delimiter=delimiter)
    required = [sample_id_column, drug_id_column, response_column]
    missing = [name for name in required if name not in fieldnames]
    if missing:
        raise ValueError(
            "Response table is missing required columns: "
            + ", ".join(missing)
            + f". Available columns: {', '.join(fieldnames)}"
        )
    records: list[ResponseRecord] = []
    n_missing_response = 0
    n_missing_ids = 0
    for row in rows:
        sample_id = str(row.get(sample_id_column, "")).strip()
        drug_id = str(row.get(drug_id_column, "")).strip()
        value = _parse_float(row.get(response_column))
        if not sample_id or not drug_id:
            n_missing_ids += 1
            continue
        if value is None:
            n_missing_response += 1
            continue
        records.append(ResponseRecord(sample_id=sample_id, drug_id=drug_id, response=float(value)))
    summary = {
        "mode": "generic",
        "path": str(path),
        "columns": {
            "sample_id": sample_id_column,
            "drug_id": drug_id_column,
            "response": response_column,
        },
        "n_rows": len(rows),
        "n_records_parsed": len(records),
        "n_missing_response": n_missing_response,
        "n_missing_ids": n_missing_ids,
        "n_samples": len({r.sample_id for r in records}),
        "n_drugs": len({r.drug_id for r in records}),
    }
    return records, summary


def read_drug_targets_tsv(
    *,
    path: str | Path,
    drug_id_column: str = "drug_id",
    gene_symbol_column: str = "gene_symbol",
    weight_column: str = "weight",
    source_column: str = "source",
    delimiter: str = "\t",
) -> tuple[list[dict[str, object]], dict[str, object]]:
    fieldnames, rows = _read_table(path, delimiter=delimiter)
    required = [drug_id_column, gene_symbol_column]
    missing = [name for name in required if name not in fieldnames]
    if missing:
        raise ValueError(
            "Drug target table is missing required columns: "
            + ", ".join(missing)
            + f". Available columns: {', '.join(fieldnames)}"
        )
    out: list[dict[str, object]] = []
    n_missing = 0
    n_weight_non_numeric = 0
    for row in rows:
        drug_id = str(row.get(drug_id_column, "")).strip()
        gene_symbol = str(row.get(gene_symbol_column, "")).strip()
        if not drug_id or not gene_symbol:
            n_missing += 1
            continue
        weight_raw = row.get(weight_column) if weight_column in row else None
        weight = _parse_float(weight_raw)
        if weight_raw is not None and str(weight_raw).strip() and weight is None:
            n_weight_non_numeric += 1
        source = str(row.get(source_column, "")).strip() if source_column in row else ""
        out.append(
            {
                "drug_id": drug_id,
                "gene_symbol": gene_symbol,
                "weight": weight,
                "source": source,
            }
        )
    summary = {
        "path": str(path),
        "n_rows": len(rows),
        "n_rows_parsed": len(out),
        "n_missing_required": n_missing,
        "n_weight_non_numeric": n_weight_non_numeric,
        "columns": {
            "drug_id": drug_id_column,
            "gene_symbol": gene_symbol_column,
            "weight": weight_column if weight_column in fieldnames else None,
            "source": source_column if source_column in fieldnames else None,
        },
    }
    return out, summary


def read_sample_metadata_tsv(
    *,
    path: str | Path,
    sample_id_column: str = "sample_id",
    delimiter: str = "\t",
) -> tuple[dict[str, dict[str, str]], dict[str, object]]:
    fieldnames, rows = _read_table(path, delimiter=delimiter)
    if sample_id_column not in fieldnames:
        raise ValueError(
            f"sample_metadata_tsv is missing sample id column '{sample_id_column}'. "
            f"Available columns: {', '.join(fieldnames)}"
        )
    out: dict[str, dict[str, str]] = {}
    n_missing = 0
    for row in rows:
        sample_id = str(row.get(sample_id_column, "")).strip()
        if not sample_id:
            n_missing += 1
            continue
        out[sample_id] = row
    return out, {
        "path": str(path),
        "n_rows": len(rows),
        "n_samples": len(out),
        "n_missing_sample_id": n_missing,
        "sample_id_column": sample_id_column,
    }


def read_groups_tsv(
    *,
    path: str | Path,
    sample_id_column: str = "sample_id",
    group_column: str = "group",
    delimiter: str = "\t",
) -> tuple[dict[str, str], dict[str, object]]:
    fieldnames, rows = _read_table(path, delimiter=delimiter)
    required = [sample_id_column, group_column]
    missing = [name for name in required if name not in fieldnames]
    if missing:
        raise ValueError(
            "groups_tsv is missing required columns: "
            + ", ".join(missing)
            + f". Available columns: {', '.join(fieldnames)}"
        )
    groups: dict[str, str] = {}
    n_missing = 0
    for row in rows:
        sample_id = str(row.get(sample_id_column, "")).strip()
        group = str(row.get(group_column, "")).strip()
        if not sample_id or not group:
            n_missing += 1
            continue
        groups[sample_id] = group
    return groups, {
        "path": str(path),
        "n_rows": len(rows),
        "n_assignments": len(groups),
        "n_missing_required": n_missing,
        "sample_id_column": sample_id_column,
        "group_column": group_column,
        "n_groups": len(set(groups.values())),
    }


def read_case_control_tsv(
    *,
    path: str | Path,
    sample_id_column: str = "sample_id",
    is_case_column: str = "is_case",
    delimiter: str = "\t",
) -> tuple[dict[str, bool], dict[str, object]]:
    fieldnames, rows = _read_table(path, delimiter=delimiter)
    required = [sample_id_column, is_case_column]
    missing = [name for name in required if name not in fieldnames]
    if missing:
        raise ValueError(
            "case_control_tsv is missing required columns: "
            + ", ".join(missing)
            + f". Available columns: {', '.join(fieldnames)}"
        )
    labels: dict[str, bool] = {}
    n_missing = 0
    n_invalid = 0
    for row in rows:
        sample_id = str(row.get(sample_id_column, "")).strip()
        value = _parse_bool(row.get(is_case_column))
        if not sample_id:
            n_missing += 1
            continue
        if value is None:
            n_invalid += 1
            continue
        labels[sample_id] = bool(value)
    return labels, {
        "path": str(path),
        "n_rows": len(rows),
        "n_labels": len(labels),
        "n_missing_sample_id": n_missing,
        "n_invalid_labels": n_invalid,
        "sample_id_column": sample_id_column,
        "is_case_column": is_case_column,
    }


def read_drug_id_list_tsv(path: str | Path, delimiter: str = "\t") -> tuple[set[str], dict[str, object]]:
    p = Path(path)
    drug_ids: set[str] = set()
    n_rows = 0
    with _open_text(p) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            n_rows += 1
            token = line.split(delimiter)[0] if delimiter in line else line.split()[0]
            drug_id = str(token).strip()
            if not drug_id:
                continue
            if drug_id.lower() in {"drug_id", "drug", "id"} and n_rows == 1:
                continue
            drug_ids.add(drug_id)
    return drug_ids, {"path": str(path), "n_rows": n_rows, "n_drug_ids": len(drug_ids)}


def read_prism_treatment_info_csv(
    *,
    path: str | Path,
    column_name_column: str = "column_name",
    broad_id_column: str = "broad_id",
    target_column: str = "target",
) -> tuple[dict[str, str], dict[str, str], dict[str, object]]:
    fieldnames, rows = _read_table(path, delimiter=",")
    if column_name_column not in fieldnames:
        raise ValueError(
            f"PRISM treatment info is missing '{column_name_column}'. "
            f"Available columns: {', '.join(fieldnames)}"
        )
    column_to_drug: dict[str, str] = {}
    drug_to_target_text: dict[str, str] = {}
    n_missing_column_name = 0
    n_missing_drug = 0
    for row in rows:
        column_name = str(row.get(column_name_column, "")).strip()
        if not column_name:
            n_missing_column_name += 1
            continue
        broad_id = str(row.get(broad_id_column, "")).strip() if broad_id_column in row else ""
        drug_id = broad_id or column_name
        if not drug_id:
            n_missing_drug += 1
            continue
        column_to_drug[column_name] = drug_id
        if target_column in row:
            target_text = str(row.get(target_column, "")).strip()
            if target_text and drug_id not in drug_to_target_text:
                drug_to_target_text[drug_id] = target_text
    summary = {
        "path": str(path),
        "n_rows": len(rows),
        "n_columns_mapped": len(column_to_drug),
        "n_drugs_with_target_text": len(drug_to_target_text),
        "n_missing_column_name": n_missing_column_name,
        "n_missing_drug": n_missing_drug,
        "columns": {
            "column_name": column_name_column,
            "broad_id": broad_id_column if broad_id_column in fieldnames else None,
            "target": target_column if target_column in fieldnames else None,
        },
    }
    return column_to_drug, drug_to_target_text, summary


def read_prism_cell_line_info_csv(
    *,
    path: str | Path,
    sample_id_column: str = "row_name",
) -> tuple[dict[str, dict[str, str]], dict[str, object]]:
    fieldnames, rows = _read_table(path, delimiter=",")
    if sample_id_column not in fieldnames:
        raise ValueError(
            f"PRISM cell-line info is missing sample id column '{sample_id_column}'. "
            f"Available columns: {', '.join(fieldnames)}"
        )
    out: dict[str, dict[str, str]] = {}
    n_missing = 0
    for row in rows:
        sample_id = str(row.get(sample_id_column, "")).strip()
        if not sample_id:
            n_missing += 1
            continue
        out[sample_id] = row
    summary = {
        "path": str(path),
        "n_rows": len(rows),
        "n_samples": len(out),
        "n_missing_sample_id": n_missing,
        "sample_id_column": sample_id_column,
        "columns": fieldnames,
    }
    return out, summary


def read_prism_matrix_csv(
    *,
    path: str | Path,
    column_to_drug: dict[str, str],
    sample_id_column: str = "row_name",
) -> tuple[list[ResponseRecord], dict[str, object]]:
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=",")
        if not reader.fieldnames:
            raise ValueError(f"PRISM matrix has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        if sample_id_column not in fieldnames:
            sample_id_column = fieldnames[0]
        records: list[ResponseRecord] = []
        n_rows = 0
        n_missing_sample = 0
        n_non_numeric = 0
        n_unmapped_columns = 0
        for row in reader:
            n_rows += 1
            sample_id = str(row.get(sample_id_column, "")).strip()
            if not sample_id:
                n_missing_sample += 1
                continue
            for column_name in fieldnames:
                if column_name == sample_id_column:
                    continue
                drug_id = column_to_drug.get(column_name, column_name)
                if not drug_id:
                    n_unmapped_columns += 1
                    continue
                value = _parse_float(row.get(column_name))
                if value is None:
                    raw = row.get(column_name)
                    if raw is not None and str(raw).strip():
                        n_non_numeric += 1
                    continue
                records.append(
                    ResponseRecord(sample_id=sample_id, drug_id=drug_id, response=float(value))
                )
    summary = {
        "mode": "prism",
        "path": str(path),
        "sample_id_column": sample_id_column,
        "n_rows": n_rows,
        "n_records_parsed": len(records),
        "n_samples": len({r.sample_id for r in records}),
        "n_drugs": len({r.drug_id for r in records}),
        "n_missing_sample_id": n_missing_sample,
        "n_non_numeric_values": n_non_numeric,
        "n_unmapped_columns": n_unmapped_columns,
    }
    return records, summary
