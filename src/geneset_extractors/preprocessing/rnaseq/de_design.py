from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from geneset_extractors.preprocessing.rnaseq.de_io import read_tsv_rows


@dataclass(frozen=True)
class ComparisonSpec:
    comparison_id: str
    comparison_kind: str
    group_column: str
    group_a: str
    group_b: str | None
    strata: dict[str, str]



def parse_csv_columns(value: str | None) -> list[str]:
    if value is None:
        return []
    out: list[str] = []
    for token in str(value).split(","):
        item = token.strip()
        if item:
            out.append(item)
    return out



def _stratum_key(row: dict[str, str], stratify_by: list[str]) -> tuple[tuple[str, str], ...]:
    return tuple((col, str(row.get(col, "")).strip()) for col in stratify_by)



def _format_strata(parts: dict[str, str]) -> str:
    if not parts:
        return ""
    return "__" + "__".join(f"{col}={parts[col]}" for col in sorted(parts))



def _rows_by_stratum(rows: list[dict[str, str]], stratify_by: list[str]) -> dict[tuple[tuple[str, str], ...], list[dict[str, str]]]:
    grouped: dict[tuple[tuple[str, str], ...], list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(_stratum_key(row, stratify_by), []).append(row)
    return grouped



def build_cli_comparisons(
    metadata_rows: list[dict[str, str]],
    *,
    comparison_mode: str,
    group_column: str,
    condition_a: str | None,
    condition_b: str | None,
    reference_level: str | None,
    stratify_by: list[str],
) -> list[ComparisonSpec]:
    if not metadata_rows:
        raise ValueError("No metadata rows were available to generate comparisons")
    specs: list[ComparisonSpec] = []
    grouped = _rows_by_stratum(metadata_rows, stratify_by)
    for stratum_key, rows in sorted(grouped.items()):
        strata = {key: value for key, value in stratum_key}
        levels = sorted({str(row.get(group_column, "")).strip() for row in rows if str(row.get(group_column, "")).strip()})
        if comparison_mode == "condition_a_vs_b":
            if not condition_a or not condition_b:
                raise ValueError("condition_a_vs_b requires --condition_a and --condition_b")
            if condition_a not in levels or condition_b not in levels:
                continue
            comparison_id = f"{group_column}={condition_a}_vs_{condition_b}{_format_strata(strata)}"
            specs.append(
                ComparisonSpec(
                    comparison_id=comparison_id,
                    comparison_kind=comparison_mode,
                    group_column=group_column,
                    group_a=condition_a,
                    group_b=condition_b,
                    strata=strata,
                )
            )
            continue
        if comparison_mode == "group_vs_rest":
            if not condition_a:
                raise ValueError("group_vs_rest requires --condition_a as the focal group")
            if condition_a not in levels or len(levels) < 2:
                continue
            comparison_id = f"{group_column}={condition_a}_vs_rest{_format_strata(strata)}"
            specs.append(
                ComparisonSpec(
                    comparison_id=comparison_id,
                    comparison_kind=comparison_mode,
                    group_column=group_column,
                    group_a=condition_a,
                    group_b=None,
                    strata=strata,
                )
            )
            continue
        if comparison_mode == "reference_level":
            if not reference_level:
                raise ValueError("reference_level requires --reference_level")
            if reference_level not in levels:
                continue
            for level in levels:
                if level == reference_level:
                    continue
                comparison_id = f"{group_column}={level}_vs_{reference_level}{_format_strata(strata)}"
                specs.append(
                    ComparisonSpec(
                        comparison_id=comparison_id,
                        comparison_kind=comparison_mode,
                        group_column=group_column,
                        group_a=level,
                        group_b=reference_level,
                        strata=strata,
                    )
                )
            continue
        raise ValueError(f"Unsupported comparison_mode: {comparison_mode}")
    if not specs:
        raise ValueError("No comparisons were generated from the supplied metadata and comparison settings")
    return specs



def read_comparisons_tsv(
    path: str | Path,
    *,
    default_group_column: str | None,
    stratify_by: list[str],
) -> list[ComparisonSpec]:
    table = read_tsv_rows(path)
    required = {"comparison_id", "group_a"}
    missing = sorted(col for col in required if col not in table.fieldnames)
    if missing:
        raise ValueError(f"comparisons_tsv is missing required columns: {', '.join(missing)}")
    specs: list[ComparisonSpec] = []
    for row in table.rows:
        comparison_id = str(row.get("comparison_id", "")).strip()
        if not comparison_id:
            continue
        group_column = str(row.get("group_column", default_group_column or "")).strip()
        if not group_column:
            raise ValueError("Each comparisons_tsv row needs group_column or a workflow-level --group_column")
        comparison_kind = str(row.get("comparison_kind", "condition_a_vs_b")).strip() or "condition_a_vs_b"
        strata = {col: str(row.get(col, "")).strip() for col in stratify_by if str(row.get(col, "")).strip()}
        specs.append(
            ComparisonSpec(
                comparison_id=comparison_id,
                comparison_kind=comparison_kind,
                group_column=group_column,
                group_a=str(row.get("group_a", "")).strip(),
                group_b=(str(row.get("group_b", "")).strip() or None),
                strata=strata,
            )
        )
    if not specs:
        raise ValueError("No valid comparison rows were parsed from comparisons_tsv")
    return specs



def comparison_manifest_rows(specs: list[ComparisonSpec]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for spec in specs:
        row: dict[str, Any] = {
            "comparison_id": spec.comparison_id,
            "comparison_kind": spec.comparison_kind,
            "group_column": spec.group_column,
            "group_a": spec.group_a,
            "group_b": spec.group_b or "",
        }
        for key, value in sorted(spec.strata.items()):
            row[key] = value
        rows.append(row)
    return rows
