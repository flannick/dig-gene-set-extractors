from __future__ import annotations

import csv
import gzip
import importlib.resources as importlib_resources
import json
from pathlib import Path
import re
from statistics import mean, median
from typing import Callable


META_PREFIX = "Metadata_"
MORPHOLOGY_TARGET_ANNOTATION_FIELDS = ["target_family", "target_class", "mechanism_label", "pathway_seed"]
DEFAULT_MORPHOLOGY_TARGET_ANNOTATIONS_RESOURCE = "morphology_target_annotations_human_v1.tsv"


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
    if text.lower() in {"na", "nan", "none", "null", "inf", "-inf"}:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    if parsed != parsed:
        return None
    return float(parsed)


def _sniff_delimiter(path: Path) -> str:
    if path.suffix.lower() == ".csv":
        return ","
    return "\t"


def read_profiles_table(
    path: str | Path,
    *,
    id_column: str = "sample_id",
    delimiter: str | None = None,
) -> tuple[dict[str, dict[str, float]], dict[str, object]]:
    p = Path(path)
    if p.suffix.lower() == ".parquet":
        try:
            import pandas as pd  # type: ignore
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                "Parquet input requires pandas/pyarrow. Use TSV/CSV instead or install optional dependencies."
            ) from exc
        df = pd.read_parquet(p)
        fieldnames = [str(c) for c in df.columns]
        if id_column not in fieldnames:
            raise ValueError(f"Profile table missing id column '{id_column}'. Available columns: {', '.join(fieldnames)}")
        profiles: dict[str, dict[str, float]] = {}
        numeric_columns: list[str] = []
        ignored_metadata: list[str] = []
        for col in fieldnames:
            if col == id_column:
                continue
            if col.startswith(META_PREFIX):
                ignored_metadata.append(col)
                continue
            series = df[col]
            if pd.api.types.is_numeric_dtype(series):
                numeric_columns.append(col)
        for _, row in df.iterrows():
            sample_id = str(row[id_column]).strip()
            if not sample_id:
                continue
            feat = {}
            for col in numeric_columns:
                val = _parse_float(row[col])
                if val is not None:
                    feat[col] = val
            profiles[sample_id] = feat
        return profiles, {
            "path": str(path),
            "format": "parquet",
            "id_column": id_column,
            "n_rows": int(len(df)),
            "n_profiles": len(profiles),
            "n_numeric_features": len(numeric_columns),
            "ignored_metadata_columns": ignored_metadata,
            "feature_names": numeric_columns,
        }

    delim = delimiter or _sniff_delimiter(p)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if not reader.fieldnames:
            raise ValueError(f"Profile table has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        if id_column not in fieldnames:
            raise ValueError(f"Profile table missing id column '{id_column}'. Available columns: {', '.join(fieldnames)}")
        rows = list(reader)
    profiles: dict[str, dict[str, float]] = {}
    numeric_columns: list[str] = []
    ignored_metadata: list[str] = []
    for col in fieldnames:
        if col == id_column:
            continue
        if col.startswith(META_PREFIX):
            ignored_metadata.append(col)
            continue
        any_numeric = False
        for row in rows[:50]:
            if _parse_float(row.get(col)) is not None:
                any_numeric = True
                break
        if any_numeric:
            numeric_columns.append(col)
    for row in rows:
        sample_id = str(row.get(id_column, "")).strip()
        if not sample_id:
            continue
        feat: dict[str, float] = {}
        for col in numeric_columns:
            val = _parse_float(row.get(col))
            if val is not None:
                feat[col] = val
        profiles[sample_id] = feat
    return profiles, {
        "path": str(path),
        "format": "delimited",
        "delimiter": delim,
        "id_column": id_column,
        "n_rows": len(rows),
        "n_profiles": len(profiles),
        "n_numeric_features": len(numeric_columns),
        "ignored_metadata_columns": ignored_metadata,
        "feature_names": numeric_columns,
    }


def read_metadata_tsv(path: str | Path, *, id_column: str, delimiter: str = "\t") -> tuple[dict[str, dict[str, str]], dict[str, object]]:
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Metadata table has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        if id_column not in fieldnames:
            raise ValueError(f"Metadata table missing id column '{id_column}'. Available columns: {', '.join(fieldnames)}")
        rows = list(reader)
    out: dict[str, dict[str, str]] = {}
    missing = 0
    for row in rows:
        sample_id = str(row.get(id_column, "")).strip()
        if not sample_id:
            missing += 1
            continue
        out[sample_id] = {str(k): "" if v is None else str(v) for k, v in row.items() if k is not None}
    return out, {"path": str(path), "id_column": id_column, "n_rows": len(rows), "n_records": len(out), "n_missing_id": missing}


def read_compound_targets_tsv(
    path: str | Path,
    *,
    compound_id_column: str = "compound_id",
    gene_symbol_column: str = "gene_symbol",
    weight_column: str = "weight",
    delimiter: str = "\t",
) -> tuple[dict[str, dict[str, float]], dict[str, object]]:
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Compound target table has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        for required in (compound_id_column, gene_symbol_column):
            if required not in fieldnames:
                raise ValueError(f"Compound target table missing column '{required}'. Available columns: {', '.join(fieldnames)}")
        rows = list(reader)
    raw: dict[str, dict[str, float]] = {}
    missing = 0
    for row in rows:
        compound_id = str(row.get(compound_id_column, "")).strip()
        gene_symbol = str(row.get(gene_symbol_column, "")).strip().upper()
        if not compound_id or not gene_symbol:
            missing += 1
            continue
        weight = _parse_float(row.get(weight_column)) if weight_column in row else None
        if weight is None or weight <= 0.0:
            weight = 1.0
        raw.setdefault(compound_id, {})
        raw[compound_id][gene_symbol] = float(raw[compound_id].get(gene_symbol, 0.0)) + float(weight)
    normalized: dict[str, dict[str, float]] = {}
    for compound_id, mapping in raw.items():
        total = sum(mapping.values())
        if total <= 0.0:
            continue
        normalized[compound_id] = {gene: float(value) / total for gene, value in mapping.items()}
    return normalized, {"path": str(path), "n_rows": len(rows), "n_compounds": len(normalized), "n_missing_required": missing}


def read_compound_targets_auto(
    path: str | Path,
    *,
    compound_id_column: str | None = None,
    gene_symbol_column: str | None = None,
    weight_column: str | None = None,
    delimiter: str = "\t",
) -> tuple[dict[str, dict[str, float]], dict[str, object]]:
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Compound target table has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        rows = list(reader)

    compound_candidates = [
        compound_id_column,
        "compound_id",
        "broad_id",
        "broad_sample",
        "perturbation_id",
        "pert_iname",
        "name",
    ]
    gene_candidates = [
        gene_symbol_column,
        "gene_symbol",
        "target",
        "targets",
        "targets_list",
        "moa_targets",
        "gene_targets",
        "target_gene",
        "target_genes",
    ]
    weight_candidates = [weight_column, "weight", "target_weight"]

    resolved_compound = next((c for c in compound_candidates if c and c in fieldnames), None)
    if resolved_compound is None:
        raise ValueError(
            f"Could not identify compound id column in {p}. Available columns: {', '.join(fieldnames)}"
        )
    resolved_gene = next((c for c in gene_candidates if c and c in fieldnames), None)
    if resolved_gene is None:
        raise ValueError(
            f"Could not identify target gene column in {p}. Available columns: {', '.join(fieldnames)}"
        )
    resolved_weight = next((c for c in weight_candidates if c and c in fieldnames), None)

    raw: dict[str, dict[str, float]] = {}
    n_missing_required = 0
    n_multi_target_rows = 0
    fallback_used = resolved_gene != "gene_symbol"
    for row in rows:
        compound_id = str(row.get(resolved_compound, "")).strip()
        target_value = str(row.get(resolved_gene, "")).strip()
        if not compound_id or not target_value:
            n_missing_required += 1
            continue
        targets = [
            token.strip().upper()
            for token in re.split(r"[;,|]", target_value)
            if token and token.strip()
        ]
        targets = [token for token in targets if token not in {"NA", "NAN", "NONE", "NULL"}]
        if not targets:
            n_missing_required += 1
            continue
        if len(targets) > 1:
            n_multi_target_rows += 1
        weight = _parse_float(row.get(resolved_weight)) if resolved_weight else None
        if weight is not None and weight > 0.0 and len(targets) == 1:
            per_target_weight = float(weight)
        else:
            per_target_weight = 1.0 / float(len(targets))
        raw.setdefault(compound_id, {})
        for target in targets:
            raw[compound_id][target] = float(raw[compound_id].get(target, 0.0)) + float(per_target_weight)

    normalized: dict[str, dict[str, float]] = {}
    zero_target_compounds = 0
    for compound_id, mapping in raw.items():
        total = sum(mapping.values())
        if total <= 0.0:
            zero_target_compounds += 1
            continue
        normalized[compound_id] = {gene: float(value) / total for gene, value in mapping.items()}
    return normalized, {
        "path": str(path),
        "n_rows": len(rows),
        "n_compounds": len(normalized),
        "n_missing_required": n_missing_required,
        "n_zero_target_compounds": zero_target_compounds,
        "resolved_compound_id_column": resolved_compound,
        "resolved_gene_symbol_column": resolved_gene,
        "resolved_weight_column": resolved_weight,
        "fallback_used": fallback_used,
        "n_multi_target_rows": n_multi_target_rows,
    }


def read_feature_schema_tsv(path: str | Path, *, feature_column: str = "feature", delimiter: str = "\t") -> tuple[list[str], dict[str, object]]:
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Feature schema has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        if feature_column not in fieldnames:
            raise ValueError(f"Feature schema missing column '{feature_column}'. Available columns: {', '.join(fieldnames)}")
        features = [str(row.get(feature_column, "")).strip() for row in reader]
    out = [f for f in features if f]
    return out, {"path": str(path), "n_features": len(out), "feature_column": feature_column}


def read_feature_stats_tsv(
    path: str | Path,
    *,
    feature_column: str = "feature",
    center_column: str = "center",
    scale_column: str = "scale",
    delimiter: str = "\t",
) -> tuple[dict[str, tuple[float, float]], dict[str, object]]:
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Feature stats table has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        for required in (feature_column, center_column, scale_column):
            if required not in fieldnames:
                raise ValueError(f"Feature stats missing column '{required}'. Available columns: {', '.join(fieldnames)}")
        rows = list(reader)
    out: dict[str, tuple[float, float]] = {}
    missing = 0
    for row in rows:
        feature = str(row.get(feature_column, "")).strip()
        center = _parse_float(row.get(center_column))
        scale = _parse_float(row.get(scale_column))
        if not feature or center is None or scale is None:
            missing += 1
            continue
        out[feature] = (float(center), float(scale if scale != 0.0 else 1.0))
    return out, {"path": str(path), "n_rows": len(rows), "n_features": len(out), "n_missing_required": missing}


def read_target_annotations_tsv(
    path: str | Path,
    *,
    gene_symbol_column: str = "gene_symbol",
    delimiter: str = "\t",
) -> tuple[dict[str, dict[str, str]], dict[str, object]]:
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Target annotation table has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        if gene_symbol_column not in fieldnames:
            raise ValueError(f"Target annotation table missing column '{gene_symbol_column}'. Available columns: {', '.join(fieldnames)}")
        rows = list(reader)
    out: dict[str, dict[str, str]] = {}
    missing = 0
    for row in rows:
        gene_symbol = str(row.get(gene_symbol_column, "")).strip().upper()
        if not gene_symbol:
            missing += 1
            continue
        out[gene_symbol] = {
            str(key): "" if value is None else str(value).strip()
            for key, value in row.items()
            if key is not None and str(key) != gene_symbol_column
        }
    return out, {
        "path": str(path),
        "n_rows": len(rows),
        "n_genes": len(out),
        "n_missing_gene_symbol": missing,
        "fields": [field for field in fieldnames if field != gene_symbol_column],
    }


def read_packaged_target_annotations_tsv() -> tuple[dict[str, dict[str, str]], dict[str, object]]:
    resource = importlib_resources.files("geneset_extractors.resources").joinpath(
        DEFAULT_MORPHOLOGY_TARGET_ANNOTATIONS_RESOURCE
    )
    with importlib_resources.as_file(resource) as resolved_path:
        payload, summary = read_target_annotations_tsv(resolved_path)
    summary = dict(summary)
    summary["source"] = "package_default"
    summary["resource_name"] = DEFAULT_MORPHOLOGY_TARGET_ANNOTATIONS_RESOURCE
    return payload, summary


def aggregate_profiles(
    profiles: dict[str, dict[str, float]],
    metadata: dict[str, dict[str, str]] | None,
    *,
    group_by: str | None,
    aggregate: str,
) -> tuple[dict[str, dict[str, float]], dict[str, object], dict[str, list[str]]]:
    if not group_by or aggregate == "none":
        membership = {sample_id: [sample_id] for sample_id in profiles}
        return dict(profiles), {"group_by": group_by, "aggregate": "none", "n_groups": len(profiles)}, membership
    if metadata is None:
        raise ValueError("--group_query_by requires --query_metadata_tsv")
    groups: dict[str, list[str]] = {}
    for sample_id in profiles:
        group = str(metadata.get(sample_id, {}).get(group_by, "")).strip()
        if not group:
            group = sample_id
        groups.setdefault(group, []).append(sample_id)
    reducer: Callable[[list[float]], float]
    if aggregate == "median":
        reducer = median
    elif aggregate == "mean":
        reducer = mean
    else:
        raise ValueError(f"Unsupported query_aggregate: {aggregate}")
    aggregated: dict[str, dict[str, float]] = {}
    for group, members in groups.items():
        features = sorted({feat for member in members for feat in profiles.get(member, {})})
        out: dict[str, float] = {}
        for feat in features:
            vals = [profiles[m][feat] for m in members if feat in profiles.get(m, {})]
            if vals:
                out[feat] = float(reducer(vals))
        aggregated[group] = out
    return aggregated, {"group_by": group_by, "aggregate": aggregate, "n_groups": len(aggregated)}, groups


def write_bundle_manifest(path: str | Path, payload: dict[str, object]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
