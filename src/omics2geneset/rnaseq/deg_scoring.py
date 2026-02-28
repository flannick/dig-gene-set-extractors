from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path
import re

from omics2geneset.io.gtf import read_genes_from_gtf


DEFAULT_GENE_SYMBOL_COLUMNS = ("gene_symbol", "gene_name")
DEFAULT_STAT_COLUMNS = ("stat",)
DEFAULT_LOGFC_COLUMNS = ("log2fc", "log2FoldChange", "logFC")
DEFAULT_PADJ_COLUMNS = ("padj", "FDR", "adj.P.Val")
DEFAULT_PVALUE_COLUMNS = ("pvalue", "PValue", "P.Value")

DEFAULT_EXCLUDE_GENE_REGEXES = (
    r"^MT-",
    r"^RPL",
    r"^RPS",
    r"^mt-",
    r"^Rpl",
    r"^Rps",
)


@dataclass(frozen=True)
class DEGRow:
    line_no: int
    values: dict[str, str]


def read_deg_tsv(path: str | Path) -> tuple[list[str], list[DEGRow]]:
    rows: list[DEGRow] = []
    with Path(path).open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"DE table has no header: {path}")
        fieldnames = [str(x) for x in reader.fieldnames]
        for idx, row in enumerate(reader, start=2):
            rows.append(DEGRow(line_no=idx, values={k: str(v) for k, v in row.items() if k is not None}))
    return fieldnames, rows


def _normalize_optional(value: str | None) -> str | None:
    if value is None:
        return None
    v = str(value).strip()
    return v or None


def resolve_column(
    fieldnames: list[str],
    explicit: str | None,
    candidates: tuple[str, ...],
    label: str,
    required: bool,
) -> str | None:
    explicit_name = _normalize_optional(explicit)
    if explicit_name is not None:
        if explicit_name not in fieldnames:
            raise ValueError(
                f"Column '{explicit_name}' for {label} not found in DE table. "
                f"Available columns: {', '.join(fieldnames)}"
            )
        return explicit_name
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    if required:
        raise ValueError(
            f"Could not resolve required column for {label}. "
            f"Tried defaults: {', '.join(candidates)}. "
            f"Available columns: {', '.join(fieldnames)}"
        )
    return None


def resolve_deg_columns(
    fieldnames: list[str],
    *,
    gene_id_column: str,
    gene_symbol_column: str | None,
    stat_column: str | None,
    logfc_column: str | None,
    padj_column: str | None,
    pvalue_column: str | None,
    score_column: str | None,
) -> dict[str, str | None]:
    return {
        "gene_id_column": resolve_column(
            fieldnames,
            gene_id_column,
            (gene_id_column,),
            "gene_id",
            required=True,
        ),
        "gene_symbol_column": resolve_column(
            fieldnames,
            gene_symbol_column,
            DEFAULT_GENE_SYMBOL_COLUMNS,
            "gene_symbol",
            required=False,
        ),
        "stat_column": resolve_column(
            fieldnames,
            stat_column,
            DEFAULT_STAT_COLUMNS,
            "stat",
            required=False,
        ),
        "logfc_column": resolve_column(
            fieldnames,
            logfc_column,
            DEFAULT_LOGFC_COLUMNS,
            "logFC",
            required=False,
        ),
        "padj_column": resolve_column(
            fieldnames,
            padj_column,
            DEFAULT_PADJ_COLUMNS,
            "adjusted p-value",
            required=False,
        ),
        "pvalue_column": resolve_column(
            fieldnames,
            pvalue_column,
            DEFAULT_PVALUE_COLUMNS,
            "p-value",
            required=False,
        ),
        "score_column": _normalize_optional(score_column),
    }


def compile_exclude_gene_regexes(
    extra_patterns: list[str] | None,
    disable_default_excludes: bool,
) -> list[re.Pattern[str]]:
    patterns: list[str] = []
    if not disable_default_excludes:
        patterns.extend(DEFAULT_EXCLUDE_GENE_REGEXES)
    for raw in extra_patterns or []:
        p = str(raw).strip()
        if not p:
            raise ValueError("exclude_gene_regex entries must be non-empty")
        patterns.append(p)
    return [re.compile(p, flags=re.IGNORECASE) for p in patterns]


def _strip_version(gene_id: str) -> str:
    return str(gene_id).split(".", 1)[0]


def _resolve_pvalue_column(rows: list[DEGRow], padj_column: str | None, pvalue_column: str | None) -> str | None:
    if padj_column and _column_has_parseable_numeric(rows, padj_column):
        return padj_column
    if pvalue_column and _column_has_parseable_numeric(rows, pvalue_column):
        return pvalue_column
    return None


def _column_has_parseable_numeric(rows: list[DEGRow], column: str | None) -> bool:
    if not column:
        return False
    for row in rows:
        parsed = _parse_float_soft(row.values.get(column))
        if parsed is not None:
            return True
    return False


def _resolve_score_mode(
    rows: list[DEGRow],
    requested_mode: str,
    stat_column: str | None,
    logfc_column: str | None,
    padj_column: str | None,
    pvalue_column: str | None,
    score_column: str | None,
) -> tuple[str, str | None]:
    if requested_mode == "custom_column":
        if not score_column:
            raise ValueError("score_mode=custom_column requires --score_column")
        return "custom_column", None

    if requested_mode == "stat":
        if not stat_column:
            raise ValueError("score_mode=stat requires --stat_column (or a table column named 'stat')")
        return "stat", None

    if requested_mode == "logfc_times_neglog10p":
        if not logfc_column:
            raise ValueError(
                "score_mode=logfc_times_neglog10p requires --logfc_column (or a default logFC column)"
            )
        p_col = _resolve_pvalue_column(rows, padj_column, pvalue_column)
        if not p_col:
            raise ValueError(
                "score_mode=logfc_times_neglog10p requires --padj_column or --pvalue_column "
                "(or default columns padj/FDR/adj.P.Val or pvalue/PValue/P.Value)."
            )
        return "logfc_times_neglog10p", p_col

    if requested_mode != "auto":
        raise ValueError(f"Unsupported score_mode: {requested_mode}")

    if stat_column and _column_has_parseable_numeric(rows, stat_column):
        return "stat", None

    if logfc_column and _column_has_parseable_numeric(rows, logfc_column):
        p_col = _resolve_pvalue_column(rows, padj_column, pvalue_column)
        if p_col:
            return "logfc_times_neglog10p", p_col

    raise ValueError(
        "score_mode=auto could not find sufficient columns. "
        "Expected either a parseable stat column, or parseable logFC plus p/padj columns. "
        "Set --score_mode explicitly and provide column flags "
        "(--stat_column / --logfc_column / --padj_column / --pvalue_column / --score_column)."
    )


def _parse_float_soft(raw: object) -> float | None:
    if raw is None:
        return None
    text = str(raw).strip()
    if not text or text.lower() in {"na", "nan", "inf", "-inf"}:
        return None
    try:
        val = float(text)
    except ValueError:
        return None
    if not math.isfinite(val):
        return None
    return val


def _row_score(
    row: DEGRow,
    mode: str,
    score_column: str | None,
    stat_column: str | None,
    logfc_column: str | None,
    pvalue_column: str | None,
    neglog10p_eps: float,
    neglog10p_cap: float,
) -> float | None:
    if mode == "custom_column":
        if not score_column:
            return None
        return _parse_float_soft(row.values.get(score_column))
    if mode == "stat":
        if not stat_column:
            return None
        return _parse_float_soft(row.values.get(stat_column))
    if mode != "logfc_times_neglog10p":
        return None
    if not logfc_column or not pvalue_column:
        return None
    logfc = _parse_float_soft(row.values.get(logfc_column))
    p_val = _parse_float_soft(row.values.get(pvalue_column))
    if logfc is None or p_val is None:
        return None
    p = min(1.0, max(float(neglog10p_eps), float(p_val)))
    score_mag = min(-math.log10(p + float(neglog10p_eps)), float(neglog10p_cap))
    return float(logfc) * float(score_mag)


def score_deg_rows(
    rows: list[DEGRow],
    *,
    gene_id_column: str,
    gene_symbol_column: str | None,
    stat_column: str | None,
    logfc_column: str | None,
    padj_column: str | None,
    pvalue_column: str | None,
    score_column: str | None,
    score_mode: str,
    neglog10p_eps: float,
    neglog10p_cap: float,
) -> tuple[str, dict[str, dict[str, object]], int]:
    resolved_mode, resolved_pvalue_column = _resolve_score_mode(
        rows=rows,
        requested_mode=score_mode,
        stat_column=stat_column,
        logfc_column=logfc_column,
        padj_column=padj_column,
        pvalue_column=pvalue_column,
        score_column=score_column,
    )

    records: dict[str, dict[str, object]] = {}
    skipped_rows = 0
    for row in rows:
        gene_id = str(row.values.get(gene_id_column, "")).strip()
        if not gene_id:
            skipped_rows += 1
            continue
        score = _row_score(
            row=row,
            mode=resolved_mode,
            score_column=score_column,
            stat_column=stat_column,
            logfc_column=logfc_column,
            pvalue_column=resolved_pvalue_column,
            neglog10p_eps=neglog10p_eps,
            neglog10p_cap=neglog10p_cap,
        )
        if score is None:
            skipped_rows += 1
            continue
        rec = records.setdefault(
            gene_id,
            {
                "gene_id": gene_id,
                "score": 0.0,
                "gene_symbol": None,
                "gene_biotype": None,
            },
        )
        rec["score"] = float(rec.get("score", 0.0)) + float(score)
        if gene_symbol_column:
            symbol = str(row.values.get(gene_symbol_column, "")).strip()
            if symbol and not rec.get("gene_symbol"):
                rec["gene_symbol"] = symbol
        row_biotype = str(row.values.get("gene_biotype", "") or row.values.get("gene_type", "")).strip()
        if row_biotype and not rec.get("gene_biotype"):
            rec["gene_biotype"] = row_biotype

    if not records:
        raise ValueError(
            "No valid gene scores were parsed from the DE table. "
            "Check column mappings and score_mode."
        )
    return resolved_mode, records, skipped_rows


def apply_exclude_filters(
    records: dict[str, dict[str, object]],
    patterns: list[re.Pattern[str]],
) -> tuple[dict[str, dict[str, object]], int]:
    if not patterns:
        return records, 0
    filtered: dict[str, dict[str, object]] = {}
    removed = 0
    for gene_id, rec in records.items():
        symbol = str(rec.get("gene_symbol", "")).strip()
        if symbol and any(p.search(symbol) for p in patterns):
            removed += 1
            continue
        filtered[gene_id] = rec
    return filtered, removed


def annotate_records_with_gtf(
    records: dict[str, dict[str, object]],
    gtf_path: str | Path,
    gtf_gene_id_field: str,
) -> dict[str, dict[str, object]]:
    genes = read_genes_from_gtf(gtf_path, gene_id_field=gtf_gene_id_field)
    exact: dict[str, dict[str, str]] = {}
    nover: dict[str, dict[str, str] | None] = {}
    for gene in genes:
        entry = {
            "gene_symbol": str(gene.gene_symbol or "").strip(),
            "gene_biotype": str(gene.gene_biotype or "").strip(),
        }
        gid = str(gene.gene_id)
        exact[gid] = entry
        short = _strip_version(gid)
        if short in nover and nover[short] is not None and nover[short] != entry:
            nover[short] = None
        elif short not in nover:
            nover[short] = entry

    annotated: dict[str, dict[str, object]] = {}
    for gene_id, rec in records.items():
        match = exact.get(gene_id)
        if match is None:
            match = nover.get(_strip_version(gene_id))
        out = dict(rec)
        if match:
            if not str(out.get("gene_symbol", "")).strip():
                out["gene_symbol"] = match.get("gene_symbol") or None
            if not str(out.get("gene_biotype", "")).strip():
                out["gene_biotype"] = match.get("gene_biotype") or None
        annotated[gene_id] = out
    return annotated
