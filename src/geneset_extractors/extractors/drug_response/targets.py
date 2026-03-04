from __future__ import annotations

import csv
import gzip
from pathlib import Path
import re


_SPLIT_RE = re.compile(r"[;,]")
_GENE_TOKEN_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")
_INVALID_TOKENS = {"", "na", "n/a", "none", "null", "nan"}


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
    try:
        parsed = float(text)
    except ValueError:
        return None
    if parsed != parsed:
        return None
    return float(parsed)


def load_target_aliases_tsv(
    path: str | Path,
    *,
    alias_column: str = "alias",
    gene_symbol_column: str = "gene_symbol",
    delimiter: str = "\t",
) -> tuple[dict[str, str], dict[str, object]]:
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Alias table has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        if alias_column not in fieldnames or gene_symbol_column not in fieldnames:
            raise ValueError(
                "Alias table is missing required columns "
                f"'{alias_column}' and/or '{gene_symbol_column}'. "
                f"Available columns: {', '.join(fieldnames)}"
            )
        aliases: dict[str, str] = {}
        n_rows = 0
        n_missing = 0
        for row in reader:
            n_rows += 1
            alias = str(row.get(alias_column, "")).strip()
            symbol = str(row.get(gene_symbol_column, "")).strip()
            if not alias or not symbol:
                n_missing += 1
                continue
            aliases[alias.upper()] = symbol.upper()
    summary = {
        "path": str(path),
        "n_rows": n_rows,
        "n_aliases": len(aliases),
        "n_missing_required": n_missing,
        "columns": {"alias": alias_column, "gene_symbol": gene_symbol_column},
    }
    return aliases, summary


def parse_target_string(
    raw: str,
    *,
    aliases: dict[str, str] | None = None,
    uppercase: bool = True,
) -> tuple[list[str], dict[str, int]]:
    tokens = [tok.strip() for tok in _SPLIT_RE.split(str(raw))]
    out: list[str] = []
    n_total = 0
    n_dropped = 0
    n_alias_mapped = 0
    for token in tokens:
        if not token:
            continue
        n_total += 1
        norm = token.upper() if uppercase else token
        if norm.lower() in _INVALID_TOKENS:
            n_dropped += 1
            continue
        if aliases is not None:
            mapped = aliases.get(norm.upper())
            if mapped:
                norm = mapped
                n_alias_mapped += 1
        if not _GENE_TOKEN_RE.match(norm):
            n_dropped += 1
            continue
        out.append(norm)
    unique = sorted(set(out))
    return unique, {
        "n_tokens_total": n_total,
        "n_tokens_kept": len(unique),
        "n_tokens_dropped": n_dropped,
        "n_alias_mapped": n_alias_mapped,
    }


def build_drug_targets_from_table_rows(
    rows: list[dict[str, object]],
    *,
    aliases: dict[str, str] | None = None,
) -> tuple[dict[str, dict[str, float]], dict[str, object]]:
    raw_by_drug: dict[str, dict[str, float]] = {}
    n_rows = 0
    n_rows_missing = 0
    n_invalid_tokens = 0
    n_tokens_total = 0
    n_alias_mapped = 0
    for row in rows:
        n_rows += 1
        drug_id = str(row.get("drug_id", "")).strip()
        gene_symbol = str(row.get("gene_symbol", "")).strip()
        if not drug_id or not gene_symbol:
            n_rows_missing += 1
            continue
        weight = _parse_float(row.get("weight"))
        if weight is None:
            weight = 1.0
        if weight <= 0:
            continue
        genes, parse_stats = parse_target_string(gene_symbol, aliases=aliases, uppercase=True)
        n_tokens_total += int(parse_stats["n_tokens_total"])
        n_invalid_tokens += int(parse_stats["n_tokens_dropped"])
        n_alias_mapped += int(parse_stats["n_alias_mapped"])
        if not genes:
            continue
        per_gene = float(weight) / float(len(genes))
        raw_by_drug.setdefault(drug_id, {})
        for gene in genes:
            raw_by_drug[drug_id][gene] = float(raw_by_drug[drug_id].get(gene, 0.0)) + per_gene
    normalized: dict[str, dict[str, float]] = {}
    for drug_id, gene_weights in raw_by_drug.items():
        total = sum(float(v) for v in gene_weights.values() if float(v) > 0.0)
        if total <= 0:
            continue
        normalized[drug_id] = {gene: float(v) / total for gene, v in gene_weights.items() if float(v) > 0.0}
    n_drugs = len(raw_by_drug)
    n_drugs_with_targets = len(normalized)
    summary = {
        "source": "drug_targets_tsv",
        "n_rows": n_rows,
        "n_rows_missing_required": n_rows_missing,
        "n_drugs_seen": n_drugs,
        "n_drugs_with_targets": n_drugs_with_targets,
        "fraction_drugs_with_targets": (float(n_drugs_with_targets) / float(n_drugs) if n_drugs else 0.0),
        "n_tokens_total": n_tokens_total,
        "n_tokens_invalid": n_invalid_tokens,
        "n_alias_mapped": n_alias_mapped,
    }
    return normalized, summary


def build_drug_targets_from_target_text(
    *,
    drug_to_target_text: dict[str, str],
    aliases: dict[str, str] | None = None,
) -> tuple[dict[str, dict[str, float]], dict[str, object]]:
    out: dict[str, dict[str, float]] = {}
    n_tokens_total = 0
    n_tokens_invalid = 0
    n_alias_mapped = 0
    n_drugs_with_text = 0
    for drug_id, raw in drug_to_target_text.items():
        text = str(raw).strip()
        if not text:
            continue
        n_drugs_with_text += 1
        genes, parse_stats = parse_target_string(text, aliases=aliases, uppercase=True)
        n_tokens_total += int(parse_stats["n_tokens_total"])
        n_tokens_invalid += int(parse_stats["n_tokens_dropped"])
        n_alias_mapped += int(parse_stats["n_alias_mapped"])
        if not genes:
            continue
        uniform = 1.0 / float(len(genes))
        out[drug_id] = {gene: uniform for gene in genes}
    summary = {
        "source": "target_text",
        "n_drugs_with_target_text": n_drugs_with_text,
        "n_drugs_with_targets": len(out),
        "fraction_drugs_with_targets": (
            float(len(out)) / float(n_drugs_with_text) if n_drugs_with_text else 0.0
        ),
        "n_tokens_total": n_tokens_total,
        "n_tokens_invalid": n_tokens_invalid,
        "n_alias_mapped": n_alias_mapped,
    }
    return out, summary


def apply_drug_blacklist(
    *,
    drug_targets: dict[str, dict[str, float]],
    blacklist: set[str],
) -> tuple[dict[str, dict[str, float]], dict[str, object]]:
    if not blacklist:
        return dict(drug_targets), {"n_blacklisted": 0, "n_remaining": len(drug_targets)}
    blocked = {str(x).strip() for x in blacklist if str(x).strip()}
    out = {drug: genes for drug, genes in drug_targets.items() if drug not in blocked}
    return out, {"n_blacklisted": len(drug_targets) - len(out), "n_remaining": len(out)}
