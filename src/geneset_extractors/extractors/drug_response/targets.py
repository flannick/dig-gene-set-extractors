from __future__ import annotations

import csv
import gzip
from pathlib import Path
import re
import statistics


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


def load_drug_alias_map_tsv(
    path: str | Path,
    *,
    input_drug_id_column: str = "input_drug_id",
    canonical_drug_id_column: str = "canonical_drug_id",
    namespace_column: str = "namespace",
    source_column: str = "source",
    delimiter: str = "\t",
) -> tuple[dict[str, str], dict[str, object]]:
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Drug alias map has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        required = [input_drug_id_column, canonical_drug_id_column]
        missing = [name for name in required if name not in fieldnames]
        if missing:
            raise ValueError(
                "Drug alias map is missing required columns: "
                + ", ".join(missing)
                + f". Available columns: {', '.join(fieldnames)}"
            )
        aliases: dict[str, str] = {}
        n_rows = 0
        n_missing = 0
        namespaces: set[str] = set()
        sources: set[str] = set()
        for row in reader:
            n_rows += 1
            raw_id = str(row.get(input_drug_id_column, "")).strip()
            canonical_id = str(row.get(canonical_drug_id_column, "")).strip()
            if not raw_id or not canonical_id:
                n_missing += 1
                continue
            namespace = str(row.get(namespace_column, "")).strip()
            source = str(row.get(source_column, "")).strip()
            if namespace:
                namespaces.add(namespace)
            if source:
                sources.add(source)
            aliases[raw_id] = canonical_id
        summary = {
            "path": str(path),
            "n_rows": n_rows,
            "n_aliases": len(aliases),
            "n_missing_required": n_missing,
            "namespaces": sorted(namespaces),
            "sources": sorted(sources),
            "columns": {
                "input_drug_id": input_drug_id_column,
                "canonical_drug_id": canonical_drug_id_column,
                "namespace": namespace_column if namespace_column in fieldnames else None,
                "source": source_column if source_column in fieldnames else None,
            },
        }
    return aliases, summary


def normalize_drug_target_ids(
    drug_targets: dict[str, dict[str, float]],
    *,
    alias_map: dict[str, str] | None,
) -> tuple[dict[str, dict[str, float]], dict[str, object]]:
    if not alias_map:
        return dict(drug_targets), {"n_input_drugs": len(drug_targets), "n_alias_hits": 0, "n_merged": 0}
    out: dict[str, dict[str, float]] = {}
    n_alias_hits = 0
    n_merged = 0
    for drug_id, targets in drug_targets.items():
        canonical = str(alias_map.get(drug_id, drug_id)).strip()
        if canonical != drug_id:
            n_alias_hits += 1
        dest = out.setdefault(canonical, {})
        if dest:
            n_merged += 1
        for gene_id, weight in targets.items():
            dest[gene_id] = float(dest.get(gene_id, 0.0)) + float(weight)
    normalized: dict[str, dict[str, float]] = {}
    for drug_id, targets in out.items():
        total = sum(float(v) for v in targets.values() if float(v) > 0.0)
        if total <= 0.0:
            continue
        normalized[drug_id] = {gene_id: float(v) / total for gene_id, v in targets.items() if float(v) > 0.0}
    return normalized, {
        "n_input_drugs": len(drug_targets),
        "n_output_drugs": len(normalized),
        "n_alias_hits": n_alias_hits,
        "n_merged": n_merged,
    }


def load_target_ubiquity_tsv(
    path: str | Path,
    *,
    gene_symbol_column: str = "gene_symbol",
    idf_column: str = "idf",
    family_column: str = "family",
    delimiter: str = "\t",
) -> tuple[dict[str, float], dict[str, object]]:
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Target ubiquity table has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        required = [gene_symbol_column, idf_column]
        missing = [name for name in required if name not in fieldnames]
        if missing:
            raise ValueError(
                "Target ubiquity table is missing required columns: "
                + ", ".join(missing)
                + f". Available columns: {', '.join(fieldnames)}"
            )
        weights: dict[str, float] = {}
        n_rows = 0
        n_missing = 0
        n_non_numeric = 0
        families: set[str] = set()
        for row in reader:
            n_rows += 1
            gene_symbol = str(row.get(gene_symbol_column, "")).strip().upper()
            weight = _parse_float(row.get(idf_column))
            if not gene_symbol:
                n_missing += 1
                continue
            if weight is None:
                n_non_numeric += 1
                continue
            family = str(row.get(family_column, "")).strip()
            if family:
                families.add(family)
            weights[gene_symbol] = float(weight)
    return weights, {
        "path": str(path),
        "n_rows": n_rows,
        "n_genes": len(weights),
        "n_missing_required": n_missing,
        "n_non_numeric_idf": n_non_numeric,
        "n_families": len(families),
        "columns": {
            "gene_symbol": gene_symbol_column,
            "idf": idf_column,
            "family": family_column if family_column in fieldnames else None,
        },
    }


def load_compound_qc_tsv(
    path: str | Path,
    *,
    drug_id_column: str = "drug_id",
    recommended_use_column: str = "recommended_use",
    pan_toxic_column: str = "pan_toxic_flag",
    polypharm_column: str = "polypharm_flag",
    reason_column: str = "blacklist_reason",
    delimiter: str = "\t",
) -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    p = Path(path)
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Compound QC table has no header: {p}")
        fieldnames = [str(f) for f in reader.fieldnames]
        required = [drug_id_column]
        missing = [name for name in required if name not in fieldnames]
        if missing:
            raise ValueError(
                "Compound QC table is missing required columns: "
                + ", ".join(missing)
                + f". Available columns: {', '.join(fieldnames)}"
            )
        rows_by_drug: dict[str, dict[str, object]] = {}
        n_rows = 0
        n_missing = 0
        recommended_counts: dict[str, int] = {}
        for row in reader:
            n_rows += 1
            drug_id = str(row.get(drug_id_column, "")).strip()
            if not drug_id:
                n_missing += 1
                continue
            recommended_use = str(row.get(recommended_use_column, "")).strip().lower() or "keep"
            pan_toxic = str(row.get(pan_toxic_column, "")).strip().lower() in {"1", "true", "t", "yes", "y"}
            polypharm = str(row.get(polypharm_column, "")).strip().lower() in {"1", "true", "t", "yes", "y"}
            reason = str(row.get(reason_column, "")).strip()
            rows_by_drug[drug_id] = {
                "recommended_use": recommended_use,
                "pan_toxic_flag": pan_toxic,
                "polypharm_flag": polypharm,
                "blacklist_reason": reason,
            }
            recommended_counts[recommended_use] = int(recommended_counts.get(recommended_use, 0)) + 1
    return rows_by_drug, {
        "path": str(path),
        "n_rows": n_rows,
        "n_drugs": len(rows_by_drug),
        "n_missing_required": n_missing,
        "recommended_use_counts": recommended_counts,
        "columns": {
            "drug_id": drug_id_column,
            "recommended_use": recommended_use_column if recommended_use_column in fieldnames else None,
            "pan_toxic_flag": pan_toxic_column if pan_toxic_column in fieldnames else None,
            "polypharm_flag": polypharm_column if polypharm_column in fieldnames else None,
            "blacklist_reason": reason_column if reason_column in fieldnames else None,
        },
    }


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


def summarize_targets_per_drug(drug_targets: dict[str, dict[str, float]]) -> dict[str, object]:
    counts = [len(genes) for genes in drug_targets.values()]
    if not counts:
        return {"n_drugs": 0}
    counts_sorted = sorted(int(x) for x in counts)
    p90_index = max(0, min(len(counts_sorted) - 1, int(round(0.9 * (len(counts_sorted) - 1)))))
    return {
        "n_drugs": len(counts_sorted),
        "min": int(counts_sorted[0]),
        "median": float(statistics.median(counts_sorted)),
        "p90": int(counts_sorted[p90_index]),
        "max": int(counts_sorted[-1]),
    }


def filter_promiscuous_targets(
    *,
    drug_targets: dict[str, dict[str, float]],
    max_targets_per_drug: int,
    policy: str,
) -> tuple[dict[str, dict[str, float]], dict[str, object], list[str]]:
    limit = max(1, int(max_targets_per_drug))
    mode = str(policy).strip().lower()
    if mode not in {"warn", "drop", "cap"}:
        raise ValueError(f"Unsupported target_promiscuity_policy: {policy}")

    out: dict[str, dict[str, float]] = {}
    n_over_limit = 0
    n_dropped = 0
    n_capped = 0
    warned: list[str] = []
    for drug_id, targets in drug_targets.items():
        n_targets = len(targets)
        if n_targets <= limit:
            out[drug_id] = dict(targets)
            continue
        n_over_limit += 1
        warned.append(
            f"warning: promiscuous target list for drug={drug_id}; n_targets={n_targets} exceeds max_targets_per_drug={limit}"
        )
        if mode == "warn":
            out[drug_id] = dict(targets)
            continue
        if mode == "drop":
            n_dropped += 1
            continue
        ranked = sorted(targets.items(), key=lambda kv: (-float(kv[1]), str(kv[0])))[:limit]
        total = sum(float(w) for _, w in ranked)
        if total > 0:
            out[drug_id] = {gene: float(weight) / float(total) for gene, weight in ranked}
        else:
            out[drug_id] = {gene: 1.0 / float(len(ranked)) for gene, _ in ranked}
        n_capped += 1
        warned.append(
            f"warning: capped targets for drug={drug_id} to top {limit} genes by weight."
        )

    summary = {
        "max_targets_per_drug": limit,
        "policy": mode,
        "n_drugs_input": len(drug_targets),
        "n_drugs_output": len(out),
        "n_over_limit": n_over_limit,
        "n_dropped": n_dropped,
        "n_capped": n_capped,
    }
    return out, summary, warned
