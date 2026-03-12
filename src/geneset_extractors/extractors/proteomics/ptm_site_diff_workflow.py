from __future__ import annotations

import csv
from dataclasses import dataclass
import gzip
import math
from pathlib import Path
import sys

from geneset_extractors.core.gmt import (
    build_gmt_sets_from_rows,
    parse_int_list_csv,
    parse_mass_list_csv,
    resolve_gmt_out_path,
    write_gmt,
)
from geneset_extractors.core.metadata import make_metadata, write_metadata
from geneset_extractors.core.qc import write_run_summary_files
from geneset_extractors.core.selection import (
    global_l1_weights,
    ranked_gene_ids,
    select_quantile,
    select_threshold,
    select_top_k,
    within_set_l1_weights,
)
from geneset_extractors.extractors.rnaseq.deg_scoring import sanitize_name_component


DEFAULT_COLUMNS = {
    "site_id": ("site_id",),
    "site_group": ("site_group_id", "site_group"),
    "gene_id": ("gene_id",),
    "gene_symbol": ("gene_symbol", "gene_name"),
    "protein_accession": ("protein_accession", "accession", "uniprot_accession"),
    "residue": ("residue",),
    "position": ("position",),
    "stat": ("stat",),
    "logfc": ("log2fc", "logFC", "log2_fold_change"),
    "padj": ("padj", "FDR", "adj.P.Val"),
    "pvalue": ("pvalue", "PValue", "P.Value"),
    "localization_prob": ("localization_prob", "site_localization_prob"),
    "peptide_count": ("peptide_count", "n_peptides", "peptides"),
    "protein_logfc": ("protein_log2fc", "protein_logFC", "protein_logfc"),
    "protein_stat": ("protein_stat",),
}


@dataclass(frozen=True)
class PTMRow:
    line_no: int
    values: dict[str, str]


@dataclass
class PTMSiteRecord:
    canonical_site_key: str
    gene_id: str
    gene_symbol: str
    protein_accession: str
    residue: str
    position: str
    ptm_type: str
    base_score: float
    adjusted_score: float
    final_score: float
    confidence_weight: float
    ubiquity_weight: float
    protein_adjustment_value: float | None
    quality_value: float
    source_line_no: int


@dataclass
class PTMWorkflowConfig:
    converter_name: str
    out_dir: Path
    organism: str
    genome_build: str
    dataset_label: str
    signature_name: str
    ptm_type: str
    site_id_column: str | None
    site_group_column: str | None
    gene_id_column: str | None
    gene_symbol_column: str | None
    protein_accession_column: str | None
    residue_column: str | None
    position_column: str | None
    score_column: str | None
    stat_column: str | None
    logfc_column: str | None
    padj_column: str | None
    pvalue_column: str | None
    localization_prob_column: str | None
    peptide_count_column: str | None
    protein_logfc_column: str | None
    protein_stat_column: str | None
    score_mode: str
    score_transform: str
    protein_adjustment: str
    protein_adjustment_lambda: float
    confidence_weight_mode: str
    min_localization_prob: float
    site_dup_policy: str
    gene_aggregation: str
    gene_topk_sites: int
    ambiguous_gene_policy: str
    select: str
    top_k: int
    quantile: float
    min_score: float
    normalize: str
    emit_full: bool
    emit_gmt: bool
    gmt_out: str | None
    gmt_prefer_symbol: bool
    gmt_require_symbol: bool
    gmt_biotype_allowlist: str
    gmt_min_genes: int
    gmt_max_genes: int
    gmt_topk_list: str
    gmt_mass_list: str
    gmt_split_signed: bool
    emit_small_gene_sets: bool
    neglog10p_cap: float
    neglog10p_eps: float


@dataclass
class AliasEntry:
    canonical_site_key: str
    gene_id: str
    gene_symbol: str
    protein_accession: str
    residue: str
    position: str
    ptm_type: str


@dataclass
class UbiquityEntry:
    canonical_site_key: str
    idf_ref: float
    df_ref: float | None
    n_samples_ref: float | None
    gene_id: str
    gene_symbol: str
    ptm_type: str



def _open_text(path: str | Path):
    p = Path(path)
    if p.suffix.lower() == ".gz":
        return gzip.open(p, "rt", encoding="utf-8")
    return p.open("r", encoding="utf-8")



def read_tsv_rows(path: str | Path) -> tuple[list[str], list[PTMRow]]:
    rows: list[PTMRow] = []
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"PTM table has no header: {path}")
        fieldnames = [str(x) for x in reader.fieldnames]
        for idx, row in enumerate(reader, start=2):
            rows.append(PTMRow(line_no=idx, values={k: str(v) for k, v in row.items() if k is not None}))
    return fieldnames, rows



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



def _resolve_column(fieldnames: list[str], explicit: str | None, logical_name: str, required: bool = False) -> str | None:
    name = _clean(explicit) or None
    if name is not None:
        if name not in fieldnames:
            raise ValueError(
                f"Column '{name}' for {logical_name} not found. Available columns: {', '.join(fieldnames)}"
            )
        return name
    for candidate in DEFAULT_COLUMNS.get(logical_name, ()):  # pragma: no branch
        if candidate in fieldnames:
            return candidate
    if required:
        raise ValueError(
            f"Could not resolve required column for {logical_name}. Available columns: {', '.join(fieldnames)}"
        )
    return None



def resolve_ptm_columns(fieldnames: list[str], cfg: PTMWorkflowConfig) -> dict[str, str | None]:
    return {
        "site_id": _resolve_column(fieldnames, cfg.site_id_column, "site_id", required=False),
        "site_group": _resolve_column(fieldnames, cfg.site_group_column, "site_group", required=False),
        "gene_id": _resolve_column(fieldnames, cfg.gene_id_column, "gene_id", required=False),
        "gene_symbol": _resolve_column(fieldnames, cfg.gene_symbol_column, "gene_symbol", required=False),
        "protein_accession": _resolve_column(fieldnames, cfg.protein_accession_column, "protein_accession", required=False),
        "residue": _resolve_column(fieldnames, cfg.residue_column, "residue", required=False),
        "position": _resolve_column(fieldnames, cfg.position_column, "position", required=False),
        "score": _clean(cfg.score_column) or None,
        "stat": _resolve_column(fieldnames, cfg.stat_column, "stat", required=False),
        "logfc": _resolve_column(fieldnames, cfg.logfc_column, "logfc", required=False),
        "padj": _resolve_column(fieldnames, cfg.padj_column, "padj", required=False),
        "pvalue": _resolve_column(fieldnames, cfg.pvalue_column, "pvalue", required=False),
        "localization_prob": _resolve_column(fieldnames, cfg.localization_prob_column, "localization_prob", required=False),
        "peptide_count": _resolve_column(fieldnames, cfg.peptide_count_column, "peptide_count", required=False),
        "protein_logfc": _resolve_column(fieldnames, cfg.protein_logfc_column, "protein_logfc", required=False),
        "protein_stat": _resolve_column(fieldnames, cfg.protein_stat_column, "protein_stat", required=False),
    }



def read_site_alias_tsv(path: str | Path) -> dict[str, AliasEntry]:
    out: dict[str, AliasEntry] = {}
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Alias TSV has no header: {path}")
        required = {"input_site_key", "canonical_site_key"}
        missing = required.difference(set(reader.fieldnames))
        if missing:
            raise ValueError(f"Alias TSV missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            input_key = _clean(row.get("input_site_key"))
            canonical = _clean(row.get("canonical_site_key"))
            if not input_key or not canonical:
                continue
            out[input_key] = AliasEntry(
                canonical_site_key=canonical,
                gene_id=_clean(row.get("gene_id")),
                gene_symbol=_clean(row.get("gene_symbol")),
                protein_accession=_clean(row.get("protein_accession")),
                residue=_clean(row.get("residue")),
                position=_clean(row.get("position")),
                ptm_type=_clean(row.get("ptm_type")) or "generic",
            )
    return out



def read_site_ubiquity_tsv(path: str | Path) -> dict[str, UbiquityEntry]:
    out: dict[str, UbiquityEntry] = {}
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Ubiquity TSV has no header: {path}")
        if "canonical_site_key" not in reader.fieldnames:
            raise ValueError("Ubiquity TSV missing canonical_site_key column")
        for row in reader:
            canonical = _clean(row.get("canonical_site_key"))
            if not canonical:
                continue
            idf = _parse_float(row.get("idf_ref"))
            if idf is None:
                df_ref = _parse_float(row.get("df_ref"))
                n_ref = _parse_float(row.get("n_samples_ref"))
                if df_ref is not None and n_ref is not None and n_ref >= 0:
                    idf = math.log((n_ref + 1.0) / (df_ref + 1.0))
                else:
                    idf = 1.0
            out[canonical] = UbiquityEntry(
                canonical_site_key=canonical,
                idf_ref=float(idf),
                df_ref=_parse_float(row.get("df_ref")),
                n_samples_ref=_parse_float(row.get("n_samples_ref")),
                gene_id=_clean(row.get("gene_id")),
                gene_symbol=_clean(row.get("gene_symbol")),
                ptm_type=_clean(row.get("ptm_type")) or "generic",
            )
    return out



def _site_candidates(row: PTMRow, columns: dict[str, str | None], ptm_type: str) -> list[str]:
    candidates: list[str] = []
    for key in ("site_id", "site_group"):
        col = columns.get(key)
        if col:
            value = _clean(row.values.get(col))
            if value:
                candidates.append(value)
    protein = _clean(row.values.get(columns.get("protein_accession", ""))) if columns.get("protein_accession") else ""
    residue = _clean(row.values.get(columns.get("residue", ""))) if columns.get("residue") else ""
    position = _clean(row.values.get(columns.get("position", ""))) if columns.get("position") else ""
    gene_symbol = _clean(row.values.get(columns.get("gene_symbol", ""))) if columns.get("gene_symbol") else ""
    if protein and residue and position:
        candidates.append(f"{protein}|{residue}|{position}|{ptm_type}")
    if gene_symbol and residue and position:
        candidates.append(f"{gene_symbol}|{residue}|{position}|{ptm_type}")
    return candidates



def _default_canonical_site_key(row: PTMRow, columns: dict[str, str | None], ptm_type: str) -> str | None:
    protein = _clean(row.values.get(columns.get("protein_accession", ""))) if columns.get("protein_accession") else ""
    residue = _clean(row.values.get(columns.get("residue", ""))) if columns.get("residue") else ""
    position = _clean(row.values.get(columns.get("position", ""))) if columns.get("position") else ""
    if protein and residue and position:
        return f"{protein}|{residue}|{position}|{ptm_type}"
    for key in ("site_id", "site_group"):
        col = columns.get(key)
        if col:
            value = _clean(row.values.get(col))
            if value:
                return f"{value}|{ptm_type}"
    return None



def _split_gene_tokens(value: str) -> list[str]:
    text = _clean(value)
    if not text:
        return []
    for delim in (";", ",", "|", "/"):
        text = text.replace(delim, ";")
    return [tok.strip() for tok in text.split(";") if tok.strip()]



def _resolve_gene_assignments(
    *,
    raw_gene_id: str,
    raw_gene_symbol: str,
    alias_gene_id: str,
    alias_gene_symbol: str,
    policy: str,
) -> tuple[list[tuple[str, str, float]], int]:
    gene_ids = _split_gene_tokens(alias_gene_id or raw_gene_id)
    gene_symbols = _split_gene_tokens(alias_gene_symbol or raw_gene_symbol)
    n = max(len(gene_ids), len(gene_symbols))
    if n == 0:
        return [], 0
    pairs: list[tuple[str, str]] = []
    for idx in range(n):
        gid = gene_ids[idx] if idx < len(gene_ids) else ""
        gsym = gene_symbols[idx] if idx < len(gene_symbols) else ""
        if not gid and gsym:
            gid = gsym
        if not gsym and gid:
            gsym = gid
        if gid:
            pairs.append((gid, gsym))
    if not pairs:
        return [], 0
    ambiguous = 1 if len(pairs) > 1 else 0
    if len(pairs) == 1:
        gid, gsym = pairs[0]
        return [(gid, gsym, 1.0)], ambiguous
    if policy == "drop":
        return [], ambiguous
    if policy == "first":
        gid, gsym = pairs[0]
        return [(gid, gsym, 1.0)], ambiguous
    if policy == "split_equal":
        frac = 1.0 / float(len(pairs))
        return [(gid, gsym, frac) for gid, gsym in pairs], ambiguous
    raise ValueError(f"Unsupported ambiguous_gene_policy: {policy}")



def _resolve_score_mode(columns: dict[str, str | None], rows: list[PTMRow], requested: str) -> str:
    if requested != "auto":
        return requested
    stat_col = columns.get("stat")
    if stat_col:
        for row in rows:
            if _parse_float(row.values.get(stat_col)) is not None:
                return "stat"
    logfc_col = columns.get("logfc")
    p_col = columns.get("padj") or columns.get("pvalue")
    if logfc_col and p_col:
        for row in rows:
            if _parse_float(row.values.get(logfc_col)) is not None and _parse_float(row.values.get(p_col)) is not None:
                return "logfc_times_neglog10p"
    score_col = columns.get("score")
    if score_col and rows and score_col in rows[0].values:
        return "custom_column"
    raise ValueError(
        "score_mode=auto could not resolve a usable score definition. Provide a parseable stat column, "
        "or logfc+pvalue/padj columns, or set --score_mode custom_column with --score_column."
    )



def _compute_base_score(
    row: PTMRow,
    columns: dict[str, str | None],
    score_mode: str,
    neglog10p_eps: float,
    neglog10p_cap: float,
) -> tuple[float | None, str]:
    if score_mode == "custom_column":
        column = columns.get("score")
        if not column:
            raise ValueError("score_mode=custom_column requires --score_column")
        return _parse_float(row.values.get(column)), "custom_column"
    if score_mode == "stat":
        column = columns.get("stat")
        if not column:
            raise ValueError("score_mode=stat requires a stat column")
        return _parse_float(row.values.get(column)), "stat"
    if score_mode == "logfc_times_neglog10p":
        logfc_col = columns.get("logfc")
        p_col = columns.get("padj") or columns.get("pvalue")
        if not logfc_col or not p_col:
            raise ValueError("score_mode=logfc_times_neglog10p requires logfc and pvalue/padj columns")
        logfc = _parse_float(row.values.get(logfc_col))
        pval = _parse_float(row.values.get(p_col))
        if logfc is None or pval is None:
            return None, "logfc_times_neglog10p"
        neglog = min(-math.log10(max(float(pval), float(neglog10p_eps))), float(neglog10p_cap))
        return float(logfc) * float(neglog), "logfc_times_neglog10p"
    raise ValueError(f"Unsupported score_mode: {score_mode}")



def _protein_value(row: PTMRow, columns: dict[str, str | None]) -> float | None:
    for key in ("protein_logfc", "protein_stat"):
        col = columns.get(key)
        if col:
            value = _parse_float(row.values.get(col))
            if value is not None:
                return float(value)
    return None



def _confidence_weight(
    row: PTMRow,
    *,
    columns: dict[str, str | None],
    mode: str,
    min_localization_prob: float,
    neglog10p_eps: float,
    neglog10p_cap: float,
    peptide_count_max: float,
) -> tuple[float, dict[str, float]]:
    p_col = columns.get("padj") or columns.get("pvalue")
    p_weight = 1.0
    if p_col:
        pval = _parse_float(row.values.get(p_col))
        if pval is not None:
            p_weight = min(-math.log10(max(float(pval), float(neglog10p_eps))), float(neglog10p_cap)) / float(neglog10p_cap)
    loc_weight = 1.0
    loc_col = columns.get("localization_prob")
    if loc_col:
        loc = _parse_float(row.values.get(loc_col))
        if loc is not None:
            if float(loc) < float(min_localization_prob):
                loc_weight = 0.0
            else:
                denom = max(1e-12, 1.0 - float(min_localization_prob))
                loc_weight = max(0.0, min(1.0, (float(loc) - float(min_localization_prob)) / denom))
    pep_weight = 1.0
    pep_col = columns.get("peptide_count")
    if pep_col and peptide_count_max > 0:
        pep = _parse_float(row.values.get(pep_col))
        if pep is not None and pep >= 0:
            pep_weight = min(1.0, math.log1p(float(pep)) / math.log1p(float(peptide_count_max)))
    if mode == "none":
        weight = 1.0
    elif mode == "pvalue":
        weight = p_weight
    elif mode == "localization":
        weight = loc_weight
    elif mode == "combined":
        weight = p_weight * loc_weight * pep_weight
    else:
        raise ValueError(f"Unsupported confidence_weight_mode: {mode}")
    return float(weight), {
        "pvalue_weight": float(p_weight),
        "localization_weight": float(loc_weight),
        "peptide_weight": float(pep_weight),
    }



def _apply_score_transform(value: float, mode: str) -> float:
    x = float(value)
    if mode == "signed":
        return x
    if mode == "abs":
        return abs(x)
    if mode == "positive":
        return max(x, 0.0)
    if mode == "negative":
        return max(-x, 0.0)
    raise ValueError(f"Unsupported score_transform: {mode}")



def _collapse_duplicate_sites(records: list[PTMSiteRecord], policy: str) -> list[PTMSiteRecord]:
    grouped: dict[tuple[str, str], list[PTMSiteRecord]] = {}
    for rec in records:
        grouped.setdefault((rec.canonical_site_key, rec.gene_id), []).append(rec)
    out: list[PTMSiteRecord] = []
    for _key, items in grouped.items():
        if len(items) == 1:
            out.append(items[0])
            continue
        if policy == "highest_confidence":
            out.append(max(items, key=lambda r: (float(r.quality_value), abs(float(r.final_score)), -int(r.source_line_no))))
            continue
        if policy == "max_abs":
            out.append(max(items, key=lambda r: (abs(float(r.final_score)), float(r.quality_value), -int(r.source_line_no))))
            continue
        if policy in {"mean", "sum"}:
            denom = float(len(items)) if policy == "mean" else 1.0
            first = items[0]
            out.append(
                PTMSiteRecord(
                    canonical_site_key=first.canonical_site_key,
                    gene_id=first.gene_id,
                    gene_symbol=first.gene_symbol,
                    protein_accession=first.protein_accession,
                    residue=first.residue,
                    position=first.position,
                    ptm_type=first.ptm_type,
                    base_score=sum(float(x.base_score) for x in items) / denom,
                    adjusted_score=sum(float(x.adjusted_score) for x in items) / denom,
                    final_score=sum(float(x.final_score) for x in items) / denom,
                    confidence_weight=sum(float(x.confidence_weight) for x in items) / float(len(items)),
                    ubiquity_weight=sum(float(x.ubiquity_weight) for x in items) / float(len(items)),
                    protein_adjustment_value=(sum(float(x.protein_adjustment_value or 0.0) for x in items) / float(len(items))),
                    quality_value=max(float(x.quality_value) for x in items),
                    source_line_no=min(int(x.source_line_no) for x in items),
                )
            )
            continue
        raise ValueError(f"Unsupported site_dup_policy: {policy}")
    return out



def _aggregate_sites_to_genes(records: list[PTMSiteRecord], method: str, topk_sites: int) -> tuple[dict[str, float], dict[str, str], dict[str, list[PTMSiteRecord]]]:
    by_gene: dict[str, list[PTMSiteRecord]] = {}
    gene_symbols: dict[str, str] = {}
    for rec in records:
        by_gene.setdefault(rec.gene_id, []).append(rec)
        if rec.gene_symbol:
            gene_symbols.setdefault(rec.gene_id, rec.gene_symbol)
    scores: dict[str, float] = {}
    for gene_id, items in by_gene.items():
        ranked = sorted(items, key=lambda r: (-abs(float(r.final_score)), str(r.canonical_site_key)))
        if method == "signed_topk_mean":
            kept = ranked[: max(1, int(topk_sites))]
            scores[gene_id] = sum(float(r.final_score) for r in kept) / float(len(kept))
        elif method == "max_abs":
            best = ranked[0]
            scores[gene_id] = float(best.final_score)
        elif method == "sum":
            scores[gene_id] = sum(float(r.final_score) for r in items)
        elif method == "mean":
            scores[gene_id] = sum(float(r.final_score) for r in items) / float(len(items))
        else:
            raise ValueError(f"Unsupported gene_aggregation: {method}")
    return scores, gene_symbols, by_gene



def _select_gene_ids(scores: dict[str, float], cfg: PTMWorkflowConfig) -> list[str]:
    if cfg.select == "none":
        return ranked_gene_ids(scores)
    if cfg.select == "top_k":
        return select_top_k(scores, int(cfg.top_k))
    if cfg.select == "quantile":
        return select_quantile(scores, float(cfg.quantile))
    if cfg.select == "threshold":
        return select_threshold(scores, float(cfg.min_score))
    raise ValueError(f"Unsupported selection method: {cfg.select}")



def _selected_weights(magnitude_scores: dict[str, float], selected_gene_ids: list[str], normalize: str) -> dict[str, float]:
    if normalize == "none":
        return {g: float(magnitude_scores.get(g, 0.0)) for g in selected_gene_ids}
    if normalize == "l1":
        global_weights = global_l1_weights(magnitude_scores)
        return {g: float(global_weights.get(g, 0.0)) for g in selected_gene_ids}
    if normalize == "within_set_l1":
        return within_set_l1_weights(magnitude_scores, selected_gene_ids)
    raise ValueError(f"Unsupported normalization method: {normalize}")



def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["gene_id", "score"]
    if any("weight" in row for row in rows):
        fieldnames.append("weight")
    if any(str(row.get("gene_symbol", "")).strip() for row in rows):
        fieldnames.append("gene_symbol")
    if any("rank" in row for row in rows):
        fieldnames.append("rank")
    if any("n_supporting_sites" in row for row in rows):
        fieldnames.append("n_supporting_sites")
    if any(str(row.get("top_sites", "")).strip() for row in rows):
        fieldnames.append("top_sites")
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)



def _warn_gmt_diagnostics(gmt_diagnostics: list[dict[str, object]]) -> None:
    for diag in gmt_diagnostics:
        code = str(diag.get("code", ""))
        base_name = str(diag.get("base_name", "unknown"))
        n_genes = int(diag.get("n_genes", 0) or 0)
        reason = str(diag.get("reason", "")).strip()
        suggestion = str(diag.get("suggestion", "")).strip()
        if code in {"small_gene_set_skipped", "no_positive_genes"}:
            print(f"warning: skipped GMT output for {base_name}; n_genes={n_genes}. {reason}", file=sys.stderr)
        elif code == "small_gene_set_emitted":
            print(f"warning: emitted small GMT output for {base_name}; n_genes={n_genes}. {reason}", file=sys.stderr)
        elif code == "marginal_gene_count":
            print(f"warning: marginal GMT signal for {base_name}; n_genes={n_genes}. {reason}", file=sys.stderr)
        else:
            continue
        if suggestion:
            print(f"warning:   {suggestion}", file=sys.stderr)



def _collect_skipped_gmt_outputs(gmt_diagnostics: list[dict[str, object]]) -> list[dict[str, object]]:
    skipped: list[dict[str, object]] = []
    for diag in gmt_diagnostics:
        code = str(diag.get("code", ""))
        if code not in {"small_gene_set_skipped", "no_positive_genes"}:
            continue
        skipped.append(
            {
                "reason": str(diag.get("reason", "")),
                "code": code,
                "n_genes": int(diag.get("n_genes", 0) or 0),
                "variant": str(diag.get("variant", "primary")),
            }
        )
    return skipped



def run_ptm_site_diff_workflow(
    *,
    cfg: PTMWorkflowConfig,
    fieldnames: list[str],
    rows: list[PTMRow],
    alias_map: dict[str, AliasEntry] | None,
    ubiquity_map: dict[str, UbiquityEntry] | None,
    input_files: list[dict[str, str]],
    resources_info: dict[str, object] | None,
) -> dict[str, object]:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    columns = resolve_ptm_columns(fieldnames, cfg)
    resolved_score_mode = _resolve_score_mode(columns, rows, cfg.score_mode)

    peptide_count_max = 0.0
    pep_col = columns.get("peptide_count")
    if pep_col:
        for row in rows:
            pep = _parse_float(row.values.get(pep_col))
            if pep is not None and pep > peptide_count_max:
                peptide_count_max = float(pep)

    base_rows: list[dict[str, object]] = []
    protein_pairs: list[tuple[float, float]] = []
    parse_summary = {
        "n_raw_rows": len(rows),
        "n_rows_dropped_missing_score": 0,
        "n_rows_unresolved_site": 0,
        "n_rows_ambiguous_gene": 0,
        "n_rows_dropped_ambiguous_gene": 0,
        "n_rows_alias_matched": 0,
        "n_sites_matched_to_ubiquity_prior": 0,
        "n_sites_with_protein_adjustment": 0,
    }

    for row in rows:
        base_score, _used_mode = _compute_base_score(
            row,
            columns,
            resolved_score_mode,
            cfg.neglog10p_eps,
            cfg.neglog10p_cap,
        )
        if base_score is None:
            parse_summary["n_rows_dropped_missing_score"] = int(parse_summary["n_rows_dropped_missing_score"]) + 1
            continue
        alias_entry = None
        for candidate in _site_candidates(row, columns, cfg.ptm_type):
            if alias_map and candidate in alias_map:
                alias_entry = alias_map[candidate]
                break
        if alias_entry is not None:
            parse_summary["n_rows_alias_matched"] = int(parse_summary["n_rows_alias_matched"]) + 1
        canonical_site_key = alias_entry.canonical_site_key if alias_entry is not None else _default_canonical_site_key(row, columns, cfg.ptm_type)
        if not canonical_site_key:
            parse_summary["n_rows_unresolved_site"] = int(parse_summary["n_rows_unresolved_site"]) + 1
            continue

        raw_gene_id = _clean(row.values.get(columns.get("gene_id", ""))) if columns.get("gene_id") else ""
        raw_gene_symbol = _clean(row.values.get(columns.get("gene_symbol", ""))) if columns.get("gene_symbol") else ""
        gene_assignments, ambiguous_flag = _resolve_gene_assignments(
            raw_gene_id=raw_gene_id,
            raw_gene_symbol=raw_gene_symbol,
            alias_gene_id=alias_entry.gene_id if alias_entry is not None else "",
            alias_gene_symbol=alias_entry.gene_symbol if alias_entry is not None else "",
            policy=cfg.ambiguous_gene_policy,
        )
        if ambiguous_flag:
            parse_summary["n_rows_ambiguous_gene"] = int(parse_summary["n_rows_ambiguous_gene"]) + 1
        if ambiguous_flag and not gene_assignments:
            parse_summary["n_rows_dropped_ambiguous_gene"] = int(parse_summary["n_rows_dropped_ambiguous_gene"]) + 1
            continue
        if not gene_assignments:
            continue

        protein_value = _protein_value(row, columns)
        if protein_value is not None:
            for gene_id, _gene_symbol, _weight_frac in gene_assignments:
                protein_pairs.append((float(base_score), float(protein_value)))
                parse_summary["n_sites_with_protein_adjustment"] = int(parse_summary["n_sites_with_protein_adjustment"]) + 1
                break

        confidence_weight, confidence_parts = _confidence_weight(
            row,
            columns=columns,
            mode=cfg.confidence_weight_mode,
            min_localization_prob=cfg.min_localization_prob,
            neglog10p_eps=cfg.neglog10p_eps,
            neglog10p_cap=cfg.neglog10p_cap,
            peptide_count_max=peptide_count_max,
        )
        base_rows.append(
            {
                "row": row,
                "canonical_site_key": canonical_site_key,
                "alias_entry": alias_entry,
                "gene_assignments": gene_assignments,
                "base_score": float(base_score),
                "protein_value": protein_value,
                "confidence_weight": float(confidence_weight),
                "confidence_parts": confidence_parts,
            }
        )

    beta = 1.0
    if cfg.protein_adjustment == "subtract":
        beta = float(cfg.protein_adjustment_lambda)
    elif cfg.protein_adjustment == "residual" and protein_pairs:
        numerator = sum(x * y for x, y in protein_pairs)
        denominator = sum(y * y for _x, y in protein_pairs)
        beta = numerator / denominator if denominator > 0 else 0.0
    elif cfg.protein_adjustment == "none":
        beta = 0.0
    elif cfg.protein_adjustment not in {"subtract", "residual", "none"}:
        raise ValueError(f"Unsupported protein_adjustment: {cfg.protein_adjustment}")

    site_records: list[PTMSiteRecord] = []
    for entry in base_rows:
        canonical_site_key = str(entry["canonical_site_key"])
        alias_entry = entry["alias_entry"]
        base_score = float(entry["base_score"])
        protein_value = entry["protein_value"]
        adjusted_score = float(base_score)
        if protein_value is not None and cfg.protein_adjustment != "none":
            adjusted_score = float(base_score) - float(beta) * float(protein_value)
        transformed = _apply_score_transform(adjusted_score, cfg.score_transform)
        ubiquity_weight = 1.0
        ubiquity_entry = ubiquity_map.get(canonical_site_key) if ubiquity_map else None
        if ubiquity_entry is not None:
            ubiquity_weight = float(ubiquity_entry.idf_ref)
            parse_summary["n_sites_matched_to_ubiquity_prior"] = int(parse_summary["n_sites_matched_to_ubiquity_prior"]) + 1
        final_score = float(transformed) * float(entry["confidence_weight"]) * float(ubiquity_weight)
        row = entry["row"]
        for gene_id, gene_symbol, gene_weight in entry["gene_assignments"]:
            final_gene_score = final_score * float(gene_weight)
            site_records.append(
                PTMSiteRecord(
                    canonical_site_key=canonical_site_key,
                    gene_id=str(gene_id),
                    gene_symbol=str(gene_symbol),
                    protein_accession=(alias_entry.protein_accession if alias_entry is not None else _clean(row.values.get(columns.get("protein_accession", ""))) if columns.get("protein_accession") else ""),
                    residue=(alias_entry.residue if alias_entry is not None else _clean(row.values.get(columns.get("residue", ""))) if columns.get("residue") else ""),
                    position=(alias_entry.position if alias_entry is not None else _clean(row.values.get(columns.get("position", ""))) if columns.get("position") else ""),
                    ptm_type=(alias_entry.ptm_type if alias_entry is not None else cfg.ptm_type),
                    base_score=float(base_score),
                    adjusted_score=float(adjusted_score),
                    final_score=float(final_gene_score),
                    confidence_weight=float(entry["confidence_weight"]),
                    ubiquity_weight=float(ubiquity_weight),
                    protein_adjustment_value=(float(protein_value) if protein_value is not None else None),
                    quality_value=float(entry["confidence_weight"]) * max(1e-6, abs(float(adjusted_score))),
                    source_line_no=int(row.line_no),
                )
            )

    collapsed_sites = _collapse_duplicate_sites(site_records, cfg.site_dup_policy)
    gene_scores, gene_symbols, by_gene = _aggregate_sites_to_genes(collapsed_sites, cfg.gene_aggregation, cfg.gene_topk_sites)
    magnitude = {gene_id: abs(float(score)) for gene_id, score in gene_scores.items()}
    selected_gene_ids = _select_gene_ids(magnitude, cfg)
    weights = _selected_weights(magnitude, selected_gene_ids, cfg.normalize)

    full_ranked_ids = sorted(gene_scores, key=lambda g: (-abs(float(gene_scores[g])), str(g)))
    full_rows: list[dict[str, object]] = []
    for idx, gene_id in enumerate(full_ranked_ids, start=1):
        support = sorted(by_gene.get(gene_id, []), key=lambda r: (-abs(float(r.final_score)), str(r.canonical_site_key)))
        top_sites = ";".join(str(r.canonical_site_key) for r in support[: min(3, len(support))])
        full_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene_symbols.get(gene_id, gene_id),
                "score": float(gene_scores[gene_id]),
                "rank": idx,
                "n_supporting_sites": len(support),
                "top_sites": top_sites,
            }
        )

    selected_rows: list[dict[str, object]] = []
    for idx, gene_id in enumerate(selected_gene_ids, start=1):
        support = sorted(by_gene.get(gene_id, []), key=lambda r: (-abs(float(r.final_score)), str(r.canonical_site_key)))
        top_sites = ";".join(str(r.canonical_site_key) for r in support[: min(3, len(support))])
        selected_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene_symbols.get(gene_id, gene_id),
                "score": float(gene_scores.get(gene_id, 0.0)),
                "weight": float(weights.get(gene_id, 0.0)),
                "rank": idx,
                "n_supporting_sites": len(support),
                "top_sites": top_sites,
            }
        )

    _write_rows(out_dir / "geneset.tsv", selected_rows)
    output_files = [{"path": str(out_dir / "geneset.tsv"), "role": "selected_gene_program"}]
    if cfg.emit_full:
        _write_rows(out_dir / "geneset.full.tsv", full_rows)
        output_files.append({"path": str(out_dir / "geneset.full.tsv"), "role": "full_scores"})

    gmt_diagnostics: list[dict[str, object]] = []
    gmt_sets: list[tuple[str, list[str]]] = []
    gmt_plans: list[dict[str, object]] = []
    if cfg.emit_gmt:
        gmt_source_rows = full_rows
        base_name = sanitize_name_component(
            f"{cfg.converter_name}__signature={cfg.signature_name}__ptm_type={cfg.ptm_type}__score_mode={resolved_score_mode}"
        )
        gmt_sets, gmt_plans = build_gmt_sets_from_rows(
            gmt_source_rows,
            base_name=base_name,
            prefer_symbol=bool(cfg.gmt_prefer_symbol),
            min_genes=int(cfg.gmt_min_genes),
            max_genes=int(cfg.gmt_max_genes),
            topk_list=parse_int_list_csv(cfg.gmt_topk_list),
            mass_list=parse_mass_list_csv(cfg.gmt_mass_list),
            split_signed=bool(cfg.gmt_split_signed),
            require_symbol=bool(cfg.gmt_require_symbol),
            allowed_biotypes=None,
            emit_small_gene_sets=bool(cfg.emit_small_gene_sets),
            diagnostics=gmt_diagnostics,
            context={"program_method": "ptm_site", "contrast_method": "input_table", "link_method": "gene_assignment"},
        )
        if gmt_sets:
            gmt_path = resolve_gmt_out_path(out_dir, cfg.gmt_out)
            write_gmt(gmt_sets, gmt_path)
            output_files.append({"path": str(gmt_path), "role": "gmt"})
        _warn_gmt_diagnostics(gmt_diagnostics)

    skipped_gmt_outputs = _collect_skipped_gmt_outputs(gmt_diagnostics)
    run_summary_payload = {
        "converter": cfg.converter_name,
        "dataset_label": cfg.dataset_label,
        "signature_name": cfg.signature_name,
        "ptm_type": cfg.ptm_type,
        "score_mode": resolved_score_mode,
        "score_transform": cfg.score_transform,
        "protein_adjustment": cfg.protein_adjustment,
        "protein_adjustment_beta": float(beta),
        "confidence_weight_mode": cfg.confidence_weight_mode,
        "site_dup_policy": cfg.site_dup_policy,
        "gene_aggregation": cfg.gene_aggregation,
        "n_raw_rows": parse_summary["n_raw_rows"],
        "n_canonical_sites": len({rec.canonical_site_key for rec in collapsed_sites}),
        "n_rows_alias_matched": parse_summary["n_rows_alias_matched"],
        "n_rows_ambiguous_gene": parse_summary["n_rows_ambiguous_gene"],
        "n_rows_dropped_ambiguous_gene": parse_summary["n_rows_dropped_ambiguous_gene"],
        "n_rows_dropped_missing_score": parse_summary["n_rows_dropped_missing_score"],
        "n_rows_unresolved_site": parse_summary["n_rows_unresolved_site"],
        "n_sites_with_protein_adjustment": parse_summary["n_sites_with_protein_adjustment"],
        "n_sites_matched_to_ubiquity_prior": parse_summary["n_sites_matched_to_ubiquity_prior"],
        "fraction_sites_matched_to_ubiquity_prior": (
            float(parse_summary["n_sites_matched_to_ubiquity_prior"]) / float(len(base_rows)) if base_rows else 0.0
        ),
        "n_genes_pre_selection": len(gene_scores),
        "n_genes_selected": len(selected_rows),
        "resources": resources_info,
        "gmt_diagnostics": gmt_diagnostics,
        "skipped_programs": skipped_gmt_outputs,
    }
    run_json_path, run_txt_path = write_run_summary_files(out_dir, run_summary_payload)
    output_files.append({"path": str(run_json_path), "role": "run_summary_json"})
    output_files.append({"path": str(run_txt_path), "role": "run_summary_text"})

    meta_parameters = {
        "dataset_label": cfg.dataset_label,
        "signature_name": cfg.signature_name,
        "ptm_type": cfg.ptm_type,
        "score_mode": resolved_score_mode,
        "score_transform": cfg.score_transform,
        "protein_adjustment": cfg.protein_adjustment,
        "protein_adjustment_lambda": cfg.protein_adjustment_lambda,
        "confidence_weight_mode": cfg.confidence_weight_mode,
        "min_localization_prob": cfg.min_localization_prob,
        "site_dup_policy": cfg.site_dup_policy,
        "gene_aggregation": cfg.gene_aggregation,
        "gene_topk_sites": cfg.gene_topk_sites,
        "ambiguous_gene_policy": cfg.ambiguous_gene_policy,
        "select": cfg.select,
        "top_k": cfg.top_k,
        "quantile": cfg.quantile,
        "min_score": cfg.min_score,
        "normalize": cfg.normalize,
        "emit_full": cfg.emit_full,
        "emit_gmt": cfg.emit_gmt,
        "gmt_split_signed": cfg.gmt_split_signed,
        "columns": columns,
        "resources": resources_info,
    }
    if resources_info is None:
        meta_parameters["resources"] = None

    gmt_path = resolve_gmt_out_path(out_dir, cfg.gmt_out) if cfg.emit_gmt else None

    meta = make_metadata(
        converter_name=cfg.converter_name,
        parameters=meta_parameters,
        data_type="proteomics_ptm",
        assay="site_level",
        organism=cfg.organism,
        genome_build=cfg.genome_build,
        files=input_files,
        gene_annotation={
            "mode": "none",
            "source": "ptm_table_and_optional_alias",
            "gene_id_field": columns.get("gene_id") or columns.get("gene_symbol") or "gene_symbol",
        },
        weights={
            "weight_type": "signed" if cfg.score_transform == "signed" else "nonnegative",
            "normalization": {
                "method": cfg.normalize,
                "target_sum": 1.0 if cfg.normalize in {"within_set_l1", "l1"} else None,
            },
            "aggregation": cfg.gene_aggregation,
        },
        summary={
            "n_input_features": parse_summary["n_raw_rows"],
            "n_raw_rows": parse_summary["n_raw_rows"],
            "n_canonical_sites": len({rec.canonical_site_key for rec in collapsed_sites}),
            "n_rows_ambiguous_gene": parse_summary["n_rows_ambiguous_gene"],
            "n_rows_dropped_missing_score": parse_summary["n_rows_dropped_missing_score"],
            "n_sites_with_protein_adjustment": parse_summary["n_sites_with_protein_adjustment"],
            "n_sites_matched_to_ubiquity_prior": parse_summary["n_sites_matched_to_ubiquity_prior"],
            "fraction_sites_matched_to_ubiquity_prior": (
                float(parse_summary["n_sites_matched_to_ubiquity_prior"]) / float(len(base_rows)) if base_rows else 0.0
            ),
            "n_genes_pre_selection": len(gene_scores),
            "n_genes": len(selected_rows),
            "n_features_assigned": len(collapsed_sites),
            "fraction_features_assigned": (float(len(collapsed_sites)) / float(parse_summary["n_raw_rows"]) if parse_summary["n_raw_rows"] else 0.0),
        },
        program_extraction={
            "selection_method": cfg.select,
            "selection_params": {
                "k": cfg.top_k,
                "quantile": cfg.quantile,
                "min_score": cfg.min_score,
            },
            "normalize": cfg.normalize,
            "n_selected_genes": len(selected_rows),
            "score_definition": "site-level PTM evidence aggregated to genes after canonical-site collapse",
        },
        output_files=output_files,
        gmt={
            "written": bool(cfg.emit_gmt),
            "path": (
                str(gmt_path.relative_to(out_dir)) if gmt_path is not None and gmt_path.is_relative_to(out_dir) else str(gmt_path)
            )
            if gmt_path is not None
            else None,
            "prefer_symbol": bool(cfg.gmt_prefer_symbol),
            "require_symbol": bool(cfg.gmt_require_symbol),
            "biotype_allowlist": [x.strip() for x in str(cfg.gmt_biotype_allowlist).split(",") if x.strip()],
            "min_genes": int(cfg.gmt_min_genes),
            "max_genes": int(cfg.gmt_max_genes),
            "emit_small_gene_sets": bool(cfg.emit_small_gene_sets),
            "requested_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
            "emitted_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
            "skipped_outputs": skipped_gmt_outputs,
            "diagnostics": gmt_diagnostics,
            "plans": gmt_plans,
        },
    )
    write_metadata(out_dir / "geneset.meta.json", meta)

    return {
        "n_input_features": parse_summary["n_raw_rows"],
        "n_genes": len(selected_rows),
        "n_canonical_sites": len({rec.canonical_site_key for rec in collapsed_sites}),
        "out_dir": str(out_dir),
        "resolved_score_mode": resolved_score_mode,
    }
