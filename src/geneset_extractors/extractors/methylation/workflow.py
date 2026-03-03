from __future__ import annotations

import csv
from dataclasses import dataclass
import gzip
import math
from pathlib import Path
import re
import sys

from geneset_extractors.core.gmt import (
    build_gmt_sets_from_rows,
    parse_int_list_csv,
    parse_mass_list_csv,
    parse_str_list_csv,
    resolve_gmt_out_path,
    write_gmt,
)
from geneset_extractors.core.metadata import make_metadata, write_metadata
from geneset_extractors.core.methylation_programs import (
    PROGRAM_DISTAL_ACTIVITY,
    PROGRAM_LINKED_ACTIVITY,
    PROGRAM_PRESET_CONNECTABLE,
    PROGRAM_PROMOTER_ACTIVITY,
    link_methods_for_program,
    normalize_program_preset,
    resolve_link_methods,
    resolve_program_methods,
)
from geneset_extractors.core.peak_to_gene import (
    link_distance_decay,
    link_external_regions,
    link_nearest_tss,
    link_promoter_overlap,
)
from geneset_extractors.core.qc import collect_emitted_method_combinations, write_run_summary_files
from geneset_extractors.core.scoring import transform_peak_weight
from geneset_extractors.core.selection import (
    global_l1_weights,
    ranked_gene_ids,
    select_quantile,
    select_threshold,
    select_top_k,
    within_set_l1_weights,
)
from geneset_extractors.io.bed import read_bed
from geneset_extractors.io.gtf import read_genes_from_gtf
from geneset_extractors.io.region_gene_links import read_region_gene_links_tsv


DEFAULT_DELTA_COLUMNS = ("delta_beta", "delta_methylation", "delta_meth", "delta")
DEFAULT_PADJ_COLUMNS = ("padj", "fdr", "adj.p.val", "qvalue")
DEFAULT_PVALUE_COLUMNS = ("pvalue", "p_value", "p.val", "p")
SEX_CHROMS = {"x", "y", "chrx", "chry"}
CHROM_MISMATCH_WARN_FRACTION = 0.2
ARTIFACT_FAMILY_PATTERNS = (
    r"^OR[0-9A-Z]+",
    r"^KRT[0-9A-Z]+",
)
ARTIFACT_FAMILY_WARN_THRESHOLD = 0.2


@dataclass(frozen=True)
class MethylationFeature:
    feature_id: str
    chrom: str
    start: int
    end: int
    delta: float
    pvalue: float | None
    score_raw: float
    probe_id: str | None = None


@dataclass
class MethylationWorkflowConfig:
    converter_name: str
    out_dir: Path
    organism: str
    genome_build: str
    dataset_label: str
    gtf: str
    gtf_source: str | None
    program_preset: str
    program_methods: str | None
    link_method: str
    region_gene_links_tsv: str | None
    promoter_upstream_bp: int
    promoter_downstream_bp: int
    max_distance_bp: int | None
    decay_length_bp: int
    max_genes_per_peak: int
    score_transform: str
    normalize: str
    aggregation: str
    select: str
    top_k: int
    quantile: float
    min_score: float
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
    resource_policy: str
    score_mode: str
    delta_orientation: str
    distal_mode: str
    enhancer_bed: str | None
    exclude_gene_symbol_regex: list[str] | None
    exclude_gene_symbols_tsv: str | None


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_tsv_rows(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    with _open_text(Path(path)) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Input table has no header: {path}")
        fieldnames = [str(x) for x in reader.fieldnames]
        for row in reader:
            rows.append({k: str(v) for k, v in row.items() if k is not None})
    return fieldnames, rows


def read_probe_blacklist(path: str | Path | None) -> set[str]:
    if not path:
        return set()
    out: set[str] = set()
    with _open_text(Path(path)) as fh:
        for line in fh:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            token = text.split("\t", 1)[0].strip()
            if token:
                out.add(token)
    return out


def read_exclude_gene_symbols(path: str | Path | None) -> set[str]:
    if not path:
        return set()
    out: set[str] = set()
    with _open_text(Path(path)) as fh:
        for line in fh:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            token = text.split("\t", 1)[0].strip()
            if token:
                out.add(token.upper())
    return out


def compile_symbol_regexes(raw_patterns: list[str] | None) -> list[re.Pattern[str]]:
    out: list[re.Pattern[str]] = []
    for raw in raw_patterns or []:
        text = str(raw).strip()
        if not text:
            continue
        for piece in text.split(","):
            pattern = piece.strip()
            if not pattern:
                continue
            out.append(re.compile(pattern, flags=re.IGNORECASE))
    return out


def _orientation_factor(delta_orientation: str) -> float:
    if delta_orientation == "activity_oriented":
        return -1.0
    if delta_orientation == "raw":
        return 1.0
    raise ValueError(
        "Unsupported delta_orientation: "
        f"{delta_orientation}. Expected one of activity_oriented, raw."
    )


def _resolve_column(
    fieldnames: list[str],
    explicit: str | None,
    defaults: tuple[str, ...],
    required: bool,
    label: str,
) -> str | None:
    if explicit is not None and str(explicit).strip():
        name = str(explicit).strip()
        if name not in fieldnames:
            if not required:
                return None
            raise ValueError(
                f"Column '{name}' for {label} not found. Available columns: {', '.join(fieldnames)}"
            )
        return name
    lower_map = {x.lower(): x for x in fieldnames}
    for candidate in defaults:
        if candidate in fieldnames:
            return candidate
        mapped = lower_map.get(candidate.lower())
        if mapped is not None:
            return mapped
    if required:
        raise ValueError(
            f"Could not resolve required column for {label}. Tried defaults: {', '.join(defaults)}. "
            f"Available columns: {', '.join(fieldnames)}"
        )
    return None


def _parse_float(raw: object) -> float | None:
    if raw is None:
        return None
    text = str(raw).strip()
    if not text:
        return None
    low = text.lower()
    if low in {"na", "nan", "none", "null", "inf", "-inf"}:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def _normalize_chrom_label(chrom: str) -> str:
    text = str(chrom).strip()
    if not text:
        return ""
    return text


def _normalize_chrom_to_gene_space(chrom: str, gene_chroms: set[str]) -> tuple[str, bool]:
    c = _normalize_chrom_label(chrom)
    if not c:
        return "", False
    if c in gene_chroms:
        return c, False
    if c.lower().startswith("chr"):
        alt = c[3:]
        if alt in gene_chroms:
            return alt, True
    else:
        alt = "chr" + c
        if alt in gene_chroms:
            return alt, True
    return c, False


def _pos_to_interval(pos: int, pos_is_0based: bool) -> tuple[int, int]:
    if pos_is_0based:
        return int(pos), int(pos) + 1
    return max(0, int(pos) - 1), int(pos)


def read_probe_manifest_tsv(
    path: str | Path,
    probe_id_column: str = "probe_id",
    chrom_column: str = "chrom",
    pos_column: str = "pos",
    start_column: str | None = None,
    end_column: str | None = None,
    pos_is_0based: bool = False,
) -> dict[str, tuple[str, int, int]]:
    fieldnames, rows = read_tsv_rows(path)
    probe_col = _resolve_column(fieldnames, probe_id_column, ("probe_id", "id", "cg_id"), True, "probe_id")
    chrom_col = _resolve_column(fieldnames, chrom_column, ("chrom", "chr"), True, "chrom")
    start_col = _resolve_column(fieldnames, start_column, ("start",), False, "start")
    end_col = _resolve_column(fieldnames, end_column, ("end",), False, "end")
    pos_col = _resolve_column(fieldnames, pos_column, ("pos", "position"), False, "pos")
    if not (start_col and end_col) and not pos_col:
        raise ValueError("Probe manifest must include either (start,end) or pos columns")

    out: dict[str, tuple[str, int, int]] = {}
    for row in rows:
        probe_id = str(row.get(probe_col or "", "")).strip()
        chrom = str(row.get(chrom_col or "", "")).strip()
        if not probe_id or not chrom:
            continue
        start: int | None = None
        end: int | None = None
        if start_col and end_col:
            sv = _parse_float(row.get(start_col))
            ev = _parse_float(row.get(end_col))
            if sv is not None and ev is not None:
                start = int(sv)
                end = int(ev)
        if (start is None or end is None) and pos_col:
            pv = _parse_float(row.get(pos_col))
            if pv is not None:
                start, end = _pos_to_interval(int(pv), pos_is_0based)
        if start is None or end is None or end <= start:
            continue
        out[probe_id] = (chrom, int(start), int(end))
    if not out:
        raise ValueError(f"Probe manifest has no usable rows: {path}")
    return out


def _resolve_pvalue_column(
    fieldnames: list[str],
    padj_column: str | None,
    pvalue_column: str | None,
) -> str | None:
    p_col = _resolve_column(fieldnames, padj_column, DEFAULT_PADJ_COLUMNS, False, "padj")
    if p_col is not None:
        return p_col
    return _resolve_column(fieldnames, pvalue_column, DEFAULT_PVALUE_COLUMNS, False, "pvalue")


def parse_cpg_diff_features(
    *,
    fieldnames: list[str],
    rows: list[dict[str, str]],
    probe_id_column: str,
    chrom_column: str,
    pos_column: str,
    start_column: str | None,
    end_column: str | None,
    delta_column: str,
    padj_column: str | None,
    pvalue_column: str | None,
    score_mode: str,
    delta_orientation: str,
    neglog10p_eps: float,
    neglog10p_cap: float,
    input_pos_is_0based: bool,
    probe_manifest: dict[str, tuple[str, int, int]] | None,
    probe_blacklist: set[str],
    drop_sex_chrom: bool,
    gene_chroms: set[str],
) -> tuple[list[MethylationFeature], dict[str, object]]:
    probe_col = _resolve_column(fieldnames, probe_id_column, ("probe_id", "id", "cg_id"), False, "probe_id")
    chrom_col = _resolve_column(fieldnames, chrom_column, ("chrom", "chr"), False, "chrom")
    pos_col = _resolve_column(fieldnames, pos_column, ("pos", "position"), False, "pos")
    start_col = _resolve_column(fieldnames, start_column, ("start",), False, "start")
    end_col = _resolve_column(fieldnames, end_column, ("end",), False, "end")
    delta_col = _resolve_column(fieldnames, delta_column, DEFAULT_DELTA_COLUMNS, True, "delta")
    p_col = _resolve_pvalue_column(fieldnames, padj_column, pvalue_column)

    features: list[MethylationFeature] = []
    unresolved_probe_ids = 0
    dropped_blacklist = 0
    dropped_sex = 0
    chrom_converted = 0
    chrom_mismatch_skipped = 0
    rows_with_probe_id = 0
    skipped_rows = 0
    missing_p_for_mode = 0
    orient = _orientation_factor(delta_orientation)

    for idx, row in enumerate(rows, start=1):
        probe_id = str(row.get(probe_col or "", "")).strip() if probe_col else ""
        if probe_id:
            rows_with_probe_id += 1
        if probe_id and probe_id in probe_blacklist:
            dropped_blacklist += 1
            continue

        delta = _parse_float(row.get(delta_col or ""))
        if delta is None:
            skipped_rows += 1
            continue

        chrom = ""
        start: int | None = None
        end: int | None = None
        if chrom_col and start_col and end_col:
            chrom = str(row.get(chrom_col, "")).strip()
            sv = _parse_float(row.get(start_col))
            ev = _parse_float(row.get(end_col))
            if chrom and sv is not None and ev is not None:
                start = int(sv)
                end = int(ev)
        if (start is None or end is None) and chrom_col and pos_col:
            chrom = str(row.get(chrom_col, "")).strip()
            pv = _parse_float(row.get(pos_col))
            if chrom and pv is not None:
                start, end = _pos_to_interval(int(pv), input_pos_is_0based)
        if (start is None or end is None or not chrom) and probe_id and probe_manifest is not None:
            mapped = probe_manifest.get(probe_id)
            if mapped is not None:
                chrom, start, end = mapped
        if start is None or end is None or not chrom or end <= start:
            if probe_id:
                unresolved_probe_ids += 1
            else:
                skipped_rows += 1
            continue

        chrom_norm, converted = _normalize_chrom_to_gene_space(chrom, gene_chroms)
        if converted:
            chrom_converted += 1
        if chrom_norm not in gene_chroms:
            chrom_mismatch_skipped += 1
            skipped_rows += 1
            continue

        if drop_sex_chrom and chrom_norm.lower() in SEX_CHROMS:
            dropped_sex += 1
            continue

        pvalue = _parse_float(row.get(p_col or "")) if p_col else None
        if score_mode == "delta_times_neglog10p":
            if pvalue is None:
                missing_p_for_mode += 1
                continue
            p = min(1.0, max(float(neglog10p_eps), float(pvalue)))
            neglog10p = min(-math.log10(p + float(neglog10p_eps)), float(neglog10p_cap))
            score_raw = float(delta) * float(orient) * float(neglog10p)
        elif score_mode == "delta_only":
            score_raw = float(delta) * float(orient)
        else:
            raise ValueError(f"Unsupported methylation score_mode: {score_mode}")

        fid = probe_id if probe_id else f"row_{idx}"
        features.append(
            MethylationFeature(
                feature_id=fid,
                chrom=chrom_norm,
                start=int(start),
                end=int(end),
                delta=float(delta),
                pvalue=(float(pvalue) if pvalue is not None else None),
                score_raw=float(score_raw),
                probe_id=probe_id or None,
            )
        )

    summary = {
        "n_rows": len(rows),
        "n_features_parsed": len(features),
        "n_probe_unresolved": unresolved_probe_ids,
        "n_rows_with_probe_id": rows_with_probe_id,
        "n_rows_skipped": skipped_rows,
        "n_rows_skipped_chrom_mismatch": chrom_mismatch_skipped,
        "n_rows_missing_pvalue_for_score_mode": missing_p_for_mode,
        "n_dropped_probe_blacklist": dropped_blacklist,
        "n_dropped_sex_chrom": dropped_sex,
        "n_chrom_name_converted": chrom_converted,
        "score_mode": score_mode,
        "delta_orientation": delta_orientation,
    }
    return features, summary


def parse_dmr_region_features(
    *,
    fieldnames: list[str],
    rows: list[dict[str, str]],
    chrom_column: str,
    start_column: str,
    end_column: str,
    delta_column: str,
    padj_column: str | None,
    pvalue_column: str | None,
    score_mode: str,
    delta_orientation: str,
    neglog10p_eps: float,
    neglog10p_cap: float,
    drop_sex_chrom: bool,
    gene_chroms: set[str],
) -> tuple[list[MethylationFeature], dict[str, object]]:
    chrom_col = _resolve_column(fieldnames, chrom_column, ("chrom", "chr"), True, "chrom")
    start_col = _resolve_column(fieldnames, start_column, ("start",), True, "start")
    end_col = _resolve_column(fieldnames, end_column, ("end",), True, "end")
    delta_col = _resolve_column(fieldnames, delta_column, DEFAULT_DELTA_COLUMNS, True, "delta")
    p_col = _resolve_pvalue_column(fieldnames, padj_column, pvalue_column)

    features: list[MethylationFeature] = []
    dropped_sex = 0
    chrom_converted = 0
    chrom_mismatch_skipped = 0
    skipped_rows = 0
    missing_p_for_mode = 0
    orient = _orientation_factor(delta_orientation)

    for idx, row in enumerate(rows, start=1):
        chrom = str(row.get(chrom_col or "", "")).strip()
        sv = _parse_float(row.get(start_col))
        ev = _parse_float(row.get(end_col))
        delta = _parse_float(row.get(delta_col or ""))
        if not chrom or sv is None or ev is None or delta is None:
            skipped_rows += 1
            continue
        start = int(sv)
        end = int(ev)
        if end <= start:
            skipped_rows += 1
            continue

        chrom_norm, converted = _normalize_chrom_to_gene_space(chrom, gene_chroms)
        if converted:
            chrom_converted += 1
        if chrom_norm not in gene_chroms:
            chrom_mismatch_skipped += 1
            skipped_rows += 1
            continue

        if drop_sex_chrom and chrom_norm.lower() in SEX_CHROMS:
            dropped_sex += 1
            continue

        pvalue = _parse_float(row.get(p_col or "")) if p_col else None
        if score_mode == "delta_times_neglog10p":
            if pvalue is None:
                missing_p_for_mode += 1
                continue
            p = min(1.0, max(float(neglog10p_eps), float(pvalue)))
            neglog10p = min(-math.log10(p + float(neglog10p_eps)), float(neglog10p_cap))
            score_raw = float(delta) * float(orient) * float(neglog10p)
        elif score_mode == "delta_only":
            score_raw = float(delta) * float(orient)
        else:
            raise ValueError(f"Unsupported methylation score_mode: {score_mode}")

        features.append(
            MethylationFeature(
                feature_id=f"dmr_{idx}",
                chrom=chrom_norm,
                start=start,
                end=end,
                delta=float(delta),
                pvalue=(float(pvalue) if pvalue is not None else None),
                score_raw=float(score_raw),
                probe_id=None,
            )
        )

    summary = {
        "n_rows": len(rows),
        "n_features_parsed": len(features),
        "n_rows_skipped": skipped_rows,
        "n_rows_skipped_chrom_mismatch": chrom_mismatch_skipped,
        "n_rows_missing_pvalue_for_score_mode": missing_p_for_mode,
        "n_dropped_sex_chrom": dropped_sex,
        "n_chrom_name_converted": chrom_converted,
        "score_mode": score_mode,
        "delta_orientation": delta_orientation,
    }
    return features, summary


def _resolve_max_distance(cfg: MethylationWorkflowConfig, link_method: str) -> int:
    if link_method == "external":
        return 0
    if cfg.max_distance_bp is not None:
        return int(cfg.max_distance_bp)
    if link_method == "distance_decay":
        return 500000
    if link_method == "nearest_tss":
        return 100000
    return 100000


def _strip_version(gene_id: str) -> str:
    return str(gene_id).split(".", 1)[0]


def _resolve_external_link_gene_ids(region_gene_links, genes) -> tuple[list[dict[str, object]], int]:
    full_map = {str(g.gene_id): str(g.gene_id) for g in genes}
    nover_map: dict[str, str | None] = {}
    for g in genes:
        k = _strip_version(str(g.gene_id))
        if k in nover_map and nover_map[k] != str(g.gene_id):
            nover_map[k] = None
        else:
            nover_map[k] = str(g.gene_id)

    resolved: list[dict[str, object]] = []
    unresolved = 0
    for row in region_gene_links:
        gid = str(row["gene_id"])
        mapped = full_map.get(gid)
        if mapped is None:
            mapped = nover_map.get(_strip_version(gid))
        if not mapped:
            unresolved += 1
            continue
        resolved.append({**row, "gene_id": mapped})
    return resolved, unresolved


def _link(
    peaks: list[dict[str, object]],
    genes,
    cfg: MethylationWorkflowConfig,
    link_method: str,
    resolved_external_links: list[dict[str, object]] | None,
) -> list[dict[str, object]]:
    if link_method == "external":
        if resolved_external_links is None:
            raise ValueError("internal error: external links requested without resolved links")
        return link_external_regions(peaks, resolved_external_links)
    max_distance_bp = _resolve_max_distance(cfg, link_method)
    if link_method == "promoter_overlap":
        return link_promoter_overlap(
            peaks,
            genes,
            cfg.promoter_upstream_bp,
            cfg.promoter_downstream_bp,
        )
    if link_method == "nearest_tss":
        return link_nearest_tss(peaks, genes, max_distance_bp)
    if link_method == "distance_decay":
        return link_distance_decay(
            peaks,
            genes,
            max_distance_bp,
            cfg.decay_length_bp,
            cfg.max_genes_per_peak,
        )
    raise ValueError(f"Unsupported link_method: {link_method}")


def _select_gene_ids(magnitude_scores: dict[str, float], cfg: MethylationWorkflowConfig) -> list[str]:
    if cfg.select == "none":
        return ranked_gene_ids(magnitude_scores)
    if cfg.select == "top_k":
        return select_top_k(magnitude_scores, int(cfg.top_k))
    if cfg.select == "quantile":
        return select_quantile(magnitude_scores, float(cfg.quantile))
    if cfg.select == "threshold":
        return select_threshold(magnitude_scores, float(cfg.min_score))
    raise ValueError(f"Unsupported selection method: {cfg.select}")


def _selected_weights(
    magnitude_scores: dict[str, float],
    selected_gene_ids: list[str],
    normalize: str,
) -> dict[str, float]:
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
    fieldnames = ["gene_id", "score", "rank"]
    if any("weight" in row for row in rows):
        fieldnames.append("weight")
    if any(str(row.get("gene_symbol", "")).strip() for row in rows):
        fieldnames.append("gene_symbol")
    if any(str(row.get("gene_biotype", "")).strip() for row in rows):
        fieldnames.append("gene_biotype")
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _aggregate_gene_scores(
    feature_values: list[float],
    links: list[dict[str, object]],
    aggregation: str,
    eps: float = 1e-12,
) -> dict[str, float]:
    if aggregation not in {"weighted_mean", "sum", "mean", "max_abs"}:
        raise ValueError(f"Unsupported aggregation method: {aggregation}")

    if aggregation == "weighted_mean":
        num: dict[str, float] = {}
        den: dict[str, float] = {}
        for link in links:
            p = int(link["peak_index"])
            g = str(link["gene_id"])
            lw = float(link["link_weight"])
            val = float(feature_values[p]) * lw
            num[g] = num.get(g, 0.0) + val
            den[g] = den.get(g, 0.0) + abs(lw)
        out: dict[str, float] = {}
        for gene_id, numerator in num.items():
            denominator = float(den.get(gene_id, 0.0))
            if denominator <= float(eps):
                continue
            out[gene_id] = float(numerator) / denominator
        return out

    if aggregation == "sum":
        out: dict[str, float] = {}
        for link in links:
            p = int(link["peak_index"])
            g = str(link["gene_id"])
            lw = float(link["link_weight"])
            out[g] = out.get(g, 0.0) + float(feature_values[p]) * lw
        return out

    if aggregation == "mean":
        num: dict[str, float] = {}
        cnt: dict[str, int] = {}
        for link in links:
            p = int(link["peak_index"])
            g = str(link["gene_id"])
            lw = float(link["link_weight"])
            num[g] = num.get(g, 0.0) + float(feature_values[p]) * lw
            cnt[g] = cnt.get(g, 0) + 1
        return {g: float(num[g]) / float(cnt[g]) for g in num if cnt[g] > 0}

    out_max: dict[str, float] = {}
    out_abs: dict[str, float] = {}
    for link in links:
        p = int(link["peak_index"])
        g = str(link["gene_id"])
        lw = float(link["link_weight"])
        val = float(feature_values[p]) * lw
        mag = abs(val)
        if g not in out_abs or mag > out_abs[g]:
            out_abs[g] = mag
            out_max[g] = val
    return out_max


def _rows_from_scores(
    scores: dict[str, float],
    genes_by_id: dict[str, dict[str, str]],
    cfg: MethylationWorkflowConfig,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    magnitude = {gid: abs(float(v)) for gid, v in scores.items()}
    selected_gene_ids = _select_gene_ids(magnitude, cfg)
    selected_gene_ids = sorted(selected_gene_ids, key=lambda gid: (-float(magnitude.get(gid, 0.0)), str(gid)))
    weights = _selected_weights(magnitude, selected_gene_ids, cfg.normalize)

    selected_rows: list[dict[str, object]] = []
    for rank, gene_id in enumerate(selected_gene_ids, start=1):
        info = genes_by_id.get(gene_id, {})
        row: dict[str, object] = {
            "gene_id": gene_id,
            "score": float(scores.get(gene_id, 0.0)),
            "weight": float(weights.get(gene_id, 0.0)),
            "rank": rank,
        }
        symbol = str(info.get("gene_symbol", "")).strip()
        biotype = str(info.get("gene_biotype", "")).strip()
        if symbol:
            row["gene_symbol"] = symbol
        if biotype:
            row["gene_biotype"] = biotype
        selected_rows.append(row)

    full_gene_ids = sorted(
        [gid for gid, score in scores.items() if float(score) != 0.0],
        key=lambda gid: (-abs(float(scores[gid])), str(gid)),
    )
    full_rows: list[dict[str, object]] = []
    for rank, gene_id in enumerate(full_gene_ids, start=1):
        info = genes_by_id.get(gene_id, {})
        row = {
            "gene_id": gene_id,
            "score": float(scores[gene_id]),
            "rank": rank,
        }
        symbol = str(info.get("gene_symbol", "")).strip()
        biotype = str(info.get("gene_biotype", "")).strip()
        if symbol:
            row["gene_symbol"] = symbol
        if biotype:
            row["gene_biotype"] = biotype
        full_rows.append(row)
    return selected_rows, full_rows


def _warn_program_methods_skipped(program_methods_skipped: dict[str, str]) -> None:
    if not program_methods_skipped:
        return
    for method in sorted(program_methods_skipped):
        print(
            f"warning: program_method={method} skipped: {program_methods_skipped[method]}",
            file=sys.stderr,
        )


def _warn_gmt_output_diagnostics(gmt_diagnostics: list[dict[str, object]]) -> None:
    for diag in gmt_diagnostics:
        code = str(diag.get("code", ""))
        base_name = str(diag.get("base_name", "unknown"))
        n_genes = diag.get("n_genes", "na")
        reason = str(diag.get("reason", "")).strip()
        suggestion = str(diag.get("suggestion", "")).strip()
        if code in {"small_gene_set_skipped", "no_positive_genes"}:
            print(
                f"warning: skipped GMT output for {base_name}; n_genes={n_genes}. {reason}",
                file=sys.stderr,
            )
        elif code == "small_gene_set_emitted":
            print(
                f"warning: emitted small GMT output for {base_name}; n_genes={n_genes}. {reason}",
                file=sys.stderr,
            )
        elif code == "marginal_gene_count":
            print(
                f"warning: marginal GMT signal for {base_name}; n_genes={n_genes}. {reason}",
                file=sys.stderr,
            )
        elif code == "require_symbol_heavy_drop":
            print(
                f"warning: GMT symbol requirement dropped many rows for {base_name}. {reason}",
                file=sys.stderr,
            )
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
                "program_method": str(diag.get("program_method", "")),
                "link_method": str(diag.get("link_method", "")),
                "reason": str(diag.get("reason", "")),
                "n_genes": int(diag.get("n_genes", 0) or 0),
                "code": code,
            }
        )
    return skipped


def _apply_gene_symbol_exclusion(
    scores: dict[str, float],
    genes_by_id: dict[str, dict[str, str]],
    patterns: list[re.Pattern[str]],
    symbol_blocklist: set[str],
) -> tuple[dict[str, float], dict[str, int]]:
    if not patterns and not symbol_blocklist:
        return scores, {"n_excluded_regex": 0, "n_excluded_symbol_list": 0, "n_excluded_total": 0}
    out: dict[str, float] = {}
    excluded_regex = 0
    excluded_symbol_list = 0
    for gene_id, score in scores.items():
        info = genes_by_id.get(gene_id, {})
        symbol = str(info.get("gene_symbol", "")).strip()
        token = symbol if symbol else str(gene_id)
        upper = token.upper()
        drop_by_list = upper in symbol_blocklist if symbol_blocklist else False
        drop_by_regex = any(p.search(token) for p in patterns) if patterns else False
        if drop_by_list:
            excluded_symbol_list += 1
            continue
        if drop_by_regex:
            excluded_regex += 1
            continue
        out[gene_id] = float(score)
    return out, {
        "n_excluded_regex": excluded_regex,
        "n_excluded_symbol_list": excluded_symbol_list,
        "n_excluded_total": excluded_regex + excluded_symbol_list,
    }


def _artifact_family_diagnostics(gmt_sets: list[tuple[str, list[str]]]) -> dict[str, object]:
    patterns = [re.compile(p, flags=re.IGNORECASE) for p in ARTIFACT_FAMILY_PATTERNS]
    warnings: list[dict[str, object]] = []
    for set_name, genes in gmt_sets:
        n = len(genes)
        if n <= 0:
            continue
        for pattern_text, pattern in zip(ARTIFACT_FAMILY_PATTERNS, patterns):
            matched = sum(1 for gene in genes if pattern.search(str(gene)))
            frac = float(matched) / float(n)
            if frac > ARTIFACT_FAMILY_WARN_THRESHOLD:
                warnings.append(
                    {
                        "set_name": set_name,
                        "pattern": pattern_text,
                        "n_matched": matched,
                        "n_genes": n,
                        "fraction": frac,
                        "threshold": ARTIFACT_FAMILY_WARN_THRESHOLD,
                    }
                )
    return {
        "pattern_threshold": ARTIFACT_FAMILY_WARN_THRESHOLD,
        "patterns": list(ARTIFACT_FAMILY_PATTERNS),
        "warnings": warnings,
    }


def _feature_overlaps_intervals(feature: MethylationFeature, intervals: list[dict[str, object]]) -> bool:
    for iv in intervals:
        if str(iv.get("chrom", "")) != feature.chrom:
            continue
        start = int(iv.get("start", 0))
        end = int(iv.get("end", 0))
        if feature.start < end and start < feature.end:
            return True
    return False


def run_methylation_workflow(
    *,
    cfg: MethylationWorkflowConfig,
    features: list[MethylationFeature],
    parse_summary: dict[str, object],
    input_files: list[dict[str, str]],
    resources_info: dict[str, object] | None = None,
) -> dict[str, object]:
    if cfg.resource_policy not in {"skip", "fail"}:
        raise ValueError(f"Unsupported resource_policy: {cfg.resource_policy}")
    if not features:
        raise ValueError("No methylation features with usable coordinates/scores were parsed")

    n_rows = int(parse_summary.get("n_rows", 0) or 0)
    n_dropped_sex = int(parse_summary.get("n_dropped_sex_chrom", 0) or 0)
    if n_dropped_sex > 0:
        print(
            "warning: drop_sex_chrom=true removed "
            f"{n_dropped_sex} features. Set --drop_sex_chrom false to keep sex chromosomes.",
            file=sys.stderr,
        )
    n_chrom_mismatch = int(parse_summary.get("n_rows_skipped_chrom_mismatch", 0) or 0)
    if n_rows > 0 and (float(n_chrom_mismatch) / float(n_rows)) > CHROM_MISMATCH_WARN_FRACTION:
        print(
            "warning: many rows could not be mapped to GTF chromosome names/genome build "
            f"({n_chrom_mismatch}/{n_rows}). Check --genome_build and chr-prefix conventions.",
            file=sys.stderr,
        )

    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    program_preset = normalize_program_preset(cfg.program_preset)
    program_methods = resolve_program_methods(program_preset, cfg.program_methods)
    link_methods = resolve_link_methods(cfg.link_method, program_preset)

    genes = read_genes_from_gtf(cfg.gtf)
    if not genes:
        raise ValueError("GTF parsing produced no genes")

    genes_by_id = {
        str(g.gene_id): {
            "gene_symbol": str(g.gene_symbol or ""),
            "gene_biotype": str(g.gene_biotype or ""),
        }
        for g in genes
    }
    exclude_patterns = compile_symbol_regexes(cfg.exclude_gene_symbol_regex)
    exclude_symbol_list = read_exclude_gene_symbols(cfg.exclude_gene_symbols_tsv)

    peaks = [
        {
            "chrom": f.chrom,
            "start": int(f.start),
            "end": int(f.end),
            "feature_id": f.feature_id,
            "probe_id": f.probe_id,
        }
        for f in features
    ]
    base_feature_scores = [float(f.score_raw) for f in features]

    resolved_external_links: list[dict[str, object]] | None = None
    n_external_links_unresolved = 0
    link_methods_effective = list(link_methods)
    if "external" in link_methods_effective:
        if not cfg.region_gene_links_tsv:
            message = "link_method includes external but --region_gene_links_tsv was not provided"
            if cfg.resource_policy == "fail":
                raise ValueError(message)
            print(f"warning: {message}; skipping external link method", file=sys.stderr)
            link_methods_effective = [m for m in link_methods_effective if m != "external"]
        else:
            raw_links = read_region_gene_links_tsv(cfg.region_gene_links_tsv)
            resolved_external_links, n_external_links_unresolved = _resolve_external_link_gene_ids(raw_links, genes)
            if not resolved_external_links:
                message = (
                    "No external region-gene links mapped to GTF gene_ids; "
                    "check gene_id namespace/version"
                )
                if cfg.resource_policy == "fail":
                    raise ValueError(message)
                print(f"warning: {message}; skipping external link method", file=sys.stderr)
                link_methods_effective = [m for m in link_methods_effective if m != "external"]

    if not link_methods_effective:
        raise ValueError("No link methods available after resolving requested configuration")

    links_by_method = {
        method: _link(peaks, genes, cfg, method, resolved_external_links)
        for method in link_methods_effective
    }

    promoter_links_mask = link_promoter_overlap(
        peaks,
        genes,
        cfg.promoter_upstream_bp,
        cfg.promoter_downstream_bp,
    )
    promoter_peak_indices = {int(link["peak_index"]) for link in promoter_links_mask}

    enhancer_intervals: list[dict[str, object]] = []
    enhancer_peak_indices: set[int] = set()
    if cfg.enhancer_bed:
        enhancer_intervals = read_bed(cfg.enhancer_bed)
        enhancer_peak_indices = {
            idx
            for idx, feature in enumerate(features)
            if _feature_overlaps_intervals(feature, enhancer_intervals)
        }

    program_methods_skipped: dict[str, str] = {}
    scored_rows_by_program_by_link: dict[str, dict[str, list[dict[str, object]]]] = {}
    selected_rows_by_program_by_link: dict[str, dict[str, list[dict[str, object]]]] = {}
    filter_counts_by_output: dict[str, dict[str, int]] = {}

    primary_program_method: str | None = None
    primary_link_method: str | None = None
    primary_selected_rows: list[dict[str, object]] | None = None
    primary_full_rows: list[dict[str, object]] | None = None

    for program_method in program_methods:
        allowed_methods = link_methods_for_program(program_preset, program_method, link_methods_effective)
        if not allowed_methods:
            program_methods_skipped[program_method] = "no compatible link_method selected for active preset"
            continue

        scored_rows_by_program_by_link.setdefault(program_method, {})
        selected_rows_by_program_by_link.setdefault(program_method, {})

        for link_method in allowed_methods:
            masked_scores = list(base_feature_scores)
            if program_method == PROGRAM_PROMOTER_ACTIVITY:
                masked_scores = [
                    (float(v) if idx in promoter_peak_indices else 0.0)
                    for idx, v in enumerate(base_feature_scores)
                ]
            elif program_method == PROGRAM_DISTAL_ACTIVITY:
                masked_scores = [
                    (float(v) if idx not in promoter_peak_indices else 0.0)
                    for idx, v in enumerate(base_feature_scores)
                ]
                if cfg.distal_mode == "enhancer_only":
                    if not cfg.enhancer_bed:
                        message = (
                            "distal_mode=enhancer_only requires --enhancer_bed or --enhancer_resource_id"
                        )
                        if cfg.resource_policy == "fail":
                            raise ValueError(message)
                        program_methods_skipped[PROGRAM_DISTAL_ACTIVITY] = message
                        continue
                    masked_scores = [
                        (float(v) if idx in enhancer_peak_indices else 0.0)
                        for idx, v in enumerate(masked_scores)
                    ]

            transformed = [transform_peak_weight(float(v), cfg.score_transform) for v in masked_scores]
            links = links_by_method[link_method]
            scores = _aggregate_gene_scores(transformed, links, cfg.aggregation)
            scores, filter_counts = _apply_gene_symbol_exclusion(
                scores,
                genes_by_id,
                exclude_patterns,
                exclude_symbol_list,
            )
            filter_counts_by_output[f"{program_method}__{link_method}"] = filter_counts
            selected_rows, full_rows = _rows_from_scores(scores, genes_by_id, cfg)
            scored_rows_by_program_by_link[program_method][link_method] = full_rows
            selected_rows_by_program_by_link[program_method][link_method] = selected_rows

            if primary_selected_rows is None:
                primary_program_method = program_method
                primary_link_method = link_method
                primary_selected_rows = selected_rows
                primary_full_rows = full_rows

    _warn_program_methods_skipped(program_methods_skipped)

    if primary_selected_rows is None or primary_full_rows is None or primary_program_method is None or primary_link_method is None:
        raise ValueError("No methylation program outputs were produced for the requested configuration")

    _write_rows(out_dir / "geneset.tsv", primary_selected_rows)
    if cfg.emit_full:
        _write_rows(out_dir / "geneset.full.tsv", primary_full_rows)

    gmt_sets: list[tuple[str, list[str]]] = []
    gmt_plans: list[dict[str, object]] = []
    gmt_diagnostics: list[dict[str, object]] = []
    gmt_path = resolve_gmt_out_path(out_dir, cfg.gmt_out)
    gmt_topk_list = parse_int_list_csv(str(cfg.gmt_topk_list))
    gmt_mass_list = parse_mass_list_csv(str(cfg.gmt_mass_list))
    gmt_biotype_allowlist = parse_str_list_csv(str(cfg.gmt_biotype_allowlist))

    if cfg.emit_gmt:
        for program_method, by_link in scored_rows_by_program_by_link.items():
            for link_method, full_rows in by_link.items():
                base_name = (
                    f"{cfg.converter_name}"
                    f"__program={program_method}"
                    f"__link_method={link_method}"
                )
                method_sets, method_plans = build_gmt_sets_from_rows(
                    rows=full_rows,
                    base_name=base_name,
                    prefer_symbol=bool(cfg.gmt_prefer_symbol),
                    min_genes=int(cfg.gmt_min_genes),
                    max_genes=int(cfg.gmt_max_genes),
                    topk_list=gmt_topk_list,
                    mass_list=gmt_mass_list,
                    split_signed=bool(cfg.gmt_split_signed),
                    require_symbol=bool(cfg.gmt_require_symbol),
                    allowed_biotypes={x.lower() for x in gmt_biotype_allowlist} if gmt_biotype_allowlist else None,
                    emit_small_gene_sets=bool(cfg.emit_small_gene_sets),
                    diagnostics=gmt_diagnostics,
                    context={
                        "program_method": program_method,
                        "link_method": link_method,
                    },
                )
                gmt_sets.extend(method_sets)
                for plan in method_plans:
                    params = dict(plan.get("parameters", {}))
                    params["program_method"] = program_method
                    params["calibration_method"] = "none"
                    params["link_method"] = link_method
                    plan["parameters"] = params
                gmt_plans.extend(method_plans)

        write_gmt(gmt_sets, gmt_path)
    artifact_diagnostics = _artifact_family_diagnostics(gmt_sets) if cfg.emit_gmt else {"warnings": []}
    for warning in artifact_diagnostics.get("warnings", []):
        if not isinstance(warning, dict):
            continue
        frac = float(warning.get("fraction", 0.0) or 0.0)
        print(
            "warning: potential gene-family dominance in "
            f"{warning.get('set_name','unknown')}: pattern={warning.get('pattern')} "
            f"fraction={frac:.3f}. "
            "Consider --exclude_gene_symbol_regex or interpret cautiously.",
            file=sys.stderr,
        )

    _warn_gmt_output_diagnostics(gmt_diagnostics)

    files = list(input_files)
    output_files = [{"path": str(out_dir / "geneset.tsv"), "role": "selected_program"}]
    if cfg.emit_full:
        output_files.append({"path": str(out_dir / "geneset.full.tsv"), "role": "full_scores"})
    if cfg.emit_gmt:
        output_files.append({"path": str(gmt_path), "role": "gmt"})

    emitted_combinations = collect_emitted_method_combinations(gmt_plans)
    requested_gmt_outputs = [{"name": p.get("name", "")} for p in gmt_plans]
    skipped_gmt_outputs = _collect_skipped_gmt_outputs(gmt_diagnostics)
    sign_semantics = {
        "delta_orientation": cfg.delta_orientation,
        "pos_set_meaning": (
            "hypomethylation-associated activity increase"
            if cfg.delta_orientation == "activity_oriented"
            else "positive raw delta signal"
        ),
        "neg_set_meaning": (
            "hypermethylation-associated activity decrease"
            if cfg.delta_orientation == "activity_oriented"
            else "negative raw delta signal"
        ),
    }

    assigned_primary = len({int(link["peak_index"]) for link in links_by_method.get(primary_link_method, [])})
    link_assignment: dict[str, dict[str, object]] = {}
    for method, links in links_by_method.items():
        assigned = len({int(link["peak_index"]) for link in links})
        link_assignment[method] = {
            "n_features_assigned": assigned,
            "fraction_features_assigned": (assigned / len(peaks) if peaks else 0.0),
        }

    run_summary_payload: dict[str, object] = {
        "converter": cfg.converter_name,
        "dataset_label": cfg.dataset_label,
        "program_preset": program_preset,
        "primary_program_method": primary_program_method,
        "primary_link_method": primary_link_method,
        "primary_calibration_method": "none",
        "selected_direction": "PRIMARY",
        "n_input_peaks": len(peaks),
        "link_assignment": link_assignment,
        "promoter_peak_count": len(promoter_peak_indices),
        "distal_peak_count": len(peaks) - len(promoter_peak_indices),
        "emitted_method_combinations": emitted_combinations,
        "requested_gmt_outputs": requested_gmt_outputs,
        "skipped_gmt_outputs": skipped_gmt_outputs,
        "program_methods_skipped": program_methods_skipped,
        "calibration_methods_skipped": {},
        "parse_summary": parse_summary,
        "sign_semantics": sign_semantics,
        "artifact_diagnostics": artifact_diagnostics,
        "symbol_filter": {
            "exclude_gene_symbol_regex": cfg.exclude_gene_symbol_regex or [],
            "exclude_gene_symbols_tsv": cfg.exclude_gene_symbols_tsv,
            "counts_by_output": filter_counts_by_output,
        },
    }
    if resources_info is not None:
        run_summary_payload["resources"] = resources_info
    run_summary_json_path, run_summary_txt_path = write_run_summary_files(out_dir, run_summary_payload)
    output_files.append({"path": str(run_summary_json_path), "role": "run_summary_json"})
    output_files.append({"path": str(run_summary_txt_path), "role": "run_summary_text"})

    params: dict[str, object] = {
        "dataset_label": cfg.dataset_label,
        "program_preset": program_preset,
        "program_methods": cfg.program_methods,
        "program_methods_active": sorted(scored_rows_by_program_by_link.keys()),
        "program_methods_skipped": program_methods_skipped,
        "link_method": cfg.link_method,
        "link_methods_evaluated": link_methods_effective,
        "primary_program_method": primary_program_method,
        "primary_link_method": primary_link_method,
        "score_mode": cfg.score_mode,
        "delta_orientation": cfg.delta_orientation,
        "score_transform": cfg.score_transform,
        "aggregation": cfg.aggregation,
        "select": cfg.select,
        "top_k": cfg.top_k,
        "quantile": cfg.quantile,
        "min_score": cfg.min_score,
        "normalize": cfg.normalize,
        "promoter_upstream_bp": cfg.promoter_upstream_bp,
        "promoter_downstream_bp": cfg.promoter_downstream_bp,
        "max_distance_bp": cfg.max_distance_bp,
        "decay_length_bp": cfg.decay_length_bp,
        "max_genes_per_peak": cfg.max_genes_per_peak,
        "emit_full": cfg.emit_full,
        "emit_gmt": cfg.emit_gmt,
        "gmt_min_genes": cfg.gmt_min_genes,
        "gmt_max_genes": cfg.gmt_max_genes,
        "gmt_topk_list": gmt_topk_list,
        "gmt_mass_list": gmt_mass_list,
        "gmt_split_signed": cfg.gmt_split_signed,
        "emit_small_gene_sets": cfg.emit_small_gene_sets,
        "distal_mode": cfg.distal_mode,
        "enhancer_bed": cfg.enhancer_bed,
        "exclude_gene_symbol_regex": cfg.exclude_gene_symbol_regex or [],
        "exclude_gene_symbols_tsv": cfg.exclude_gene_symbols_tsv,
        "resource_policy": cfg.resource_policy,
        "n_external_links_unresolved_gene_id": n_external_links_unresolved,
    }

    oriented_delta_expr = "-delta" if cfg.delta_orientation == "activity_oriented" else "delta"
    score_definition = (
        f"gene_score = aggregate(link_weight * transform(({oriented_delta_expr}) * neglog10p))"
        if cfg.score_mode == "delta_times_neglog10p"
        else f"gene_score = aggregate(link_weight * transform({oriented_delta_expr}))"
    )

    meta = make_metadata(
        converter_name=cfg.converter_name,
        parameters=params,
        data_type="dna_methylation",
        assay="bulk",
        organism=cfg.organism,
        genome_build=cfg.genome_build,
        files=files,
        gene_annotation={
            "mode": "gtf",
            "gtf_path": str(cfg.gtf),
            "source": cfg.gtf_source or "user",
            "gene_id_field": "gene_id",
        },
        weights={
            "weight_type": "signed" if cfg.score_transform == "signed" else "nonnegative",
            "normalization": {
                "method": cfg.normalize,
                "target_sum": 1.0 if cfg.normalize in {"within_set_l1", "l1"} else None,
            },
            "aggregation": cfg.aggregation,
        },
        summary={
            "n_input_features": int(len(peaks)),
            "n_genes": int(len(primary_selected_rows)),
            "n_features_assigned": int(assigned_primary),
            "fraction_features_assigned": (assigned_primary / len(peaks) if peaks else 0.0),
            "parse_summary": parse_summary,
            "n_program_methods": len(scored_rows_by_program_by_link),
            "n_program_methods_skipped": len(program_methods_skipped),
            "n_external_links_unresolved_gene_id": n_external_links_unresolved,
            "symbol_filter_counts_by_output": filter_counts_by_output,
        },
        program_extraction={
            "selection_method": cfg.select,
            "selection_params": {
                "k": cfg.top_k,
                "quantile": cfg.quantile,
                "min_score": cfg.min_score,
            },
            "normalize": cfg.normalize,
            "n_selected_genes": len(primary_selected_rows),
            "score_definition": score_definition,
            "program_methods": sorted(scored_rows_by_program_by_link.keys()),
            "program_methods_skipped": program_methods_skipped,
            "primary_program_method": primary_program_method,
            "primary_link_method": primary_link_method,
            "delta_orientation": cfg.delta_orientation,
        },
        output_files=output_files,
        gmt={
            "written": cfg.emit_gmt,
            "path": (
                str(gmt_path.relative_to(out_dir)) if gmt_path.is_relative_to(out_dir) else str(gmt_path)
            )
            if cfg.emit_gmt
            else None,
            "prefer_symbol": cfg.gmt_prefer_symbol,
            "require_symbol": cfg.gmt_require_symbol,
            "biotype_allowlist": gmt_biotype_allowlist,
            "min_genes": cfg.gmt_min_genes,
            "max_genes": cfg.gmt_max_genes,
            "emit_small_gene_sets": cfg.emit_small_gene_sets,
            "requested_outputs": requested_gmt_outputs,
            "emitted_outputs": emitted_combinations,
            "skipped_outputs": skipped_gmt_outputs,
            "artifact_diagnostics": artifact_diagnostics,
            "diagnostics": gmt_diagnostics,
            "plans": gmt_plans,
        },
    )

    if resources_info:
        meta["resources"] = resources_info
    write_metadata(out_dir / "geneset.meta.json", meta)

    return {
        "n_peaks": len(peaks),
        "n_genes": len(primary_selected_rows),
        "out_dir": str(out_dir),
        "program_methods": sorted(scored_rows_by_program_by_link.keys()),
        "program_methods_skipped": program_methods_skipped,
        "calibration_methods": ["none"],
        "calibration_methods_skipped": {},
    }
