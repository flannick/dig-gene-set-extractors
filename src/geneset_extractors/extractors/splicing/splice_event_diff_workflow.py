from __future__ import annotations

import csv
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import sys
from typing import Iterable

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


GLOBAL_DEFAULT_COLUMNS = {
    "event_id": ("event_id", "event", "splice_event_id", "spliceseq_event_id", "as_id", "Event"),
    "event_group": ("event_group", "cluster_id", "cluster", "lsv_id", "event_cluster"),
    "event_type": ("event_type", "event_class", "splice_type", "Type"),
    "gene_id": ("gene_id", "ensembl_gene_id", "GeneID"),
    "gene_symbol": ("gene_symbol", "gene_name", "gene", "Gene", "symbol"),
    "chrom": ("chrom", "chr", "Chromosome"),
    "start": ("start", "donor_start", "Start"),
    "end": ("end", "acceptor_end", "End"),
    "strand": ("strand", "Strand"),
    "score": ("score", "custom_score"),
    "stat": ("stat", "t", "z", "statistic", "dPSI_stat"),
    "delta_psi": ("delta_psi", "dpsi", "DeltaPsi", "dPSI", "deltaPsi"),
    "psi": ("psi", "PSI", "mean_psi"),
    "padj": ("padj", "fdr", "qvalue", "adj_pvalue", "QValue"),
    "pvalue": ("pvalue", "p_value", "pval", "PValue", "probability_pvalue"),
    "probability": ("probability", "posterior", "prob", "Probability", "confidence"),
    "read_support": ("read_support", "junction_reads", "coverage", "counts"),
    "novel_flag": ("novel_flag", "novel", "is_novel"),
    "annotation_status": ("annotation_status", "annotation", "status"),
}

TOOL_DEFAULT_COLUMNS = {
    "generic": {},
    "leafcutter": {
        "event_id": ("intron", "event_id"),
        "event_group": ("cluster", "cluster_id", "event_group"),
        "gene_symbol": ("gene_symbol", "gene", "Gene"),
        "gene_id": ("gene_id",),
        "chrom": ("chrom", "chr"),
        "start": ("start",),
        "end": ("end",),
        "strand": ("strand",),
        "delta_psi": ("delta_psi", "deltapsi", "dpsi"),
        "padj": ("padj", "fdr", "p.adjust"),
        "pvalue": ("pvalue", "pval"),
        "read_support": ("read_support", "counts", "junction_reads"),
        "event_type": ("event_type",),
    },
    "majiq": {
        "event_id": ("event_id", "junction_id", "lsv_id"),
        "event_group": ("lsv_id", "event_group"),
        "gene_symbol": ("gene_symbol", "gene_name", "Gene"),
        "delta_psi": ("delta_psi", "dpsi", "DeltaPsi"),
        "probability": ("probability", "Probability", "posterior"),
        "event_type": ("event_type", "splice_type"),
    },
    "whippet": {
        "event_id": ("Event", "event_id"),
        "gene_symbol": ("Gene", "gene_symbol"),
        "event_type": ("Type", "event_type"),
        "delta_psi": ("DeltaPsi", "delta_psi"),
        "probability": ("Probability", "probability"),
        "psi": ("Psi", "psi"),
        "annotation_status": ("Complexity", "annotation_status"),
    },
    "tcga_spliceseq": {
        "event_id": ("event_id", "splice_event_id", "Event"),
        "gene_symbol": ("gene_symbol", "gene_name", "Gene"),
        "event_type": ("event_type", "splice_type"),
        "chrom": ("chrom", "chr"),
        "start": ("start",),
        "end": ("end",),
        "strand": ("strand",),
        "psi": ("psi", "PSI"),
        "annotation_status": ("annotation_status", "annotation"),
    },
}


@dataclass(frozen=True)
class SpliceRow:
    line_no: int
    values: dict[str, str]


@dataclass
class SpliceEventRecord:
    canonical_event_key: str
    event_group: str
    canonicalization_status: str
    canonicalization_confidence: str
    event_key_namespace: str
    event_type: str
    gene_id: str
    gene_symbol: str
    chrom: str
    start: str
    end: str
    strand: str
    base_score: float
    transformed_score: float
    final_score: float
    confidence_weight: float
    impact_weight: float
    ubiquity_weight: float
    delta_psi_effect_weight: float
    prior_confidence_weight: float
    quality_value: float
    source_line_no: int


@dataclass
class SplicingWorkflowConfig:
    converter_name: str
    out_dir: Path
    organism: str
    genome_build: str
    dataset_label: str
    signature_name: str
    tool_family: str
    event_id_column: str | None
    event_group_column: str | None
    event_type_column: str | None
    gene_id_column: str | None
    gene_symbol_column: str | None
    chrom_column: str | None
    start_column: str | None
    end_column: str | None
    strand_column: str | None
    score_column: str | None
    stat_column: str | None
    delta_psi_column: str | None
    psi_column: str | None
    padj_column: str | None
    pvalue_column: str | None
    probability_column: str | None
    read_support_column: str | None
    novel_flag_column: str | None
    annotation_status_column: str | None
    score_mode: str
    score_transform: str
    confidence_weight_mode: str
    min_probability: float
    min_read_support: float
    neglog10p_cap: float
    neglog10p_eps: float
    delta_psi_soft_floor: float
    delta_psi_soft_floor_mode: str
    event_dup_policy: str
    gene_aggregation: str
    gene_topk_events: int
    gene_burden_penalty_mode: str
    min_gene_burden_penalty: float
    gene_support_penalty_mode: str
    locus_density_penalty_mode: str
    locus_density_window_bp: int
    locus_density_top_n: int
    source_dataset: str | None
    bundle_same_dataset_policy: str
    ambiguous_gene_policy: str
    impact_mode: str
    impact_min: float
    impact_max: float
    select: str
    top_k: int
    quantile: float
    min_score: float
    normalize: str
    emit_full: bool
    emit_gmt: bool
    gmt_out: str | None
    gmt_format: str
    gmt_prefer_symbol: bool
    gmt_require_symbol: bool
    gmt_biotype_allowlist: str
    gmt_min_genes: int
    gmt_max_genes: int
    gmt_topk_list: str
    gmt_mass_list: str
    gmt_split_signed: bool
    emit_small_gene_sets: bool


@dataclass
class AliasEntry:
    canonical_event_key: str
    canonicalization_status: str
    canonicalization_confidence: str
    event_key_namespace: str
    gene_id: str
    gene_symbol: str
    event_type: str
    chrom: str
    start: str
    end: str
    strand: str


@dataclass
class UbiquityEntry:
    canonical_event_key: str
    canonicalization_status: str
    canonicalization_confidence: str
    event_key_namespace: str
    gene_id: str
    gene_symbol: str
    event_type: str
    idf_ref: float
    df_ref: float | None
    n_samples_ref: float | None
    n_datasets_ref: float | None


@dataclass
class ImpactEntry:
    canonical_event_key: str
    canonicalization_status: str
    canonicalization_confidence: str
    event_key_namespace: str
    gene_id: str
    gene_symbol: str
    event_type: str
    impact_weight_raw: float
    impact_evidence: float
    annotation_status: str
    n_datasets_ref: float | None


@dataclass
class GeneBurdenEntry:
    gene_symbol: str
    n_canonical_events_ref: float
    n_high_confidence_events_ref: float | None
    n_medium_confidence_events_ref: float | None
    n_low_confidence_events_ref: float | None
    n_unique_event_groups_ref: float | None
    n_studies_ref: float | None
    n_studies_high_confidence_ref: float | None
    fraction_low_confidence_events_ref: float | None
    median_unique_groups_per_study: float | None


@dataclass
class UbiquityByDatasetEntry:
    canonical_event_key: str
    source_dataset: str
    canonicalization_status: str
    canonicalization_confidence: str
    event_key_namespace: str
    gene_id: str
    gene_symbol: str
    event_type: str
    df_ref: float | None
    n_samples_ref: float | None


@dataclass
class GeneBurdenByDatasetEntry:
    gene_symbol: str
    source_dataset: str
    n_canonical_events_ref: float
    n_high_confidence_events_ref: float | None
    n_medium_confidence_events_ref: float | None
    n_low_confidence_events_ref: float | None
    n_unique_event_groups_ref: float | None


@dataclass
class GeneLocusEntry:
    gene_id: str
    gene_symbol: str
    chrom: str
    start: float | None
    end: float | None
    chromosome_arm: str


EVENT_TYPE_ALIASES = {
    "ri": "retained_intron",
    "retained_intron": "retained_intron",
    "retained intron": "retained_intron",
    "intron_retention": "retained_intron",
    "ir": "retained_intron",
    "se": "exon_skip",
    "skipped_exon": "exon_skip",
    "cassette_exon": "exon_skip",
    "cassette exon": "exon_skip",
    "exon_skip": "exon_skip",
    "es": "exon_skip",
    "a5ss": "alt_donor",
    "alt_donor": "alt_donor",
    "alternative_5prime": "alt_donor",
    "alternative_5'_splice_site": "alt_donor",
    "a3ss": "alt_acceptor",
    "alt_acceptor": "alt_acceptor",
    "alternative_3prime": "alt_acceptor",
    "alternative_3'_splice_site": "alt_acceptor",
    "mxe": "mutually_exclusive_exon",
    "mutually_exclusive_exon": "mutually_exclusive_exon",
    "afe": "alt_first_exon",
    "ale": "alt_last_exon",
}


def _open_text(path: str | Path):
    p = Path(path)
    if p.suffix.lower() == ".gz":
        return gzip.open(p, "rt", encoding="utf-8")
    return p.open("r", encoding="utf-8")


def read_tsv_rows(path: str | Path) -> tuple[list[str], list[SpliceRow]]:
    rows: list[SpliceRow] = []
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Splicing table has no header: {path}")
        fieldnames = [str(x) for x in reader.fieldnames]
        for idx, row in enumerate(reader, start=2):
            rows.append(SpliceRow(line_no=idx, values={k: str(v) for k, v in row.items() if k is not None}))
    return fieldnames, rows


def _clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _parse_float(value: object) -> float | None:
    text = _clean(value)
    if not text or text.lower() in {"na", "nan", "none", "null", "inf", "-inf"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _parse_boolish(value: object) -> bool:
    text = _clean(value).lower()
    return text in {"1", "true", "yes", "y", "novel"}


def _normalize_token(value: str) -> str:
    out = _clean(value).lower().replace("-", "_").replace("/", "_").replace(" ", "_")
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_")


def _normalize_event_type(value: str) -> str:
    token = _normalize_token(value)
    return EVENT_TYPE_ALIASES.get(token, token or "unknown")


def _tool_default_candidates(tool_family: str, logical_name: str) -> tuple[str, ...]:
    defaults = TOOL_DEFAULT_COLUMNS.get(tool_family, {})
    candidates = list(defaults.get(logical_name, ()))
    for candidate in GLOBAL_DEFAULT_COLUMNS.get(logical_name, ()):
        if candidate not in candidates:
            candidates.append(candidate)
    return tuple(candidates)


def _resolve_column(fieldnames: list[str], explicit: str | None, logical_name: str, tool_family: str, required: bool = False) -> str | None:
    name = _clean(explicit) or None
    if name is not None:
        if name not in fieldnames:
            raise ValueError(
                f"Column '{name}' for {logical_name} not found. Available columns: {', '.join(fieldnames)}"
            )
        return name
    for candidate in _tool_default_candidates(tool_family, logical_name):
        if candidate in fieldnames:
            return candidate
    if required:
        raise ValueError(
            f"Could not resolve required column for {logical_name}. Available columns: {', '.join(fieldnames)}"
        )
    return None


def detect_tool_family(fieldnames: list[str]) -> str:
    names = set(fieldnames)
    if {"Event", "Gene", "DeltaPsi", "Probability"}.issubset(names):
        return "whippet"
    if {"lsv_id", "dpsi"}.issubset(names) or {"event_id", "probability", "delta_psi"}.issubset(names):
        return "majiq"
    if {"cluster", "deltapsi"}.issubset(names) or {"cluster_id", "delta_psi"}.issubset(names):
        return "leafcutter"
    if {"splice_event_id", "PSI"}.issubset(names) or {"event_id", "psi"}.issubset(names):
        return "tcga_spliceseq"
    return "generic"


def resolve_splicing_columns(fieldnames: list[str], cfg: SplicingWorkflowConfig) -> dict[str, str | None]:
    tool_family = cfg.tool_family if cfg.tool_family != "auto" else detect_tool_family(fieldnames)
    return {
        "tool_family": tool_family,
        "event_id": _resolve_column(fieldnames, cfg.event_id_column, "event_id", tool_family, required=False),
        "event_group": _resolve_column(fieldnames, cfg.event_group_column, "event_group", tool_family, required=False),
        "event_type": _resolve_column(fieldnames, cfg.event_type_column, "event_type", tool_family, required=False),
        "gene_id": _resolve_column(fieldnames, cfg.gene_id_column, "gene_id", tool_family, required=False),
        "gene_symbol": _resolve_column(fieldnames, cfg.gene_symbol_column, "gene_symbol", tool_family, required=False),
        "chrom": _resolve_column(fieldnames, cfg.chrom_column, "chrom", tool_family, required=False),
        "start": _resolve_column(fieldnames, cfg.start_column, "start", tool_family, required=False),
        "end": _resolve_column(fieldnames, cfg.end_column, "end", tool_family, required=False),
        "strand": _resolve_column(fieldnames, cfg.strand_column, "strand", tool_family, required=False),
        "score": _clean(cfg.score_column) or None,
        "stat": _resolve_column(fieldnames, cfg.stat_column, "stat", tool_family, required=False),
        "delta_psi": _resolve_column(fieldnames, cfg.delta_psi_column, "delta_psi", tool_family, required=False),
        "psi": _resolve_column(fieldnames, cfg.psi_column, "psi", tool_family, required=False),
        "padj": _resolve_column(fieldnames, cfg.padj_column, "padj", tool_family, required=False),
        "pvalue": _resolve_column(fieldnames, cfg.pvalue_column, "pvalue", tool_family, required=False),
        "probability": _resolve_column(fieldnames, cfg.probability_column, "probability", tool_family, required=False),
        "read_support": _resolve_column(fieldnames, cfg.read_support_column, "read_support", tool_family, required=False),
        "novel_flag": _resolve_column(fieldnames, cfg.novel_flag_column, "novel_flag", tool_family, required=False),
        "annotation_status": _resolve_column(fieldnames, cfg.annotation_status_column, "annotation_status", tool_family, required=False),
    }


def read_event_alias_tsv(path: str | Path) -> dict[str, AliasEntry]:
    out: dict[str, AliasEntry] = {}
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Alias TSV has no header: {path}")
        if "input_event_key" not in reader.fieldnames or "canonical_event_key" not in reader.fieldnames:
            raise ValueError("Alias TSV must include input_event_key and canonical_event_key")
        for row in reader:
            input_key = _clean(row.get("input_event_key"))
            canonical = _clean(row.get("canonical_event_key"))
            if not input_key or not canonical:
                continue
            out[input_key] = AliasEntry(
                canonical_event_key=canonical,
                canonicalization_status=_clean(row.get("canonicalization_status")) or ("coordinate_canonical" if canonical.startswith("chr") else "raw_id_fallback"),
                canonicalization_confidence=_clean(row.get("canonicalization_confidence")).lower() or ("high" if canonical.startswith("chr") else "low"),
                event_key_namespace=_clean(row.get("event_key_namespace")) or ("global_coordinate" if canonical.startswith("chr") else "raw_id_fallback"),
                gene_id=_clean(row.get("gene_id")),
                gene_symbol=_clean(row.get("gene_symbol")),
                event_type=_normalize_event_type(_clean(row.get("event_type"))),
                chrom=_clean(row.get("chrom")),
                start=_clean(row.get("start")),
                end=_clean(row.get("end")),
                strand=_clean(row.get("strand")),
            )
    return out


def read_event_ubiquity_tsv(path: str | Path) -> dict[str, UbiquityEntry]:
    out: dict[str, UbiquityEntry] = {}
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Ubiquity TSV has no header: {path}")
        if "canonical_event_key" not in reader.fieldnames:
            raise ValueError("Ubiquity TSV must include canonical_event_key")
        for row in reader:
            canonical = _clean(row.get("canonical_event_key"))
            if not canonical:
                continue
            idf = _parse_float(row.get("idf_ref"))
            df_ref = _parse_float(row.get("df_ref"))
            n_ref = _parse_float(row.get("n_samples_ref"))
            if idf is None:
                if df_ref is not None and n_ref is not None:
                    idf = math.log1p((n_ref + 1.0) / (df_ref + 1.0))
                else:
                    idf = 1.0
            out[canonical] = UbiquityEntry(
                canonical_event_key=canonical,
                canonicalization_status=_clean(row.get("canonicalization_status")) or ("coordinate_canonical" if canonical.startswith("chr") else "raw_id_fallback"),
                canonicalization_confidence=_clean(row.get("canonicalization_confidence")).lower() or ("high" if canonical.startswith("chr") else "low"),
                event_key_namespace=_clean(row.get("event_key_namespace")) or ("global_coordinate" if canonical.startswith("chr") else "raw_id_fallback"),
                gene_id=_clean(row.get("gene_id")),
                gene_symbol=_clean(row.get("gene_symbol")),
                event_type=_normalize_event_type(_clean(row.get("event_type"))),
                idf_ref=float(idf),
                df_ref=df_ref,
                n_samples_ref=n_ref,
                n_datasets_ref=_parse_float(row.get("n_datasets_ref")),
            )
    return out


def read_event_impact_tsv(path: str | Path) -> dict[str, ImpactEntry]:
    out: dict[str, ImpactEntry] = {}
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Impact TSV has no header: {path}")
        if "canonical_event_key" not in reader.fieldnames:
            raise ValueError("Impact TSV must include canonical_event_key")
        for row in reader:
            canonical = _clean(row.get("canonical_event_key"))
            if not canonical:
                continue
            out[canonical] = ImpactEntry(
                canonical_event_key=canonical,
                canonicalization_status=_clean(row.get("canonicalization_status")) or ("coordinate_canonical" if canonical.startswith("chr") else "raw_id_fallback"),
                canonicalization_confidence=_clean(row.get("canonicalization_confidence")).lower() or ("high" if canonical.startswith("chr") else "low"),
                event_key_namespace=_clean(row.get("event_key_namespace")) or ("global_coordinate" if canonical.startswith("chr") else "raw_id_fallback"),
                gene_id=_clean(row.get("gene_id")),
                gene_symbol=_clean(row.get("gene_symbol")),
                event_type=_normalize_event_type(_clean(row.get("event_type"))),
                impact_weight_raw=float(_parse_float(row.get("impact_weight_raw")) or 1.0),
                impact_evidence=max(0.0, min(1.0, float(_parse_float(row.get("impact_evidence")) or 0.0))),
                annotation_status=_clean(row.get("annotation_status")),
                n_datasets_ref=_parse_float(row.get("n_datasets_ref")),
            )
    return out


def read_gene_burden_tsv(path: str | Path) -> dict[str, GeneBurdenEntry]:
    out: dict[str, GeneBurdenEntry] = {}
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Gene burden TSV has no header: {path}")
        if "gene_symbol" not in reader.fieldnames:
            raise ValueError("Gene burden TSV must include gene_symbol")
        for row in reader:
            gene_symbol = _clean(row.get("gene_symbol"))
            if not gene_symbol:
                continue
            out[gene_symbol] = GeneBurdenEntry(
                gene_symbol=gene_symbol,
                n_canonical_events_ref=float(_parse_float(row.get("n_canonical_events_ref")) or 0.0),
                n_high_confidence_events_ref=_parse_float(row.get("n_high_confidence_events_ref")),
                n_medium_confidence_events_ref=_parse_float(row.get("n_medium_confidence_events_ref")),
                n_low_confidence_events_ref=_parse_float(row.get("n_low_confidence_events_ref")),
                n_unique_event_groups_ref=_parse_float(row.get("n_unique_event_groups_ref")),
                n_studies_ref=_parse_float(row.get("n_studies_ref")),
                n_studies_high_confidence_ref=_parse_float(row.get("n_studies_high_confidence_ref")),
                fraction_low_confidence_events_ref=_parse_float(row.get("fraction_low_confidence_events_ref")),
                median_unique_groups_per_study=_parse_float(row.get("median_unique_groups_per_study")),
            )
    return out


def read_event_ubiquity_by_dataset_tsv(path: str | Path) -> dict[str, list[UbiquityByDatasetEntry]]:
    out: dict[str, list[UbiquityByDatasetEntry]] = {}
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Dataset ubiquity TSV has no header: {path}")
        if "canonical_event_key" not in reader.fieldnames or "source_dataset" not in reader.fieldnames:
            raise ValueError("Dataset ubiquity TSV must include canonical_event_key and source_dataset")
        for row in reader:
            canonical = _clean(row.get("canonical_event_key"))
            source_dataset = _clean(row.get("source_dataset"))
            if not canonical or not source_dataset:
                continue
            out.setdefault(canonical, []).append(
                UbiquityByDatasetEntry(
                    canonical_event_key=canonical,
                    source_dataset=source_dataset,
                    canonicalization_status=_clean(row.get("canonicalization_status")) or ("coordinate_canonical" if canonical.startswith("chr") else "raw_id_fallback"),
                    canonicalization_confidence=_clean(row.get("canonicalization_confidence")).lower() or ("high" if canonical.startswith("chr") else "low"),
                    event_key_namespace=_clean(row.get("event_key_namespace")) or ("global_coordinate" if canonical.startswith("chr") else "raw_id_fallback"),
                    gene_id=_clean(row.get("gene_id")),
                    gene_symbol=_clean(row.get("gene_symbol")),
                    event_type=_normalize_event_type(_clean(row.get("event_type"))),
                    df_ref=_parse_float(row.get("df_ref")),
                    n_samples_ref=_parse_float(row.get("n_samples_ref")),
                )
            )
    return out


def read_gene_burden_by_dataset_tsv(path: str | Path) -> dict[str, list[GeneBurdenByDatasetEntry]]:
    out: dict[str, list[GeneBurdenByDatasetEntry]] = {}
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Dataset burden TSV has no header: {path}")
        if "gene_symbol" not in reader.fieldnames or "source_dataset" not in reader.fieldnames:
            raise ValueError("Dataset burden TSV must include gene_symbol and source_dataset")
        for row in reader:
            gene_symbol = _clean(row.get("gene_symbol"))
            source_dataset = _clean(row.get("source_dataset"))
            if not gene_symbol or not source_dataset:
                continue
            out.setdefault(gene_symbol, []).append(
                GeneBurdenByDatasetEntry(
                    gene_symbol=gene_symbol,
                    source_dataset=source_dataset,
                    n_canonical_events_ref=float(_parse_float(row.get("n_canonical_events_ref")) or 0.0),
                    n_high_confidence_events_ref=_parse_float(row.get("n_high_confidence_events_ref")),
                    n_medium_confidence_events_ref=_parse_float(row.get("n_medium_confidence_events_ref")),
                    n_low_confidence_events_ref=_parse_float(row.get("n_low_confidence_events_ref")),
                    n_unique_event_groups_ref=_parse_float(row.get("n_unique_event_groups_ref")),
                )
            )
    return out


def read_gene_locus_tsv(path: str | Path) -> dict[str, GeneLocusEntry]:
    out: dict[str, GeneLocusEntry] = {}
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Gene locus TSV has no header: {path}")
        for row in reader:
            gene_symbol = _clean(row.get("gene_symbol"))
            gene_id = _clean(row.get("gene_id")) or gene_symbol
            chrom = _clean(row.get("chrom"))
            if not gene_symbol or not chrom:
                continue
            entry = GeneLocusEntry(
                gene_id=gene_id,
                gene_symbol=gene_symbol,
                chrom=chrom,
                start=_parse_float(row.get("start")),
                end=_parse_float(row.get("end")),
                chromosome_arm=_clean(row.get("chromosome_arm")),
            )
            out[gene_symbol] = entry
            if gene_id:
                out.setdefault(gene_id, entry)
    return out


def _split_gene_tokens(value: str) -> list[str]:
    text = _clean(value)
    if not text:
        return []
    for delim in (";", ",", "|", "/"):
        text = text.replace(delim, ";")
    return [tok.strip() for tok in text.split(";") if tok.strip()]


def _resolve_gene_assignments(*, raw_gene_id: str, raw_gene_symbol: str, alias_gene_id: str, alias_gene_symbol: str, policy: str) -> tuple[list[tuple[str, str, float]], int]:
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
            pairs.append((gid, gsym or gid))
    ambiguous = 1 if len(pairs) > 1 else 0
    if len(pairs) == 1:
        gid, gsym = pairs[0]
        return [(gid, gsym, 1.0)], ambiguous
    if policy == "drop":
        return [], ambiguous
    if policy == "first":
        gid, gsym = sorted(pairs)[0]
        return [(gid, gsym, 1.0)], ambiguous
    if policy == "split_equal":
        frac = 1.0 / float(len(pairs))
        return [(gid, gsym, frac) for gid, gsym in sorted(pairs)], ambiguous
    raise ValueError(f"Unsupported ambiguous_gene_policy: {policy}")


def _string_value(row: SpliceRow, columns: dict[str, str | None], logical_name: str) -> str:
    column = columns.get(logical_name)
    if not column:
        return ""
    return _clean(row.values.get(column))


def _coordinate_event_key(row: SpliceRow, columns: dict[str, str | None]) -> str | None:
    chrom = _string_value(row, columns, "chrom")
    start = _string_value(row, columns, "start")
    end = _string_value(row, columns, "end")
    strand = _string_value(row, columns, "strand") or "."
    event_type = _normalize_event_type(_string_value(row, columns, "event_type"))
    event_group = _string_value(row, columns, "event_group")
    if chrom and start and end:
        tail = event_group or _string_value(row, columns, "event_id")
        if tail:
            return f"{chrom}:{start}-{end}:{strand}:{event_type}:{tail}"
        return f"{chrom}:{start}-{end}:{strand}:{event_type}"
    return None


def _canonical_event_candidates(row: SpliceRow, columns: dict[str, str | None]) -> list[str]:
    out: list[str] = []
    event_id = _string_value(row, columns, "event_id")
    event_group = _string_value(row, columns, "event_group")
    if event_id:
        out.append(event_id)
    if event_group:
        out.append(event_group)
    coord_key = _coordinate_event_key(row, columns)
    if coord_key:
        out.append(coord_key)
    return out


def _default_canonical_event_key(row: SpliceRow, columns: dict[str, str | None]) -> str | None:
    coord_key = _coordinate_event_key(row, columns)
    if coord_key:
        return coord_key
    tool_family = str(columns.get("tool_family", "generic"))
    event_id = _string_value(row, columns, "event_id")
    gene_symbol = _string_value(row, columns, "gene_symbol") or _string_value(row, columns, "gene_id")
    event_type = _normalize_event_type(_string_value(row, columns, "event_type"))
    if tool_family == "tcga_spliceseq" and event_id and gene_symbol and event_type and event_type != "unknown":
        return f"tcga_spliceseq_asid::{gene_symbol}::{event_type}::{event_id}"
    event_id = _string_value(row, columns, "event_id") or _string_value(row, columns, "event_group")
    return event_id or None


def _default_canonicalization_status(row: SpliceRow, columns: dict[str, str | None]) -> tuple[str, str, str]:
    coord_key = _coordinate_event_key(row, columns)
    if coord_key:
        return "coordinate_canonical", "high", "global_coordinate"
    tool_family = str(columns.get("tool_family", "generic"))
    event_id = _string_value(row, columns, "event_id")
    gene_symbol = _string_value(row, columns, "gene_symbol") or _string_value(row, columns, "gene_id")
    event_type = _normalize_event_type(_string_value(row, columns, "event_type"))
    if tool_family == "tcga_spliceseq" and event_id and gene_symbol and event_type and event_type != "unknown":
        return "source_family_stable_id", "medium", "tcga_spliceseq_asid"
    if _string_value(row, columns, "event_id") or _string_value(row, columns, "event_group"):
        return "raw_id_fallback", "low", "raw_id_fallback"
    return "unresolved", "low", "raw_id_fallback"


def _prior_confidence_weight(canonicalization_confidence: str) -> float:
    confidence = str(canonicalization_confidence).lower()
    if confidence == "high":
        return 1.0
    if confidence == "medium":
        return 0.6
    return 0.25


def _shrink_prior_toward_neutral(value: float, canonicalization_confidence: str) -> float:
    weight = _prior_confidence_weight(canonicalization_confidence)
    return 1.0 + weight * (float(value) - 1.0)


def _event_group_key(record: SpliceEventRecord) -> str:
    return record.event_group or record.canonical_event_key


def _namespace_matches(runtime_namespace: str, prior_namespace: str) -> bool:
    runtime_value = _clean(runtime_namespace)
    prior_value = _clean(prior_namespace)
    if not runtime_value or not prior_value:
        return runtime_value == prior_value
    return runtime_value == prior_value


def _build_effective_ubiquity_map(
    *,
    pooled_map: dict[str, UbiquityEntry] | None,
    by_dataset_map: dict[str, list[UbiquityByDatasetEntry]] | None,
    source_dataset: str | None,
    policy: str,
) -> tuple[dict[str, UbiquityEntry] | None, dict[str, float], dict[str, object]]:
    if pooled_map is None:
        return None, {}, {"policy": policy, "source_dataset": source_dataset, "dataset_level_available": False}
    effective = dict(pooled_map)
    recurrence: dict[str, float] = {key: float(entry.n_datasets_ref or 0.0) for key, entry in pooled_map.items()}
    summary: dict[str, object] = {
        "policy": policy,
        "source_dataset": source_dataset,
        "dataset_level_available": bool(by_dataset_map),
        "bundle_source_datasets": sorted({entry.source_dataset for entries in (by_dataset_map or {}).values() for entry in entries}),
        "n_events_excluded": 0,
        "n_events_neutralized": 0,
        "warning": "",
    }
    target = _clean(source_dataset)
    if not target or not by_dataset_map or policy in {"ignore", "warn"}:
        if target and by_dataset_map and target in set(summary["bundle_source_datasets"]) and policy == "warn":
            summary["warning"] = f"bundle_same_dataset_match_detected source_dataset={target} action=warn"
        return effective, recurrence, summary
    bundle_datasets = set(summary["bundle_source_datasets"])
    if target not in bundle_datasets:
        return effective, recurrence, summary
    if policy == "fail":
        raise ValueError(
            f"Bundle contains matching source_dataset={target} and --bundle_same_dataset_policy=fail. "
            "Use a leave-one-dataset-out bundle or change the policy."
        )
    n_total_excluding_self = max(1, len(bundle_datasets - {target}))
    for canonical, pooled_entry in pooled_map.items():
        entries = list(by_dataset_map.get(canonical, []))
        if not any(entry.source_dataset == target for entry in entries):
            continue
        remaining = [entry for entry in entries if entry.source_dataset != target]
        if not remaining:
            effective[canonical] = UbiquityEntry(
                canonical_event_key=pooled_entry.canonical_event_key,
                canonicalization_status=pooled_entry.canonicalization_status,
                canonicalization_confidence=pooled_entry.canonicalization_confidence,
                event_key_namespace=pooled_entry.event_key_namespace,
                gene_id=pooled_entry.gene_id,
                gene_symbol=pooled_entry.gene_symbol,
                event_type=pooled_entry.event_type,
                idf_ref=1.0,
                df_ref=None,
                n_samples_ref=None,
                n_datasets_ref=0.0,
            )
            recurrence[canonical] = 0.0
            summary["n_events_neutralized"] = int(summary["n_events_neutralized"]) + 1
            continue
        n_event_datasets = len({entry.source_dataset for entry in remaining})
        df_ref = sum(float(entry.df_ref or 0.0) for entry in remaining)
        n_samples_ref = sum(float(entry.n_samples_ref or 0.0) for entry in remaining)
        effective[canonical] = UbiquityEntry(
            canonical_event_key=pooled_entry.canonical_event_key,
            canonicalization_status=pooled_entry.canonicalization_status,
            canonicalization_confidence=pooled_entry.canonicalization_confidence,
            event_key_namespace=pooled_entry.event_key_namespace,
            gene_id=pooled_entry.gene_id,
            gene_symbol=pooled_entry.gene_symbol,
            event_type=pooled_entry.event_type,
            idf_ref=math.log1p((float(n_total_excluding_self) + 1.0) / (float(n_event_datasets) + 1.0)),
            df_ref=df_ref,
            n_samples_ref=n_samples_ref,
            n_datasets_ref=float(n_event_datasets),
        )
        recurrence[canonical] = float(n_event_datasets)
        summary["n_events_excluded"] = int(summary["n_events_excluded"]) + 1
    if bundle_datasets == {target}:
        summary["warning"] = f"bundle_same_dataset_only source_dataset={target} action=neutralize"
    return effective, recurrence, summary


def _build_effective_gene_burden_map(
    *,
    pooled_map: dict[str, GeneBurdenEntry] | None,
    by_dataset_map: dict[str, list[GeneBurdenByDatasetEntry]] | None,
    source_dataset: str | None,
    policy: str,
) -> tuple[dict[str, GeneBurdenEntry] | None, dict[str, object]]:
    if pooled_map is None:
        return None, {"policy": policy, "source_dataset": source_dataset, "dataset_level_available": False}
    effective = dict(pooled_map)
    summary: dict[str, object] = {
        "policy": policy,
        "source_dataset": source_dataset,
        "dataset_level_available": bool(by_dataset_map),
        "n_genes_excluded": 0,
        "n_genes_neutralized": 0,
    }
    target = _clean(source_dataset)
    if not target or not by_dataset_map or policy in {"ignore", "warn"}:
        return effective, summary
    if policy == "fail":
        bundle_datasets = {entry.source_dataset for entries in by_dataset_map.values() for entry in entries}
        if target in bundle_datasets:
            raise ValueError(
                f"Bundle contains matching source_dataset={target} and --bundle_same_dataset_policy=fail. "
                "Use a leave-one-dataset-out bundle or change the policy."
            )
        return effective, summary
    for gene_symbol, pooled_entry in pooled_map.items():
        entries = list(by_dataset_map.get(gene_symbol, []))
        if not any(entry.source_dataset == target for entry in entries):
            continue
        remaining = [entry for entry in entries if entry.source_dataset != target]
        if not remaining:
            effective[gene_symbol] = GeneBurdenEntry(
                gene_symbol=pooled_entry.gene_symbol,
                n_canonical_events_ref=0.0,
                n_high_confidence_events_ref=0.0,
                n_medium_confidence_events_ref=0.0,
                n_low_confidence_events_ref=0.0,
                n_unique_event_groups_ref=0.0,
                n_studies_ref=0.0,
                n_studies_high_confidence_ref=0.0,
                fraction_low_confidence_events_ref=0.0,
                median_unique_groups_per_study=0.0,
            )
            summary["n_genes_neutralized"] = int(summary["n_genes_neutralized"]) + 1
            continue
        unique_groups = [float(entry.n_unique_event_groups_ref or 0.0) for entry in remaining]
        total_events = sum(float(entry.n_canonical_events_ref or 0.0) for entry in remaining)
        low_events = sum(float(entry.n_low_confidence_events_ref or 0.0) for entry in remaining)
        effective[gene_symbol] = GeneBurdenEntry(
            gene_symbol=pooled_entry.gene_symbol,
            n_canonical_events_ref=total_events,
            n_high_confidence_events_ref=sum(float(entry.n_high_confidence_events_ref or 0.0) for entry in remaining),
            n_medium_confidence_events_ref=sum(float(entry.n_medium_confidence_events_ref or 0.0) for entry in remaining),
            n_low_confidence_events_ref=low_events,
            n_unique_event_groups_ref=sum(unique_groups),
            n_studies_ref=float(len(remaining)),
            n_studies_high_confidence_ref=float(sum(1 for entry in remaining if float(entry.n_high_confidence_events_ref or 0.0) > 0.0)),
            fraction_low_confidence_events_ref=(low_events / float(max(1.0, total_events))),
            median_unique_groups_per_study=float(median(unique_groups)) if unique_groups else 0.0,
        )
        summary["n_genes_excluded"] = int(summary["n_genes_excluded"]) + 1
    return effective, summary


def _delta_psi_effect_weight(row: SpliceRow, columns: dict[str, str | None], cfg: SplicingWorkflowConfig) -> tuple[float, bool]:
    if cfg.delta_psi_soft_floor_mode in {"off", "none"}:
        return 1.0, False
    delta_col = columns.get("delta_psi")
    if not delta_col:
        return 1.0, False
    delta = _parse_float(row.values.get(delta_col))
    if delta is None:
        return 1.0, False
    floor = max(0.0, float(cfg.delta_psi_soft_floor))
    if floor <= 0.0:
        return 1.0, False
    return min(1.0, abs(float(delta)) / floor), True


def _resolve_score_mode(columns: dict[str, str | None], rows: list[SpliceRow], requested: str) -> str:
    if requested != "auto":
        return requested
    stat_col = columns.get("stat")
    if stat_col:
        for row in rows:
            if _parse_float(row.values.get(stat_col)) is not None:
                return "stat"
    delta_col = columns.get("delta_psi")
    p_col = columns.get("padj") or columns.get("pvalue")
    if delta_col and p_col:
        for row in rows:
            if _parse_float(row.values.get(delta_col)) is not None and _parse_float(row.values.get(p_col)) is not None:
                return "delta_psi_times_neglog10p"
    if delta_col and (columns.get("probability") or columns.get("read_support")):
        for row in rows:
            if _parse_float(row.values.get(delta_col)) is not None:
                return "delta_psi_times_confidence"
    score_col = columns.get("score")
    if score_col and rows and score_col in rows[0].values:
        return "custom_column"
    raise ValueError(
        "score_mode=auto could not resolve a usable score definition. Provide a stat column, delta_psi with pvalue/padj, "
        "delta_psi with probability/read_support, or use --score_mode custom_column with --score_column."
    )


def _confidence_proxy_from_row(row: SpliceRow, columns: dict[str, str | None], max_read_support: float) -> float:
    prob = _parse_float(row.values.get(columns.get("probability", ""))) if columns.get("probability") else None
    if prob is not None:
        return max(0.0, min(1.0, float(prob)))
    read_support = _parse_float(row.values.get(columns.get("read_support", ""))) if columns.get("read_support") else None
    if read_support is not None and max_read_support > 0:
        return min(1.0, math.log1p(float(read_support)) / math.log1p(float(max_read_support)))
    return 1.0


def _compute_base_score(row: SpliceRow, columns: dict[str, str | None], score_mode: str, neglog10p_eps: float, neglog10p_cap: float, max_read_support: float) -> tuple[float | None, str]:
    if score_mode == "custom_column":
        col = columns.get("score")
        if not col:
            raise ValueError("score_mode=custom_column requires --score_column")
        return _parse_float(row.values.get(col)), "custom_column"
    if score_mode == "stat":
        col = columns.get("stat")
        if not col:
            raise ValueError("score_mode=stat requires a stat column")
        return _parse_float(row.values.get(col)), "stat"
    if score_mode == "delta_psi_times_neglog10p":
        delta_col = columns.get("delta_psi")
        p_col = columns.get("padj") or columns.get("pvalue")
        if not delta_col or not p_col:
            raise ValueError("score_mode=delta_psi_times_neglog10p requires delta_psi and pvalue/padj")
        delta = _parse_float(row.values.get(delta_col))
        pval = _parse_float(row.values.get(p_col))
        if delta is None or pval is None:
            return None, "delta_psi_times_neglog10p"
        neglog = min(-math.log10(max(float(pval), float(neglog10p_eps))), float(neglog10p_cap))
        return float(delta) * float(neglog), "delta_psi_times_neglog10p"
    if score_mode == "delta_psi_times_confidence":
        delta_col = columns.get("delta_psi")
        if not delta_col:
            raise ValueError("score_mode=delta_psi_times_confidence requires delta_psi")
        delta = _parse_float(row.values.get(delta_col))
        if delta is None:
            return None, "delta_psi_times_confidence"
        return float(delta) * _confidence_proxy_from_row(row, columns, max_read_support), "delta_psi_times_confidence"
    raise ValueError(f"Unsupported score_mode: {score_mode}")


def _confidence_weight(row: SpliceRow, *, columns: dict[str, str | None], mode: str, min_probability: float, min_read_support: float, max_read_support: float, neglog10p_eps: float, neglog10p_cap: float) -> tuple[float, dict[str, float]]:
    p_weight = 1.0
    p_col = columns.get("padj") or columns.get("pvalue")
    if p_col:
        pval = _parse_float(row.values.get(p_col))
        if pval is not None:
            p_weight = min(-math.log10(max(float(pval), float(neglog10p_eps))), float(neglog10p_cap)) / float(neglog10p_cap)

    prob_weight = 1.0
    prob_col = columns.get("probability")
    if prob_col:
        prob = _parse_float(row.values.get(prob_col))
        if prob is not None:
            if float(prob) < float(min_probability):
                prob_weight = 0.0
            else:
                denom = max(1e-12, 1.0 - float(min_probability))
                prob_weight = max(0.0, min(1.0, (float(prob) - float(min_probability)) / denom))

    read_weight = 1.0
    read_col = columns.get("read_support")
    if read_col:
        reads = _parse_float(row.values.get(read_col))
        if reads is not None:
            if float(reads) < float(min_read_support):
                read_weight = 0.0
            elif max_read_support > 0:
                read_weight = min(1.0, math.log1p(float(reads)) / math.log1p(float(max_read_support)))

    if mode == "none":
        weight = 1.0
    elif mode == "pvalue":
        weight = p_weight
    elif mode == "probability":
        weight = prob_weight
    elif mode == "read_support":
        weight = read_weight
    elif mode == "combined":
        weight = p_weight * prob_weight * read_weight
    else:
        raise ValueError(f"Unsupported confidence_weight_mode: {mode}")
    return float(weight), {
        "pvalue_weight": float(p_weight),
        "probability_weight": float(prob_weight),
        "read_support_weight": float(read_weight),
    }


def _impact_weight(
    impact_entry: ImpactEntry | None,
    cfg: SplicingWorkflowConfig,
    event_type: str,
    annotation_status: str,
    novel_flag: bool,
    canonicalization_confidence: str,
) -> float:
    if cfg.impact_mode == "none":
        return 1.0
    if impact_entry is not None:
        if (
            str(impact_entry.canonicalization_confidence).lower() != "high"
            and float(impact_entry.n_datasets_ref or 0.0) < 2.0
        ):
            return 1.0
        raw = float(impact_entry.impact_weight_raw)
        evidence = float(impact_entry.impact_evidence)
        base = 1.0 + evidence * (raw - 1.0)
        shrunk = _shrink_prior_toward_neutral(base, impact_entry.canonicalization_confidence or canonicalization_confidence)
        return max(float(cfg.impact_min), min(float(cfg.impact_max), shrunk))
    if cfg.impact_mode in {"conservative", "custom_bundle"}:
        raw = 1.0
        evidence = 0.2
        if event_type in {"exon_skip", "alt_donor", "alt_acceptor", "mutually_exclusive_exon"}:
            raw = 1.1
            evidence = 0.4
        elif event_type == "retained_intron":
            raw = 0.92
            evidence = 0.4
        if "coding" in annotation_status.lower():
            raw = max(raw, 1.12)
            evidence = max(evidence, 0.5)
        if novel_flag:
            raw = min(raw, 0.95)
            evidence = max(evidence, 0.3)
        base = 1.0 + evidence * (raw - 1.0)
        shrunk = _shrink_prior_toward_neutral(base, canonicalization_confidence)
        return max(float(cfg.impact_min), min(float(cfg.impact_max), shrunk))
    raise ValueError(f"Unsupported impact_mode: {cfg.impact_mode}")


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


def _collapse_duplicate_events(records: list[SpliceEventRecord], policy: str) -> list[SpliceEventRecord]:
    grouped: dict[tuple[str, str], list[SpliceEventRecord]] = {}
    for rec in records:
        grouped.setdefault((rec.canonical_event_key, rec.gene_id), []).append(rec)
    out: list[SpliceEventRecord] = []
    for items in grouped.values():
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
                SpliceEventRecord(
                    canonical_event_key=first.canonical_event_key,
                    event_group=first.event_group,
                    canonicalization_status=first.canonicalization_status,
                    canonicalization_confidence=first.canonicalization_confidence,
                    event_type=first.event_type,
                    gene_id=first.gene_id,
                    gene_symbol=first.gene_symbol,
                    chrom=first.chrom,
                    start=first.start,
                    end=first.end,
                    strand=first.strand,
                    base_score=sum(float(x.base_score) for x in items) / denom,
                    transformed_score=sum(float(x.transformed_score) for x in items) / denom,
                    final_score=sum(float(x.final_score) for x in items) / denom,
                    confidence_weight=sum(float(x.confidence_weight) for x in items) / float(len(items)),
                    impact_weight=sum(float(x.impact_weight) for x in items) / float(len(items)),
                    ubiquity_weight=sum(float(x.ubiquity_weight) for x in items) / float(len(items)),
                    delta_psi_effect_weight=sum(float(x.delta_psi_effect_weight) for x in items) / float(len(items)),
                    prior_confidence_weight=sum(float(x.prior_confidence_weight) for x in items) / float(len(items)),
                    quality_value=max(float(x.quality_value) for x in items),
                    source_line_no=min(int(x.source_line_no) for x in items),
                )
            )
            continue
        raise ValueError(f"Unsupported event_dup_policy: {policy}")
    return out


def _collapse_event_groups(records: list[SpliceEventRecord]) -> list[SpliceEventRecord]:
    grouped: dict[tuple[str, str], list[SpliceEventRecord]] = {}
    for rec in records:
        grouped.setdefault((rec.gene_id, _event_group_key(rec)), []).append(rec)
    out: list[SpliceEventRecord] = []
    for items in grouped.values():
        out.append(max(items, key=lambda r: (abs(float(r.final_score)), float(r.quality_value), -int(r.source_line_no))))
    return out


def _aggregate_events_to_genes(records: list[SpliceEventRecord], method: str, topk_events: int) -> tuple[dict[str, float], dict[str, str], dict[str, list[SpliceEventRecord]]]:
    by_gene: dict[str, list[SpliceEventRecord]] = {}
    gene_symbols: dict[str, str] = {}
    for rec in records:
        by_gene.setdefault(rec.gene_id, []).append(rec)
        if rec.gene_symbol:
            gene_symbols.setdefault(rec.gene_id, rec.gene_symbol)
    scores: dict[str, float] = {}
    for gene_id, items in by_gene.items():
        ranked = sorted(items, key=lambda r: (-abs(float(r.final_score)), str(r.canonical_event_key)))
        if method == "signed_topk_mean":
            kept = ranked[: max(1, int(topk_events))]
            scores[gene_id] = sum(float(r.final_score) for r in kept) / float(len(kept))
        elif method == "max_abs":
            scores[gene_id] = float(ranked[0].final_score)
        elif method == "sum":
            scores[gene_id] = sum(float(r.final_score) for r in items)
        elif method == "mean":
            scores[gene_id] = sum(float(r.final_score) for r in items) / float(len(items))
        else:
            raise ValueError(f"Unsupported gene_aggregation: {method}")
    return scores, gene_symbols, by_gene


def _support_subset(records: list[SpliceEventRecord], method: str, topk_events: int) -> list[SpliceEventRecord]:
    ranked = sorted(records, key=lambda r: (-abs(float(r.final_score)), str(r.canonical_event_key)))
    if method == "signed_topk_mean":
        return ranked[: max(1, int(topk_events))]
    if method == "max_abs":
        return ranked[:1]
    if method in {"sum", "mean"}:
        return ranked
    raise ValueError(f"Unsupported gene_aggregation: {method}")


def _apply_gene_support_penalty(
    *,
    gene_scores: dict[str, float],
    by_gene: dict[str, list[SpliceEventRecord]],
    cfg: SplicingWorkflowConfig,
) -> tuple[dict[str, float], dict[str, float], dict[str, int], dict[str, int], dict[str, float]]:
    penalized: dict[str, float] = {}
    penalties: dict[str, float] = {}
    group_counts: dict[str, int] = {}
    event_counts: dict[str, int] = {}
    sign_coherence: dict[str, float] = {}
    for gene_id, score in gene_scores.items():
        support = _support_subset(by_gene.get(gene_id, []), cfg.gene_aggregation, cfg.gene_topk_events)
        n_groups = len({_event_group_key(rec) for rec in support})
        n_events = len(support)
        abs_sum = sum(abs(float(rec.final_score)) for rec in support)
        signed_sum = sum(float(rec.final_score) for rec in support)
        coherence = abs(signed_sum) / abs_sum if abs_sum > 0.0 else 0.0
        if cfg.gene_support_penalty_mode == "none":
            penalty = 1.0
        elif cfg.gene_support_penalty_mode == "independent_groups":
            penalty = min(1.0, math.sqrt(float(max(1, n_groups)) / 2.0))
        elif cfg.gene_support_penalty_mode == "auto":
            group_factor = 1.0 if n_groups >= 2 else 0.85
            coherence_factor = 0.75 + 0.25 * coherence
            penalty = max(0.55, min(1.0, group_factor * coherence_factor))
        else:
            raise ValueError(f"Unsupported gene_support_penalty_mode: {cfg.gene_support_penalty_mode}")
        group_counts[gene_id] = n_groups
        event_counts[gene_id] = n_events
        sign_coherence[gene_id] = coherence
        penalties[gene_id] = float(penalty)
        penalized[gene_id] = float(score) * float(penalty)
    return penalized, penalties, group_counts, event_counts, sign_coherence


def _gene_primary_locus(
    by_gene: dict[str, list[SpliceEventRecord]],
    window_bp: int,
    gene_symbols: dict[str, str],
    gene_locus_map: dict[str, GeneLocusEntry] | None,
) -> dict[str, dict[str, object]]:
    loci: dict[str, dict[str, object]] = {}
    window = max(1, int(window_bp))
    for gene_id, items in by_gene.items():
        if not items:
            continue
        locus_entry = None
        if gene_locus_map:
            locus_entry = gene_locus_map.get(gene_symbols.get(gene_id, gene_id)) or gene_locus_map.get(gene_id)
        if locus_entry is not None:
            midpoint = None
            if locus_entry.start is not None and locus_entry.end is not None:
                midpoint = 0.5 * (float(locus_entry.start) + float(locus_entry.end))
            elif locus_entry.start is not None:
                midpoint = float(locus_entry.start)
            elif locus_entry.end is not None:
                midpoint = float(locus_entry.end)
            locus_key = ""
            if locus_entry.chrom and midpoint is not None:
                locus_key = f"{locus_entry.chrom}:{int(midpoint // float(window))}"
            elif locus_entry.chrom:
                locus_key = locus_entry.chrom
            loci[gene_id] = {
                "chrom": locus_entry.chrom,
                "midpoint": midpoint,
                "locus_key": locus_key,
                "chromosome_arm": locus_entry.chromosome_arm,
            }
            continue
        top = max(items, key=lambda r: (abs(float(r.final_score)), float(r.quality_value), -int(r.source_line_no)))
        chrom = str(top.chrom or "").strip()
        midpoint = None
        start = _parse_float(top.start)
        end = _parse_float(top.end)
        if start is not None and end is not None:
            midpoint = 0.5 * (float(start) + float(end))
        elif start is not None:
            midpoint = float(start)
        elif end is not None:
            midpoint = float(end)
        if chrom and midpoint is not None:
            window_index = int(midpoint // float(window))
            locus_key = f"{chrom}:{window_index}"
        elif chrom:
            locus_key = chrom
        else:
            locus_key = ""
        loci[gene_id] = {
            "chrom": chrom,
            "midpoint": midpoint,
            "locus_key": locus_key,
            "chromosome_arm": "",
        }
    return loci


def _apply_locus_density_penalty(
    *,
    gene_scores: dict[str, float],
    by_gene: dict[str, list[SpliceEventRecord]],
    gene_symbols: dict[str, str],
    cfg: SplicingWorkflowConfig,
) -> tuple[dict[str, float], dict[str, float], dict[str, object]]:
    gene_locus_map = getattr(cfg, "_gene_locus_map", None)
    loci = _gene_primary_locus(by_gene, cfg.locus_density_window_bp, gene_symbols, gene_locus_map)
    ranked = sorted(gene_scores, key=lambda gid: (-abs(float(gene_scores[gid])), str(gid)))
    top_ids = ranked[: max(1, int(cfg.locus_density_top_n))]
    chrom_counts: dict[str, int] = {}
    arm_counts: dict[str, int] = {}
    locus_counts: dict[str, int] = {}
    for gene_id in top_ids:
        chrom = str(loci.get(gene_id, {}).get("chrom", ""))
        arm = str(loci.get(gene_id, {}).get("chromosome_arm", ""))
        locus_key = str(loci.get(gene_id, {}).get("locus_key", ""))
        if chrom:
            chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        if arm:
            arm_counts[arm] = arm_counts.get(arm, 0) + 1
        if locus_key:
            locus_counts[locus_key] = locus_counts.get(locus_key, 0) + 1
    dominant_chrom, dominant_chrom_count = ("", 0)
    if chrom_counts:
        dominant_chrom, dominant_chrom_count = max(chrom_counts.items(), key=lambda item: (item[1], item[0]))
    dominant_locus, dominant_locus_count = ("", 0)
    if locus_counts:
        dominant_locus, dominant_locus_count = max(locus_counts.items(), key=lambda item: (item[1], item[0]))
    dominant_arm, dominant_arm_count = ("", 0)
    if arm_counts:
        dominant_arm, dominant_arm_count = max(arm_counts.items(), key=lambda item: (item[1], item[0]))
    top_n = float(max(1, len(top_ids)))
    dominant_chrom_fraction = float(dominant_chrom_count) / top_n
    dominant_locus_fraction = float(dominant_locus_count) / top_n
    dominant_arm_fraction = float(dominant_arm_count) / top_n
    likely_artifact = bool(
        len(top_ids) >= 5
        and (
            dominant_locus_fraction >= 0.5
            or dominant_arm_fraction >= 0.6
            or dominant_chrom_fraction >= 0.7
        )
    )
    penalized = dict(gene_scores)
    penalties: dict[str, float] = {gene_id: 1.0 for gene_id in gene_scores}
    applied = False
    if cfg.locus_density_penalty_mode == "window_diversity":
        locus_seen: dict[str, int] = {}
        for gene_id in ranked:
            locus_key = str(loci.get(gene_id, {}).get("locus_key", ""))
            if not locus_key:
                continue
            count = locus_seen.get(locus_key, 0)
            locus_seen[locus_key] = count + 1
            if count == 0:
                continue
            penalties[gene_id] = 1.0 / math.sqrt(float(count + 1))
            penalized[gene_id] = float(penalized[gene_id]) * penalties[gene_id]
            applied = True
    elif cfg.locus_density_penalty_mode == "chromosome_diversity":
        chrom_seen: dict[str, int] = {}
        for gene_id in ranked:
            chrom = str(loci.get(gene_id, {}).get("chrom", ""))
            if not chrom:
                continue
            count = chrom_seen.get(chrom, 0)
            chrom_seen[chrom] = count + 1
            if count == 0:
                continue
            penalties[gene_id] = 1.0 / math.sqrt(float(count + 1))
            penalized[gene_id] = float(penalized[gene_id]) * penalties[gene_id]
            applied = True
    elif cfg.locus_density_penalty_mode == "auto" and likely_artifact:
        window_seen: dict[str, int] = {}
        arm_seen: dict[str, int] = {}
        chrom_seen: dict[str, int] = {}
        window_trigger = dominant_locus_fraction >= 0.5
        arm_trigger = dominant_arm_fraction >= 0.6 and bool(dominant_arm)
        chrom_trigger = dominant_chrom_fraction >= 0.7
        for gene_id in ranked:
            locus = loci.get(gene_id, {})
            local_penalties = [1.0]
            locus_key = str(locus.get("locus_key", ""))
            chrom = str(locus.get("chrom", ""))
            arm = str(locus.get("chromosome_arm", ""))
            if window_trigger and locus_key:
                count = window_seen.get(locus_key, 0)
                window_seen[locus_key] = count + 1
                if count > 0:
                    local_penalties.append(1.0 / math.sqrt(float(count + 1)))
            if arm_trigger and arm:
                count = arm_seen.get(arm, 0)
                arm_seen[arm] = count + 1
                if count > 0:
                    local_penalties.append(1.0 / math.sqrt(float(count + 1)))
            if chrom_trigger and chrom:
                count = chrom_seen.get(chrom, 0)
                chrom_seen[chrom] = count + 1
                if count > 0:
                    local_penalties.append(1.0 / (1.0 + 0.3 * float(count)))
            penalty = min(local_penalties)
            penalties[gene_id] = penalty
            if penalty < 1.0:
                penalized[gene_id] = float(penalized[gene_id]) * penalty
                applied = True
    n_penalized = sum(1 for value in penalties.values() if float(value) < 1.0)
    mode_used = cfg.locus_density_penalty_mode
    if cfg.locus_density_penalty_mode == "auto" and not likely_artifact:
        mode_used = "auto_inactive"
    return penalized, penalties, {
        "mode": cfg.locus_density_penalty_mode,
        "mode_used": mode_used,
        "window_bp": int(cfg.locus_density_window_bp),
        "top_n": int(cfg.locus_density_top_n),
        "dominant_chrom": dominant_chrom,
        "dominant_chrom_fraction": dominant_chrom_fraction,
        "dominant_chromosome_arm": dominant_arm,
        "dominant_chromosome_arm_fraction": dominant_arm_fraction,
        "dominant_locus": dominant_locus,
        "dominant_locus_fraction": dominant_locus_fraction,
        "top_n_fraction_same_chromosome": dominant_chrom_fraction,
        "top_n_fraction_same_chromosome_arm": dominant_arm_fraction,
        "top_n_fraction_same_window": dominant_locus_fraction,
        "likely_artifact": likely_artifact,
        "penalty_applied": applied,
        "n_genes_penalized_for_locus_density": n_penalized,
    }


def _apply_gene_burden_penalty(
    *,
    gene_scores: dict[str, float],
    gene_symbols: dict[str, str],
    by_gene: dict[str, list[SpliceEventRecord]],
    cfg: SplicingWorkflowConfig,
    gene_burden_map: dict[str, GeneBurdenEntry] | None,
) -> tuple[dict[str, float], dict[str, int], dict[str, str], dict[str, float]]:
    penalized: dict[str, float] = {}
    counts: dict[str, int] = {}
    sources: dict[str, str] = {}
    penalties: dict[str, float] = {}
    for gene_id, score in gene_scores.items():
        symbol = gene_symbols.get(gene_id, gene_id)
        ref_entry = gene_burden_map.get(symbol) if gene_burden_map else None
        current_group_count = len({_event_group_key(rec) for rec in by_gene.get(gene_id, [])})
        if cfg.gene_burden_penalty_mode == "none":
            burden_count = current_group_count
            source = "none"
            penalty = 1.0
        elif cfg.gene_burden_penalty_mode == "reference_bundle" or (
            cfg.gene_burden_penalty_mode == "auto" and ref_entry is not None
        ):
            if ref_entry is not None:
                unique_groups = float(ref_entry.n_unique_event_groups_ref or ref_entry.n_canonical_events_ref or 0.0)
                median_per_study = float(ref_entry.median_unique_groups_per_study or 0.0)
                high_conf_studies = float(ref_entry.n_studies_high_confidence_ref or ref_entry.n_studies_ref or 0.0)
                effective = max(unique_groups, median_per_study * max(1.0, math.sqrt(max(1.0, high_conf_studies))))
                burden_count = int(max(1.0, round(effective))) if effective > 0 else 0
            else:
                burden_count = current_group_count
            source = "reference_bundle" if ref_entry is not None else "current_input"
            penalty = math.sqrt(float(cfg.gene_topk_events) / float(max(int(cfg.gene_topk_events), burden_count))) if burden_count > 0 else 1.0
        else:
            burden_count = current_group_count
            source = "current_input"
            penalty = math.sqrt(float(cfg.gene_topk_events) / float(max(int(cfg.gene_topk_events), burden_count))) if burden_count > 0 else 1.0
        penalty = max(float(cfg.min_gene_burden_penalty), min(1.0, float(penalty)))
        counts[gene_id] = int(burden_count)
        sources[gene_id] = source
        penalties[gene_id] = float(penalty)
        penalized[gene_id] = float(score) * float(penalty)
    return penalized, counts, sources, penalties


def _select_gene_ids(scores: dict[str, float], cfg: SplicingWorkflowConfig) -> list[str]:
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
    fieldnames = [
        "gene_id",
        "gene_symbol",
        "score",
        "weight",
        "rank",
        "n_events_total",
        "n_events_used",
        "n_independent_event_groups_used",
        "sign_coherence",
        "gene_support_penalty",
        "top_event_key",
        "top_event_type",
        "max_abs_event_score",
        "mean_confidence_weight",
        "mean_impact_weight",
        "mean_ubiquity_weight",
        "fraction_high_confidence_events",
        "gene_event_burden_count",
        "gene_event_burden_source",
        "gene_burden_penalty",
        "locus_density_penalty",
    ]
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


def _warning_flags(
    parse_summary: dict[str, object],
    collapsed_events: list[SpliceEventRecord],
    selected_rows: list[dict[str, object]],
    by_gene: dict[str, list[SpliceEventRecord]],
    locus_density_summary: dict[str, object] | None = None,
) -> list[str]:
    warnings: list[str] = []
    n_events = len(collapsed_events)
    if n_events == 0:
        return warnings
    retained_introns = sum(1 for rec in collapsed_events if rec.event_type == "retained_intron")
    if (retained_introns / float(n_events)) >= 0.5:
        warnings.append("retained_intron_dominance")
    if float(parse_summary.get("fraction_novel", 0.0) or 0.0) >= 0.5:
        warnings.append("novel_event_dominance")
    if float(parse_summary.get("fraction_low_support", 0.0) or 0.0) >= 0.5:
        warnings.append("low_confidence_dominance")
    multi_event_heavy = sum(1 for row in selected_rows if int(row.get("n_events_total", 0) or 0) >= 5)
    if selected_rows and (multi_event_heavy / float(len(selected_rows))) >= 0.4:
        warnings.append("likely_event_multiplicity_bias")
    if locus_density_summary and bool(locus_density_summary.get("likely_artifact")):
        warnings.append("likely_locus_density_artifact")
    return warnings


def _warn_summary_flags(flags: Iterable[str]) -> None:
    messages = {
        "retained_intron_dominance": "warning: retained introns dominate the retained event set; interpret selected genes cautiously.",
        "novel_event_dominance": "warning: novel or weakly annotated events dominate the retained event set.",
        "low_confidence_dominance": "warning: low-confidence events dominate the retained event set; consider stronger confidence thresholds.",
        "likely_event_multiplicity_bias": "warning: selected genes appear enriched for high event multiplicity rather than isolated strong events.",
        "likely_locus_density_artifact": "warning: top genes are unusually concentrated within one genomic locus or chromosome; this may reflect a local splicing artifact.",
        "bundle_same_dataset_match": "warning: bundle contains matching source-dataset priors; same-dataset contributions were excluded or neutralized.",
    }
    for flag in flags:
        if flag in messages:
            print(messages[flag], file=sys.stderr)


def run_splice_event_diff_workflow(
    *,
    cfg: SplicingWorkflowConfig,
    fieldnames: list[str],
    rows: list[SpliceRow],
    alias_map: dict[str, AliasEntry] | None,
    ubiquity_map: dict[str, UbiquityEntry] | None,
    event_ubiquity_by_dataset_map: dict[str, list[UbiquityByDatasetEntry]] | None,
    impact_map: dict[str, ImpactEntry] | None,
    gene_burden_map: dict[str, GeneBurdenEntry] | None,
    gene_burden_by_dataset_map: dict[str, list[GeneBurdenByDatasetEntry]] | None,
    gene_locus_map: dict[str, GeneLocusEntry] | None,
    input_files: list[dict[str, str]],
    resources_info: dict[str, object] | None,
) -> dict[str, object]:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    columns = resolve_splicing_columns(fieldnames, cfg)
    resolved_tool_family = str(columns.get("tool_family", cfg.tool_family))
    resolved_score_mode = _resolve_score_mode(columns, rows, cfg.score_mode)
    setattr(cfg, "_gene_locus_map", gene_locus_map)

    effective_ubiquity_map, effective_event_dataset_counts, same_dataset_event_summary = _build_effective_ubiquity_map(
        pooled_map=ubiquity_map,
        by_dataset_map=event_ubiquity_by_dataset_map,
        source_dataset=cfg.source_dataset,
        policy=cfg.bundle_same_dataset_policy,
    )
    effective_gene_burden_map, same_dataset_gene_summary = _build_effective_gene_burden_map(
        pooled_map=gene_burden_map,
        by_dataset_map=gene_burden_by_dataset_map,
        source_dataset=cfg.source_dataset,
        policy=cfg.bundle_same_dataset_policy,
    )

    max_read_support = 0.0
    read_col = columns.get("read_support")
    if read_col:
        for row in rows:
            value = _parse_float(row.values.get(read_col))
            if value is not None and value > max_read_support:
                max_read_support = float(value)

    parse_summary: dict[str, object] = {
        "n_input_rows": len(rows),
        "n_rows_scored": 0,
        "n_rows_dropped_missing_score": 0,
        "n_rows_unresolved_event": 0,
        "n_rows_ambiguous_gene": 0,
        "n_rows_dropped_ambiguous_gene": 0,
        "n_events_matched_to_alias_prior": 0,
        "n_events_matched_to_ubiquity_prior": 0,
        "n_events_matched_to_impact_prior": 0,
        "n_events_with_low_confidence_prior_match": 0,
        "n_low_confidence_ubiquity_priors_neutralized": 0,
        "n_low_confidence_impact_priors_neutralized": 0,
        "n_medium_confidence_namespace_mismatches_neutralized": 0,
        "n_low_support_rows": 0,
        "n_novel_rows": 0,
        "n_delta_psi_soft_floor_downweighted": 0,
        "n_delta_psi_soft_floor_full_weight": 0,
        "n_delta_psi_soft_floor_inactive": 0,
        "n_same_dataset_ubiquity_priors_excluded": int(same_dataset_event_summary.get("n_events_excluded", 0)),
        "n_same_dataset_ubiquity_priors_neutralized": int(same_dataset_event_summary.get("n_events_neutralized", 0)),
        "n_same_dataset_gene_burden_excluded": int(same_dataset_gene_summary.get("n_genes_excluded", 0)),
        "n_same_dataset_gene_burden_neutralized": int(same_dataset_gene_summary.get("n_genes_neutralized", 0)),
        "event_type_counts": {},
        "canonicalization_confidence_counts": {"high": 0, "medium": 0, "low": 0},
    }

    event_records: list[SpliceEventRecord] = []
    for row in rows:
        base_score, _used_mode = _compute_base_score(row, columns, resolved_score_mode, cfg.neglog10p_eps, cfg.neglog10p_cap, max_read_support)
        if base_score is None:
            parse_summary["n_rows_dropped_missing_score"] = int(parse_summary["n_rows_dropped_missing_score"]) + 1
            continue

        default_canonical_event_key = _default_canonical_event_key(row, columns)
        default_canonicalization_status, default_canonicalization_confidence, default_event_key_namespace = _default_canonicalization_status(row, columns)
        alias_entry = None
        for candidate in _canonical_event_candidates(row, columns):
            if alias_map and candidate in alias_map:
                candidate_entry = alias_map[candidate]
                if (
                    candidate_entry.canonicalization_confidence == "medium"
                    and not _namespace_matches(default_event_key_namespace, candidate_entry.event_key_namespace)
                ):
                    parse_summary["n_medium_confidence_namespace_mismatches_neutralized"] = int(
                        parse_summary["n_medium_confidence_namespace_mismatches_neutralized"]
                    ) + 1
                    continue
                alias_entry = candidate_entry
                break
        if alias_entry is not None:
            parse_summary["n_events_matched_to_alias_prior"] = int(parse_summary["n_events_matched_to_alias_prior"]) + 1
        canonical_event_key = alias_entry.canonical_event_key if alias_entry is not None else default_canonical_event_key
        if not canonical_event_key:
            parse_summary["n_rows_unresolved_event"] = int(parse_summary["n_rows_unresolved_event"]) + 1
            continue
        canonicalization_status, canonicalization_confidence, event_key_namespace = (
            (alias_entry.canonicalization_status, alias_entry.canonicalization_confidence, alias_entry.event_key_namespace)
            if alias_entry is not None
            else (default_canonicalization_status, default_canonicalization_confidence, default_event_key_namespace)
        )
        parse_summary["canonicalization_confidence_counts"][canonicalization_confidence] = int(
            parse_summary["canonicalization_confidence_counts"].get(canonicalization_confidence, 0)
        ) + 1

        raw_gene_id = _string_value(row, columns, "gene_id")
        raw_gene_symbol = _string_value(row, columns, "gene_symbol")
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

        confidence_weight, confidence_parts = _confidence_weight(
            row,
            columns=columns,
            mode=cfg.confidence_weight_mode,
            min_probability=cfg.min_probability,
            min_read_support=cfg.min_read_support,
            max_read_support=max_read_support,
            neglog10p_eps=cfg.neglog10p_eps,
            neglog10p_cap=cfg.neglog10p_cap,
        )
        if confidence_weight <= 0.0:
            parse_summary["n_low_support_rows"] = int(parse_summary["n_low_support_rows"]) + 1
        novel_flag = _parse_boolish(_string_value(row, columns, "novel_flag"))
        annotation_status = _string_value(row, columns, "annotation_status")
        if novel_flag or "novel" in annotation_status.lower():
            parse_summary["n_novel_rows"] = int(parse_summary["n_novel_rows"]) + 1
        event_type = alias_entry.event_type if alias_entry is not None and alias_entry.event_type else _normalize_event_type(_string_value(row, columns, "event_type"))
        if not event_type:
            event_type = "unknown"
        parse_summary["event_type_counts"].setdefault(event_type, 0)
        parse_summary["event_type_counts"][event_type] += 1

        prior_confidence_weight = _prior_confidence_weight(canonicalization_confidence)
        low_conf_prior_matched = False
        ubiquity_entry = effective_ubiquity_map.get(canonical_event_key) if effective_ubiquity_map else None
        if ubiquity_entry is None:
            ubiquity_weight = 1.0
        elif ubiquity_entry.canonicalization_confidence == "high":
            ubiquity_weight = float(ubiquity_entry.idf_ref)
        elif ubiquity_entry.canonicalization_confidence == "medium":
            if _namespace_matches(event_key_namespace, ubiquity_entry.event_key_namespace):
                ubiquity_weight = _shrink_prior_toward_neutral(float(ubiquity_entry.idf_ref), "medium")
            else:
                ubiquity_weight = 1.0
                parse_summary["n_medium_confidence_namespace_mismatches_neutralized"] = int(
                    parse_summary["n_medium_confidence_namespace_mismatches_neutralized"]
                ) + 1
        else:
            ubiquity_weight = 1.0
            parse_summary["n_low_confidence_ubiquity_priors_neutralized"] = int(
                parse_summary["n_low_confidence_ubiquity_priors_neutralized"]
            ) + 1
        if ubiquity_entry is not None:
            parse_summary["n_events_matched_to_ubiquity_prior"] = int(parse_summary["n_events_matched_to_ubiquity_prior"]) + 1
            if ubiquity_entry.canonicalization_confidence != "high":
                low_conf_prior_matched = True
        impact_entry = impact_map.get(canonical_event_key) if impact_map else None
        effective_impact_entry = impact_entry
        if impact_entry is not None:
            effective_impact_entry = ImpactEntry(
                canonical_event_key=impact_entry.canonical_event_key,
                canonicalization_status=impact_entry.canonicalization_status,
                canonicalization_confidence=impact_entry.canonicalization_confidence,
                event_key_namespace=impact_entry.event_key_namespace,
                gene_id=impact_entry.gene_id,
                gene_symbol=impact_entry.gene_symbol,
                event_type=impact_entry.event_type,
                impact_weight_raw=impact_entry.impact_weight_raw,
                impact_evidence=impact_entry.impact_evidence,
                annotation_status=impact_entry.annotation_status,
                n_datasets_ref=effective_event_dataset_counts.get(canonical_event_key, impact_entry.n_datasets_ref or 0.0),
            )
            if (
                effective_impact_entry.canonicalization_confidence == "medium"
                and not _namespace_matches(event_key_namespace, effective_impact_entry.event_key_namespace)
            ):
                effective_impact_entry = None
                parse_summary["n_medium_confidence_namespace_mismatches_neutralized"] = int(
                    parse_summary["n_medium_confidence_namespace_mismatches_neutralized"]
                ) + 1
        impact_weight = _impact_weight(effective_impact_entry, cfg, event_type, annotation_status, novel_flag, canonicalization_confidence)
        if impact_entry is not None:
            parse_summary["n_events_matched_to_impact_prior"] = int(parse_summary["n_events_matched_to_impact_prior"]) + 1
            if impact_entry.canonicalization_confidence == "low":
                low_conf_prior_matched = True
                if abs(float(impact_weight) - 1.0) <= 1e-9:
                    parse_summary["n_low_confidence_impact_priors_neutralized"] = int(
                        parse_summary["n_low_confidence_impact_priors_neutralized"]
                    ) + 1
        if low_conf_prior_matched:
            parse_summary["n_events_with_low_confidence_prior_match"] = int(parse_summary["n_events_with_low_confidence_prior_match"]) + 1

        transformed = _apply_score_transform(float(base_score), cfg.score_transform)
        delta_psi_effect_weight, effect_active = _delta_psi_effect_weight(row, columns, cfg)
        if effect_active and delta_psi_effect_weight < 1.0:
            parse_summary["n_delta_psi_soft_floor_downweighted"] = int(parse_summary["n_delta_psi_soft_floor_downweighted"]) + 1
        elif effect_active:
            parse_summary["n_delta_psi_soft_floor_full_weight"] = int(parse_summary["n_delta_psi_soft_floor_full_weight"]) + 1
        else:
            parse_summary["n_delta_psi_soft_floor_inactive"] = int(parse_summary["n_delta_psi_soft_floor_inactive"]) + 1
        final_score = float(transformed) * float(delta_psi_effect_weight) * float(confidence_weight) * float(impact_weight) * float(ubiquity_weight)
        chrom = alias_entry.chrom if alias_entry is not None and alias_entry.chrom else _string_value(row, columns, "chrom")
        start = alias_entry.start if alias_entry is not None and alias_entry.start else _string_value(row, columns, "start")
        end = alias_entry.end if alias_entry is not None and alias_entry.end else _string_value(row, columns, "end")
        strand = alias_entry.strand if alias_entry is not None and alias_entry.strand else _string_value(row, columns, "strand")

        for gene_id, gene_symbol, gene_weight in gene_assignments:
            event_records.append(
                SpliceEventRecord(
                    canonical_event_key=canonical_event_key,
                    event_group=_string_value(row, columns, "event_group"),
                    canonicalization_status=canonicalization_status,
                    canonicalization_confidence=canonicalization_confidence,
                    event_key_namespace=event_key_namespace,
                    event_type=event_type,
                    gene_id=str(gene_id),
                    gene_symbol=str(gene_symbol),
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    base_score=float(base_score),
                    transformed_score=float(transformed),
                    final_score=float(final_score) * float(gene_weight),
                    confidence_weight=float(confidence_weight),
                    impact_weight=float(impact_weight),
                    ubiquity_weight=float(ubiquity_weight),
                    delta_psi_effect_weight=float(delta_psi_effect_weight),
                    prior_confidence_weight=float(prior_confidence_weight),
                    quality_value=float(confidence_parts.get("pvalue_weight", 1.0)) * max(1e-6, abs(float(transformed))),
                    source_line_no=int(row.line_no),
                )
            )
        parse_summary["n_rows_scored"] = int(parse_summary["n_rows_scored"]) + 1

    collapsed_events = _collapse_duplicate_events(event_records, cfg.event_dup_policy)
    grouped_events = _collapse_event_groups(collapsed_events)
    gene_scores, gene_symbols, by_gene = _aggregate_events_to_genes(grouped_events, cfg.gene_aggregation, cfg.gene_topk_events)
    gene_scores, gene_support_penalties, gene_support_group_counts, gene_support_event_counts, gene_sign_coherence = _apply_gene_support_penalty(
        gene_scores=gene_scores,
        by_gene=by_gene,
        cfg=cfg,
    )
    gene_scores, gene_burden_counts, gene_burden_sources, gene_burden_penalties = _apply_gene_burden_penalty(
        gene_scores=gene_scores,
        gene_symbols=gene_symbols,
        by_gene=by_gene,
        cfg=cfg,
        gene_burden_map=effective_gene_burden_map,
    )
    gene_scores, locus_density_penalties, locus_density_summary = _apply_locus_density_penalty(
        gene_scores=gene_scores,
        by_gene=by_gene,
        gene_symbols=gene_symbols,
        cfg=cfg,
    )
    magnitude = {gene_id: abs(float(score)) for gene_id, score in gene_scores.items()}
    selected_gene_ids = _select_gene_ids(magnitude, cfg)
    weights = _selected_weights(magnitude, selected_gene_ids, cfg.normalize)

    parse_summary["fraction_retained_intron"] = (
        sum(1 for rec in grouped_events if rec.event_type == "retained_intron") / float(len(grouped_events)) if grouped_events else 0.0
    )
    parse_summary["fraction_novel"] = (
        float(parse_summary["n_novel_rows"]) / float(parse_summary["n_rows_scored"]) if parse_summary["n_rows_scored"] else 0.0
    )
    parse_summary["fraction_low_support"] = (
        float(parse_summary["n_low_support_rows"]) / float(parse_summary["n_rows_scored"]) if parse_summary["n_rows_scored"] else 0.0
    )

    full_ranked_ids = sorted(gene_scores, key=lambda g: (-abs(float(gene_scores[g])), str(g)))
    full_rows: list[dict[str, object]] = []
    for idx, gene_id in enumerate(full_ranked_ids, start=1):
        support = sorted(by_gene.get(gene_id, []), key=lambda r: (-abs(float(r.final_score)), str(r.canonical_event_key)))
        used = support[: max(1, int(cfg.gene_topk_events))] if cfg.gene_aggregation == "signed_topk_mean" else support[:1]
        full_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene_symbols.get(gene_id, gene_id),
                "score": float(gene_scores[gene_id]),
                "weight": float(weights.get(gene_id, 0.0)),
                "rank": idx,
                "n_events_total": len(support),
                "n_events_used": int(gene_support_event_counts.get(gene_id, len(used))),
                "n_independent_event_groups_used": int(gene_support_group_counts.get(gene_id, len(used))),
                "sign_coherence": float(gene_sign_coherence.get(gene_id, 0.0)),
                "gene_support_penalty": float(gene_support_penalties.get(gene_id, 1.0)),
                "top_event_key": support[0].canonical_event_key if support else "",
                "top_event_type": support[0].event_type if support else "",
                "max_abs_event_score": abs(float(support[0].final_score)) if support else 0.0,
                "mean_confidence_weight": (sum(float(r.confidence_weight) for r in support) / float(len(support))) if support else 0.0,
                "mean_impact_weight": (sum(float(r.impact_weight) for r in support) / float(len(support))) if support else 0.0,
                "mean_ubiquity_weight": (sum(float(r.ubiquity_weight) for r in support) / float(len(support))) if support else 0.0,
                "fraction_high_confidence_events": (
                    sum(1 for r in support if r.canonicalization_confidence == "high") / float(len(support))
                ) if support else 0.0,
                "gene_event_burden_count": int(gene_burden_counts.get(gene_id, 0)),
                "gene_event_burden_source": gene_burden_sources.get(gene_id, "current_input"),
                "gene_burden_penalty": float(gene_burden_penalties.get(gene_id, 1.0)),
                "locus_density_penalty": float(locus_density_penalties.get(gene_id, 1.0)),
            }
        )

    selected_rows = [row for row in full_rows if row["gene_id"] in set(selected_gene_ids)]
    selected_rows.sort(key=lambda row: int(row["rank"]))
    rank = 1
    for row in selected_rows:
        row["rank"] = rank
        row["weight"] = float(weights.get(str(row["gene_id"]), 0.0))
        rank += 1

    parse_summary["fraction_single_event_genes"] = (
        sum(1 for row in selected_rows if int(row.get("n_events_total", 0) or 0) <= 1) / float(len(selected_rows)) if selected_rows else 0.0
    )
    selected_support_events = [
        rec
        for gene_id in selected_gene_ids
        for rec in by_gene.get(gene_id, [])
    ]
    selected_high_conf = sum(1 for rec in selected_support_events if rec.canonicalization_confidence == "high")
    selected_low_conf = sum(1 for rec in selected_support_events if rec.canonicalization_confidence != "high")
    selected_low_conf_prior = sum(
        1
        for rec in selected_support_events
        if rec.canonicalization_confidence != "high" and (abs(rec.impact_weight - 1.0) > 1e-9 or abs(rec.ubiquity_weight - 1.0) > 1e-9)
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
        base_name = sanitize_name_component(
            f"{cfg.converter_name}__signature={cfg.signature_name}__tool_family={resolved_tool_family}__score_mode={resolved_score_mode}"
        )
        gmt_sets, gmt_plans = build_gmt_sets_from_rows(
            full_rows,
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
            context={"program_method": "splice_event", "contrast_method": "input_table", "link_method": "event_to_gene"},
        )
        if gmt_sets:
            gmt_path = resolve_gmt_out_path(out_dir, cfg.gmt_out)
            write_gmt(gmt_sets, gmt_path, gmt_format=cfg.gmt_format)
            output_files.append({"path": str(gmt_path), "role": "gmt"})
        _warn_gmt_diagnostics(gmt_diagnostics)

    skipped_gmt_outputs = _collect_skipped_gmt_outputs(gmt_diagnostics)
    warning_flags = _warning_flags(parse_summary, grouped_events, selected_rows, by_gene, locus_density_summary)
    if _clean(same_dataset_event_summary.get("warning")):
        warning_flags.append("bundle_same_dataset_match")
        print(f"warning: {_clean(same_dataset_event_summary.get('warning'))}", file=sys.stderr)
    _warn_summary_flags(warning_flags)

    run_summary_payload = {
        "converter": cfg.converter_name,
        "dataset_label": cfg.dataset_label,
        "signature_name": cfg.signature_name,
        "tool_family": resolved_tool_family,
        "score_mode": resolved_score_mode,
        "score_transform": cfg.score_transform,
        "confidence_weight_mode": cfg.confidence_weight_mode,
        "impact_mode": cfg.impact_mode,
        "event_dup_policy": cfg.event_dup_policy,
        "gene_aggregation": cfg.gene_aggregation,
        "n_input_rows": parse_summary["n_input_rows"],
        "n_rows_scored": parse_summary["n_rows_scored"],
        "n_canonical_events": len({rec.canonical_event_key for rec in collapsed_events}),
        "n_independent_event_groups": len(grouped_events),
        "n_output_genes": len(selected_rows),
        "n_events_matched_to_alias_prior": parse_summary["n_events_matched_to_alias_prior"],
        "n_events_matched_to_ubiquity_prior": parse_summary["n_events_matched_to_ubiquity_prior"],
        "n_events_matched_to_impact_prior": parse_summary["n_events_matched_to_impact_prior"],
        "n_events_with_low_confidence_prior_match": parse_summary["n_events_with_low_confidence_prior_match"],
        "n_low_confidence_ubiquity_priors_neutralized": parse_summary["n_low_confidence_ubiquity_priors_neutralized"],
        "n_low_confidence_impact_priors_neutralized": parse_summary["n_low_confidence_impact_priors_neutralized"],
        "n_medium_confidence_namespace_mismatches_neutralized": parse_summary["n_medium_confidence_namespace_mismatches_neutralized"],
        "n_same_dataset_ubiquity_priors_excluded": parse_summary["n_same_dataset_ubiquity_priors_excluded"],
        "n_same_dataset_ubiquity_priors_neutralized": parse_summary["n_same_dataset_ubiquity_priors_neutralized"],
        "n_same_dataset_gene_burden_excluded": parse_summary["n_same_dataset_gene_burden_excluded"],
        "n_same_dataset_gene_burden_neutralized": parse_summary["n_same_dataset_gene_burden_neutralized"],
        "selected_events_high_confidence": selected_high_conf,
        "selected_events_low_confidence": selected_low_conf,
        "selected_events_fraction_high_confidence": (float(selected_high_conf) / float(len(selected_support_events))) if selected_support_events else 0.0,
        "selected_events_fraction_low_confidence": (float(selected_low_conf) / float(len(selected_support_events))) if selected_support_events else 0.0,
        "selected_events_low_confidence_prior_applied": selected_low_conf_prior,
        "fraction_retained_intron": parse_summary["fraction_retained_intron"],
        "fraction_novel": parse_summary["fraction_novel"],
        "fraction_low_support": parse_summary["fraction_low_support"],
        "fraction_single_event_genes": parse_summary["fraction_single_event_genes"],
        "delta_psi_soft_floor_mode": cfg.delta_psi_soft_floor_mode,
        "delta_psi_soft_floor": cfg.delta_psi_soft_floor,
        "source_dataset": cfg.source_dataset,
        "bundle_same_dataset_policy": cfg.bundle_same_dataset_policy,
        "same_dataset_event_prior_summary": same_dataset_event_summary,
        "same_dataset_gene_burden_summary": same_dataset_gene_summary,
        "n_delta_psi_soft_floor_downweighted": parse_summary["n_delta_psi_soft_floor_downweighted"],
        "n_delta_psi_soft_floor_full_weight": parse_summary["n_delta_psi_soft_floor_full_weight"],
        "n_delta_psi_soft_floor_inactive": parse_summary["n_delta_psi_soft_floor_inactive"],
        "canonicalization_confidence_counts": parse_summary["canonicalization_confidence_counts"],
        "gene_burden_penalty_mode": cfg.gene_burden_penalty_mode,
        "min_gene_burden_penalty": cfg.min_gene_burden_penalty,
        "gene_support_penalty_mode": cfg.gene_support_penalty_mode,
        "locus_density_penalty_mode": cfg.locus_density_penalty_mode,
        "locus_density_window_bp": cfg.locus_density_window_bp,
        "locus_density_top_n": cfg.locus_density_top_n,
        "locus_density": locus_density_summary,
        "event_type_counts": parse_summary["event_type_counts"],
        "warnings": warning_flags,
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
        "tool_family": resolved_tool_family,
        "score_mode": resolved_score_mode,
        "score_transform": cfg.score_transform,
        "confidence_weight_mode": cfg.confidence_weight_mode,
        "min_probability": cfg.min_probability,
        "min_read_support": cfg.min_read_support,
        "impact_mode": cfg.impact_mode,
        "impact_min": cfg.impact_min,
        "impact_max": cfg.impact_max,
        "delta_psi_soft_floor": cfg.delta_psi_soft_floor,
        "delta_psi_soft_floor_mode": cfg.delta_psi_soft_floor_mode,
        "event_dup_policy": cfg.event_dup_policy,
        "gene_aggregation": cfg.gene_aggregation,
        "gene_topk_events": cfg.gene_topk_events,
        "gene_burden_penalty_mode": cfg.gene_burden_penalty_mode,
        "min_gene_burden_penalty": cfg.min_gene_burden_penalty,
        "gene_support_penalty_mode": cfg.gene_support_penalty_mode,
        "source_dataset": cfg.source_dataset,
        "bundle_same_dataset_policy": cfg.bundle_same_dataset_policy,
        "locus_density_penalty_mode": cfg.locus_density_penalty_mode,
        "locus_density_window_bp": cfg.locus_density_window_bp,
        "locus_density_top_n": cfg.locus_density_top_n,
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
        data_type="splicing_event",
        assay="event_level",
        organism=cfg.organism,
        genome_build=cfg.genome_build,
        files=input_files,
        gene_annotation={
            "mode": "none",
            "source": "splicing_table_and_optional_bundle",
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
            "n_input_features": parse_summary["n_input_rows"],
            "n_input_rows": parse_summary["n_input_rows"],
            "n_rows_scored": parse_summary["n_rows_scored"],
            "n_canonical_events": len({rec.canonical_event_key for rec in collapsed_events}),
            "n_independent_event_groups": len(grouped_events),
            "n_events_matched_to_alias_prior": parse_summary["n_events_matched_to_alias_prior"],
            "n_events_matched_to_ubiquity_prior": parse_summary["n_events_matched_to_ubiquity_prior"],
            "n_events_matched_to_impact_prior": parse_summary["n_events_matched_to_impact_prior"],
            "n_events_with_low_confidence_prior_match": parse_summary["n_events_with_low_confidence_prior_match"],
            "n_low_confidence_ubiquity_priors_neutralized": parse_summary["n_low_confidence_ubiquity_priors_neutralized"],
            "n_low_confidence_impact_priors_neutralized": parse_summary["n_low_confidence_impact_priors_neutralized"],
            "n_medium_confidence_namespace_mismatches_neutralized": parse_summary["n_medium_confidence_namespace_mismatches_neutralized"],
            "n_same_dataset_ubiquity_priors_excluded": parse_summary["n_same_dataset_ubiquity_priors_excluded"],
            "n_same_dataset_ubiquity_priors_neutralized": parse_summary["n_same_dataset_ubiquity_priors_neutralized"],
            "n_same_dataset_gene_burden_excluded": parse_summary["n_same_dataset_gene_burden_excluded"],
            "n_same_dataset_gene_burden_neutralized": parse_summary["n_same_dataset_gene_burden_neutralized"],
            "canonicalization_confidence_counts": parse_summary["canonicalization_confidence_counts"],
            "fraction_retained_intron": parse_summary["fraction_retained_intron"],
            "fraction_novel": parse_summary["fraction_novel"],
            "fraction_low_support": parse_summary["fraction_low_support"],
            "fraction_single_event_genes": parse_summary["fraction_single_event_genes"],
            "same_dataset_event_prior_summary": same_dataset_event_summary,
            "same_dataset_gene_burden_summary": same_dataset_gene_summary,
            "n_delta_psi_soft_floor_downweighted": parse_summary["n_delta_psi_soft_floor_downweighted"],
            "n_delta_psi_soft_floor_full_weight": parse_summary["n_delta_psi_soft_floor_full_weight"],
            "n_delta_psi_soft_floor_inactive": parse_summary["n_delta_psi_soft_floor_inactive"],
            "event_type_counts": parse_summary["event_type_counts"],
            "locus_density": locus_density_summary,
            "warnings": warning_flags,
            "n_genes_pre_selection": len(gene_scores),
            "n_genes": len(selected_rows),
            "n_features_assigned": len(grouped_events),
            "fraction_features_assigned": (float(len(grouped_events)) / float(parse_summary["n_input_rows"]) if parse_summary["n_input_rows"] else 0.0),
        },
        program_extraction={
            "selection_method": cfg.select,
            "selection_params": {"k": cfg.top_k, "quantile": cfg.quantile, "min_score": cfg.min_score},
            "normalize": cfg.normalize,
            "n_selected_genes": len(selected_rows),
            "score_definition": "event-level splicing evidence aggregated to genes after canonical-event and event-group collapse",
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
        "n_input_rows": parse_summary["n_input_rows"],
        "n_genes": len(selected_rows),
        "n_canonical_events": len({rec.canonical_event_key for rec in collapsed_events}),
        "out_dir": str(out_dir),
        "resolved_score_mode": resolved_score_mode,
        "resolved_tool_family": resolved_tool_family,
    }
