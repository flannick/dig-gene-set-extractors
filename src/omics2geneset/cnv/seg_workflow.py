from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from pathlib import Path
import re
import statistics
import sys

from omics2geneset.cnv.seg_io import CNVSegment
from omics2geneset.cnv.seg_scoring import (
    focal_penalty,
    gene_count_penalty,
    purity_corrected_amplitude,
)
from omics2geneset.core.gmt import (
    build_gmt_sets_from_rows,
    parse_int_list_csv,
    parse_mass_list_csv,
    parse_str_list_csv,
    resolve_gmt_out_path,
    write_gmt,
)
from omics2geneset.core.metadata import make_metadata, write_metadata
from omics2geneset.core.qc import write_run_summary_files
from omics2geneset.core.selection import (
    global_l1_weights,
    ranked_gene_ids,
    select_quantile,
    select_threshold,
    select_top_k,
    within_set_l1_weights,
)
from omics2geneset.io.gtf import read_genes_from_gtf


_SAFE_COMPONENT_RE = re.compile(r"[^A-Za-z0-9._-]+")
CHROM_MISMATCH_WARN_FRACTION = 0.2
BREADTH_MEDIAN_BP_WARN = 50_000_000
BREADTH_MEAN_FOCAL_WARN = 0.10


@dataclass
class CNVWorkflowConfig:
    converter_name: str
    out_dir: Path
    organism: str
    genome_build: str
    gtf: str
    gtf_source: str | None
    gtf_gene_id_field: str
    dataset_label: str
    program_preset: str
    program_methods: str | None
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
    coord_system: str
    chrom_prefix_mode: str
    min_abs_amplitude: float
    focal_length_scale_bp: float
    focal_length_alpha: float
    gene_count_penalty_mode: str
    aggregation: str
    use_purity_correction: bool
    purity_floor: float
    max_abs_amplitude: float
    emit_cohort_sets: bool
    cohort_score_threshold: float
    cohort_min_fraction: float
    cohort_min_samples: int


def _safe_component(value: str, fallback: str) -> str:
    out = _SAFE_COMPONENT_RE.sub("_", str(value)).strip("_")
    return out or fallback


def _parse_program_methods(raw: str | None, program_preset: str) -> list[str]:
    allowed = {"amp", "del", "abs"}
    if raw and str(raw).strip():
        out: list[str] = []
        for token in str(raw).split(","):
            item = token.strip().lower()
            if not item:
                continue
            if item not in allowed:
                raise ValueError(
                    f"Unsupported CNV program method '{item}'. Expected one of: {', '.join(sorted(allowed))}"
                )
            if item not in out:
                out.append(item)
        return out
    preset = str(program_preset).strip().lower()
    if preset in {"default", "connectable", "broad"}:
        return ["amp", "del"]
    if preset in {"qc", "all"}:
        return ["amp", "del", "abs"]
    if preset == "none":
        return []
    raise ValueError(f"Unsupported program_preset: {program_preset}")


def _select_gene_ids(scores: dict[str, float], cfg: CNVWorkflowConfig) -> list[str]:
    if cfg.select == "none":
        return ranked_gene_ids(scores)
    if cfg.select == "top_k":
        return select_top_k(scores, int(cfg.top_k))
    if cfg.select == "quantile":
        return select_quantile(scores, float(cfg.quantile))
    if cfg.select == "threshold":
        return select_threshold(scores, float(cfg.min_score))
    raise ValueError(f"Unsupported selection method: {cfg.select}")


def _selected_weights(basis_scores: dict[str, float], selected_gene_ids: list[str], normalize: str) -> dict[str, float]:
    if normalize == "none":
        return {g: float(basis_scores.get(g, 0.0)) for g in selected_gene_ids}
    if normalize == "l1":
        global_weights = global_l1_weights(basis_scores)
        return {g: float(global_weights.get(g, 0.0)) for g in selected_gene_ids}
    if normalize == "within_set_l1":
        return within_set_l1_weights(basis_scores, selected_gene_ids)
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


def _warn_gmt_diagnostics(gmt_diagnostics: list[dict[str, object]], *, sample_id: str, program: str) -> int:
    n_warnings = 0
    for diag in gmt_diagnostics:
        code = str(diag.get("code", ""))
        if code not in {"small_gene_set_skipped", "no_positive_genes", "small_gene_set_emitted", "marginal_gene_count"}:
            continue
        n_warnings += 1
        n_genes = int(diag.get("n_genes", 0) or 0)
        reason = str(diag.get("reason", "")).strip()
        suggestion = str(diag.get("suggestion", "")).strip()
        if code == "small_gene_set_emitted":
            print(
                f"warning: sample={sample_id} program={program} emitted small GMT set (n_genes={n_genes}). {reason}",
                file=sys.stderr,
            )
        elif code == "marginal_gene_count":
            print(
                f"warning: sample={sample_id} program={program} has marginal gene count (n_genes={n_genes}). {reason}",
                file=sys.stderr,
            )
        elif code == "no_positive_genes":
            print(
                f"warning: sample={sample_id} program={program} has no positive genes after filtering; GMT not emitted.",
                file=sys.stderr,
            )
        else:
            print(
                f"warning: sample={sample_id} program={program} skipped GMT set (n_genes={n_genes}). {reason}",
                file=sys.stderr,
            )
        if suggestion:
            print(f"warning:   {suggestion}", file=sys.stderr)
    return n_warnings


def _map_chrom_to_gene_space(chrom: str, gene_chroms: set[str], chrom_prefix_mode: str) -> str | None:
    c = str(chrom).strip()
    if not c:
        return None
    if c in gene_chroms:
        return c
    if chrom_prefix_mode == "auto":
        if c.lower().startswith("chr"):
            alt = c[3:]
            if alt in gene_chroms:
                return alt
        else:
            alt = f"chr{c}"
            if alt in gene_chroms:
                return alt
    return c if c in gene_chroms else None


def _segment_length(seg: CNVSegment) -> int:
    return max(1, int(seg.end) - int(seg.start))


def _filter_and_correct_segments(
    *,
    segments: list[CNVSegment],
    gene_chroms: set[str],
    cfg: CNVWorkflowConfig,
    purity_by_sample: dict[str, float] | None,
) -> tuple[list[CNVSegment], dict[str, object]]:
    out: list[CNVSegment] = []
    n_skipped_chrom = 0
    n_below_amp = 0
    n_purity_missing = 0
    for seg in segments:
        mapped_chrom = _map_chrom_to_gene_space(seg.chrom, gene_chroms, cfg.chrom_prefix_mode)
        if mapped_chrom is None:
            n_skipped_chrom += 1
            continue
        purity = None if purity_by_sample is None else purity_by_sample.get(seg.sample_id)
        if cfg.use_purity_correction and purity_by_sample is not None and purity is None:
            n_purity_missing += 1
        corrected = purity_corrected_amplitude(
            seg.amplitude,
            purity,
            use_purity_correction=cfg.use_purity_correction,
            purity_floor=cfg.purity_floor,
            max_abs_amplitude=cfg.max_abs_amplitude,
        )
        if abs(corrected) < float(cfg.min_abs_amplitude):
            n_below_amp += 1
            continue
        out.append(
            replace(
                seg,
                chrom=mapped_chrom,
                amplitude=float(corrected),
            )
        )
    summary: dict[str, object] = {
        "n_segments_input": len(segments),
        "n_segments_used": len(out),
        "n_segments_skipped_chrom_mismatch": n_skipped_chrom,
        "n_segments_below_min_abs_amplitude": n_below_amp,
        "n_segments_missing_purity": n_purity_missing,
    }
    return out, summary


def _score_sample_segments(
    *,
    sample_segments: list[CNVSegment],
    genes_by_chrom,
    cfg: CNVWorkflowConfig,
) -> tuple[dict[str, float], dict[str, object], dict[str, list[tuple[int, float]]]]:
    numerator: dict[str, float] = {}
    link_denom: dict[str, float] = {}
    contrib_count: dict[str, int] = {}
    max_abs_contrib: dict[str, float] = {}
    gene_segment_stats: dict[str, list[tuple[int, float]]] = {}

    n_segments_with_overlap = 0
    n_segments_without_overlap = 0
    n_segment_gene_links = 0
    segment_gene_counts: list[int] = []
    focal_values: list[float] = []
    segment_lengths: list[int] = []

    for chrom, gene_list in genes_by_chrom.items():
        chrom_segments = [s for s in sample_segments if s.chrom == chrom]
        if not chrom_segments:
            continue
        segs = sorted(chrom_segments, key=lambda s: (s.start, s.end))
        genes = gene_list
        gi = 0
        for seg in segs:
            while gi < len(genes) and genes[gi].gene_end <= seg.start:
                gi += 1
            gj = gi
            overlaps: list[tuple[str, float]] = []
            while gj < len(genes) and genes[gj].gene_start < seg.end:
                gene = genes[gj]
                ov_start = max(seg.start, gene.gene_start)
                ov_end = min(seg.end, gene.gene_end)
                if ov_end > ov_start:
                    gene_len = max(1, gene.gene_end - gene.gene_start)
                    overlap_frac = float(ov_end - ov_start) / float(gene_len)
                    overlaps.append((gene.gene_id, overlap_frac))
                gj += 1

            if not overlaps:
                n_segments_without_overlap += 1
                continue

            n_segments_with_overlap += 1
            n_genes_in_seg = len(overlaps)
            segment_gene_counts.append(n_genes_in_seg)
            seg_len = _segment_length(seg)
            segment_lengths.append(seg_len)
            focal = focal_penalty(seg_len, cfg.focal_length_scale_bp, cfg.focal_length_alpha)
            focal_values.append(focal)
            gcount_pen = gene_count_penalty(n_genes_in_seg, cfg.gene_count_penalty_mode)
            w_seg = float(focal) * float(gcount_pen)

            for gene_id, overlap_frac in overlaps:
                contrib = float(seg.amplitude) * w_seg * float(overlap_frac)
                numerator[gene_id] = float(numerator.get(gene_id, 0.0)) + contrib
                link_denom[gene_id] = float(link_denom.get(gene_id, 0.0)) + abs(float(overlap_frac))
                contrib_count[gene_id] = int(contrib_count.get(gene_id, 0)) + 1
                prev = float(max_abs_contrib.get(gene_id, 0.0))
                if abs(contrib) >= abs(prev):
                    max_abs_contrib[gene_id] = contrib
                gene_segment_stats.setdefault(gene_id, []).append((seg_len, focal))
                n_segment_gene_links += 1

    scores: dict[str, float] = {}
    for gene_id, numer in numerator.items():
        if cfg.aggregation == "sum":
            scores[gene_id] = float(numer)
        elif cfg.aggregation == "mean":
            denom = max(1, int(contrib_count.get(gene_id, 1)))
            scores[gene_id] = float(numer) / float(denom)
        elif cfg.aggregation == "max_abs":
            scores[gene_id] = float(max_abs_contrib.get(gene_id, 0.0))
        elif cfg.aggregation == "weighted_mean":
            denom = float(link_denom.get(gene_id, 0.0))
            scores[gene_id] = float(numer) / denom if denom > 0 else 0.0
        else:
            raise ValueError(f"Unsupported aggregation mode: {cfg.aggregation}")

    stats: dict[str, object] = {
        "n_segments_total": len(sample_segments),
        "n_segments_with_overlap": n_segments_with_overlap,
        "n_segments_without_overlap": n_segments_without_overlap,
        "n_segment_gene_links": n_segment_gene_links,
        "mean_genes_per_segment_with_overlap": (
            sum(segment_gene_counts) / float(len(segment_gene_counts))
            if segment_gene_counts
            else 0.0
        ),
        "mean_focal_penalty": (sum(focal_values) / float(len(focal_values)) if focal_values else 0.0),
        "median_segment_length_bp": (
            float(statistics.median(segment_lengths))
            if segment_lengths
            else 0.0
        ),
    }
    return scores, stats, gene_segment_stats


def _selected_breadth_stats(
    selected_gene_ids: list[str],
    gene_segment_stats: dict[str, list[tuple[int, float]]],
) -> tuple[float, float]:
    lengths: list[int] = []
    focals: list[float] = []
    for gid in selected_gene_ids:
        for seg_len, focal in gene_segment_stats.get(gid, []):
            lengths.append(int(seg_len))
            focals.append(float(focal))
    if not lengths or not focals:
        return 0.0, 0.0
    return float(statistics.median(lengths)), float(sum(focals) / float(len(focals)))


def _compute_program_scores(program: str, signed_scores: dict[str, float]) -> dict[str, float]:
    if program == "amp":
        return {g: max(float(s), 0.0) for g, s in signed_scores.items() if float(s) > 0.0}
    if program == "del":
        return {g: max(-float(s), 0.0) for g, s in signed_scores.items() if float(s) < 0.0}
    if program == "abs":
        return {g: abs(float(s)) for g, s in signed_scores.items() if float(s) != 0.0}
    raise ValueError(f"Unsupported CNV program: {program}")


def _apply_preset_overrides(cfg: CNVWorkflowConfig) -> CNVWorkflowConfig:
    preset = str(cfg.program_preset).strip().lower()
    if preset != "broad":
        return cfg
    scale = cfg.focal_length_scale_bp
    alpha = cfg.focal_length_alpha
    penalty_mode = cfg.gene_count_penalty_mode
    if float(scale) == 10_000_000:
        scale = 100_000_000
    if float(alpha) == 1.0:
        alpha = 0.25
    if penalty_mode == "inv_sqrt":
        penalty_mode = "none"
    return replace(
        cfg,
        focal_length_scale_bp=scale,
        focal_length_alpha=alpha,
        gene_count_penalty_mode=penalty_mode,
    )


def run_cnv_workflow(
    *,
    cfg: CNVWorkflowConfig,
    segments: list[CNVSegment],
    parse_summary: dict[str, object],
    purity_by_sample: dict[str, float] | None,
    purity_summary: dict[str, object] | None,
    input_files: list[dict[str, str]],
) -> dict[str, object]:
    cfg = _apply_preset_overrides(cfg)
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    genes = read_genes_from_gtf(cfg.gtf, cfg.gtf_gene_id_field)
    if not genes:
        raise ValueError("No genes parsed from GTF.")
    gene_chroms = {g.chrom for g in genes}
    genes_by_chrom = {
        chrom: sorted([g for g in genes if g.chrom == chrom], key=lambda x: (x.gene_start, x.gene_end))
        for chrom in sorted(gene_chroms)
    }
    gene_meta = {
        g.gene_id: {"gene_symbol": g.gene_symbol, "gene_biotype": g.gene_biotype}
        for g in genes
    }

    usable_segments, segment_filter_summary = _filter_and_correct_segments(
        segments=segments,
        gene_chroms=gene_chroms,
        cfg=cfg,
        purity_by_sample=purity_by_sample,
    )
    if not usable_segments:
        raise ValueError("No usable segments after coordinate/chrom/amplitude filtering.")

    if int(parse_summary.get("n_rows_invalid_coords", 0) or 0) > 0:
        print(
            "warning: segments with invalid coordinates were skipped "
            f"(n={int(parse_summary.get('n_rows_invalid_coords', 0) or 0)}).",
            file=sys.stderr,
        )
    if int(parse_summary.get("n_rows_negative_coords", 0) or 0) > 0:
        print(
            "warning: segments with negative coordinates were skipped "
            f"(n={int(parse_summary.get('n_rows_negative_coords', 0) or 0)}).",
            file=sys.stderr,
        )
    mismatch = int(segment_filter_summary.get("n_segments_skipped_chrom_mismatch", 0) or 0)
    if segments and (float(mismatch) / float(len(segments))) > CHROM_MISMATCH_WARN_FRACTION:
        print(
            "warning: many segments could not be mapped to GTF chromosome names/genome build; "
            f"skipped={mismatch}/{len(segments)}. Check --genome_build and --chrom_prefix_mode.",
            file=sys.stderr,
        )
    if cfg.use_purity_correction and purity_by_sample is not None:
        missing = int(segment_filter_summary.get("n_segments_missing_purity", 0) or 0)
        if missing > 0:
            print(
                "warning: purity correction requested but purity values were missing for some segments; "
                f"segments_without_purity={missing}.",
                file=sys.stderr,
            )

    program_methods = _parse_program_methods(cfg.program_methods, cfg.program_preset)
    if not program_methods:
        raise ValueError("No CNV program methods selected; adjust --program_preset or --program_methods.")

    gmt_topk_list = parse_int_list_csv(str(cfg.gmt_topk_list))
    gmt_mass_list = parse_mass_list_csv(str(cfg.gmt_mass_list))
    gmt_biotype_allowlist = parse_str_list_csv(str(cfg.gmt_biotype_allowlist))
    allow_biotypes = {x.lower() for x in gmt_biotype_allowlist} if gmt_biotype_allowlist else None

    sample_ids = sorted({s.sample_id for s in usable_segments})
    manifest_rows: list[dict[str, object]] = []
    combined_signed_scores: dict[str, dict[str, float]] = {}
    combined_gmt_sets: list[tuple[str, list[str]]] = []

    for sample_id in sample_ids:
        sample_segments = [s for s in usable_segments if s.sample_id == sample_id]
        signed_scores, sample_stats, gene_segment_stats = _score_sample_segments(
            sample_segments=sample_segments,
            genes_by_chrom=genes_by_chrom,
            cfg=cfg,
        )
        combined_signed_scores[sample_id] = signed_scores
        for program in program_methods:
            program_scores = _compute_program_scores(program, signed_scores)
            if not program_scores:
                print(
                    f"warning: sample={sample_id} program={program} has no positive signal after scoring; skipped.",
                    file=sys.stderr,
                )
                continue
            selected_gene_ids = _select_gene_ids(program_scores, cfg)
            selected_gene_ids = sorted(
                [gid for gid in selected_gene_ids if gid in program_scores],
                key=lambda gid: (-float(program_scores.get(gid, 0.0)), str(gid)),
            )
            weights = _selected_weights(program_scores, selected_gene_ids, cfg.normalize)

            selected_rows: list[dict[str, object]] = []
            for rank, gene_id in enumerate(selected_gene_ids, start=1):
                meta = gene_meta.get(gene_id, {})
                row: dict[str, object] = {
                    "gene_id": gene_id,
                    "score": float(program_scores.get(gene_id, 0.0)),
                    "weight": float(weights.get(gene_id, 0.0)),
                    "rank": rank,
                }
                symbol = str(meta.get("gene_symbol", "") or "").strip()
                biotype = str(meta.get("gene_biotype", "") or "").strip()
                if symbol:
                    row["gene_symbol"] = symbol
                if biotype:
                    row["gene_biotype"] = biotype
                selected_rows.append(row)

            full_gene_ids = sorted(
                [gid for gid, v in program_scores.items() if float(v) > 0.0],
                key=lambda gid: (-float(program_scores[gid]), str(gid)),
            )
            full_rows: list[dict[str, object]] = []
            for rank, gene_id in enumerate(full_gene_ids, start=1):
                meta = gene_meta.get(gene_id, {})
                row = {
                    "gene_id": gene_id,
                    "score": float(program_scores.get(gene_id, 0.0)),
                    "rank": rank,
                }
                symbol = str(meta.get("gene_symbol", "") or "").strip()
                biotype = str(meta.get("gene_biotype", "") or "").strip()
                if symbol:
                    row["gene_symbol"] = symbol
                if biotype:
                    row["gene_biotype"] = biotype
                full_rows.append(row)

            safe_sample = _safe_component(sample_id, "sample")
            safe_program = _safe_component(program, "program")
            out_subdir = out_dir / f"sample={safe_sample}" / f"program={safe_program}"
            out_subdir.mkdir(parents=True, exist_ok=True)
            _write_rows(out_subdir / "geneset.tsv", selected_rows)
            if cfg.emit_full:
                _write_rows(out_subdir / "geneset.full.tsv", full_rows)

            gmt_sets: list[tuple[str, list[str]]] = []
            gmt_plans: list[dict[str, object]] = []
            gmt_diagnostics: list[dict[str, object]] = []
            gmt_path = resolve_gmt_out_path(out_subdir, cfg.gmt_out)
            if cfg.emit_gmt:
                base_name = (
                    f"{cfg.converter_name}"
                    f"__sample={safe_sample}"
                    f"__program={safe_program}"
                )
                gmt_sets, gmt_plans = build_gmt_sets_from_rows(
                    rows=full_rows,
                    base_name=base_name,
                    prefer_symbol=bool(cfg.gmt_prefer_symbol),
                    min_genes=int(cfg.gmt_min_genes),
                    max_genes=int(cfg.gmt_max_genes),
                    topk_list=gmt_topk_list,
                    mass_list=gmt_mass_list,
                    split_signed=bool(cfg.gmt_split_signed),
                    require_symbol=bool(cfg.gmt_require_symbol),
                    allowed_biotypes=allow_biotypes,
                    emit_small_gene_sets=bool(cfg.emit_small_gene_sets),
                    diagnostics=gmt_diagnostics,
                    context={"sample_id": sample_id, "program": program},
                )
                if gmt_sets:
                    write_gmt(gmt_sets, gmt_path)
                    combined_gmt_sets.extend(gmt_sets)

            warning_count = _warn_gmt_diagnostics(gmt_diagnostics, sample_id=sample_id, program=program)

            median_sel_len, mean_sel_focal = _selected_breadth_stats(selected_gene_ids, gene_segment_stats)
            if median_sel_len > BREADTH_MEDIAN_BP_WARN or (mean_sel_focal > 0 and mean_sel_focal < BREADTH_MEAN_FOCAL_WARN):
                warning_count += 1
                print(
                    "warning: sample="
                    f"{sample_id} program={program} appears breadth-dominated "
                    f"(median_segment_length_bp={median_sel_len:.0f}, mean_focal_penalty={mean_sel_focal:.4f}). "
                    "Consider --program_preset broad or weaker focal penalties for arm-level CNVs.",
                    file=sys.stderr,
                )

            run_summary_payload: dict[str, object] = {
                "converter": cfg.converter_name,
                "dataset_label": cfg.dataset_label,
                "sample_id": sample_id,
                "program": program,
                "program_preset": cfg.program_preset,
                "n_segments_input": sample_stats["n_segments_total"],
                "n_segments_with_overlap": sample_stats["n_segments_with_overlap"],
                "n_segments_without_overlap": sample_stats["n_segments_without_overlap"],
                "n_segments_used": segment_filter_summary["n_segments_used"],
                "n_segments_below_min_abs_amplitude": segment_filter_summary["n_segments_below_min_abs_amplitude"],
                "n_genes_scored": len(program_scores),
                "n_genes_selected": len(selected_rows),
                "median_segment_length_bp_selected": median_sel_len,
                "mean_focal_penalty_selected": mean_sel_focal,
                "requested_gmt_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
                "skipped_gmt_outputs": [
                    d for d in gmt_diagnostics if str(d.get("code")) in {"small_gene_set_skipped", "no_positive_genes"}
                ],
                "parse_summary": parse_summary,
            }
            run_summary_json_path, run_summary_txt_path = write_run_summary_files(out_subdir, run_summary_payload)

            output_files: list[dict[str, object]] = [
                {"path": str(out_subdir / "geneset.tsv"), "role": "selected_program"},
                {"path": str(run_summary_json_path), "role": "run_summary_json"},
                {"path": str(run_summary_txt_path), "role": "run_summary_text"},
            ]
            if cfg.emit_full:
                output_files.append({"path": str(out_subdir / "geneset.full.tsv"), "role": "full_scores"})
            if cfg.emit_gmt and gmt_sets:
                output_files.append({"path": str(gmt_path), "role": "gmt"})

            params: dict[str, object] = {
                "dataset_label": cfg.dataset_label,
                "sample_id": sample_id,
                "program": program,
                "program_preset": cfg.program_preset,
                "program_methods": program_methods,
                "coord_system": cfg.coord_system,
                "chrom_prefix_mode": cfg.chrom_prefix_mode,
                "min_abs_amplitude": cfg.min_abs_amplitude,
                "focal_length_scale_bp": cfg.focal_length_scale_bp,
                "focal_length_alpha": cfg.focal_length_alpha,
                "gene_count_penalty": cfg.gene_count_penalty_mode,
                "aggregation": cfg.aggregation,
                "use_purity_correction": cfg.use_purity_correction,
                "purity_floor": cfg.purity_floor,
                "max_abs_amplitude": cfg.max_abs_amplitude,
                "selection_method": cfg.select,
                "top_k": cfg.top_k,
                "quantile": cfg.quantile,
                "min_score": cfg.min_score,
                "normalize": cfg.normalize,
                "emit_full": cfg.emit_full,
                "emit_gmt": cfg.emit_gmt,
                "gmt_split_signed": cfg.gmt_split_signed,
                "gmt_topk_list": gmt_topk_list,
                "gmt_mass_list": gmt_mass_list,
                "gmt_min_genes": cfg.gmt_min_genes,
                "gmt_max_genes": cfg.gmt_max_genes,
                "emit_small_gene_sets": cfg.emit_small_gene_sets,
            }
            meta = make_metadata(
                converter_name=cfg.converter_name,
                parameters=params,
                data_type="cnv_segments",
                assay="cnv",
                organism=cfg.organism,
                genome_build=cfg.genome_build,
                files=input_files,
                gene_annotation={
                    "mode": "gtf",
                    "gtf_path": str(cfg.gtf),
                    "source": cfg.gtf_source or "user",
                    "gene_id_field": cfg.gtf_gene_id_field,
                },
                weights={
                    "weight_type": "nonnegative",
                    "normalization": {
                        "method": cfg.normalize,
                        "target_sum": 1.0 if cfg.normalize in {"within_set_l1", "l1"} else None,
                    },
                    "aggregation": cfg.aggregation,
                },
                summary={
                    "n_input_features": int(sample_stats["n_segments_total"]),
                    "n_genes": int(len(selected_rows)),
                    "n_features_assigned": int(sample_stats["n_segments_with_overlap"]),
                    "fraction_features_assigned": (
                        float(sample_stats["n_segments_with_overlap"]) / float(sample_stats["n_segments_total"])
                        if int(sample_stats["n_segments_total"]) > 0
                        else 0.0
                    ),
                    "n_segments_used": int(segment_filter_summary["n_segments_used"]),
                    "n_genes_scored": int(len(program_scores)),
                    "n_genes_selected": int(len(selected_rows)),
                    "segment_stats": sample_stats,
                    "segment_filter_summary": segment_filter_summary,
                    "parse_summary": parse_summary,
                    "warnings_count": warning_count,
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
                    "score_definition": (
                        "segment overlap scoring: score_g = sum_i(amplitude_i * focal_penalty_i * "
                        "gene_count_penalty_i * overlap_frac_ig), with per-gene aggregation="
                        f"{cfg.aggregation}"
                    ),
                },
                output_files=output_files,
                gmt={
                    "written": bool(cfg.emit_gmt and gmt_sets),
                    "path": (
                        str(gmt_path.relative_to(out_subdir))
                        if cfg.emit_gmt and gmt_sets and gmt_path.is_relative_to(out_subdir)
                        else str(gmt_path)
                    )
                    if cfg.emit_gmt and gmt_sets
                    else None,
                    "prefer_symbol": bool(cfg.gmt_prefer_symbol),
                    "min_genes": int(cfg.gmt_min_genes),
                    "max_genes": int(cfg.gmt_max_genes),
                    "plans": gmt_plans,
                },
            )
            write_metadata(out_subdir / "geneset.meta.json", meta)

            manifest_rows.append(
                {
                    "sample_id": sample_id,
                    "program": program,
                    "path": str(out_subdir.relative_to(out_dir)),
                    "n_segments": int(sample_stats["n_segments_total"]),
                    "n_segments_used": int(segment_filter_summary["n_segments_used"]),
                    "n_genes_scored": int(len(program_scores)),
                    "warnings_count": int(warning_count),
                }
            )

    if not manifest_rows:
        raise ValueError("No CNV outputs were emitted after filtering/scoring.")

    manifest_path = out_dir / "manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "sample_id",
                "program",
                "path",
                "n_segments",
                "n_segments_used",
                "n_genes_scored",
                "warnings_count",
            ],
        )
        writer.writeheader()
        for row in manifest_rows:
            writer.writerow(row)

    if cfg.emit_cohort_sets:
        _emit_cohort_gmt(
            cfg=cfg,
            out_dir=out_dir,
            combined_signed_scores=combined_signed_scores,
            gene_meta=gene_meta,
            gmt_topk_list=gmt_topk_list,
            gmt_mass_list=gmt_mass_list,
            allow_biotypes=allow_biotypes,
        )

    return {
        "n_groups": len(manifest_rows),
        "n_peaks": len(segments),
        "n_genes": len({gid for scores in combined_signed_scores.values() for gid in scores}),
        "out_dir": str(out_dir),
        "program_methods": program_methods,
    }


def _emit_cohort_gmt(
    *,
    cfg: CNVWorkflowConfig,
    out_dir: Path,
    combined_signed_scores: dict[str, dict[str, float]],
    gene_meta: dict[str, dict[str, object]],
    gmt_topk_list: list[int],
    gmt_mass_list: list[float],
    allow_biotypes: set[str] | None,
) -> None:
    n_samples = len(combined_signed_scores)
    if n_samples < int(cfg.cohort_min_samples):
        print(
            "warning: cohort set emission skipped; "
            f"n_samples={n_samples} < cohort_min_samples={cfg.cohort_min_samples}.",
            file=sys.stderr,
        )
        return

    amp_counts: dict[str, int] = {}
    del_counts: dict[str, int] = {}
    for _sample, scores in combined_signed_scores.items():
        for gene_id, value in scores.items():
            v = float(value)
            if max(v, 0.0) > float(cfg.cohort_score_threshold):
                amp_counts[gene_id] = int(amp_counts.get(gene_id, 0)) + 1
            if max(-v, 0.0) > float(cfg.cohort_score_threshold):
                del_counts[gene_id] = int(del_counts.get(gene_id, 0)) + 1

    amp_rows: list[dict[str, object]] = []
    del_rows: list[dict[str, object]] = []
    for gene_id, n in amp_counts.items():
        frac = float(n) / float(n_samples)
        if frac < float(cfg.cohort_min_fraction):
            continue
        meta = gene_meta.get(gene_id, {})
        row: dict[str, object] = {"gene_id": gene_id, "score": frac, "rank": 0}
        sym = str(meta.get("gene_symbol", "") or "").strip()
        bio = str(meta.get("gene_biotype", "") or "").strip()
        if sym:
            row["gene_symbol"] = sym
        if bio:
            row["gene_biotype"] = bio
        amp_rows.append(row)
    for gene_id, n in del_counts.items():
        frac = float(n) / float(n_samples)
        if frac < float(cfg.cohort_min_fraction):
            continue
        meta = gene_meta.get(gene_id, {})
        row = {"gene_id": gene_id, "score": frac, "rank": 0}
        sym = str(meta.get("gene_symbol", "") or "").strip()
        bio = str(meta.get("gene_biotype", "") or "").strip()
        if sym:
            row["gene_symbol"] = sym
        if bio:
            row["gene_biotype"] = bio
        del_rows.append(row)

    for rows in (amp_rows, del_rows):
        rows.sort(key=lambda r: (-float(r["score"]), str(r["gene_id"])))
        for i, row in enumerate(rows, start=1):
            row["rank"] = i

    cohort_dir = out_dir / "cohort"
    cohort_dir.mkdir(parents=True, exist_ok=True)
    combined_gmt: list[tuple[str, list[str]]] = []
    for program, rows in (("recurrent_amp", amp_rows), ("recurrent_del", del_rows)):
        if not rows:
            continue
        diags: list[dict[str, object]] = []
        sets, _plans = build_gmt_sets_from_rows(
            rows=rows,
            base_name=f"cnv_gene_extractor__program={program}",
            prefer_symbol=bool(cfg.gmt_prefer_symbol),
            min_genes=int(cfg.gmt_min_genes),
            max_genes=int(cfg.gmt_max_genes),
            topk_list=gmt_topk_list,
            mass_list=gmt_mass_list,
            split_signed=False,
            require_symbol=bool(cfg.gmt_require_symbol),
            allowed_biotypes=allow_biotypes,
            emit_small_gene_sets=bool(cfg.emit_small_gene_sets),
            diagnostics=diags,
            context={"program": program, "scope": "cohort"},
        )
        _warn_gmt_diagnostics(diags, sample_id="cohort", program=program)
        combined_gmt.extend(sets)
    if combined_gmt:
        write_gmt(combined_gmt, cohort_dir / "genesets.gmt")
        run_summary_payload: dict[str, object] = {
            "converter": cfg.converter_name,
            "scope": "cohort",
            "n_samples": n_samples,
            "cohort_score_threshold": cfg.cohort_score_threshold,
            "cohort_min_fraction": cfg.cohort_min_fraction,
            "n_amp_genes": len(amp_rows),
            "n_del_genes": len(del_rows),
        }
        write_run_summary_files(cohort_dir, run_summary_payload)

