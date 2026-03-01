from __future__ import annotations

import csv
from dataclasses import dataclass
import gzip
from pathlib import Path
import re
import sys

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
from omics2geneset.rnaseq.deg_scoring import (
    annotate_records_with_gtf,
    apply_exclude_filters,
    compile_exclude_gene_regexes,
    maybe_promote_gene_id_to_symbol,
    sanitize_name_component,
)


KNOWN_LOADINGS_FORMATS = {
    "auto",
    "wide_genes_by_program",
    "wide_programs_by_gene",
    "long_tidy",
    "cnmf_gene_spectra_tpm",
    "cnmf_gene_spectra_score",
    "schpf_gene_scores",
}
WIDE_GENES_BY_PROGRAM_FORMATS = {
    "wide_genes_by_program",
    "cnmf_gene_spectra_tpm",
    "cnmf_gene_spectra_score",
    "schpf_gene_scores",
}

_SAFE_COMPONENT_RE = re.compile(r"[^A-Za-z0-9._-]+")
NEGATIVE_WARNING_FRACTION = 0.2
LARGE_PROGRAM_COUNT_WARNING_THRESHOLD = 200
LARGE_PARSED_VALUES_WARNING_THRESHOLD = 5_000_000
LARGE_INPUT_GUIDE_HINT = (
    "See docs/rna-seq-to-geneset.md#best-practices-for-large-scrna-maps-recommended-upstream-factorization-workflow"
)


@dataclass
class SCRNAProgramsWorkflowConfig:
    converter_name: str
    out_dir: Path
    organism: str
    genome_build: str
    dataset_label: str
    signature_name: str
    source_path: str
    source_kind: str
    loadings_format: str
    score_transform: str
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
    exclude_gene_regex: list[str] | None
    disable_default_excludes: bool
    gtf: str | None
    gtf_gene_id_field: str
    gtf_source: str | None
    program_id_prefix: str | None


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


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
    if value != value:
        return None
    return float(value)


def _safe_path_component(value: str, fallback: str) -> str:
    out = _SAFE_COMPONENT_RE.sub("_", str(value)).strip("_")
    return out or fallback


def _resolve_format_auto(
    path: Path,
    fieldnames: list[str],
    gene_id_column: str,
    program_id_column: str,
    loading_column: str,
) -> str:
    name = path.name.lower()
    if "gene_spectra_tpm" in name:
        return "cnmf_gene_spectra_tpm"
    if "gene_spectra_score" in name:
        return "cnmf_gene_spectra_score"
    if "schpf" in name:
        return "schpf_gene_scores"
    field_set = set(fieldnames)
    if {gene_id_column, program_id_column, loading_column}.issubset(field_set):
        return "long_tidy"
    if gene_id_column in field_set:
        return "wide_genes_by_program"
    if program_id_column in field_set:
        return "wide_programs_by_gene"
    return "wide_genes_by_program"


def _resolve_gene_column(fieldnames: list[str], gene_id_column: str) -> str:
    preferred = [
        gene_id_column,
        "gene_id",
        "gene",
        "gene_name",
        "gene_symbol",
    ]
    for candidate in preferred:
        if candidate in fieldnames:
            return candidate
    return fieldnames[0]


def _resolve_program_column(fieldnames: list[str], program_id_column: str) -> str:
    preferred = [
        program_id_column,
        "program_id",
        "program",
        "factor",
        "topic",
    ]
    for candidate in preferred:
        if candidate in fieldnames:
            return candidate
    return fieldnames[0]


def read_program_loadings(
    *,
    path: str | Path,
    loadings_format: str,
    gene_id_column: str,
    program_id_column: str,
    loading_column: str,
    transpose: bool,
) -> tuple[dict[str, dict[str, float]], dict[str, object]]:
    p = Path(path)
    requested = str(loadings_format).strip()
    if requested not in KNOWN_LOADINGS_FORMATS:
        raise ValueError(
            "Unsupported loadings_format: "
            f"{requested}. Expected one of {', '.join(sorted(KNOWN_LOADINGS_FORMATS))}"
        )
    programs: dict[str, dict[str, float]] = {}
    values_parsed = 0
    values_non_numeric = 0
    n_rows = 0
    used_gene_column: str | None = None
    used_program_column: str | None = None
    resolved_format: str
    effective_format: str

    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Program loadings table has no header: {p}")
        fieldnames = [str(x) for x in reader.fieldnames]

        resolved_format = requested
        if requested == "auto":
            resolved_format = _resolve_format_auto(
                p,
                fieldnames,
                gene_id_column=gene_id_column,
                program_id_column=program_id_column,
                loading_column=loading_column,
            )

        effective_format = resolved_format
        if transpose and resolved_format in {"wide_genes_by_program", "wide_programs_by_gene"}:
            effective_format = (
                "wide_programs_by_gene" if resolved_format == "wide_genes_by_program" else "wide_genes_by_program"
            )

        if effective_format == "long_tidy":
            if gene_id_column not in fieldnames:
                raise ValueError(
                    f"long_tidy format requires gene_id_column '{gene_id_column}'. "
                    f"Available columns: {', '.join(fieldnames)}"
                )
            if program_id_column not in fieldnames:
                raise ValueError(
                    f"long_tidy format requires program_id_column '{program_id_column}'. "
                    f"Available columns: {', '.join(fieldnames)}"
                )
            if loading_column not in fieldnames:
                raise ValueError(
                    f"long_tidy format requires loading_column '{loading_column}'. "
                    f"Available columns: {', '.join(fieldnames)}"
                )
            used_gene_column = gene_id_column
            used_program_column = program_id_column
            for row in reader:
                n_rows += 1
                gene_id = str(row.get(gene_id_column, "")).strip()
                program_id = str(row.get(program_id_column, "")).strip()
                raw_loading = str(row.get(loading_column, "")).strip()
                if not gene_id or not program_id:
                    continue
                val = _parse_float(raw_loading)
                if val is None:
                    if raw_loading:
                        values_non_numeric += 1
                    continue
                programs.setdefault(program_id, {})
                programs[program_id][gene_id] = float(programs[program_id].get(gene_id, 0.0)) + float(val)
                values_parsed += 1
        elif effective_format in WIDE_GENES_BY_PROGRAM_FORMATS:
            gene_col = _resolve_gene_column(fieldnames, gene_id_column)
            used_gene_column = gene_col
            program_cols = [name for name in fieldnames if name != gene_col]
            if not program_cols:
                raise ValueError("wide_genes_by_program format requires at least one program column")
            for row in reader:
                n_rows += 1
                gene_id = str(row.get(gene_col, "")).strip()
                if not gene_id:
                    continue
                for prog_col in program_cols:
                    raw_loading = str(row.get(prog_col, "")).strip()
                    if not raw_loading:
                        continue
                    val = _parse_float(raw_loading)
                    if val is None:
                        values_non_numeric += 1
                        continue
                    programs.setdefault(str(prog_col), {})
                    programs[str(prog_col)][gene_id] = float(programs[str(prog_col)].get(gene_id, 0.0)) + float(val)
                    values_parsed += 1
        elif effective_format == "wide_programs_by_gene":
            program_col = _resolve_program_column(fieldnames, program_id_column)
            used_program_column = program_col
            gene_cols = [name for name in fieldnames if name != program_col]
            if not gene_cols:
                raise ValueError("wide_programs_by_gene format requires at least one gene column")
            for row in reader:
                n_rows += 1
                program_id = str(row.get(program_col, "")).strip()
                if not program_id:
                    continue
                programs.setdefault(program_id, {})
                for gene_col in gene_cols:
                    raw_loading = str(row.get(gene_col, "")).strip()
                    if not raw_loading:
                        continue
                    val = _parse_float(raw_loading)
                    if val is None:
                        values_non_numeric += 1
                        continue
                    programs[program_id][str(gene_col)] = float(programs[program_id].get(str(gene_col), 0.0)) + float(val)
                    values_parsed += 1
        else:
            raise ValueError(f"Unsupported effective loadings format: {effective_format}")

    if not programs:
        raise ValueError("No program loadings were parsed from input table")

    n_genes_unique = len({gid for by_gene in programs.values() for gid in by_gene})
    parse_summary: dict[str, object] = {
        "input_path": str(p),
        "requested_loadings_format": requested,
        "resolved_loadings_format": resolved_format,
        "effective_loadings_format": effective_format,
        "transpose": bool(transpose),
        "n_rows": n_rows,
        "n_programs": len(programs),
        "n_genes_unique": n_genes_unique,
        "n_values_parsed": values_parsed,
        "n_values_non_numeric": values_non_numeric,
        "used_gene_column": used_gene_column,
        "used_program_column": used_program_column,
    }
    return programs, parse_summary


def _select_gene_ids(scores: dict[str, float], cfg: SCRNAProgramsWorkflowConfig) -> list[str]:
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


def _transform_score(raw_value: float, transform: str) -> float:
    if transform == "signed":
        return float(raw_value)
    if transform == "abs":
        return abs(float(raw_value))
    if transform == "positive":
        return max(float(raw_value), 0.0)
    if transform == "negative":
        return max(-float(raw_value), 0.0)
    raise ValueError(f"Unsupported score_transform: {transform}")


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


def run_sc_programs_workflow(
    *,
    cfg: SCRNAProgramsWorkflowConfig,
    programs: dict[str, dict[str, float]],
    parse_summary: dict[str, object],
    input_files: list[dict[str, str]],
) -> dict[str, object]:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    exclude_patterns = compile_exclude_gene_regexes(
        cfg.exclude_gene_regex,
        cfg.disable_default_excludes,
    )

    combined_gmt_sets: list[tuple[str, list[str]]] = []
    manifest_rows: list[tuple[str, str]] = []
    safe_ids_seen: set[str] = set()

    total_negative_raw = 0
    total_raw = 0
    programs_processed = 0

    safe_signature = sanitize_name_component(cfg.signature_name)
    safe_dataset = sanitize_name_component(cfg.dataset_label)
    gmt_topk_list = parse_int_list_csv(str(cfg.gmt_topk_list))
    gmt_mass_list = parse_mass_list_csv(str(cfg.gmt_mass_list))
    gmt_biotype_allowlist = parse_str_list_csv(str(cfg.gmt_biotype_allowlist))

    n_programs_summary = int(parse_summary.get("n_programs", 0) or 0)
    n_values_summary = int(parse_summary.get("n_values_parsed", 0) or 0)
    if n_programs_summary > LARGE_PROGRAM_COUNT_WARNING_THRESHOLD:
        print(
            "warning: large program-loading input detected "
            f"(n_programs={n_programs_summary} > {LARGE_PROGRAM_COUNT_WARNING_THRESHOLD}). "
            "This usually means too many factors were converted at once. "
            "Best practice: run factorization within cell types or broad compartments, "
            "cap cells per donor, reduce K, and/or convert per cell type separately. "
            f"{LARGE_INPUT_GUIDE_HINT}",
            file=sys.stderr,
        )
    if n_values_summary > LARGE_PARSED_VALUES_WARNING_THRESHOLD:
        print(
            "warning: very large parsed loading table detected "
            f"(n_values_parsed={n_values_summary} > {LARGE_PARSED_VALUES_WARNING_THRESHOLD}). "
            "Best practice: run factorization within cell types or broad compartments, "
            "cap cells per donor, reduce K, and/or split conversion by cell type. "
            f"{LARGE_INPUT_GUIDE_HINT}",
            file=sys.stderr,
        )

    for program_id in sorted(programs):
        raw_scores = {str(gid): float(val) for gid, val in programs[program_id].items()}
        if not raw_scores:
            continue
        programs_processed += 1
        n_raw = len(raw_scores)
        n_negative_raw = sum(1 for val in raw_scores.values() if float(val) < 0.0)
        total_raw += n_raw
        total_negative_raw += n_negative_raw

        if cfg.score_transform == "positive" and n_raw > 0:
            negative_fraction = float(n_negative_raw) / float(n_raw)
            if negative_fraction >= NEGATIVE_WARNING_FRACTION:
                print(
                    "warning: program="
                    f"{program_id} has {n_negative_raw}/{n_raw} negative raw loadings; "
                    "score_transform=positive drops negatives. Use --score_transform signed "
                    "and --gmt_split_signed true to keep directional sets.",
                    file=sys.stderr,
                )

        records: dict[str, dict[str, object]] = {
            str(gene_id): {
                "gene_id": str(gene_id),
                "gene_symbol": None,
                "gene_biotype": None,
            }
            for gene_id in raw_scores
        }
        if cfg.gtf:
            records = annotate_records_with_gtf(records, cfg.gtf, cfg.gtf_gene_id_field)
        symbol_promotion = maybe_promote_gene_id_to_symbol(records)
        records = symbol_promotion["records"]

        records, n_filtered = apply_exclude_filters(records, exclude_patterns)
        transformed_scores: dict[str, float] = {}
        for gene_id in records:
            transformed_scores[gene_id] = _transform_score(float(raw_scores.get(gene_id, 0.0)), cfg.score_transform)

        if cfg.score_transform == "signed":
            selection_basis = {gid: abs(float(score)) for gid, score in transformed_scores.items() if float(score) != 0.0}
        else:
            selection_basis = {gid: float(score) for gid, score in transformed_scores.items() if float(score) > 0.0}

        selected_gene_ids = _select_gene_ids(selection_basis, cfg)
        selected_gene_ids = sorted(
            [gid for gid in selected_gene_ids if gid in selection_basis],
            key=lambda gid: (-float(selection_basis.get(gid, 0.0)), str(gid)),
        )
        weights = _selected_weights(selection_basis, selected_gene_ids, cfg.normalize)

        selected_rows: list[dict[str, object]] = []
        for rank, gene_id in enumerate(selected_gene_ids, start=1):
            rec = records[gene_id]
            row: dict[str, object] = {
                "gene_id": gene_id,
                "score": float(transformed_scores.get(gene_id, 0.0)),
                "weight": float(weights.get(gene_id, 0.0)),
                "rank": rank,
            }
            symbol_obj = rec.get("gene_symbol")
            symbol = "" if symbol_obj is None else str(symbol_obj).strip()
            if symbol:
                row["gene_symbol"] = symbol
            biotype_obj = rec.get("gene_biotype")
            biotype = "" if biotype_obj is None else str(biotype_obj).strip()
            if biotype:
                row["gene_biotype"] = biotype
            selected_rows.append(row)

        full_gene_ids = sorted(
            [gid for gid, score in transformed_scores.items() if float(score) != 0.0],
            key=lambda gid: (-abs(float(transformed_scores[gid])), str(gid)),
        )
        full_rows: list[dict[str, object]] = []
        for rank, gene_id in enumerate(full_gene_ids, start=1):
            rec = records[gene_id]
            row: dict[str, object] = {
                "gene_id": gene_id,
                "score": float(transformed_scores.get(gene_id, 0.0)),
                "rank": rank,
            }
            symbol_obj = rec.get("gene_symbol")
            symbol = "" if symbol_obj is None else str(symbol_obj).strip()
            if symbol:
                row["gene_symbol"] = symbol
            biotype_obj = rec.get("gene_biotype")
            biotype = "" if biotype_obj is None else str(biotype_obj).strip()
            if biotype:
                row["gene_biotype"] = biotype
            full_rows.append(row)

        base_program_id = str(program_id)
        if cfg.program_id_prefix:
            base_program_id = f"{str(cfg.program_id_prefix).strip()}{base_program_id}"
        safe_program_id = _safe_path_component(base_program_id, "program")
        dedup_safe = safe_program_id
        suffix = 2
        while dedup_safe in safe_ids_seen:
            dedup_safe = f"{safe_program_id}_{suffix}"
            suffix += 1
        safe_ids_seen.add(dedup_safe)

        program_dir = out_dir / f"program={dedup_safe}"
        program_dir.mkdir(parents=True, exist_ok=True)
        manifest_rows.append((base_program_id, str(program_dir.relative_to(out_dir))))

        _write_rows(program_dir / "geneset.tsv", selected_rows)
        if cfg.emit_full:
            _write_rows(program_dir / "geneset.full.tsv", full_rows)

        gmt_sets: list[tuple[str, list[str]]] = []
        gmt_plans: list[dict[str, object]] = []
        gmt_diagnostics: list[dict[str, object]] = []
        gmt_path = resolve_gmt_out_path(program_dir, cfg.gmt_out)
        if cfg.emit_gmt:
            base_name = (
                f"{cfg.converter_name}"
                f"__dataset={safe_dataset}"
                f"__signature={safe_signature}"
                f"__program={sanitize_name_component(base_program_id)}"
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
                allowed_biotypes={x.lower() for x in gmt_biotype_allowlist} if gmt_biotype_allowlist else None,
                emit_small_gene_sets=bool(cfg.emit_small_gene_sets),
                diagnostics=gmt_diagnostics,
                context={
                    "converter": cfg.converter_name,
                    "program_id": base_program_id,
                },
            )
            write_gmt(gmt_sets, gmt_path)
            combined_gmt_sets.extend(gmt_sets)
        _warn_gmt_diagnostics(gmt_diagnostics)

        output_files: list[dict[str, object]] = [
            {"path": str(program_dir / "geneset.tsv"), "role": "selected_program"},
        ]
        if cfg.emit_full:
            output_files.append({"path": str(program_dir / "geneset.full.tsv"), "role": "full_scores"})
        if cfg.emit_gmt:
            output_files.append({"path": str(gmt_path), "role": "gmt"})

        run_summary_payload: dict[str, object] = {
            "converter": cfg.converter_name,
            "dataset_label": cfg.dataset_label,
            "signature_name": cfg.signature_name,
            "program_id": base_program_id,
            "score_transform": cfg.score_transform,
            "n_input_features": n_raw,
            "n_negative_raw_loadings": n_negative_raw,
            "fraction_negative_raw_loadings": (float(n_negative_raw) / float(n_raw) if n_raw else 0.0),
            "n_genes_after_filter": len(records),
            "n_genes_selected": len(selected_rows),
            "n_genes_filtered_by_symbol_regex": n_filtered,
            "symbol_promotion": {
                "promoted": bool(symbol_promotion.get("promoted", False)),
                "reason": str(symbol_promotion.get("reason", "")),
                "n_promoted": int(symbol_promotion.get("n_promoted", 0) or 0),
                "ensembl_like_fraction": float(symbol_promotion.get("ensembl_like_fraction", 0.0) or 0.0),
            },
            "requested_gmt_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
            "skipped_gmt_outputs": _collect_skipped_gmt_outputs(gmt_diagnostics),
        }
        run_summary_json_path, run_summary_txt_path = write_run_summary_files(program_dir, run_summary_payload)
        output_files.append({"path": str(run_summary_json_path), "role": "run_summary_json"})
        output_files.append({"path": str(run_summary_txt_path), "role": "run_summary_text"})

        params: dict[str, object] = {
            "dataset_label": cfg.dataset_label,
            "signature_name": cfg.signature_name,
            "program_id": base_program_id,
            "source_kind": cfg.source_kind,
            "source_path": cfg.source_path,
            "loadings_format": cfg.loadings_format,
            "score_transform": cfg.score_transform,
            "selection_method": cfg.select,
            "top_k": cfg.top_k,
            "quantile": cfg.quantile,
            "min_score": cfg.min_score,
            "normalize": cfg.normalize,
            "exclude_gene_regex": cfg.exclude_gene_regex or [],
            "disable_default_excludes": cfg.disable_default_excludes,
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
            data_type="rna_seq",
            assay="single_cell",
            organism=cfg.organism,
            genome_build=cfg.genome_build,
            files=input_files,
            gene_annotation=(
                {
                    "mode": "gtf",
                    "gtf_path": str(cfg.gtf),
                    "source": cfg.gtf_source or "user",
                    "gene_id_field": cfg.gtf_gene_id_field,
                }
                if cfg.gtf
                else {
                    "mode": "none",
                    "source": "none",
                    "gene_id_field": "gene_id",
                }
            ),
            weights={
                "weight_type": ("signed" if cfg.score_transform == "signed" else "nonnegative"),
                "normalization": {
                    "method": cfg.normalize,
                    "target_sum": 1.0 if cfg.normalize in {"within_set_l1", "l1"} else None,
                },
                "aggregation": "none",
            },
            summary={
                "n_input_features": int(n_raw),
                "n_genes": int(len(selected_rows)),
                "n_features_assigned": int(n_raw),
                "fraction_features_assigned": (1.0 if n_raw else 0.0),
                "n_negative_raw_loadings": int(n_negative_raw),
                "n_genes_after_filter": int(len(records)),
                "n_genes_filtered_by_symbol_regex": int(n_filtered),
                "symbol_promotion": {
                    "promoted": bool(symbol_promotion.get("promoted", False)),
                    "reason": str(symbol_promotion.get("reason", "")),
                    "n_promoted": int(symbol_promotion.get("n_promoted", 0) or 0),
                    "ensembl_like_fraction": float(symbol_promotion.get("ensembl_like_fraction", 0.0) or 0.0),
                },
                "parse_summary": parse_summary,
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
                    "program loading score from external factorization output; "
                    f"score_transform={cfg.score_transform}; selection per program"
                ),
            },
            output_files=output_files,
            gmt={
                "written": bool(cfg.emit_gmt),
                "path": (
                    str(gmt_path.relative_to(program_dir))
                    if cfg.emit_gmt and gmt_path.is_relative_to(program_dir)
                    else str(gmt_path)
                )
                if cfg.emit_gmt
                else None,
                "prefer_symbol": bool(cfg.gmt_prefer_symbol),
                "require_symbol": bool(cfg.gmt_require_symbol),
                "biotype_allowlist": gmt_biotype_allowlist,
                "min_genes": int(cfg.gmt_min_genes),
                "max_genes": int(cfg.gmt_max_genes),
                "emit_small_gene_sets": bool(cfg.emit_small_gene_sets),
                "requested_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
                "emitted_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
                "skipped_outputs": _collect_skipped_gmt_outputs(gmt_diagnostics),
                "diagnostics": gmt_diagnostics,
                "plans": gmt_plans,
            },
        )
        write_metadata(program_dir / "geneset.meta.json", meta)

    if not manifest_rows:
        raise ValueError("No valid programs were parsed after filtering.")

    manifest_path = out_dir / "manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["program_id", "path"])
        writer.writerows(manifest_rows)

    if cfg.emit_gmt and combined_gmt_sets:
        write_gmt(combined_gmt_sets, out_dir / "genesets.gmt")

    if cfg.score_transform == "positive" and total_raw > 0:
        overall_fraction = float(total_negative_raw) / float(total_raw)
        if overall_fraction >= NEGATIVE_WARNING_FRACTION:
            print(
                "warning: input loadings contain substantial negative values overall "
                f"({total_negative_raw}/{total_raw}); score_transform=positive drops negatives. "
                "Use --score_transform signed and optionally --gmt_split_signed true "
                "if directional programs are expected.",
                file=sys.stderr,
            )

    return {
        "n_groups": len(manifest_rows),
        "out_dir": str(out_dir),
        "n_programs_processed": programs_processed,
        "parse_summary": parse_summary,
    }
