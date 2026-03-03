from __future__ import annotations

import csv
from pathlib import Path
import re

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import input_file_record
from geneset_extractors.extractors.rnaseq.deg_scoring import DEGRow, read_deg_tsv, sanitize_name_component
from geneset_extractors.extractors.rnaseq.deg_workflow import DEGWorkflowConfig, run_deg_workflow


def _safe_name(value: str) -> str:
    out = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_")
    return out or "comparison"


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fieldnames, rows = read_deg_tsv(args.deg_tsv)
    if args.comparison_column not in fieldnames:
        raise ValueError(
            f"comparison_column '{args.comparison_column}' not found in DE table. "
            f"Available columns: {', '.join(fieldnames)}"
        )

    grouped: dict[str, list[DEGRow]] = {}
    for row in rows:
        comparison = str(row.values.get(args.comparison_column, "")).strip()
        if not comparison:
            continue
        grouped.setdefault(comparison, []).append(row)
    if not grouped:
        raise ValueError("No non-empty comparison labels found in comparison_column.")

    signature_name = str(args.signature_name or "").strip()
    if not signature_name or signature_name == "contrast":
        signature_name = sanitize_name_component(Path(args.deg_tsv).stem)

    files = [input_file_record(args.deg_tsv, "deg_tsv")]
    if args.gtf:
        files.append(input_file_record(args.gtf, "gtf"))

    used_paths: set[str] = set()
    manifest_rows: list[tuple[str, str]] = []
    combined_gmt_sets: list[tuple[str, list[str]]] = []
    biotype_warning_seen = False

    for comparison in sorted(grouped):
        base = _safe_name(comparison)
        safe = base
        suffix = 2
        while safe in used_paths:
            safe = f"{base}_{suffix}"
            suffix += 1
        used_paths.add(safe)

        group_dir = out_dir / f"comparison={safe}"
        cfg = DEGWorkflowConfig(
            converter_name="rna_deg_multi",
            out_dir=group_dir,
            organism=args.organism,
            genome_build=args.genome_build,
            signature_name=signature_name,
            deg_tsv_label=Path(args.deg_tsv).name,
            comparison_label=comparison,
            gene_id_column=args.gene_id_column,
            gene_symbol_column=args.gene_symbol_column,
            stat_column=args.stat_column,
            logfc_column=args.logfc_column,
            padj_column=args.padj_column,
            pvalue_column=args.pvalue_column,
            score_column=args.score_column,
            score_mode=args.score_mode,
            neglog10p_cap=args.neglog10p_cap,
            neglog10p_eps=args.neglog10p_eps,
            duplicate_gene_policy=args.duplicate_gene_policy,
            exclude_gene_regex=args.exclude_gene_regex,
            disable_default_excludes=args.disable_default_excludes,
            gtf=args.gtf,
            gtf_gene_id_field=args.gtf_gene_id_field,
            gtf_source=args.gtf_source,
            select=args.select,
            top_k=args.top_k,
            quantile=args.quantile,
            min_score=args.min_score,
            normalize=args.normalize,
            emit_full=args.emit_full,
            emit_gmt=args.emit_gmt,
            gmt_out=args.gmt_out,
            gmt_prefer_symbol=args.gmt_prefer_symbol,
            gmt_require_symbol=args.gmt_require_symbol,
            gmt_biotype_allowlist=args.gmt_biotype_allowlist,
            gmt_min_genes=args.gmt_min_genes,
            gmt_max_genes=args.gmt_max_genes,
            gmt_topk_list=args.gmt_topk_list,
            gmt_mass_list=args.gmt_mass_list,
            gmt_split_signed=args.gmt_split_signed,
            gmt_emit_abs=args.gmt_emit_abs,
            gmt_source=args.gmt_source,
            emit_small_gene_sets=args.emit_small_gene_sets,
            warn_biotype_missing=not biotype_warning_seen,
        )
        result = run_deg_workflow(
            cfg=cfg,
            fieldnames=fieldnames,
            rows=grouped[comparison],
            input_files=files,
        )
        if bool(result.get("biotype_warning_emitted", False)):
            biotype_warning_seen = True
        manifest_rows.append((comparison, str(group_dir.relative_to(out_dir))))
        if args.emit_gmt:
            combined_gmt_sets.extend(result.get("gmt_sets", []))

    with (out_dir / "manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["comparison", "path"])
        writer.writerows(manifest_rows)

    if args.emit_gmt and combined_gmt_sets:
        write_gmt(combined_gmt_sets, out_dir / "genesets.gmt")

    return {
        "n_groups": len(manifest_rows),
        "out_dir": str(out_dir),
    }
