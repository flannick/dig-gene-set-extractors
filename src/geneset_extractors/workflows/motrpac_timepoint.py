from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path

from geneset_extractors.workflows.gtex_runtime_common import write_tsv, write_workflow_provenance_graph
from geneset_extractors.workflows.motrpac_common import prepare_tissue_inputs, write_prepared_tissue_inputs


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def _resolve_rscript_bin(rscript_bin: str) -> str:
    resolved = shutil.which(rscript_bin) if not Path(rscript_bin).is_absolute() else rscript_bin
    if not resolved or not Path(resolved).exists():
        raise ValueError(f"Rscript not found: {rscript_bin}")
    return resolved


def _slugify_tissue_id(tissue_id: str) -> str:
    return str(tissue_id).strip().lower().replace("_", "-")


def _build_comparison_rows(
    *,
    metadata_rows: list[dict[str, str]],
    tissue_slug: str,
    min_samples_per_group: int,
) -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    counts_by_timepoint: dict[tuple[str, str], dict[str, int]] = {}
    for row in metadata_rows:
        timepoint_label = str(row.get("timepoint_label", "")).strip()
        tissue_code_no = str(row.get("tissue_code_no", "")).strip().lower()
        intervention = str(row.get("intervention", "")).strip().lower()
        if not timepoint_label or not tissue_code_no or intervention not in {"control", "training"}:
            continue
        key = (timepoint_label, tissue_code_no)
        counts_by_timepoint.setdefault(key, {"control": 0, "training": 0})
        counts_by_timepoint[key][intervention] += 1

    comparisons: list[dict[str, str]] = []
    summary_rows: list[dict[str, object]] = []
    for timepoint_label, tissue_code_no in sorted(counts_by_timepoint):
        counts = counts_by_timepoint[(timepoint_label, tissue_code_no)]
        comparison_id = f"{tissue_code_no}-{tissue_slug}_{timepoint_label}"
        summary_rows.append(
            {
                "comparison_id": comparison_id,
                "timepoint_label": timepoint_label,
                "n_control": counts["control"],
                "n_training": counts["training"],
            }
        )
        if counts["control"] < int(min_samples_per_group) or counts["training"] < int(min_samples_per_group):
            continue
        comparisons.append(
            {
                "comparison_id": comparison_id,
                "comparison_kind": "condition_a_vs_b",
                "group_column": "intervention",
                "group_a": "training",
                "group_b": "control",
                "timepoint_label": timepoint_label,
                "tissue_code_no": tissue_code_no,
                "tissue_slug": tissue_slug,
            }
        )
    return comparisons, summary_rows


def _write_timepoint_r_script(
    *,
    script_path: Path,
    counts_tsv: Path,
    metadata_tsv: Path,
    comparisons_tsv: Path,
    output_tsv: Path,
) -> None:
    script = f'''suppressPackageStartupMessages({{
  library(edgeR)
  library(limma)
}})

counts <- read.delim("{counts_tsv}", check.names=FALSE)
meta <- read.delim("{metadata_tsv}", check.names=FALSE)
comparisons <- read.delim("{comparisons_tsv}", check.names=FALSE)
count_cols <- setdiff(colnames(counts), c("gene_id", "gene_symbol"))
count_mat <- as.matrix(counts[, count_cols, drop=FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- counts$gene_id
gene_ids <- counts$gene_id
gene_symbols <- as.character(counts$gene_symbol)
meta$sample_id <- as.character(meta$sample_id)
meta$sex <- factor(as.character(meta$sex), levels=c("M", "F"))
meta$intervention <- factor(as.character(meta$intervention), levels=c("control", "training"))
meta$timepoint_label <- as.character(meta$timepoint_label)

all_results <- list()
for (i in seq_len(nrow(comparisons))) {{
  comp <- comparisons[i, , drop=FALSE]
  tp <- as.character(comp$timepoint_label[1])
  comparison_id <- as.character(comp$comparison_id[1])
  meta_tp <- meta[meta$timepoint_label == tp, , drop=FALSE]
  count_mat_tp <- count_mat[, meta_tp$sample_id, drop=FALSE]
  y <- DGEList(counts=count_mat_tp)
  design_formula <- if (length(unique(as.character(meta_tp$sex))) > 1) ~ intervention + sex else ~ intervention
  design <- model.matrix(design_formula, data=meta_tp)
  keep_genes <- filterByExpr(y, design=design)
  y <- y[keep_genes, , keep.lib.sizes=FALSE]
  gene_ids_tp <- gene_ids[keep_genes]
  gene_symbols_tp <- gene_symbols[keep_genes]
  y <- calcNormFactors(y)
  v <- voom(y, design, plot=FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef="interventiontraining", number=Inf, sort.by="none")
  tt$comparison_id <- comparison_id
  tt$gene_id <- gene_ids_tp
  tt$gene_symbol <- gene_symbols_tp
  tt$group_a <- "training"
  tt$group_b <- "control"
  tt$stratum <- paste("timepoint=", tp, sep="")
  tt$backend <- "r_limma_voom_motrpac_timepoint"
  tt$n_group_a <- sum(meta_tp$intervention == "training")
  tt$n_group_b <- sum(meta_tp$intervention == "control")
  tt$mean_expr <- tt$AveExpr
  tt$model_formula <- if (length(unique(as.character(meta_tp$sex))) > 1) "intervention + sex" else "intervention"
  keep_cols <- c("comparison_id", "gene_id", "gene_symbol", "logFC", "t", "P.Value", "adj.P.Val", "group_a", "group_b", "stratum", "backend", "n_group_a", "n_group_b", "mean_expr", "model_formula")
  tt <- tt[, keep_cols, drop=FALSE]
  colnames(tt)[colnames(tt) == "t"] <- "stat"
  colnames(tt)[colnames(tt) == "P.Value"] <- "pvalue"
  colnames(tt)[colnames(tt) == "adj.P.Val"] <- "padj"
  all_results[[length(all_results) + 1]] <- tt
}}

result <- do.call(rbind, all_results)
write.table(result, file="{output_tsv}", sep="\\t", row.names=FALSE, quote=FALSE)
'''
    _write_text(script_path, script)


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    source_counts_tsv = Path(args.counts_tsv).resolve()

    workflow_outputs: list[tuple[Path, str]] = []
    provenance_inputs: list[tuple[Path, str]]
    if getattr(args, "sample_metadata_tsv", None):
        counts_tsv = source_counts_tsv
        sample_metadata_tsv = Path(args.sample_metadata_tsv).resolve()
        provenance_inputs = [
            (counts_tsv, "counts_tsv"),
            (sample_metadata_tsv, "sample_metadata_tsv"),
        ]
    else:
        required = [
            ("transcript_metadata_tsv", getattr(args, "transcript_metadata_tsv", None)),
            ("phenotype_metadata_tsv", getattr(args, "phenotype_metadata_tsv", None)),
            ("feature_to_gene_tsv", getattr(args, "feature_to_gene_tsv", None)),
            ("rat_to_human_tsv", getattr(args, "rat_to_human_tsv", None)),
            ("tissue_label", getattr(args, "tissue_label", None)),
            ("transcript_tissue_label", getattr(args, "transcript_tissue_label", None)),
        ]
        missing = [name for name, value in required if not value]
        if missing:
            raise ValueError(
                "motrpac_timepoint requires either --sample_metadata_tsv or the raw-input prep arguments: "
                + ", ".join(missing)
            )
        prepared = prepare_tissue_inputs(
            counts_tsv=source_counts_tsv,
            transcript_metadata_tsv=Path(args.transcript_metadata_tsv).resolve(),
            phenotype_metadata_tsv=Path(args.phenotype_metadata_tsv).resolve(),
            feature_to_gene_tsv=Path(args.feature_to_gene_tsv).resolve(),
            rat_to_human_tsv=Path(args.rat_to_human_tsv).resolve(),
            tissue_label=str(args.tissue_label),
            transcript_tissue_label=str(args.transcript_tissue_label),
        )
        prepared_paths = write_prepared_tissue_inputs(out_dir=out_dir, prepared=prepared)
        counts_tsv = prepared_paths["counts_tsv"]
        sample_metadata_tsv = prepared_paths["sample_metadata_tsv"]
        workflow_outputs.extend(
            [
                (prepared_paths["counts_tsv"], "prepared_counts_tsv"),
                (prepared_paths["sample_metadata_tsv"], "prepared_sample_metadata_tsv"),
                (prepared_paths["prepare_summary_json"], "prepare_summary_json"),
                (prepared_paths["prepare_log"], "prepare_log"),
            ]
        )
        provenance_inputs = [
            (source_counts_tsv, "raw_counts_tsv"),
            (Path(args.transcript_metadata_tsv).resolve(), "transcript_metadata_tsv"),
            (Path(args.phenotype_metadata_tsv).resolve(), "phenotype_metadata_tsv"),
            (Path(args.feature_to_gene_tsv).resolve(), "feature_to_gene_tsv"),
            (Path(args.rat_to_human_tsv).resolve(), "rat_to_human_tsv"),
        ]

    metadata_rows = _read_tsv_rows(sample_metadata_tsv)
    tissue_slug = _slugify_tissue_id(str(args.tissue_id))
    comparisons, comparison_summary = _build_comparison_rows(
        metadata_rows=metadata_rows,
        tissue_slug=tissue_slug,
        min_samples_per_group=int(args.min_samples_per_group),
    )
    if not comparisons:
        raise ValueError("No runnable MoTrPAC timepoint comparisons were identified from the prepared sample metadata.")

    comparisons_path = out_dir / "comparisons.tsv"
    comparison_summary_path = out_dir / "comparison_summary.tsv"
    r_script_path = out_dir / "run_motrpac_timepoint_limma_voom.R"
    deg_long_path = out_dir / "deg_long.tsv"
    write_tsv(
        comparisons_path,
        comparisons,
        ["comparison_id", "comparison_kind", "group_column", "group_a", "group_b", "timepoint_label", "tissue_code_no", "tissue_slug"],
    )
    write_tsv(
        comparison_summary_path,
        comparison_summary,
        ["comparison_id", "timepoint_label", "n_control", "n_training"],
    )
    _write_timepoint_r_script(
        script_path=r_script_path,
        counts_tsv=counts_tsv,
        metadata_tsv=sample_metadata_tsv,
        comparisons_tsv=comparisons_path,
        output_tsv=deg_long_path,
    )
    rscript_bin = _resolve_rscript_bin(str(args.rscript_bin))
    completed = subprocess.run(
        [rscript_bin, str(r_script_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"motrpac_timepoint limma/voom workflow failed with exit code {completed.returncode}. "
            f"stderr: {completed.stdout.strip()}"
        )

    write_workflow_provenance_graph(
        workflow_name="motrpac_timepoint",
        module_name="geneset_extractors.workflows.motrpac_timepoint",
        output_dir=out_dir,
        focus_output_path=deg_long_path,
        output_paths=workflow_outputs + [
            (deg_long_path, "deg_tsv"),
            (comparisons_path, "comparisons_tsv"),
            (comparison_summary_path, "comparison_summary"),
            (r_script_path, "workflow_r_script"),
        ],
        input_paths=provenance_inputs,
        parameters={
            "group_column": "intervention",
            "group_a": "training",
            "group_b": "control",
            "stratify_by": ["timepoint_label", "tissue_code_no", "tissue_slug"],
            "covariates": ["sex"],
            "min_samples_per_group": int(args.min_samples_per_group),
            "n_comparisons": len(comparisons),
        },
    )
    return {"deg_long_path": str(deg_long_path), "n_comparisons": len(comparisons), "out_dir": str(out_dir)}
