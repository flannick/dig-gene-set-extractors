from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path
from typing import Any

from geneset_extractors.workflows.gtex_runtime_common import (
    build_sample_metadata_rows,
    parse_gct_header,
    read_tsv,
    write_filtered_counts,
    write_text,
    write_tsv,
    write_workflow_provenance_graph,
)


def _parse_age_midpoint(age_bin: str) -> float:
    left, right = [part.strip() for part in str(age_bin).split("-", 1)]
    return (float(left) + float(right)) / 2.0


def _resolve_rscript_bin(rscript_bin: str) -> str:
    resolved = shutil.which(rscript_bin) if not Path(rscript_bin).is_absolute() else rscript_bin
    if not resolved or not Path(resolved).exists():
        raise ValueError(f"Rscript not found: {rscript_bin}")
    return resolved


def _write_continuous_age_r_script(
    *,
    script_path: Path,
    counts_tsv: Path,
    metadata_tsv: Path,
    output_tsv: Path,
    include_sex: bool,
) -> None:
    formula_terms = ["age_mid"]
    if include_sex:
        formula_terms.append("SEX")
    formula_expr = " + ".join(formula_terms)
    script = f'''suppressPackageStartupMessages({{
  library(edgeR)
  library(limma)
}})

counts <- read.delim("{counts_tsv}", check.names=FALSE)
meta <- read.delim("{metadata_tsv}", check.names=FALSE)
feature_ids <- counts[[1]]
gene_symbols <- if ("gene_symbol" %in% colnames(counts)) as.character(counts[["gene_symbol"]]) else as.character(feature_ids)
count_cols <- setdiff(colnames(counts), c(colnames(counts)[1], "gene_symbol"))
count_mat <- as.matrix(counts[, count_cols, drop=FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- feature_ids
meta$sample_id <- as.character(meta$sample_id)
meta$SEX <- factor(as.character(meta$SEX))
meta$age_mid <- as.numeric(meta$age_mid)
count_mat <- count_mat[, meta$sample_id, drop=FALSE]
y <- DGEList(counts=count_mat)
keep_genes <- filterByExpr(y)
y <- y[keep_genes, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(as.formula("~ {formula_expr}"), data=meta)
v <- voom(y, design, plot=FALSE)
fit <- lmFit(v, design)
coef_name <- "age_mid"
fit <- eBayes(fit)
tt <- topTable(fit, coef=coef_name, number=Inf, sort.by="none")
tt$comparison_id <- "continuous_age"
tt$gene_id <- rownames(tt)
tt$gene_symbol <- gene_symbols[match(rownames(tt), feature_ids)]
tt$group_a <- "older"
tt$group_b <- "younger"
tt$stratum <- ""
tt$backend <- "r_limma_voom_continuous_age"
tt$n_group_a <- nrow(meta)
tt$n_group_b <- nrow(meta)
tt$mean_expr <- tt$AveExpr
tt$model_formula <- "{formula_expr}"
keep_cols <- c("comparison_id", "gene_id", "gene_symbol", "logFC", "t", "P.Value", "adj.P.Val", "group_a", "group_b", "stratum", "backend", "n_group_a", "n_group_b", "mean_expr", "model_formula")
tt <- tt[, keep_cols, drop=FALSE]
colnames(tt)[colnames(tt) == "t"] <- "stat"
colnames(tt)[colnames(tt) == "P.Value"] <- "pvalue"
colnames(tt)[colnames(tt) == "adj.P.Val"] <- "padj"
write.table(tt, file="{output_tsv}", sep="\\t", row.names=FALSE, quote=FALSE)
'''
    write_text(script_path, script)


def run(args) -> dict[str, object]:
    expression_gct = Path(args.expression_gct).resolve()
    sample_attributes_tsv = Path(args.sample_attributes_tsv).resolve()
    subject_phenotypes_tsv = Path(args.subject_phenotypes_tsv).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_rows = read_tsv(sample_attributes_tsv)
    subject_rows = read_tsv(subject_phenotypes_tsv)
    _header, gct_sample_ids = parse_gct_header(expression_gct)
    tissue_label = str(args.tissue_label).strip()
    tissue_id = str(args.tissue_id).strip()
    tissue_column = str(getattr(args, "tissue_column", "") or "").strip() or None
    tissue_value = str(getattr(args, "tissue_value", "") or "").strip() or None

    prepared_meta = build_sample_metadata_rows(
        sample_rows=sample_rows,
        subject_rows=subject_rows,
        gct_sample_ids=gct_sample_ids,
        tissue_label=tissue_label,
        tissue_id=tissue_id,
        tissue_column=tissue_column,
        tissue_value=tissue_value,
    )
    retained_sample_ids = [row["sample_id"] for row in prepared_meta]
    retained_sample_id_set = set(retained_sample_ids)
    sample_index = [idx for idx, sample_id in enumerate(gct_sample_ids) if sample_id in retained_sample_id_set]
    sample_columns = [gct_sample_ids[idx] for idx in sample_index]

    counts_path = out_dir / "tissue_counts.tsv"
    sample_metadata_path = out_dir / "sample_metadata.tsv"
    continuous_metadata_path = out_dir / "continuous_sample_metadata.tsv"
    r_script_path = out_dir / "run_continuous_age_limma_voom.R"
    deg_long_path = out_dir / "deg_long.tsv"

    n_genes_retained = write_filtered_counts(
        counts_gct=expression_gct,
        sample_columns=sample_columns,
        sample_index=sample_index,
        out_path=counts_path,
    )
    write_tsv(
        sample_metadata_path,
        prepared_meta,
        ["sample_id", "subject_id", "age_bin", "SEX", "primary_tissue", "detailed_tissue", "tissue_id", "tissue_label"],
    )

    continuous_rows: list[dict[str, str]] = []
    age_values: list[float] = []
    for row in prepared_meta:
        age_mid = _parse_age_midpoint(row["age_bin"])
        age_values.append(age_mid)
        continuous_rows.append(
            {
                "sample_id": row["sample_id"],
                "subject_id": row["subject_id"],
                "age_bin": row["age_bin"],
                "age_mid": f"{age_mid:.1f}",
                "SEX": row["SEX"],
                "primary_tissue": row["primary_tissue"],
                "detailed_tissue": row["detailed_tissue"],
                "tissue_id": row["tissue_id"],
                "tissue_label": row["tissue_label"],
            }
        )
    write_tsv(
        continuous_metadata_path,
        continuous_rows,
        ["sample_id", "subject_id", "age_bin", "age_mid", "SEX", "primary_tissue", "detailed_tissue", "tissue_id", "tissue_label"],
    )

    include_sex = str(getattr(args, "covariates", "") or "").strip().lower() != "none" and bool(str(getattr(args, "covariates", "") or "").strip())
    _write_continuous_age_r_script(
        script_path=r_script_path,
        counts_tsv=counts_path,
        metadata_tsv=continuous_metadata_path,
        output_tsv=deg_long_path,
        include_sex=include_sex,
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
            f"gtex_continuous_age limma/voom workflow failed with exit code {completed.returncode}. "
            f"stderr: {completed.stdout.strip()}"
        )

    write_workflow_provenance_graph(
        workflow_name="gtex_continuous_age",
        module_name="geneset_extractors.workflows.gtex_continuous_age",
        output_dir=out_dir,
        focus_output_path=deg_long_path,
        output_paths=[
            (deg_long_path, "deg_tsv"),
            (counts_path, "counts_tsv"),
            (sample_metadata_path, "sample_metadata_tsv"),
            (continuous_metadata_path, "continuous_sample_metadata_tsv"),
            (r_script_path, "workflow_r_script"),
        ],
        input_paths=[
            (expression_gct, "expression_gct"),
            (sample_attributes_tsv, "sample_attributes_tsv"),
            (subject_phenotypes_tsv, "subject_phenotypes_tsv"),
        ],
        parameters={
            "tissue_id": tissue_id,
            "tissue_label": tissue_label,
            "tissue_column": tissue_column,
            "tissue_value": tissue_value,
            "covariates": str(getattr(args, "covariates", "") or ""),
            "n_samples_retained": len(prepared_meta),
            "n_genes_retained": n_genes_retained,
            "age_mid_min": min(age_values) if age_values else None,
            "age_mid_max": max(age_values) if age_values else None,
        },
    )
    return {
        "deg_long_path": str(deg_long_path),
        "n_samples": len(prepared_meta),
        "n_genes_retained": n_genes_retained,
    }
