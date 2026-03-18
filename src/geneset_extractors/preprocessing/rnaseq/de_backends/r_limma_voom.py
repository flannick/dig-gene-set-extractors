from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


REQUIRED_PACKAGES = ("edgeR", "limma")



def backend_available() -> tuple[bool, str]:
    rscript = shutil.which("Rscript") or shutil.which("/opt/homebrew/bin/Rscript")
    if not rscript:
        return False, "Rscript not found on PATH"
    expr = " && ".join([f"requireNamespace('{pkg}', quietly=TRUE)" for pkg in REQUIRED_PACKAGES])
    cmd = [rscript, "--vanilla", "-e", f"quit(status=if ({expr}) 0 else 1)"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except Exception as exc:
        return False, f"Failed to probe R backend availability: {exc}"
    if result.returncode != 0:
        return False, f"Missing required R packages: {', '.join(REQUIRED_PACKAGES)}"
    return True, rscript



def write_script(
    *,
    script_path: str | Path,
    counts_tsv: str | Path,
    metadata_tsv: str | Path,
    comparisons_tsv: str | Path,
    selected_samples_tsv: str | Path,
    output_tsv: str | Path,
    covariates_csv: str,
    batch_columns_csv: str,
    gene_filter_scope: str,
) -> Path:
    path = Path(script_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    script = f'''suppressPackageStartupMessages({{
  library(edgeR)
  library(limma)
}})

counts <- read.delim("{Path(counts_tsv)}", check.names=FALSE)
meta <- read.delim("{Path(metadata_tsv)}", check.names=FALSE)
comps <- read.delim("{Path(comparisons_tsv)}", check.names=FALSE)
selected <- read.delim("{Path(selected_samples_tsv)}", check.names=FALSE)
feature_ids <- counts[[1]]
gene_symbols <- if ("gene_symbol" %in% colnames(counts)) as.character(counts[["gene_symbol"]]) else as.character(feature_ids)
count_cols <- setdiff(colnames(counts), c(colnames(counts)[1], "gene_symbol"))
count_mat <- as.matrix(counts[, count_cols, drop=FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- feature_ids
meta$sample_id <- as.character(meta$sample_id)
selected$sample_id <- as.character(selected$sample_id)
selected$comparison_id <- as.character(selected$comparison_id)
selected$group_label <- as.character(selected$group_label)
count_mat <- count_mat[, meta$sample_id, drop=FALSE]
all_rows <- list()
extra_cols <- unique(c(strsplit("{covariates_csv}", ",", fixed=TRUE)[[1]], strsplit("{batch_columns_csv}", ",", fixed=TRUE)[[1]]))
extra_cols <- extra_cols[nzchar(extra_cols)]
for (i in seq_len(nrow(comps))) {{
  comp <- comps[i, , drop=FALSE]
  gcol <- as.character(comp$group_column[1])
  ga <- as.character(comp$group_a[1])
  gb <- as.character(comp$group_b[1])
  sel <- selected[selected$comparison_id == as.character(comp$comparison_id[1]), , drop=FALSE]
  if (!nrow(sel)) next
  sub_meta <- merge(sel, meta, by="sample_id", all.x=TRUE, sort=FALSE)
  sub_meta <- sub_meta[match(sel$sample_id, sub_meta$sample_id), , drop=FALSE]
  sub_meta$.__group <- factor(as.character(sub_meta$group_label), levels=c(gb, ga))
  sub_meta <- sub_meta[!is.na(sub_meta$.__group), , drop=FALSE]
  if (sum(sub_meta$.__group == ga) < 2 || sum(sub_meta$.__group == gb) < 2) next
  sub_counts <- count_mat[, sub_meta$sample_id, drop=FALSE]
  y <- DGEList(counts=sub_counts)
  if ("{gene_filter_scope}" == "stratum") {{
    strata_cols <- setdiff(colnames(comp), c("comparison_id", "comparison_kind", "group_column", "group_a", "group_b"))
    scope_meta <- meta
    for (col_name in strata_cols) {{
      scope_meta <- scope_meta[as.character(scope_meta[[col_name]]) == as.character(comp[[col_name]][1]), , drop=FALSE]
    }}
    scope_meta <- scope_meta[nzchar(as.character(scope_meta[[gcol]])), , drop=FALSE]
    scope_counts <- count_mat[, scope_meta$sample_id, drop=FALSE]
    keep_genes <- filterByExpr(DGEList(counts=scope_counts))
  }} else {{
    keep_genes <- filterByExpr(y, group=sub_meta$.__group)
  }}
  y <- y[keep_genes, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  formula_terms <- c(".__group", extra_cols[extra_cols %in% colnames(sub_meta)])
  design <- model.matrix(as.formula(paste("~", paste(formula_terms, collapse=" + "))), data=sub_meta)
  v <- voom(y, design, plot=FALSE)
  fit <- lmFit(v, design)
  coef_name <- grep("^.__group", colnames(design), value=TRUE)[1]
  fit <- eBayes(fit)
  tt <- topTable(fit, coef=coef_name, number=Inf, sort.by="none")
  tt$comparison_id <- as.character(comp$comparison_id[1])
  tt$gene_id <- rownames(tt)
  tt$gene_symbol <- gene_symbols[match(rownames(tt), feature_ids)]
  tt$group_a <- ga
  tt$group_b <- gb
  tt$stratum <- paste(paste(setdiff(colnames(comp), c("comparison_id", "comparison_kind", "group_column", "group_a", "group_b")), as.character(comp[1, setdiff(colnames(comp), c("comparison_id", "comparison_kind", "group_column", "group_a", "group_b"))]), sep="="), collapse="|")
  tt$backend <- "r_limma_voom"
  tt$n_group_a <- sum(sub_meta$.__group == ga)
  tt$n_group_b <- sum(sub_meta$.__group == gb)
  tt$mean_expr <- tt$AveExpr
  tt$model_formula <- paste(formula_terms, collapse=" + ")
  keep_cols <- c("comparison_id", "gene_id", "gene_symbol", "logFC", "t", "P.Value", "adj.P.Val", "group_a", "group_b", "stratum", "backend", "n_group_a", "n_group_b", "mean_expr", "model_formula")
  tt <- tt[, keep_cols, drop=FALSE]
  colnames(tt)[colnames(tt) == "t"] <- "stat"
  colnames(tt)[colnames(tt) == "P.Value"] <- "pvalue"
  colnames(tt)[colnames(tt) == "adj.P.Val"] <- "padj"
  all_rows[[length(all_rows) + 1]] <- tt
}}
if (!length(all_rows)) stop("No contrasts produced output")
out <- do.call(rbind, all_rows)
write.table(out, file="{Path(output_tsv)}", sep="\t", row.names=FALSE, quote=FALSE)
'''
    path.write_text(script, encoding="utf-8")
    return path
