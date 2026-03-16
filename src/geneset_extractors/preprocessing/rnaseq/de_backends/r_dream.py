from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


REQUIRED_PACKAGES = ("edgeR", "limma", "variancePartition")



def backend_available() -> tuple[bool, str]:
    rscript = shutil.which("Rscript") or shutil.which("/opt/homebrew/bin/Rscript")
    if not rscript:
        return False, "Rscript not found on PATH"
    expr = " && ".join([f"requireNamespace('{pkg}', quietly=TRUE)" for pkg in REQUIRED_PACKAGES])
    cmd = [rscript, "--vanilla", "-e", f"quit(status=if ({expr}) 0 else 1)"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except Exception as exc:
        return False, f"Failed to probe dream backend availability: {exc}"
    if result.returncode != 0:
        return False, f"Missing required R packages: {', '.join(REQUIRED_PACKAGES)}"
    return True, rscript



def write_script(
    *,
    script_path: str | Path,
    counts_tsv: str | Path,
    metadata_tsv: str | Path,
    comparisons_tsv: str | Path,
    output_tsv: str | Path,
    covariates_csv: str,
    batch_columns_csv: str,
    random_effect_column: str,
) -> Path:
    path = Path(script_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    script = f'''suppressPackageStartupMessages({{
  library(edgeR)
  library(limma)
  library(variancePartition)
}})

counts <- read.delim("{Path(counts_tsv)}", check.names=FALSE)
meta <- read.delim("{Path(metadata_tsv)}", check.names=FALSE)
comps <- read.delim("{Path(comparisons_tsv)}", check.names=FALSE)
feature_ids <- counts[[1]]
count_mat <- as.matrix(counts[, -1, drop=FALSE])
rownames(count_mat) <- feature_ids
meta$sample_id <- as.character(meta$sample_id)
count_mat <- count_mat[, meta$sample_id, drop=FALSE]
all_rows <- list()
extra_cols <- unique(c(strsplit("{covariates_csv}", ",", fixed=TRUE)[[1]], strsplit("{batch_columns_csv}", ",", fixed=TRUE)[[1]]))
extra_cols <- extra_cols[nzchar(extra_cols)]
for (i in seq_len(nrow(comps))) {{
  comp <- comps[i, , drop=FALSE]
  keep <- rep(TRUE, nrow(meta))
  for (col in setdiff(colnames(comp), c("comparison_id", "comparison_kind", "group_column", "group_a", "group_b"))) {{
    val <- as.character(comp[[col]][1])
    if (!nzchar(val)) next
    keep <- keep & as.character(meta[[col]]) == val
  }}
  sub_meta <- meta[keep, , drop=FALSE]
  gcol <- as.character(comp$group_column[1])
  ga <- as.character(comp$group_a[1])
  gb <- as.character(comp$group_b[1])
  if (!nzchar(gb) || as.character(comp$comparison_kind[1]) == "group_vs_rest") {{
    sub_meta <- sub_meta[sub_meta[[gcol]] %in% c(ga, setdiff(unique(as.character(sub_meta[[gcol]])), ga)), , drop=FALSE]
    sub_meta$.__group <- ifelse(as.character(sub_meta[[gcol]]) == ga, ga, "rest")
    gb <- "rest"
  }} else {{
    sub_meta <- sub_meta[sub_meta[[gcol]] %in% c(ga, gb), , drop=FALSE]
    sub_meta$.__group <- as.character(sub_meta[[gcol]])
  }}
  if (sum(sub_meta$.__group == ga) < 2 || sum(sub_meta$.__group == gb) < 2) next
  sub_counts <- count_mat[, sub_meta$sample_id, drop=FALSE]
  y <- DGEList(counts=sub_counts)
  keep_genes <- filterByExpr(y, group=sub_meta$.__group)
  y <- y[keep_genes, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  form_terms <- c(".__group", extra_cols[extra_cols %in% colnames(sub_meta)])
  form <- as.formula(paste("~", paste(form_terms, collapse=" + "), "+ (1|{random_effect_column})"))
  vobj <- voomWithDreamWeights(y, form, sub_meta)
  fit <- dream(vobj, form, sub_meta)
  tt <- topTable(fit, coef=grep("^.__group", colnames(coef(fit)), value=TRUE)[1], number=Inf, sort.by="none")
  tt$comparison_id <- as.character(comp$comparison_id[1])
  tt$gene_id <- rownames(tt)
  tt$gene_symbol <- rownames(tt)
  tt$group_a <- ga
  tt$group_b <- gb
  tt$stratum <- paste(paste(setdiff(colnames(comp), c("comparison_id", "comparison_kind", "group_column", "group_a", "group_b")), as.character(comp[1, setdiff(colnames(comp), c("comparison_id", "comparison_kind", "group_column", "group_a", "group_b"))]), sep="="), collapse="|")
  tt$backend <- "r_dream"
  tt$n_group_a <- sum(sub_meta$.__group == ga)
  tt$n_group_b <- sum(sub_meta$.__group == gb)
  tt$mean_expr <- tt$AveExpr
  tt$model_formula <- paste(form_terms, collapse=" + ")
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
