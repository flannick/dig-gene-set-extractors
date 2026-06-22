library(rliger)
library(Matrix)
library(Seurat)
library(dplyr)
library(purrr)
library(clue)
library(proxy)
library(reticulate)
library(anndata)

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

parse_args <- function(args) {
  if (length(args) < 13) {
    stop(
      paste(
        "Usage: Rscript liger_inmf.R",
        "<input_mode> <input_path> <output_dir> <dataset_column>",
        "<cell_type_column> <max_cells_total> <min_cells_per_cell_type>",
        "<seed> <top_n_genes> <meta_path> <cell_type_label>",
        "<k_grid_csv> <n_reps> [fixed_k] [min_cells_per_dataset] [min_features] [min_umi] [max_mito]"
      ),
      call. = FALSE
    )
  }
  list(
    input_mode = args[[1]],
    input_path = args[[2]],
    output_dir = args[[3]],
    dataset_column = args[[4]],
    cell_type_column = args[[5]],
    max_cells_total = as.integer(args[[6]]),
    min_cells_per_cell_type = as.integer(args[[7]]),
    seed = as.integer(args[[8]]),
    top_n_genes = as.integer(args[[9]]),
    meta_path = args[[10]],
    cell_type_label = args[[11]],
    k_grid = as.integer(strsplit(args[[12]], ",", fixed = TRUE)[[1]]),
    n_reps = as.integer(args[[13]]),
    fixed_k = if (length(args) >= 14 && nzchar(args[[14]])) as.integer(args[[14]]) else NA_integer_,
    min_cells_per_dataset = if (length(args) >= 15 && nzchar(args[[15]])) as.integer(args[[15]]) else 30L,
    min_features = if (length(args) >= 16 && nzchar(args[[16]])) as.integer(args[[16]]) else 200L,
    min_umi = if (length(args) >= 17 && nzchar(args[[17]])) as.numeric(args[[17]]) else 500,
    max_mito = if (length(args) >= 18 && nzchar(args[[18]])) as.numeric(args[[18]]) else 5
  )
}

sanitize_component <- function(value, fallback = "unknown") {
  out <- gsub("[^A-Za-z0-9._=-]+", "_", as.character(value))
  out <- gsub("^_+|_+$", "", out)
  if (!nzchar(out)) {
    return(fallback)
  }
  out
}

load_matrix_mode <- function(input_path, meta_path) {
  expr <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(expr) < 2) {
    stop("matrix_tsv input must contain a cell_id column and at least one gene column")
  }
  meta <- read.delim(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"cell_id" %in% colnames(meta)) {
    stop("matrix_tsv mode requires meta_tsv with a cell_id column")
  }
  cell_ids <- as.character(expr[[1]])
  gene_names <- colnames(expr)[-1]
  mat <- as.matrix(expr[, -1, drop = FALSE])
  rownames(mat) <- cell_ids
  mode(mat) <- "numeric"
  meta <- meta[match(cell_ids, as.character(meta$cell_id)), , drop = FALSE]
  rownames(meta) <- cell_ids
  CreateSeuratObject(counts = t(mat), meta.data = meta)
}

load_h5ad_mode <- function(input_path) {
  adata <- anndata::read_h5ad(input_path)
  raw_x <- adata$X
  if (inherits(raw_x, "Matrix")) {
    counts <- Matrix::t(Matrix::Matrix(raw_x, sparse = TRUE))
    counts <- methods::as(counts, "dgCMatrix")
  } else {
    counts <- t(as.matrix(raw_x))
    counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  }
  meta <- py_to_r(adata$obs)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  if (nrow(meta) == 0) {
    meta <- data.frame(row.names = colnames(counts))
  }
  rownames(meta) <- colnames(counts)
  CreateSeuratObject(counts = counts, meta.data = meta)
}

load_seurat_rds_mode <- function(input_path) {
  obj <- readRDS(input_path)
  if (!inherits(obj, "Seurat")) {
    stop("seurat_rds input did not contain a Seurat object")
  }
  obj
}

load_mtx_mode <- function(input_path, meta_path) {
  counts <- Seurat::Read10X(data.dir = input_path)
  if (is.list(counts)) {
    counts <- counts[[1]]
  }
  meta <- data.frame(row.names = colnames(counts))
  if (nzchar(meta_path)) {
    meta <- read.delim(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!"cell_id" %in% colnames(meta)) {
      stop("mtx_dir mode with meta_tsv requires a cell_id column")
    }
    rownames(meta) <- as.character(meta$cell_id)
    meta <- meta[colnames(counts), , drop = FALSE]
  }
  CreateSeuratObject(counts = counts, meta.data = meta)
}

load_input <- function(cfg) {
  if (cfg$input_mode == "matrix_tsv") {
    return(load_matrix_mode(cfg$input_path, cfg$meta_path))
  }
  if (cfg$input_mode == "h5ad") {
    return(load_h5ad_mode(cfg$input_path))
  }
  if (cfg$input_mode == "seurat_rds") {
    return(load_seurat_rds_mode(cfg$input_path))
  }
  if (cfg$input_mode == "mtx_dir") {
    return(load_mtx_mode(cfg$input_path, cfg$meta_path))
  }
  stop(paste("Unsupported input_mode:", cfg$input_mode))
}

subset_for_cell_type <- function(seurat_obj, cfg, cell_type_value) {
  if (nzchar(cfg$cell_type_column) && cfg$cell_type_column %in% colnames(seurat_obj@meta.data)) {
    subset(seurat_obj, cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data[[cfg$cell_type_column]] == cell_type_value, , drop = FALSE]))
  } else {
    seurat_obj
  }
}

downsample_cells <- function(seurat_obj, max_cells_total, seed) {
  if (ncol(seurat_obj) <= max_cells_total) {
    return(seurat_obj)
  }
  set.seed(seed)
  keep <- sample(colnames(seurat_obj), size = max_cells_total)
  subset(seurat_obj, cells = keep)
}

exclude_small_datasets_liger <- function(counts_list, min_cells = 30) {
  counts_list[sapply(counts_list, ncol) >= min_cells]
}

make_liger_object <- function(seurat_sub, batch, min_cells_per_dataset, min_features, min_umi, max_mito) {
  if (nzchar(batch) && batch %in% colnames(seurat_sub@meta.data)) {
    counts_list <- SplitObject(seurat_sub, split.by = batch) |>
      lapply(function(x) GetAssayData(x, layer = "counts"))
  } else {
    counts_list <- list(all = GetAssayData(seurat_sub, layer = "counts"))
  }
  counts_list <- exclude_small_datasets_liger(counts_list, min_cells = min_cells_per_dataset)
  if (length(counts_list) == 0) {
    stop("No datasets retained after min_cells_per_dataset filtering")
  }
  liger_obj <- createLiger(counts_list)
  liger_obj <- removeMissing(liger_obj, minFeatures = min_features)
  if ("mito" %in% colnames(cellMeta(liger_obj))) {
    liger_obj <- liger_obj[, liger_obj$nUMI > min_umi & liger_obj$mito < max_mito]
  } else {
    liger_obj <- liger_obj[, liger_obj$nUMI > min_umi]
  }
  cell_counts <- sapply(liger_obj@datasets, function(d) ncol(d@rawData))
  keep_datasets <- names(cell_counts[cell_counts >= min_cells_per_dataset])
  liger_obj@datasets <- liger_obj@datasets[keep_datasets]
  liger_obj <- normalize(liger_obj)
  liger_obj <- selectGenes(liger_obj)
  liger_obj <- scaleNotCenter(liger_obj)
  liger_obj
}

run_inmf <- function(liger_obj, k, seed = 1) {
  set.seed(seed)
  runIntegration(
    liger_obj,
    k = k,
    method = "iNMF",
    seed = seed
  )
}

match_factors <- function(W1, W2) {
  sim <- cor(W1, W2)
  sim[is.na(sim)] <- 0
  sim[is.infinite(sim)] <- 0
  sim[sim < 0] <- 0
  sim[sim > 1] <- 1
  cost <- 1 - sim
  cost[cost < 0] <- 0
  assignment <- clue::solve_LSAP(cost)
  matched_sim <- sim[cbind(seq_len(ncol(W1)), assignment)]
  mean(matched_sim)
}

stability_for_k <- function(liger_obj, k, n_reps = 5L) {
  Ws <- vector("list", n_reps)
  for (i in seq_len(n_reps)) {
    tmp <- run_inmf(liger_obj, k, seed = i)
    Ws[[i]] <- tmp@W
  }
  pairwise <- combn(n_reps, 2)
  sims <- apply(pairwise, 2, function(idx) {
    match_factors(Ws[[idx[[1]]]], Ws[[idx[[2]]]])
  })
  mean(sims)
}

find_best_k <- function(liger_obj, k_grid, n_reps) {
  scores <- map_dbl(k_grid, ~stability_for_k(liger_obj, .x, n_reps = n_reps))
  data.frame(k = k_grid, stability = scores)
}

write_liger_outputs <- function(liger_obj, output_dir, cell_type_label, top_n_genes, best_k, stability_df) {
  ct_dir <- file.path(output_dir, sanitize_component(cell_type_label, "all"))
  ensure_dir(ct_dir)

  W <- liger_obj@W
  write.table(W, file = file.path(ct_dir, "gene_loadings.tsv"), sep = "\t", quote = FALSE)

  liger_obj_norm <- quantileNorm(liger_obj)
  Hnorm <- liger_obj_norm@H.norm
  write.table(Hnorm, file = file.path(ct_dir, "cell_scores.tsv"), sep = "\t", quote = FALSE)

  programs <- apply(W, 2, function(w) names(sort(w, decreasing = TRUE))[seq_len(min(length(w), top_n_genes))])
  write.table(
    programs,
    file = file.path(ct_dir, "gene_programs.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = paste0("Factor_", seq_len(ncol(programs)))
  )

  meta <- data.frame(
    cell_type = cell_type_label,
    k = best_k,
    method = "LIGER_iNMF",
    timestamp = as.character(Sys.time())
  )
  write.table(meta, file = file.path(ct_dir, "metadata.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(stability_df, file = file.path(ct_dir, "k_stability.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  factor_importance_norm <- sapply(seq_len(ncol(W)), function(i) {
    sqrt(sum(W[, i]^2)) * sqrt(sum(Hnorm[i, ]^2))
  })
  names(factor_importance_norm) <- colnames(W)
  write.table(
    factor_importance_norm,
    file = file.path(ct_dir, "factor_importance.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  set.seed(cfg$seed)
  ensure_dir(cfg$output_dir)

  seurat_obj <- load_input(cfg)
  meta_cols <- colnames(seurat_obj@meta.data)
  if (!nzchar(cfg$dataset_column) || !(cfg$dataset_column %in% meta_cols)) {
    seurat_obj@meta.data$dataset_auto <- "all"
    cfg$dataset_column <- "dataset_auto"
  }

  cell_types <- if (nzchar(cfg$cell_type_label)) {
    cfg$cell_type_label
  } else if (nzchar(cfg$cell_type_column) && cfg$cell_type_column %in% meta_cols) {
    unique(as.character(seurat_obj@meta.data[[cfg$cell_type_column]]))
  } else {
    "all"
  }

  for (ct in cell_types) {
    seurat_sub <- if (ct == "all" && (!nzchar(cfg$cell_type_column) || !(cfg$cell_type_column %in% meta_cols))) {
      seurat_obj
    } else {
      subset_for_cell_type(seurat_obj, cfg, ct)
    }
    seurat_sub <- downsample_cells(seurat_sub, cfg$max_cells_total, cfg$seed)
    if (ncol(seurat_sub) < cfg$min_cells_per_cell_type) {
      message("Skipping ", ct, ": too few cells after filtering/downsampling")
      next
    }

    liger_obj <- make_liger_object(
      seurat_sub,
      batch = cfg$dataset_column,
      min_cells_per_dataset = cfg$min_cells_per_dataset,
      min_features = cfg$min_features,
      min_umi = cfg$min_umi,
      max_mito = cfg$max_mito
    )

    stability_df <- if (!is.na(cfg$fixed_k)) {
      data.frame(k = cfg$fixed_k, stability = NA_real_)
    } else {
      find_best_k(liger_obj, cfg$k_grid, cfg$n_reps)
    }
    best_k <- if (!is.na(cfg$fixed_k)) cfg$fixed_k else stability_df$k[which.max(stability_df$stability)]
    final <- run_inmf(liger_obj, best_k, seed = cfg$seed)
    write_liger_outputs(final, cfg$output_dir, ct, cfg$top_n_genes, best_k, stability_df)
  }
}

main()
