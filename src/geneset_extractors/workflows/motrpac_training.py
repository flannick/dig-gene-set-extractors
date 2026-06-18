from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from geneset_extractors.workflows.gtex_runtime_common import write_workflow_provenance_graph
from geneset_extractors.workflows.motrpac_common import prepare_tissue_inputs, write_prepared_tissue_inputs


def _resolve_rscript_bin(rscript_bin: str) -> str:
    resolved = shutil.which(rscript_bin) if not Path(rscript_bin).is_absolute() else rscript_bin
    if not resolved or not Path(resolved).exists():
        raise ValueError(f"Rscript not found: {rscript_bin}")
    return resolved


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_training_r_script(
    *,
    script_path: Path,
    counts_tsv: Path,
    metadata_tsv: Path,
    output_tsv: Path,
    include_sex: bool,
) -> None:
    formula_expr = "~ intervention + sex" if include_sex else "~ intervention"
    model_formula_label = "intervention + sex" if include_sex else "intervention"
    script = f'''suppressPackageStartupMessages({{
  library(edgeR)
  library(limma)
}})

counts <- read.delim("{counts_tsv}", check.names=FALSE)
meta <- read.delim("{metadata_tsv}", check.names=FALSE)
count_cols <- setdiff(colnames(counts), c("gene_id", "gene_symbol"))
count_mat <- as.matrix(counts[, count_cols, drop=FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- counts$gene_id
gene_ids <- counts$gene_id
gene_symbols <- as.character(counts$gene_symbol)

meta$sample_id <- as.character(meta$sample_id)
meta$sex <- factor(as.character(meta$sex), levels=c("M", "F"))
meta$intervention <- factor(as.character(meta$intervention), levels=c("control", "training"))
count_mat <- count_mat[, meta$sample_id, drop=FALSE]
y <- DGEList(counts=count_mat)
design_formula <- if ({str(include_sex).upper()} && length(unique(as.character(meta$sex))) > 1) ~ intervention + sex else ~ intervention
design <- model.matrix(design_formula, data=meta)
keep_genes <- filterByExpr(y, design=design)
y <- y[keep_genes, , keep.lib.sizes=FALSE]
gene_ids <- gene_ids[keep_genes]
gene_symbols <- gene_symbols[keep_genes]
y <- calcNormFactors(y)
v <- voom(y, design, plot=FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef="interventiontraining", number=Inf, sort.by="none")
tt$comparison_id <- "training_vs_control"
tt$gene_id <- gene_ids
tt$gene_symbol <- gene_symbols
tt$group_a <- "training"
tt$group_b <- "control"
tt$stratum <- ""
tt$backend <- "r_limma_voom_motrpac_training"
tt$n_group_a <- sum(meta$intervention == "training")
tt$n_group_b <- sum(meta$intervention == "control")
tt$mean_expr <- tt$AveExpr
tt$model_formula <- if ({str(include_sex).upper()} && length(unique(as.character(meta$sex))) > 1) "intervention + sex" else "intervention"
keep_cols <- c("comparison_id", "gene_id", "gene_symbol", "logFC", "t", "P.Value", "adj.P.Val", "group_a", "group_b", "stratum", "backend", "n_group_a", "n_group_b", "mean_expr", "model_formula")
tt <- tt[, keep_cols, drop=FALSE]
colnames(tt)[colnames(tt) == "t"] <- "stat"
colnames(tt)[colnames(tt) == "P.Value"] <- "pvalue"
colnames(tt)[colnames(tt) == "adj.P.Val"] <- "padj"
write.table(tt, file="{output_tsv}", sep="\\t", row.names=FALSE, quote=FALSE)
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
                "motrpac_training requires either --sample_metadata_tsv or the raw-input prep arguments: "
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

    include_sex = str(getattr(args, "covariates", "sex") or "sex").strip().lower() != "none"
    r_script_path = out_dir / "run_motrpac_training_limma_voom.R"
    deg_path = out_dir / "training_deg.tsv"
    _write_training_r_script(
        script_path=r_script_path,
        counts_tsv=counts_tsv,
        metadata_tsv=sample_metadata_tsv,
        output_tsv=deg_path,
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
            f"motrpac_training limma/voom workflow failed with exit code {completed.returncode}. "
            f"stderr: {completed.stdout.strip()}"
        )

    write_workflow_provenance_graph(
        workflow_name="motrpac_training",
        module_name="geneset_extractors.workflows.motrpac_training",
        output_dir=out_dir,
        focus_output_path=deg_path,
        output_paths=workflow_outputs + [
            (deg_path, "deg_tsv"),
            (r_script_path, "workflow_r_script"),
        ],
        input_paths=provenance_inputs,
        parameters={
            "covariates": "sex" if include_sex else "none",
            "model_formula": "intervention + sex" if include_sex else "intervention",
        },
    )
    return {"deg_tsv_path": str(deg_path), "out_dir": str(out_dir)}
