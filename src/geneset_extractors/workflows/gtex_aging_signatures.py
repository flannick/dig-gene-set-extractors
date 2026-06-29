from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from geneset_extractors.core.metadata import _resolve_git_commit, current_invocation_context
from geneset_extractors.core.provenance import (
    REPO_URL,
    build_analysis_node,
    build_edges,
    build_file_node,
    build_output_file_record,
    get_runtime_context,
    mirror_graph_payload,
    stable_operation_id,
    write_canonical_json,
)


AGE_ORDER = ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"]
_WS_RE = re.compile(r"[^A-Za-z0-9]+")


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with _open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    _write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _choose_column(fieldnames: Iterable[str], candidates: list[str], label: str) -> str:
    present = list(fieldnames)
    for candidate in candidates:
        if candidate in present:
            return candidate
    raise ValueError(
        f"Could not resolve {label} column. Tried: {', '.join(candidates)}. "
        f"Available columns begin: {', '.join(present[:30])}"
    )


def _sample_to_subject_id(sample_id: str) -> str:
    parts = str(sample_id).split("-")
    if len(parts) >= 2:
        return "-".join(parts[:2])
    return str(sample_id)


def _parse_age_group(value: object) -> str | None:
    text = "" if value is None else str(value).strip()
    if text in AGE_ORDER:
        return text
    match = re.search(r"(\d+)", text)
    if not match:
        return None
    age = int(match.group(1))
    if 20 <= age <= 29:
        return "20-29"
    if 30 <= age <= 39:
        return "30-39"
    if 40 <= age <= 49:
        return "40-49"
    if 50 <= age <= 59:
        return "50-59"
    if 60 <= age <= 69:
        return "60-69"
    if 70 <= age <= 79:
        return "70-79"
    return None


def _sanitize_tissue_label(value: str) -> str:
    parts = [part for part in _WS_RE.sub(" ", str(value).strip()).split() if part]
    return "".join(parts)


def _strip_ensembl_version(gene_id: str) -> str:
    return str(gene_id).strip().split(".", 1)[0]


def _sample_variance(values: list[float]) -> float:
    if len(values) <= 1:
        return 0.0
    mean_value = sum(values) / float(len(values))
    return sum((value - mean_value) ** 2 for value in values) / float(len(values) - 1)


def _balanced_sample_ids(sample_ids: list[str], n: int, random_state: int) -> list[str]:
    if n <= 0:
        return []
    return (
        pd.Series(list(sample_ids))
        .sample(n=n, random_state=random_state, replace=False)
        .tolist()
    )


def _read_human_gene_info_mapping(path: Path) -> dict[str, str]:
    rows = _read_tsv_rows(path)
    if not rows:
        raise ValueError(f"No rows found in human_gene_info file: {path}")
    fieldnames = list(rows[0].keys())
    if "#tax_id" not in fieldnames or "Symbol" not in fieldnames or "dbXrefs" not in fieldnames:
        raise ValueError(
            "human_gene_info is missing required columns. "
            "Expected #tax_id, Symbol, and dbXrefs."
        )
    ensembl_re = re.compile(r"Ensembl:([^|]+)")
    mapping: dict[str, str] = {}
    for row in rows:
        if str(row.get("#tax_id", "")).strip() != "9606":
            continue
        symbol = str(row.get("Symbol", "")).strip()
        if not symbol or symbol == "-":
            continue
        match = ensembl_re.search(str(row.get("dbXrefs", "")))
        if not match:
            continue
        ensembl = _strip_ensembl_version(match.group(1).split("|", 1)[0])
        if ensembl:
            mapping.setdefault(ensembl, symbol)
    if not mapping:
        raise ValueError(f"No Ensembl-to-symbol mappings parsed from {path}")
    return mapping


def _build_sample_metadata(
    sample_attr_rows: list[dict[str, str]],
    subject_rows: list[dict[str, str]],
    *,
    tissue_column: str,
    tissue_value: str,
    tissue_label: str,
    tissue_id: str,
) -> list[dict[str, str]]:
    if not sample_attr_rows:
        raise ValueError("Sample attributes table is empty")
    if not subject_rows:
        raise ValueError("Subject phenotypes table is empty")
    sample_fieldnames = list(sample_attr_rows[0].keys())
    subject_fieldnames = list(subject_rows[0].keys())
    sample_id_col = _choose_column(sample_fieldnames, ["SAMPID", "sample_id", "SampleID"], "sample id")
    smts_col = _choose_column(sample_fieldnames, ["SMTS", "smts", "tissue", "TISSUE"], "broad tissue")
    smtsd_col = "SMTSD" if "SMTSD" in sample_fieldnames else None
    subject_id_col = _choose_column(subject_fieldnames, ["SUBJID", "subject_id", "SubjectID"], "subject id")
    age_col = _choose_column(subject_fieldnames, ["AGE", "age", "Age"], "age")
    sex_col = next((candidate for candidate in ["SEX", "sex", "Sex"] if candidate in subject_fieldnames), None)

    subject_by_id: dict[str, dict[str, str]] = {}
    for row in subject_rows:
        subject_id = str(row.get(subject_id_col, "")).strip()
        if subject_id:
            subject_by_id[subject_id] = row

    out_rows: list[dict[str, str]] = []
    tissue_compact = _sanitize_tissue_label(tissue_label)
    for row in sample_attr_rows:
        sample_id = str(row.get(sample_id_col, "")).strip()
        if not sample_id:
            continue
        broad_tissue = str(row.get(smts_col, "")).strip()
        if tissue_column == smts_col and broad_tissue != tissue_value:
            continue
        if tissue_column in row and str(row.get(tissue_column, "")).strip() != tissue_value:
            continue
        subject_id = _sample_to_subject_id(sample_id)
        subject_row = subject_by_id.get(subject_id)
        if subject_row is None:
            continue
        age_group = _parse_age_group(subject_row.get(age_col))
        if age_group is None:
            continue
        out_rows.append(
            {
                "sample_id": sample_id,
                "subject_id": subject_id,
                "tissue_id": tissue_id,
                "tissue_label": tissue_label,
                "tissue_compact": tissue_compact,
                "SMTS": broad_tissue,
                "SMTSD": str(row.get(smtsd_col, "")).strip() if smtsd_col else broad_tissue,
                "sex": str(subject_row.get(sex_col, "")).strip() if sex_col else "",
                "age_group": age_group,
            }
        )
    if not out_rows:
        raise ValueError(
            f"No samples remained after filtering sample metadata for {tissue_column}={tissue_value!r}"
        )
    return out_rows


def _gct_sample_columns(path: Path) -> list[str]:
    with _open_text(path) as handle:
        _ = handle.readline()
        _ = handle.readline()
        header = handle.readline().rstrip("\n").split("\t")
    return header[2:]


def _prepare_counts_matrix(
    *,
    expression_gct: Path,
    sample_ids: list[str],
    ensembl_to_symbol: dict[str, str],
    out_path: Path,
) -> dict[str, object]:
    requested = set(sample_ids)
    header = _gct_sample_columns(expression_gct)
    selected_sample_ids = [sample_id for sample_id in header if sample_id in requested]
    if not selected_sample_ids:
        raise ValueError("No selected samples were present in the expression GCT")

    usecols = ["Name", "Description", *selected_sample_ids]
    chunks: list[pd.DataFrame] = []
    rows_seen = 0

    reader = pd.read_csv(
        expression_gct,
        sep="\t",
        skiprows=2,
        usecols=usecols,
        chunksize=1000,
        low_memory=False,
    )

    for chunk in reader:
        rows_seen += int(chunk.shape[0])
        work = chunk[["Name", *selected_sample_ids]].copy()
        work["Ens"] = work["Name"].map(_strip_ensembl_version)
        work = work[work["Ens"].isin(ensembl_to_symbol)].copy()
        if work.empty:
            continue

        counts_only = work.loc[:, selected_sample_ids].apply(pd.to_numeric, errors="coerce").fillna(0)
        work["_Var"] = counts_only.var(axis=1)
        keep_idx = (
            work[["Ens", "_Var"]]
            .sort_values(by=["Ens", "_Var"], ascending=True)
            .drop_duplicates(subset=["Ens"], keep="last")
            .index
        )

        kept = work.loc[keep_idx, ["Ens", "_Var"]].copy()
        counts = counts_only.loc[keep_idx, :].copy()
        counts.insert(0, "_Var", kept["_Var"].to_numpy())
        counts.insert(0, "Ens", kept["Ens"].to_numpy())
        chunks.append(counts)

    if not chunks:
        raise ValueError("No mapped Ensembl genes were retained from the expression GCT")

    work_all = pd.concat(chunks, axis=0, ignore_index=True)
    keep_idx = (
        work_all[["Ens", "_Var"]]
        .sort_values(by=["Ens", "_Var"], ascending=True)
        .drop_duplicates(subset=["Ens"], keep="last")
        .index
    )
    work_all = work_all.loc[keep_idx, :].copy()

    counts_with_meta = work_all.loc[:, ["Ens", "_Var", *selected_sample_ids]].copy()
    counts_with_meta["gene_symbol"] = counts_with_meta["Ens"].map(
        lambda value: ensembl_to_symbol.get(value, value)
    ).astype(str)
    keep_symbol_idx = (
        counts_with_meta[["gene_symbol", "_Var"]]
        .sort_values(by=["gene_symbol", "_Var"], ascending=True)
        .drop_duplicates(subset=["gene_symbol"], keep="last")
        .index
    )
    counts_with_meta = counts_with_meta.loc[keep_symbol_idx, :].copy()

    count_rows: list[dict[str, Any]] = []
    for _, row in counts_with_meta.iterrows():
        out_row: dict[str, Any] = {
            "gene_id": str(row["Ens"]),
            "gene_symbol": str(row["gene_symbol"]),
        }
        for sample_id in selected_sample_ids:
            out_row[sample_id] = row[sample_id]
        count_rows.append(out_row)
    _write_tsv(out_path, count_rows, ["gene_id", "gene_symbol", *selected_sample_ids])
    return {
        "selected_sample_ids": selected_sample_ids,
        "n_selected_samples": len(selected_sample_ids),
        "n_requested_samples": len(sample_ids),
        "n_rows_after_ensembl_dedup": int(work_all["Ens"].nunique()),
        "n_rows_after_symbol_dedup": len(count_rows),
        "n_rows_seen": rows_seen,
    }


def _apply_tissue_filter(*, counts_tsv: Path, sample_metadata_tsv: Path, output_tsv: Path, rscript_bin: str) -> dict[str, int]:
    script_path = output_tsv.parent / "run_filter_by_expr.R"
    script = f'''suppressPackageStartupMessages({{
  library(edgeR)
}})

counts <- read.delim("{counts_tsv}", check.names=FALSE)
meta <- read.delim("{sample_metadata_tsv}", check.names=FALSE)
count_cols <- setdiff(colnames(counts), c("gene_id", "gene_symbol"))
count_mat <- as.matrix(counts[, count_cols, drop=FALSE])
storage.mode(count_mat) <- "numeric"
group <- factor(as.character(meta$age_group))
y <- DGEList(counts=count_mat)
keep <- filterByExpr(y, group=group)
filtered <- counts[keep, , drop=FALSE]
write.table(filtered, file="{output_tsv}", sep="\\t", row.names=FALSE, quote=FALSE)
cat(paste(sum(keep), nrow(counts), sep="\\t"))
'''
    _write_text(script_path, script)
    result = subprocess.run(
        [rscript_bin, "--vanilla", str(script_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise ValueError(
            f"filterByExpr tissue prefilter failed with exit code {result.returncode}. "
            f"stderr: {result.stderr.strip() or '[empty]'}"
        )
    parts = (result.stdout or "").strip().split("\t")
    if len(parts) != 2:
        return {"n_genes_after_filter": 0, "n_genes_before_filter": 0}
    return {"n_genes_after_filter": int(parts[0]), "n_genes_before_filter": int(parts[1])}


def _comparison_id(age_group: str, reference_age_group: str) -> str:
    return f"age{age_group.split('-', 1)[0]}_{reference_age_group.split('-', 1)[0]}"


def _build_comparisons(
    metadata_rows: list[dict[str, str]],
    *,
    tissue_label: str,
    reference_age_group: str,
    comparison_age_groups: list[str],
    random_state: int,
    min_samples_per_group: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    age_to_samples: dict[str, list[str]] = {age_group: [] for age_group in AGE_ORDER}
    for row in metadata_rows:
        age_group = str(row.get("age_group", "")).strip()
        sample_id = str(row.get("sample_id", "")).strip()
        if age_group in age_to_samples and sample_id:
            age_to_samples[age_group].append(sample_id)

    comparisons: list[dict[str, Any]] = []
    selected_samples: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    tissue_compact = _sanitize_tissue_label(tissue_label)
    for age_group in comparison_age_groups:
        if age_group == reference_age_group:
            continue
        reference_all = age_to_samples.get(reference_age_group, [])
        case_all = age_to_samples.get(age_group, [])
        n_balanced = min(len(reference_all), len(case_all))
        comparison_id = _comparison_id(age_group, reference_age_group)
        aging_signature = f"GTEx_aging_{tissue_compact}_{reference_age_group}_{age_group}"
        status = "ok"
        if n_balanced < min_samples_per_group:
            status = f"skipped: n per group < {min_samples_per_group}"
        comparison_row = {
            "comparison_id": comparison_id,
            "comparison_kind": "reference_level",
            "aging_signature": aging_signature,
            "group_column": "age_group",
            "group_a": age_group,
            "group_b": reference_age_group,
            "tissue_label": tissue_label,
        }
        comparisons.append(comparison_row)
        if status != "ok":
            audit_rows.append(
                {
                    "comparison_id": comparison_id,
                    "aging_signature": aging_signature,
                    "tissue_label": tissue_label,
                    "group_a": age_group,
                    "group_b": reference_age_group,
                    "n_group_a_available": len(case_all),
                    "n_group_b_available": len(reference_all),
                    "n_group_a_selected": 0,
                    "n_group_b_selected": 0,
                    "status": status,
                }
            )
            continue
        reference_selected = _balanced_sample_ids(reference_all, n_balanced, random_state)
        case_selected = _balanced_sample_ids(case_all, n_balanced, random_state)
        for rank, sample_id in enumerate(reference_selected, start=1):
            selected_samples.append(
                {
                    "comparison_id": comparison_id,
                    "sample_id": sample_id,
                    "group_label": reference_age_group,
                    "group_value": reference_age_group,
                    "selected_rank": rank,
                    "balance_requested": "true",
                    "balance_applied": "true",
                    "balance_seed": random_state,
                    "balance_seed_effective": random_state,
                    "age_group": reference_age_group,
                }
            )
        for rank, sample_id in enumerate(case_selected, start=1):
            selected_samples.append(
                {
                    "comparison_id": comparison_id,
                    "sample_id": sample_id,
                    "group_label": age_group,
                    "group_value": age_group,
                    "selected_rank": rank,
                    "balance_requested": "true",
                    "balance_applied": "true",
                    "balance_seed": random_state,
                    "balance_seed_effective": random_state,
                    "age_group": age_group,
                }
            )
        audit_rows.append(
            {
                "comparison_id": comparison_id,
                "aging_signature": aging_signature,
                "tissue_label": tissue_label,
                "group_a": age_group,
                "group_b": reference_age_group,
                "n_group_a_available": len(case_all),
                "n_group_b_available": len(reference_all),
                "n_group_a_selected": len(case_selected),
                "n_group_b_selected": len(reference_selected),
                "status": status,
            }
        )

    emitted_comparisons = [row for row in comparisons if any(sel["comparison_id"] == row["comparison_id"] for sel in selected_samples)]
    return emitted_comparisons, selected_samples, audit_rows


def _write_limma_script(
    *,
    script_path: Path,
    counts_tsv: Path,
    sample_metadata_tsv: Path,
    comparisons_tsv: Path,
    selected_samples_tsv: Path,
    deg_long_tsv: Path,
    comparisons_dir: Path,
    tissue_label: str,
) -> None:
    script = f'''suppressPackageStartupMessages({{
  library(edgeR)
  library(limma)
}})

counts <- read.delim("{counts_tsv}", check.names=FALSE)
meta <- read.delim("{sample_metadata_tsv}", check.names=FALSE)
comps <- read.delim("{comparisons_tsv}", check.names=FALSE)
selected <- read.delim("{selected_samples_tsv}", check.names=FALSE)
count_cols <- setdiff(colnames(counts), c("gene_id", "gene_symbol"))
count_mat <- as.matrix(counts[, count_cols, drop=FALSE])
storage.mode(count_mat) <- "numeric"
gene_ids <- as.character(counts$gene_id)
gene_symbols <- as.character(counts$gene_symbol)
rownames(count_mat) <- gene_symbols
dir.create("{comparisons_dir}", recursive=TRUE, showWarnings=FALSE)
all_rows <- list()

for (i in seq_len(nrow(comps))) {{
  comp <- comps[i, , drop=FALSE]
  comparison_id <- as.character(comp$comparison_id[1])
  group_a <- as.character(comp$group_a[1])
  group_b <- as.character(comp$group_b[1])
  aging_signature <- as.character(comp$aging_signature[1])
  sel <- selected[selected$comparison_id == comparison_id, , drop=FALSE]
  if (!nrow(sel)) next
  sel$sample_id <- as.character(sel$sample_id)
  sub_meta <- merge(sel, meta, by="sample_id", all.x=TRUE, sort=FALSE)
  sub_meta <- sub_meta[match(sel$sample_id, sub_meta$sample_id), , drop=FALSE]
  ref_samples <- as.character(sel$sample_id[as.character(sel$group_label) == group_b])
  case_samples <- as.character(sel$sample_id[as.character(sel$group_label) == group_a])
  if (length(ref_samples) < 2 || length(case_samples) < 2) next
  selected_samples <- c(ref_samples, case_samples)
  expression <- count_mat[, selected_samples, drop=FALSE]
  colnames(expression) <- paste0("s", seq_len(ncol(expression)) - 1L)
  design <- data.frame(
    A = as.integer(selected_samples %in% ref_samples),
    B = as.integer(selected_samples %in% case_samples),
    row.names = colnames(expression),
    check.names = FALSE
  )
  dge <- DGEList(counts=expression)
  dge <- calcNormFactors(dge)
  v <- voom(dge, plot=FALSE)
  fit <- lmFit(v, as.matrix(design))
  cont.matrix <- makeContrasts(de=B-A, levels=as.matrix(design))
  fit <- contrasts.fit(fit, cont.matrix)
  fit <- eBayes(fit)
  tt <- topTable(fit, adjust="BH", number=nrow(expression))
  tt$gene_symbol <- rownames(tt)
  tt$gene_id <- gene_ids[match(tt$gene_symbol, gene_symbols)]
  tt$comparison_id <- comparison_id
  tt$aging_signature <- aging_signature
  tt$group_a <- group_a
  tt$group_b <- group_b
  tt$stratum <- paste("tissue={tissue_label}", sep="")
  tt$backend <- "maayan_limma_voom_compatible"
  tt$n_group_a <- length(case_samples)
  tt$n_group_b <- length(ref_samples)
  tt$mean_expr <- tt$AveExpr
  tt$model_formula <- "B-A"
  keep_cols <- c("comparison_id", "aging_signature", "gene_id", "gene_symbol", "logFC", "t", "P.Value", "adj.P.Val", "group_a", "group_b", "stratum", "backend", "n_group_a", "n_group_b", "mean_expr", "model_formula")
  tt <- tt[, keep_cols, drop=FALSE]
  write.table(tt, file=file.path("{comparisons_dir}", paste0(comparison_id, ".tsv")), sep="\\t", row.names=FALSE, quote=FALSE)
  colnames(tt)[colnames(tt) == "t"] <- "stat"
  colnames(tt)[colnames(tt) == "P.Value"] <- "pvalue"
  colnames(tt)[colnames(tt) == "adj.P.Val"] <- "padj"
  all_rows[[length(all_rows) + 1]] <- tt
}}

if (length(all_rows)) {{
  out <- do.call(rbind, all_rows)
}} else {{
  out <- data.frame(
    comparison_id=character(),
    aging_signature=character(),
    gene_id=character(),
    gene_symbol=character(),
    logFC=numeric(),
    stat=numeric(),
    pvalue=numeric(),
    padj=numeric(),
    group_a=character(),
    group_b=character(),
    stratum=character(),
    backend=character(),
    n_group_a=integer(),
    n_group_b=integer(),
    mean_expr=numeric(),
    model_formula=character(),
    check.names=FALSE
  )
}}
write.table(out, file="{deg_long_tsv}", sep="\\t", row.names=FALSE, quote=FALSE)
'''
    _write_text(script_path, script)


def _write_deg_long_provenance_graph(
    *,
    output_dir: Path,
    deg_long_path: Path,
    comparison_manifest_path: Path,
    comparison_audit_path: Path,
    selected_samples_path: Path,
    sample_metadata_path: Path,
    counts_path: Path,
    input_paths: list[tuple[Path, str]],
    tissue_label: str,
    filter_mode: str,
    random_state: int,
    min_samples_per_group: int,
) -> Path:
    runtime_ctx = get_runtime_context()
    mirror_local_prefix = runtime_ctx.provenance_mirror_local_prefix if runtime_ctx is not None else None
    mirror_remote_prefix = runtime_ctx.provenance_mirror_remote_prefix if runtime_ctx is not None else None
    input_records = [{"path": str(path), "role": role} for path, role in input_paths]
    input_nodes = [
        build_file_node(
            record,
            {},
            mirror_local_prefix=mirror_local_prefix,
            mirror_remote_prefix=mirror_remote_prefix,
        )
        for record in input_records
    ]
    output_records = [
        build_output_file_record(output_dir, {"path": str(deg_long_path), "role": "deg_tsv"}),
        build_output_file_record(output_dir, {"path": str(comparison_manifest_path), "role": "comparison_manifest"}),
        build_output_file_record(output_dir, {"path": str(comparison_audit_path), "role": "comparison_audit"}),
        build_output_file_record(output_dir, {"path": str(selected_samples_path), "role": "comparison_selected_samples"}),
        build_output_file_record(output_dir, {"path": str(sample_metadata_path), "role": "sample_metadata"}),
        build_output_file_record(output_dir, {"path": str(counts_path), "role": "counts_tsv"}),
    ]
    output_nodes = [
        build_file_node(
            record,
            {},
            mirror_local_prefix=mirror_local_prefix,
            mirror_remote_prefix=mirror_remote_prefix,
        )
        for record in output_records
    ]
    deg_long_node = next(node for node in output_nodes if node["name"] == deg_long_path.name)
    extra_output_nodes = [node for node in output_nodes if node["id"] != deg_long_node["id"]]
    invocation = current_invocation_context()
    command = invocation.get("argv") if invocation else list(sys.argv)
    entrypoint = "geneset-extractors workflows gtex_aging_signatures"
    operation_id = stable_operation_id(
        "gtex_aging_signatures",
        str(deg_long_path),
        [str(node["id"]) for node in input_nodes],
    )
    operation = build_analysis_node(
        analysis_id=operation_id,
        method="gtex_aging_signatures",
        name=f"prepare_{deg_long_path.stem}",
        description="Analysis step that prepares notebook-style GTEx aging-signature differential expression results and emits deg_long.tsv.",
        parameters={
            "tissue_label": tissue_label,
            "filter_mode": filter_mode,
            "random_state": int(random_state),
            "min_samples_per_group": int(min_samples_per_group),
        },
        command=command,
        entrypoint=entrypoint,
        repo_url=REPO_URL,
        module="geneset_extractors.workflows.gtex_aging_signatures",
        script_url=REPO_URL,
        version=_resolve_git_commit(),
        dcc_url=REPO_URL,
        drc_url=REPO_URL,
    )
    graph_path = output_dir / "deg_long.provenance_graph.json"
    payload = mirror_graph_payload(
        {
            "deg_long": {
                "nodes": input_nodes + [operation] + output_nodes,
                "edges": build_edges(input_nodes, str(operation["id"]), str(deg_long_node["id"]), extra_output_nodes),
            }
        },
        mirror_local_prefix,
        mirror_remote_prefix,
    )
    write_canonical_json(graph_path, payload)
    return graph_path


def run(args) -> dict[str, object]:
    expression_gct = Path(args.expression_gct).resolve()
    sample_attributes_tsv = Path(args.sample_attributes_tsv).resolve()
    subject_phenotypes_tsv = Path(args.subject_phenotypes_tsv).resolve()
    human_gene_info = Path(args.human_gene_info).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_attr_rows = _read_tsv_rows(sample_attributes_tsv)
    subject_rows = _read_tsv_rows(subject_phenotypes_tsv)
    ensembl_to_symbol = _read_human_gene_info_mapping(human_gene_info)
    tissue_value = str(args.tissue_value).strip()
    tissue_label = str(getattr(args, "tissue_label", "") or tissue_value).strip()
    tissue_id = str(getattr(args, "tissue_id", "") or _WS_RE.sub("_", tissue_value.strip().lower()).strip("_")).strip()
    metadata_rows = _build_sample_metadata(
        sample_attr_rows,
        subject_rows,
        tissue_column=str(args.tissue_column),
        tissue_value=tissue_value,
        tissue_label=tissue_label,
        tissue_id=tissue_id,
    )
    sample_ids = [row["sample_id"] for row in metadata_rows]

    raw_counts_path = out_dir / "counts_gene_by_sample.raw.tsv"
    matrix_summary = _prepare_counts_matrix(
        expression_gct=expression_gct,
        sample_ids=sample_ids,
        ensembl_to_symbol=ensembl_to_symbol,
        out_path=raw_counts_path,
    )
    selected_sample_ids = {
        str(sample_id)
        for sample_id in matrix_summary.get("selected_sample_ids", [])
        if str(sample_id).strip()
    }
    metadata_rows = [row for row in metadata_rows if row["sample_id"] in selected_sample_ids]
    if not metadata_rows:
        raise ValueError(
            "No overlapping samples remained after intersecting tissue metadata with the expression GCT"
        )
    sample_metadata_path = out_dir / "sample_metadata.tsv"
    _write_tsv(
        sample_metadata_path,
        metadata_rows,
        ["sample_id", "subject_id", "tissue_id", "tissue_label", "tissue_compact", "SMTS", "SMTSD", "sex", "age_group"],
    )

    counts_path = out_dir / "counts_gene_by_sample.tsv"
    filter_summary = {
        "filter_mode": str(args.filter_mode),
        "n_genes_before_filter": matrix_summary["n_rows_after_symbol_dedup"],
        "n_genes_after_filter": matrix_summary["n_rows_after_symbol_dedup"],
    }
    if str(args.filter_mode) == "tissue":
        filter_summary = {
            "filter_mode": str(args.filter_mode),
            **_apply_tissue_filter(
                counts_tsv=raw_counts_path,
                sample_metadata_tsv=sample_metadata_path,
                output_tsv=counts_path,
                rscript_bin=str(args.rscript_bin),
            ),
        }
    else:
        shutil.copyfile(raw_counts_path, counts_path)

    comparison_age_groups = [token.strip() for token in str(args.comparison_age_groups).split(",") if token.strip()]
    comparisons, selected_samples, audit_rows = _build_comparisons(
        metadata_rows,
        tissue_label=tissue_label,
        reference_age_group=str(args.reference_age_group),
        comparison_age_groups=comparison_age_groups,
        random_state=int(args.random_state),
        min_samples_per_group=int(args.min_samples_per_group),
    )
    comparison_manifest_path = out_dir / "comparison_manifest.tsv"
    comparison_audit_path = out_dir / "comparison_audit.tsv"
    selected_samples_path = out_dir / "comparison_selected_samples.tsv"
    _write_tsv(
        comparison_manifest_path,
        comparisons,
        ["comparison_id", "comparison_kind", "aging_signature", "group_column", "group_a", "group_b", "tissue_label"],
    )
    _write_tsv(
        comparison_audit_path,
        audit_rows,
        ["comparison_id", "aging_signature", "tissue_label", "group_a", "group_b", "n_group_a_available", "n_group_b_available", "n_group_a_selected", "n_group_b_selected", "status"],
    )
    _write_tsv(
        selected_samples_path,
        selected_samples,
        ["comparison_id", "sample_id", "group_label", "group_value", "selected_rank", "balance_requested", "balance_applied", "balance_seed", "balance_seed_effective", "age_group"],
    )

    comparisons_dir = out_dir / "comparisons"
    deg_long_path = out_dir / "deg_long.tsv"
    workflow_script = out_dir / "run_limma_voom.R"
    _write_limma_script(
        script_path=workflow_script,
        counts_tsv=counts_path,
        sample_metadata_tsv=sample_metadata_path,
        comparisons_tsv=comparison_manifest_path,
        selected_samples_tsv=selected_samples_path,
        deg_long_tsv=deg_long_path,
        comparisons_dir=comparisons_dir,
        tissue_label=tissue_label,
    )
    result = subprocess.run(
        [str(args.rscript_bin), "--vanilla", str(workflow_script)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise ValueError(
            f"gtex_aging_signatures limma/voom workflow failed with exit code {result.returncode}. "
            f"stderr: {result.stderr.strip() or '[empty]'}"
        )

    summary_path = out_dir / "prepare_summary.json"
    graph_path = _write_deg_long_provenance_graph(
        output_dir=out_dir,
        deg_long_path=deg_long_path,
        comparison_manifest_path=comparison_manifest_path,
        comparison_audit_path=comparison_audit_path,
        selected_samples_path=selected_samples_path,
        sample_metadata_path=sample_metadata_path,
        counts_path=counts_path,
        input_paths=[
            (expression_gct, "expression_gct"),
            (sample_attributes_tsv, "sample_attributes_tsv"),
            (subject_phenotypes_tsv, "subject_phenotypes_tsv"),
            (human_gene_info, "human_gene_info"),
        ],
        tissue_label=tissue_label,
        filter_mode=str(args.filter_mode),
        random_state=int(args.random_state),
        min_samples_per_group=int(args.min_samples_per_group),
    )
    summary = {
        "workflow": "gtex_aging_signatures",
        "out_dir": str(out_dir),
        "expression_gct": str(expression_gct),
        "sample_attributes_tsv": str(sample_attributes_tsv),
        "subject_phenotypes_tsv": str(subject_phenotypes_tsv),
        "human_gene_info": str(human_gene_info),
        "tissue_id": tissue_id,
        "tissue_label": tissue_label,
        "tissue_column": str(args.tissue_column),
        "tissue_value": tissue_value,
        "reference_age_group": str(args.reference_age_group),
        "comparison_age_groups": comparison_age_groups,
        "random_state": int(args.random_state),
        "min_samples_per_group": int(args.min_samples_per_group),
        "filter_summary": filter_summary,
        "count_input_summary": matrix_summary,
        "n_samples": len(metadata_rows),
        "n_comparisons_requested": len(comparison_age_groups),
        "n_comparisons_emitted": len(comparisons),
        "deg_long_tsv": str(deg_long_path),
        "comparison_manifest_tsv": str(comparison_manifest_path),
        "comparison_audit_tsv": str(comparison_audit_path),
        "comparison_selected_samples_tsv": str(selected_samples_path),
        "deg_long_provenance_graph_path": str(graph_path),
        "backend_stdout": result.stdout,
        "backend_stderr": result.stderr,
    }
    _write_json(summary_path, summary)
    return {
        "out_dir": str(out_dir),
        "deg_long_tsv": str(deg_long_path),
        "comparison_manifest_tsv": str(comparison_manifest_path),
        "comparison_audit_tsv": str(comparison_audit_path),
        "comparison_selected_samples_tsv": str(selected_samples_path),
        "n_comparisons": len(comparisons),
        "n_deg_rows": max(0, sum(1 for _ in deg_long_path.open("r", encoding="utf-8")) - 1) if deg_long_path.exists() else 0,
    }
