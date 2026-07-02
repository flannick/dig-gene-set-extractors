from __future__ import annotations

from pathlib import Path

import pandas as pd

from geneset_extractors.workflows.gtex_runtime_common import write_tsv, write_workflow_provenance_graph


# DER-16_Disorder_Gene_Modules.csv columns (PsychENCODE released cross-disorder layer,
# Gandal et al. 2018, Science aat8127). One row per gene; each gene is assigned to a single
# WGCNA co-expression module via the Module column.
GENE_ID_COLUMN = "ensembl_gene_id"
GENE_SYMBOL_COLUMN = "gene_name"
MODULE_COLUMN = "Module"


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _load_modules(modules_csv: Path) -> pd.DataFrame:
    frame = pd.read_csv(modules_csv, dtype=str, usecols=lambda c: c in {GENE_ID_COLUMN, GENE_SYMBOL_COLUMN, MODULE_COLUMN})
    missing = [
        column
        for column in (GENE_ID_COLUMN, GENE_SYMBOL_COLUMN, MODULE_COLUMN)
        if column not in frame.columns
    ]
    if missing:
        raise ValueError(f"{modules_csv} is missing expected columns: {missing}")
    return frame


def run(args) -> dict[str, object]:
    modules_csv = Path(args.modules_csv).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    _require_file(modules_csv, "PsychENCODE gene modules CSV")

    # The WGCNA grey/unassigned bucket (default geneM0) is not a real co-expression
    # module and is excluded from the emitted gene sets.
    excluded = {
        token.strip()
        for token in str(getattr(args, "exclude_modules", "geneM0") or "").split(",")
        if token.strip()
    }

    frame = _load_modules(modules_csv)

    unsigned_rows: list[dict[str, str]] = []
    for _, row in frame.iterrows():
        gene_id = str(row[GENE_ID_COLUMN]).strip()
        gene_symbol = str(row[GENE_SYMBOL_COLUMN]).strip()
        module = str(row[MODULE_COLUMN]).strip()
        if not gene_id or gene_id.lower() == "nan":
            continue
        if not module or module.lower() == "nan" or module in excluded:
            continue
        unsigned_rows.append(
            {
                "term": module,
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "score": "1",
            }
        )

    unsigned_rows.sort(key=lambda r: (r["term"], r["gene_symbol"]))
    unsigned_path = out_dir / "psychencode_unsigned_term_gene.tsv"
    write_tsv(unsigned_path, unsigned_rows, ["term", "gene_id", "gene_symbol", "score"])

    n_terms = len({row["term"] for row in unsigned_rows})
    n_genes = len({row["gene_id"] for row in unsigned_rows})
    write_workflow_provenance_graph(
        workflow_name="psychencode_modules",
        module_name="geneset_extractors.workflows.psychencode_modules",
        output_dir=out_dir,
        focus_output_path=unsigned_path,
        output_paths=[(unsigned_path, "unsigned_term_gene_tsv")],
        input_paths=[(modules_csv, "gene_modules_csv")],
        parameters={
            "n_rows": len(unsigned_rows),
            "n_terms": n_terms,
            "n_genes": n_genes,
            "excluded_modules": sorted(excluded),
        },
    )
    return {"n_rows": len(unsigned_rows), "out_dir": str(out_dir)}
