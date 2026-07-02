from __future__ import annotations

from pathlib import Path

import pandas as pd

from geneset_extractors.workflows.gtex_runtime_common import write_tsv, write_workflow_provenance_graph


# DER-13_Disorder_DEX_Genes.csv columns (PsychENCODE released cross-disorder layer,
# Gandal et al. 2018, Science aat8127).
GENE_SYMBOL_COLUMN = "Gene_Name"
GENE_ID_COLUMN = "Ensembl_Name"
DISORDER_DIRECTION_COLUMN = "Disorder.DGE_RegulationDirection"

# The disorder/direction field is encoded like "ASD.DGE_down" / "SCZ.DGE_up".
_DIRECTION_SIGN = {"up": "1", "down": "-1"}


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _parse_disorder_direction(value: str) -> tuple[str, str]:
    """Split "ASD.DGE_down" into ("ASD", "down")."""
    text = str(value or "").strip()
    if ".DGE_" not in text:
        raise ValueError(f"Unexpected disorder/direction value: {value!r}")
    disorder, direction = text.split(".DGE_", 1)
    disorder = disorder.strip()
    direction = direction.strip().lower()
    if not disorder or direction not in _DIRECTION_SIGN:
        raise ValueError(f"Unexpected disorder/direction value: {value!r}")
    return disorder, direction


def _load_dex(dex_csv: Path) -> pd.DataFrame:
    frame = pd.read_csv(dex_csv, dtype=str)
    missing = [
        column
        for column in (GENE_SYMBOL_COLUMN, GENE_ID_COLUMN, DISORDER_DIRECTION_COLUMN)
        if column not in frame.columns
    ]
    if missing:
        raise ValueError(f"{dex_csv} is missing expected columns: {missing}")
    return frame


def run(args) -> dict[str, object]:
    dex_csv = Path(args.dex_csv).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    _require_file(dex_csv, "PsychENCODE disorder DEX CSV")

    frame = _load_dex(dex_csv)

    signed_rows: list[dict[str, str]] = []
    for _, row in frame.iterrows():
        gene_id = str(row[GENE_ID_COLUMN]).strip()
        gene_symbol = str(row[GENE_SYMBOL_COLUMN]).strip()
        if not gene_id or gene_id.lower() == "nan":
            continue
        disorder, direction = _parse_disorder_direction(row[DISORDER_DIRECTION_COLUMN])
        signed_rows.append(
            {
                "term": disorder,
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "score": "1",
                "sign": _DIRECTION_SIGN[direction],
            }
        )

    signed_rows.sort(key=lambda r: (r["term"], -int(r["sign"]), r["gene_symbol"]))
    signed_path = out_dir / "psychencode_signed_term_gene.tsv"
    write_tsv(signed_path, signed_rows, ["term", "gene_id", "gene_symbol", "score", "sign"])

    n_terms = len({row["term"] for row in signed_rows})
    n_genes = len({row["gene_id"] for row in signed_rows})
    write_workflow_provenance_graph(
        workflow_name="psychencode_dex",
        module_name="geneset_extractors.workflows.psychencode_dex",
        output_dir=out_dir,
        focus_output_path=signed_path,
        output_paths=[(signed_path, "signed_term_gene_tsv")],
        input_paths=[(dex_csv, "disorder_dex_csv")],
        parameters={
            "n_rows": len(signed_rows),
            "n_terms": n_terms,
            "n_genes": n_genes,
        },
    )
    return {"n_rows": len(signed_rows), "out_dir": str(out_dir)}
