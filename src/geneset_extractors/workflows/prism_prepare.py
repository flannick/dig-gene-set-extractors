from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil
import time
from urllib.request import Request, urlopen

from geneset_extractors.extractors.drug_response.io import (
    read_prism_cell_line_info_csv,
    read_prism_matrix_csv,
    read_prism_treatment_info_csv,
)
from geneset_extractors.extractors.drug_response.targets import build_drug_targets_from_target_text


def _download_with_retries(
    *,
    source: str,
    out_path: Path,
    max_retries: int,
    retry_backoff_sec: float,
    user_agent: str,
) -> dict[str, object]:
    src = str(source).strip()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    local_source = Path(src)
    if local_source.exists():
        shutil.copyfile(local_source, out_path)
        return {"source": str(local_source), "mode": "local_copy", "path": str(out_path)}

    last_exc: Exception | None = None
    retries = max(1, int(max_retries))
    backoff = max(0.1, float(retry_backoff_sec))
    for attempt in range(1, retries + 1):
        try:
            req = Request(src, headers={"User-Agent": user_agent})
            with urlopen(req, timeout=120) as resp:
                data = resp.read()
            out_path.write_bytes(data)
            return {"source": src, "mode": "download", "attempt": attempt, "path": str(out_path)}
        except Exception as exc:  # pragma: no cover - network errors vary by platform
            last_exc = exc
            if attempt == retries:
                break
            time.sleep(backoff * (2 ** (attempt - 1)))
    raise RuntimeError(f"Failed to fetch {src} after {retries} attempts: {last_exc}")


def _resolve_source_url(*, explicit_url: str | None, file_id: str, file_id_template: str) -> str:
    if explicit_url and str(explicit_url).strip():
        return str(explicit_url).strip()
    fid = str(file_id).strip()
    if not fid:
        raise ValueError("Missing file id and explicit URL")
    return str(file_id_template).format(file_id=fid)


def _write_response_long(path: Path, records) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=["sample_id", "drug_id", "response"])
        writer.writeheader()
        n_rows = 0
        for rec in records:
            writer.writerow({"sample_id": rec.sample_id, "drug_id": rec.drug_id, "response": rec.response})
            n_rows += 1
    return n_rows


def _write_groups(path: Path, sample_meta: dict[str, dict[str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=["sample_id", "group"])
        writer.writeheader()
        n_rows = 0
        for sample_id in sorted(sample_meta):
            group = str(sample_meta[sample_id].get("primary_tissue", "")).strip()
            if not group:
                continue
            writer.writerow({"sample_id": sample_id, "group": group})
            n_rows += 1
    return n_rows


def _write_drug_targets(path: Path, drug_targets: dict[str, dict[str, float]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=["drug_id", "gene_symbol", "weight"])
        writer.writeheader()
        n_rows = 0
        for drug_id in sorted(drug_targets):
            for gene_symbol, weight in sorted(drug_targets[drug_id].items()):
                writer.writerow(
                    {
                        "drug_id": drug_id,
                        "gene_symbol": gene_symbol,
                        "weight": float(weight),
                    }
                )
                n_rows += 1
    return n_rows


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    raw_dir = out_dir / "raw"
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    matrix_source = _resolve_source_url(
        explicit_url=args.matrix_url,
        file_id=args.matrix_file_id,
        file_id_template=args.depmap_file_id_url_template,
    )
    treatment_source = _resolve_source_url(
        explicit_url=args.treatment_info_url,
        file_id=args.treatment_info_file_id,
        file_id_template=args.depmap_file_id_url_template,
    )
    cell_line_source = _resolve_source_url(
        explicit_url=args.cell_line_info_url,
        file_id=args.cell_line_info_file_id,
        file_id_template=args.depmap_file_id_url_template,
    )

    matrix_csv = raw_dir / "prism_matrix.csv"
    treatment_csv = raw_dir / "prism_treatment_info.csv"
    cell_line_csv = raw_dir / "prism_cell_line_info.csv"
    fetch_matrix = _download_with_retries(
        source=matrix_source,
        out_path=matrix_csv,
        max_retries=int(args.max_retries),
        retry_backoff_sec=float(args.retry_backoff_sec),
        user_agent=str(args.user_agent),
    )
    fetch_treatment = _download_with_retries(
        source=treatment_source,
        out_path=treatment_csv,
        max_retries=int(args.max_retries),
        retry_backoff_sec=float(args.retry_backoff_sec),
        user_agent=str(args.user_agent),
    )
    fetch_cell_line = _download_with_retries(
        source=cell_line_source,
        out_path=cell_line_csv,
        max_retries=int(args.max_retries),
        retry_backoff_sec=float(args.retry_backoff_sec),
        user_agent=str(args.user_agent),
    )

    column_to_drug, drug_to_target_text, treatment_summary = read_prism_treatment_info_csv(path=treatment_csv)
    response_records, matrix_summary = read_prism_matrix_csv(
        path=matrix_csv,
        column_to_drug=column_to_drug,
    )
    sample_meta, cell_line_summary = read_prism_cell_line_info_csv(path=cell_line_csv)
    drug_targets, target_summary = build_drug_targets_from_target_text(drug_to_target_text=drug_to_target_text)

    response_long_tsv = out_dir / "response_long.tsv"
    groups_tsv = out_dir / "groups.tsv"
    targets_tsv = out_dir / "drug_targets.tsv"

    n_response_rows = _write_response_long(response_long_tsv, response_records)
    n_group_rows = _write_groups(groups_tsv, sample_meta)
    n_target_rows = _write_drug_targets(targets_tsv, drug_targets)

    summary = {
        "workflow": "prism_prepare",
        "release": str(args.release),
        "fetch": {
            "matrix": fetch_matrix,
            "treatment_info": fetch_treatment,
            "cell_line_info": fetch_cell_line,
        },
        "parse_summary": {
            "treatment": treatment_summary,
            "matrix": matrix_summary,
            "cell_line": cell_line_summary,
            "targets": target_summary,
        },
        "outputs": {
            "response_long_tsv": str(response_long_tsv),
            "groups_tsv": str(groups_tsv),
            "drug_targets_tsv": str(targets_tsv),
            "n_response_rows": n_response_rows,
            "n_group_rows": n_group_rows,
            "n_target_rows": n_target_rows,
        },
    }
    (out_dir / "prepare_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"out_dir": str(out_dir), "n_response_rows": int(n_response_rows)}

