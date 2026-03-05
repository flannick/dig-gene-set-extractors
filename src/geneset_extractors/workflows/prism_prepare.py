from __future__ import annotations

import csv
import json
from pathlib import Path
import random
import shutil
import sys
import time
from urllib.request import Request, urlopen

from geneset_extractors.extractors.drug_response.io import (
    ResponseRecord,
    read_prism_cell_line_info_csv,
    read_prism_matrix_csv,
    read_prism_treatment_info_csv,
)
from geneset_extractors.extractors.drug_response.targets import build_drug_targets_from_target_text


_FIGSHARE_TEMPLATE = "https://ndownloader.figshare.com/files/{file_id}"


def _preview_text(raw: bytes, n: int = 200) -> str:
    return raw[:n].decode("utf-8", errors="replace").replace("\n", "\\n")


def sniff_delimited_payload(raw: bytes) -> dict[str, str]:
    text = raw.decode("utf-8", errors="replace")
    stripped = text.lstrip()
    preview = _preview_text(raw)
    if not stripped:
        return {"kind": "empty", "preview": preview}
    if stripped.startswith("{") or stripped.startswith("["):
        return {"kind": "json_like", "preview": preview}
    if stripped.startswith("<"):
        return {"kind": "html_like", "preview": preview}
    first_line = ""
    for line in text.splitlines():
        if line.strip():
            first_line = line
            break
    if "\t" in first_line or "," in first_line:
        return {"kind": "delimited_text", "preview": preview}
    return {"kind": "unknown_text", "preview": preview}


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
        raw = out_path.read_bytes()
        sniff = sniff_delimited_payload(raw)
        return {
            "source": str(local_source),
            "mode": "local_copy",
            "path": str(out_path),
            "sniff": sniff,
        }

    last_exc: Exception | None = None
    retries = max(1, int(max_retries))
    backoff = max(0.1, float(retry_backoff_sec))
    for attempt in range(1, retries + 1):
        try:
            req = Request(src, headers={"User-Agent": user_agent})
            with urlopen(req, timeout=120) as resp:
                data = resp.read()
            out_path.write_bytes(data)
            sniff = sniff_delimited_payload(data)
            return {
                "source": src,
                "mode": "download",
                "attempt": attempt,
                "path": str(out_path),
                "sniff": sniff,
            }
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


def _download_prism_file(
    *,
    kind: str,
    explicit_url: str | None,
    file_id: str,
    template: str,
    out_path: Path,
    max_retries: int,
    retry_backoff_sec: float,
    user_agent: str,
) -> dict[str, object]:
    source = _resolve_source_url(explicit_url=explicit_url, file_id=file_id, file_id_template=template)
    fetch = _download_with_retries(
        source=source,
        out_path=out_path,
        max_retries=max_retries,
        retry_backoff_sec=retry_backoff_sec,
        user_agent=user_agent,
    )
    sniff = dict(fetch.get("sniff", {}))
    fallback_used = False
    template_used = template

    # If content does not look like tabular data and this came from file-id template mode,
    # retry once with the Figshare template.
    if sniff.get("kind") != "delimited_text" and explicit_url is None and template != _FIGSHARE_TEMPLATE:
        print(
            "warning: downloaded content does not look like PRISM tabular data "
            f"for {kind} (kind={sniff.get('kind')}); retrying with figshare template.",
            file=sys.stderr,
        )
        fallback_source = _resolve_source_url(
            explicit_url=None,
            file_id=file_id,
            file_id_template=_FIGSHARE_TEMPLATE,
        )
        fetch = _download_with_retries(
            source=fallback_source,
            out_path=out_path,
            max_retries=max_retries,
            retry_backoff_sec=retry_backoff_sec,
            user_agent=user_agent,
        )
        sniff = dict(fetch.get("sniff", {}))
        fallback_used = True
        template_used = _FIGSHARE_TEMPLATE

    if sniff.get("kind") != "delimited_text":
        preview = str(sniff.get("preview", ""))
        raise RuntimeError(
            f"{kind} download does not look like delimited text (kind={sniff.get('kind')}). "
            "Likely portal metadata/HTML. Try --depmap_file_id_url_template "
            f"'{_FIGSHARE_TEMPLATE}'. Preview: {preview}"
        )

    fetch["template_used"] = template_used
    fetch["fallback_to_figshare"] = fallback_used
    fetch["kind"] = kind
    return fetch


def _write_response_long(path: Path, records: list[ResponseRecord]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=["sample_id", "drug_id", "response"])
        writer.writeheader()
        n_rows = 0
        for rec in records:
            writer.writerow({"sample_id": rec.sample_id, "drug_id": rec.drug_id, "response": rec.response})
            n_rows += 1
    return n_rows


def _write_groups(path: Path, sample_meta: dict[str, dict[str, str]], *, group_column: str) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=["sample_id", "group"])
        writer.writeheader()
        n_rows = 0
        for sample_id in sorted(sample_meta):
            group = str(sample_meta[sample_id].get(group_column, "")).strip()
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
                writer.writerow({"drug_id": drug_id, "gene_symbol": gene_symbol, "weight": float(weight)})
                n_rows += 1
    return n_rows


def _sample_without_replacement(items: list[str], n: int, rng: random.Random) -> list[str]:
    if n >= len(items):
        return list(items)
    return sorted(rng.sample(items, n))


def _subset_response_records(
    *,
    records: list[ResponseRecord],
    sample_meta: dict[str, dict[str, str]],
    subset_seed: int,
    max_cell_lines_total: int | None,
    max_compounds_total: int | None,
    max_cell_lines_per_group: int | None,
    balance_by: str | None,
    min_per_balance_bin: int,
) -> tuple[list[ResponseRecord], dict[str, dict[str, str]], dict[str, object], list[str]]:
    warnings: list[str] = []
    rng = random.Random(int(subset_seed))
    all_samples = sorted({str(r.sample_id) for r in records})
    samples_with_meta = [s for s in all_samples if s in sample_meta]
    dropped_no_meta = sorted(set(all_samples) - set(samples_with_meta))
    if dropped_no_meta:
        warnings.append(f"dropped {len(dropped_no_meta)} samples with no metadata rows")

    balance_col = str(balance_by or "").strip()
    bins: dict[str, list[str]] = {}
    if balance_col:
        for sample_id in samples_with_meta:
            key = str(sample_meta.get(sample_id, {}).get(balance_col, "")).strip()
            if key:
                bins.setdefault(key, []).append(sample_id)
        if not bins:
            warnings.append(
                f"balance_by={balance_col} had no populated bins; falling back to unstratified sampling"
            )
    if not bins:
        # Optionally treat primary_tissue as group bin when per-group cap requested.
        if max_cell_lines_per_group and any("primary_tissue" in sample_meta.get(s, {}) for s in samples_with_meta):
            for sample_id in samples_with_meta:
                key = str(sample_meta.get(sample_id, {}).get("primary_tissue", "")).strip() or "unknown"
                bins.setdefault(key, []).append(sample_id)
        else:
            bins = {"all": list(samples_with_meta)}

    selected_samples: list[str] = []
    dropped_bins: dict[str, int] = {}
    min_bin = max(1, int(min_per_balance_bin))
    for bin_key in sorted(bins):
        members = sorted(bins[bin_key])
        if len(members) < min_bin and len(bins) > 1:
            dropped_bins[bin_key] = len(members)
            continue
        cap = int(max_cell_lines_per_group) if max_cell_lines_per_group else len(members)
        selected_samples.extend(_sample_without_replacement(members, cap, rng))

    if dropped_bins:
        warnings.append(
            "dropped undersized bins: "
            + ", ".join(f"{k} (n={v})" for k, v in sorted(dropped_bins.items()))
        )

    selected_samples = sorted(set(selected_samples))
    if max_cell_lines_total and len(selected_samples) > int(max_cell_lines_total):
        selected_samples = _sample_without_replacement(selected_samples, int(max_cell_lines_total), rng)

    selected_sample_set = set(selected_samples)
    selected_records = [r for r in records if str(r.sample_id) in selected_sample_set]
    all_drugs = sorted({str(r.drug_id) for r in selected_records})
    selected_drugs = list(all_drugs)
    if max_compounds_total and len(all_drugs) > int(max_compounds_total):
        selected_drugs = _sample_without_replacement(all_drugs, int(max_compounds_total), rng)
    selected_drug_set = set(selected_drugs)
    selected_records = [r for r in selected_records if str(r.drug_id) in selected_drug_set]

    selected_meta = {sample_id: sample_meta[sample_id] for sample_id in selected_samples if sample_id in sample_meta}
    summary = {
        "subset_seed": int(subset_seed),
        "balance_by": balance_col or None,
        "min_per_balance_bin": min_bin,
        "max_cell_lines_per_group": (int(max_cell_lines_per_group) if max_cell_lines_per_group else None),
        "max_cell_lines_total": (int(max_cell_lines_total) if max_cell_lines_total else None),
        "max_compounds_total": (int(max_compounds_total) if max_compounds_total else None),
        "n_samples_input": len(all_samples),
        "n_samples_with_metadata": len(samples_with_meta),
        "n_samples_selected": len(selected_samples),
        "n_compounds_selected": len(selected_drugs),
        "n_response_rows_selected": len(selected_records),
        "n_bins_total": len(bins),
        "n_bins_dropped": len(dropped_bins),
        "dropped_bins": [{"bin": k, "n": v} for k, v in sorted(dropped_bins.items())],
    }
    return selected_records, selected_meta, summary, warnings


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    raw_dir = out_dir / "raw"
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    template = str(args.depmap_file_id_url_template)
    matrix_csv = raw_dir / "prism_matrix.csv"
    treatment_csv = raw_dir / "prism_treatment_info.csv"
    cell_line_csv = raw_dir / "prism_cell_line_info.csv"

    fetch_matrix = _download_prism_file(
        kind="matrix",
        explicit_url=args.matrix_url,
        file_id=str(args.matrix_file_id),
        template=template,
        out_path=matrix_csv,
        max_retries=int(args.max_retries),
        retry_backoff_sec=float(args.retry_backoff_sec),
        user_agent=str(args.user_agent),
    )
    fetch_treatment = _download_prism_file(
        kind="treatment_info",
        explicit_url=args.treatment_info_url,
        file_id=str(args.treatment_info_file_id),
        template=template,
        out_path=treatment_csv,
        max_retries=int(args.max_retries),
        retry_backoff_sec=float(args.retry_backoff_sec),
        user_agent=str(args.user_agent),
    )
    fetch_cell_line = _download_prism_file(
        kind="cell_line_info",
        explicit_url=args.cell_line_info_url,
        file_id=str(args.cell_line_info_file_id),
        template=template,
        out_path=cell_line_csv,
        max_retries=int(args.max_retries),
        retry_backoff_sec=float(args.retry_backoff_sec),
        user_agent=str(args.user_agent),
    )

    column_to_drug, drug_to_target_text, treatment_summary = read_prism_treatment_info_csv(path=treatment_csv)
    response_records, matrix_summary = read_prism_matrix_csv(path=matrix_csv, column_to_drug=column_to_drug)
    sample_meta, cell_line_summary = read_prism_cell_line_info_csv(path=cell_line_csv)
    drug_targets, target_summary = build_drug_targets_from_target_text(drug_to_target_text=drug_to_target_text)

    selected_records, selected_meta, subset_summary, subset_warnings = _subset_response_records(
        records=response_records,
        sample_meta=sample_meta,
        subset_seed=int(args.subset_seed),
        max_cell_lines_total=(None if args.max_cell_lines_total in {None, 0} else int(args.max_cell_lines_total)),
        max_compounds_total=(None if args.max_compounds_total in {None, 0} else int(args.max_compounds_total)),
        max_cell_lines_per_group=(
            None if args.max_cell_lines_per_group in {None, 0} else int(args.max_cell_lines_per_group)
        ),
        balance_by=(str(args.balance_by).strip() if args.balance_by else None),
        min_per_balance_bin=int(args.min_per_balance_bin),
    )
    selected_drugs = sorted({str(r.drug_id) for r in selected_records})
    drug_targets = {drug_id: dict(drug_targets[drug_id]) for drug_id in selected_drugs if drug_id in drug_targets}

    group_column = str(args.balance_by).strip() if args.balance_by else "primary_tissue"
    response_long_tsv = out_dir / "response_long.tsv"
    groups_tsv = out_dir / "groups.tsv"
    targets_tsv = out_dir / "drug_targets.tsv"

    n_response_rows = _write_response_long(response_long_tsv, selected_records)
    n_group_rows = _write_groups(groups_tsv, selected_meta, group_column=group_column)
    n_target_rows = _write_drug_targets(targets_tsv, drug_targets)

    if not bool(args.keep_raw_downloads):
        shutil.rmtree(raw_dir, ignore_errors=True)

    summary = {
        "workflow": "prism_prepare",
        "release": str(args.release),
        "fetch": {
            "matrix": fetch_matrix,
            "treatment_info": fetch_treatment,
            "cell_line_info": fetch_cell_line,
            "default_template": _FIGSHARE_TEMPLATE,
        },
        "parse_summary": {
            "treatment": treatment_summary,
            "matrix": matrix_summary,
            "cell_line": cell_line_summary,
            "targets": target_summary,
        },
        "subset": subset_summary,
        "warnings": subset_warnings,
        "outputs": {
            "response_long_tsv": str(response_long_tsv),
            "groups_tsv": str(groups_tsv),
            "group_column": group_column,
            "drug_targets_tsv": str(targets_tsv),
            "n_response_rows": n_response_rows,
            "n_group_rows": n_group_rows,
            "n_target_rows": n_target_rows,
            "n_cell_lines": len({str(r.sample_id) for r in selected_records}),
            "n_compounds": len({str(r.drug_id) for r in selected_records}),
        },
    }
    (out_dir / "prepare_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"out_dir": str(out_dir), "n_response_rows": int(n_response_rows)}
