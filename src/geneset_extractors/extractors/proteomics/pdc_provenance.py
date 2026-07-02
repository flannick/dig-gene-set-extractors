"""Build the engine provenance overlay + source map from a PDC file manifest (stdlib only)."""
from __future__ import annotations

import csv
import json
from pathlib import Path

CRDC_DRC_URL = "https://datacommons.cancer.gov/repository/proteomic-data-commons"


def pdc_study_url(pdc_study_id: str) -> str:
    return f"https://pdc.cancer.gov/pdc/study/{pdc_study_id}"


def _by_role(manifest_rows: list[dict]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for row in manifest_rows:
        out[row["role"]] = row
    if "phospho" not in out:
        raise ValueError("manifest is missing a 'phospho' role row")
    return out


def _file_overlay(rec: dict, *, provider_suffix: str) -> dict:
    return {
        "persistent_id": rec["file_id"],
        "local_id": rec["drs_uri"],
        "dcc_url": pdc_study_url(rec["pdc_study_id"]),
        "drc_url": CRDC_DRC_URL,
        "provider": f"NCI Proteomic Data Commons (CPTAC) — {provider_suffix}",
        "version": rec["pdc_study_id"],
        "source": Path(rec["local_path"]).name,
    }


def build_overlay(*, manifest_rows: list[dict], prepared_dir: str, operation_meta: dict) -> dict:
    roles = _by_role(manifest_rows)
    phospho = roles["phospho"]
    prepared = Path(prepared_dir)

    inputs: dict[str, dict] = {
        str(prepared / "ptm_matrix.tsv"): _file_overlay(phospho, provider_suffix="derived from phosphosite report"),
    }
    if "proteome" in roles:
        inputs[str(prepared / "protein_matrix.tsv")] = _file_overlay(
            roles["proteome"], provider_suffix="derived from proteome report"
        )
    # sample_metadata derives from the phospho study's biospecimen records.
    meta_overlay = _file_overlay(phospho, provider_suffix="biospecimen-derived sample metadata")
    inputs[str(prepared / "sample_metadata.tsv")] = meta_overlay

    study_url = pdc_study_url(phospho["pdc_study_id"])
    operation = dict(operation_meta)
    operation.setdefault("dcc_url", study_url)
    operation.setdefault("drc_url", CRDC_DRC_URL)
    operation.setdefault("description", "CPTAC tumor-vs-normal phosphoregulation extraction via ptm_site_matrix")

    return {
        "inputs": inputs,
        "operation": operation,
        "gene_set": {"dcc_url": study_url, "drc_url": CRDC_DRC_URL},
    }


def write_overlay(*, manifest_rows: list[dict], prepared_dir: str, operation_meta: dict, out_dir: str | Path) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    overlay = build_overlay(manifest_rows=manifest_rows, prepared_dir=prepared_dir, operation_meta=operation_meta)
    overlay_path = out_dir / "provenance_overlay.json"
    overlay_path.write_text(json.dumps(overlay, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    source_map_path = out_dir / "local_input_source_map.tsv"
    with source_map_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["local_path", "source_uri"])
        for row in manifest_rows:
            writer.writerow([row["local_path"], row["drs_uri"]])
    return {"overlay_json": overlay_path, "source_map_tsv": source_map_path}


def read_manifest_tsv(path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
