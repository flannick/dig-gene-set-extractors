from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from statistics import mean

from geneset_extractors.hashing import sha256_file
from geneset_extractors.extractors.proteomics.ptm_site_diff_workflow import read_tsv_rows, resolve_ptm_columns


class _CfgShim:
    def __init__(self, ptm_type: str):
        self.site_id_column = None
        self.site_group_column = None
        self.gene_id_column = None
        self.gene_symbol_column = None
        self.protein_accession_column = None
        self.residue_column = None
        self.position_column = None
        self.score_column = None
        self.stat_column = None
        self.logfc_column = None
        self.padj_column = None
        self.pvalue_column = None
        self.localization_prob_column = None
        self.peptide_count_column = None
        self.protein_logfc_column = None
        self.protein_stat_column = None
        self.ptm_type = ptm_type



def _clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()



def _split_tokens(value: str) -> list[str]:
    text = _clean(value)
    if not text:
        return []
    for delim in (";", ",", "|", "/"):
        text = text.replace(delim, ";")
    return [tok.strip() for tok in text.split(";") if tok.strip()]



def _site_candidates(row_values: dict[str, str], columns: dict[str, str | None], ptm_type: str) -> list[str]:
    candidates: list[str] = []
    for key in ("site_id", "site_group"):
        col = columns.get(key)
        if col:
            value = _clean(row_values.get(col))
            if value:
                candidates.append(value)
    protein = _clean(row_values.get(columns.get("protein_accession", ""))) if columns.get("protein_accession") else ""
    residue = _clean(row_values.get(columns.get("residue", ""))) if columns.get("residue") else ""
    position = _clean(row_values.get(columns.get("position", ""))) if columns.get("position") else ""
    gene_symbol = _clean(row_values.get(columns.get("gene_symbol", ""))) if columns.get("gene_symbol") else ""
    if protein and residue and position:
        candidates.append(f"{protein}|{residue}|{position}|{ptm_type}")
    if gene_symbol and residue and position:
        candidates.append(f"{gene_symbol}|{residue}|{position}|{ptm_type}")
    return candidates



def _canonical_site_key(row_values: dict[str, str], columns: dict[str, str | None], ptm_type: str) -> str | None:
    protein = _clean(row_values.get(columns.get("protein_accession", ""))) if columns.get("protein_accession") else ""
    residue = _clean(row_values.get(columns.get("residue", ""))) if columns.get("residue") else ""
    position = _clean(row_values.get(columns.get("position", ""))) if columns.get("position") else ""
    if protein and residue and position:
        return f"{protein}|{residue}|{position}|{ptm_type}"
    for key in ("site_id", "site_group"):
        col = columns.get(key)
        if col:
            value = _clean(row_values.get(col))
            if value:
                return f"{value}|{ptm_type}"
    return None



def _write_tsv_gz(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)



def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sources_path = Path(args.sources_tsv)
    with sources_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        source_rows = list(reader)
    if not source_rows:
        raise ValueError("sources_tsv must contain at least one row")
    if "path" not in reader.fieldnames:
        raise ValueError("sources_tsv must contain a path column")

    ptm_type = str(args.ptm_type).strip()
    organism = str(args.organism).strip().lower()
    cfg = _CfgShim(ptm_type)
    alias_records: dict[str, dict[str, object]] = {}
    site_counts: dict[str, dict[str, object]] = {}
    source_file_records: list[dict[str, object]] = []

    for source_row in source_rows:
        source_file = Path(_clean(source_row.get("path")))
        if not source_file.is_absolute() and not source_file.exists():
            source_file = (sources_path.parent / source_file).resolve()
        source_dataset = _clean(source_row.get("source_dataset")) or source_file.stem
        fieldnames, rows = read_tsv_rows(source_file)
        columns = resolve_ptm_columns(fieldnames, cfg)
        present_in_source: set[str] = set()
        for row in rows:
            canonical = _canonical_site_key(row.values, columns, ptm_type)
            if not canonical:
                continue
            gene_ids = _split_tokens(_clean(row.values.get(columns.get("gene_id", "")))) if columns.get("gene_id") else []
            gene_symbols = _split_tokens(_clean(row.values.get(columns.get("gene_symbol", "")))) if columns.get("gene_symbol") else []
            gene_id = gene_ids[0] if gene_ids else (gene_symbols[0] if gene_symbols else "")
            gene_symbol = gene_symbols[0] if gene_symbols else gene_id
            protein_accession = _clean(row.values.get(columns.get("protein_accession", ""))) if columns.get("protein_accession") else ""
            residue = _clean(row.values.get(columns.get("residue", ""))) if columns.get("residue") else ""
            position = _clean(row.values.get(columns.get("position", ""))) if columns.get("position") else ""
            for input_site_key in _site_candidates(row.values, columns, ptm_type):
                alias_records[input_site_key] = {
                    "input_site_key": input_site_key,
                    "canonical_site_key": canonical,
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "protein_accession": protein_accession,
                    "residue": residue,
                    "position": position,
                    "ptm_type": ptm_type,
                    "source_dataset": source_dataset,
                }
            if canonical in present_in_source:
                continue
            present_in_source.add(canonical)
            stats = site_counts.setdefault(
                canonical,
                {
                    "canonical_site_key": canonical,
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "ptm_type": ptm_type,
                    "df_ref": 0,
                    "datasets": set(),
                },
            )
            stats["df_ref"] = int(stats["df_ref"]) + 1
            stats["datasets"].add(source_dataset)
        source_file_records.append(
            {
                "path": str(source_file),
                "sha256": sha256_file(source_file),
                "source_dataset": source_dataset,
                "n_rows": len(rows),
                "n_unique_sites": len(present_in_source),
            }
        )

    n_ref = len(source_rows)
    alias_rows = [alias_records[key] for key in sorted(alias_records)]
    ubiquity_rows: list[dict[str, object]] = []
    idf_values: list[float] = []
    for canonical in sorted(site_counts):
        stats = site_counts[canonical]
        df_ref = int(stats["df_ref"])
        idf_ref = float(__import__("math").log((n_ref + 1.0) / (df_ref + 1.0)))
        idf_values.append(idf_ref)
        ubiquity_rows.append(
            {
                "canonical_site_key": canonical,
                "gene_id": stats["gene_id"],
                "gene_symbol": stats["gene_symbol"],
                "ptm_type": ptm_type,
                "n_samples_ref": n_ref,
                "df_ref": df_ref,
                "fraction_ref": float(df_ref) / float(n_ref) if n_ref else 0.0,
                "idf_ref": idf_ref,
                "n_datasets_ref": len(stats["datasets"]),
            }
        )

    prefix = "phosphosite" if ptm_type == "phospho" and organism == "human" else f"{ptm_type}_site"
    alias_filename = f"{prefix}_aliases_{organism}_v1.tsv.gz"
    ubiquity_filename = f"{prefix}_ubiquity_{organism}_v1.tsv.gz"
    alias_path = out_dir / alias_filename
    ubiquity_path = out_dir / ubiquity_filename
    provenance_path = out_dir / "bundle_provenance.json"
    local_manifest_path = out_dir / "local_resources_manifest.json"

    _write_tsv_gz(
        alias_path,
        [
            "input_site_key",
            "canonical_site_key",
            "gene_id",
            "gene_symbol",
            "protein_accession",
            "residue",
            "position",
            "ptm_type",
            "source_dataset",
        ],
        alias_rows,
    )
    _write_tsv_gz(
        ubiquity_path,
        [
            "canonical_site_key",
            "gene_id",
            "gene_symbol",
            "ptm_type",
            "n_samples_ref",
            "df_ref",
            "fraction_ref",
            "idf_ref",
            "n_datasets_ref",
        ],
        ubiquity_rows,
    )

    provenance = {
        "bundle_id": args.bundle_id,
        "organism": organism,
        "ptm_type": ptm_type,
        "n_sources": len(source_rows),
        "n_alias_rows": len(alias_rows),
        "n_canonical_sites": len(ubiquity_rows),
        "idf_ref_summary": {
            "min": min(idf_values) if idf_values else 0.0,
            "mean": mean(idf_values) if idf_values else 0.0,
            "max": max(idf_values) if idf_values else 0.0,
        },
        "source_files": source_file_records,
        "outputs": {
            "alias_table": alias_filename,
            "ubiquity_table": ubiquity_filename,
        },
    }
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    manifest = {
        "resources": [
            {
                "id": f"{prefix}_aliases_{organism}_v1",
                "description": f"Canonical {ptm_type} site alias table for {organism}",
                "provider": "local_bundle",
                "stable_id": f"local:{prefix}_aliases_{organism}_v1",
                "version": "v1",
                "genome_build": organism,
                "filename": alias_filename,
                "url": "",
                "sha256": sha256_file(alias_path),
                "license": "Derived from local public source tables; see bundle_provenance.json",
            },
            {
                "id": f"{prefix}_ubiquity_{organism}_v1",
                "description": f"Canonical {ptm_type} site ubiquity prior for {organism}",
                "provider": "local_bundle",
                "stable_id": f"local:{prefix}_ubiquity_{organism}_v1",
                "version": "v1",
                "genome_build": organism,
                "filename": ubiquity_filename,
                "url": "",
                "sha256": sha256_file(ubiquity_path),
                "license": "Derived from local public source tables; see bundle_provenance.json",
            },
        ],
        "presets": {
            "phosphoproteomics_default_optional_human": [
                "phosphosite_aliases_human_v1",
                "phosphosite_ubiquity_human_v1"
            ] if organism == "human" and ptm_type == "phospho" else []
        },
        "bundle": {
            "bundle_id": args.bundle_id,
            "provenance": provenance_path.name,
        },
    }
    local_manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "bundle_id": args.bundle_id,
        "n_sources": len(source_rows),
        "n_alias_rows": len(alias_rows),
        "n_canonical_sites": len(ubiquity_rows),
        "out_dir": str(out_dir),
    }
