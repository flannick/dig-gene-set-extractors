from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import re
import shutil
import tarfile
import xml.etree.ElementTree as ET

from geneset_extractors.core.metadata import _resolve_git_commit, input_file_record
from geneset_extractors.core.provenance import (
    REPO_URL,
    build_analysis_node,
    build_edges,
    build_file_node,
    build_output_file_record,
    stable_operation_id,
    write_canonical_json,
)


MINIML_NS = {"m": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}


def _split_csv(value: str | None) -> set[str]:
    return {item.strip() for item in str(value or "").split(",") if item.strip()}


def _safe_column(value: str) -> str:
    cleaned = re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")
    return cleaned or "characteristic"


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _read_miniml_root(path: Path) -> ET.Element:
    if path.name.endswith((".tgz", ".tar.gz")):
        with tarfile.open(path, "r:gz") as archive:
            members = [member for member in archive.getmembers() if member.isfile() and member.name.endswith(".xml")]
            if len(members) != 1:
                raise ValueError(f"Expected exactly one MINiML XML file in {path}, found {len(members)}")
            handle = archive.extractfile(members[0])
            if handle is None:
                raise ValueError(f"Unable to read MINiML XML from {path}")
            return ET.fromstring(handle.read())
    return ET.parse(path).getroot()


def _parse_samples(
    miniml_path: Path,
    *,
    sample_id_field: str,
    group_characteristic: str,
    condition_a_values: set[str],
    condition_b_values: set[str],
    condition_a_label: str,
    condition_b_label: str,
) -> list[dict[str, str]]:
    root = _read_miniml_root(miniml_path)
    rows: list[dict[str, str]] = []
    seen_sample_ids: set[str] = set()
    for sample in root.findall("m:Sample", MINIML_NS):
        title = (sample.findtext("m:Title", default="", namespaces=MINIML_NS) or "").strip()
        accession = (sample.findtext("m:Accession", default="", namespaces=MINIML_NS) or "").strip()
        channel = sample.find("m:Channel", MINIML_NS)
        if channel is None or not title or not accession:
            continue
        characteristics: dict[str, str] = {}
        for element in channel.findall("m:Characteristics", MINIML_NS):
            tag = str(element.attrib.get("tag", "")).strip()
            if tag:
                characteristics[tag] = (element.text or "").strip()
        source = (channel.findtext("m:Source", default="", namespaces=MINIML_NS) or "").strip()
        reserved_group_values = {
            "__title__": title,
            "__accession__": accession,
            "__source__": source,
        }
        group_value = reserved_group_values.get(
            group_characteristic,
            characteristics.get(group_characteristic, ""),
        )
        if group_value in condition_a_values:
            condition = condition_a_label
        elif group_value in condition_b_values:
            condition = condition_b_label
        else:
            continue
        sample_id = accession if sample_id_field == "accession" else title
        if sample_id in seen_sample_ids:
            raise ValueError(f"Duplicate GEO sample identifier in MINiML metadata: {sample_id}")
        seen_sample_ids.add(sample_id)
        row = {
            "sample_id": sample_id,
            "gsm_id": accession,
            "geo_title": title,
            "condition": condition,
            "geo_group_value": group_value,
            "source": source,
            "organism": (channel.findtext("m:Organism", default="", namespaces=MINIML_NS) or "").strip(),
        }
        for tag, value in characteristics.items():
            row.setdefault(_safe_column(tag), value)
        rows.append(row)
    if not rows:
        raise ValueError("No MINiML samples matched the configured GEO comparison groups")
    counts = {condition_a_label: 0, condition_b_label: 0}
    for row in rows:
        counts[row["condition"]] = counts.get(row["condition"], 0) + 1
    if counts.get(condition_a_label, 0) < 2 or counts.get(condition_b_label, 0) < 2:
        raise ValueError(f"GEO comparison requires at least two samples per group; observed {counts}")
    return rows


def _write_counts(counts_path: Path, out_path: Path, selected_sample_ids: set[str]) -> dict[str, object]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with _open_text(counts_path) as source, out_path.open("w", encoding="utf-8", newline="") as target:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Empty GEO counts file: {counts_path}") from exc
        if len(header) < 3:
            raise ValueError("GEO count matrix must contain a feature column and at least two sample columns")
        missing = sorted(selected_sample_ids.difference(header[1:]))
        if missing:
            raise ValueError(f"Selected MINiML samples are absent from the count matrix: {missing[:10]}")
        keep_indices = [0] + [idx for idx, name in enumerate(header[1:], start=1) if name in selected_sample_ids]
        writer.writerow([header[idx] for idx in keep_indices])
        n_features = 0
        for row in reader:
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(f"Count row has {len(row)} columns; expected {len(header)}")
            writer.writerow([row[idx] for idx in keep_indices])
            n_features += 1
    return {"feature_id_column": header[0], "n_features": n_features, "n_samples": len(keep_indices) - 1}


def _write_metadata(rows: list[dict[str, str]], out_path: Path) -> None:
    preferred = ["sample_id", "gsm_id", "geo_title", "condition", "geo_group_value", "source", "organism"]
    extra = sorted({key for row in rows for key in row}.difference(preferred))
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=preferred + extra, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_feature_mapping(
    annotation_path: Path,
    out_path: Path,
    *,
    source_column: str,
    target_column: str,
) -> int:
    with _open_text(annotation_path) as source:
        reader = csv.DictReader(source, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if source_column not in fields or target_column not in fields:
            raise ValueError(
                f"Annotation table requires {source_column!r} and {target_column!r}; found {sorted(fields)}"
            )
        mappings: dict[str, str] = {}
        for row in reader:
            symbol = str(row.get(target_column, "")).strip()
            raw_ids = str(row.get(source_column, "")).strip()
            if not symbol or not raw_ids:
                continue
            for feature_id in re.split(r"(?:///|[|,;])", raw_ids):
                clean_id = feature_id.strip().split(".", 1)[0]
                if clean_id and clean_id not in mappings:
                    mappings[clean_id] = symbol
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["source_feature_id", "gene_symbol"],
            lineterminator="\n",
        )
        writer.writeheader()
        for feature_id in sorted(mappings):
            writer.writerow({"source_feature_id": feature_id, "gene_symbol": mappings[feature_id]})
    return len(mappings)


def _source_record(path: Path, role: str, *, source_url: str | None, landing_page_url: str | None) -> dict[str, object]:
    record: dict[str, object] = dict(input_file_record(path, role))
    if source_url:
        record.update(
            {
                "canonical_uri": source_url,
                "download_url": source_url,
                "provider": "NCBI GEO",
                "access_level": "public",
                "landing_page_url": landing_page_url,
            }
        )
    return record


def _write_provenance(
    *,
    out_dir: Path,
    study_id: str,
    raw_records: list[dict[str, object]],
    output_paths: list[tuple[Path, str]],
    parameters: dict[str, object],
) -> Path:
    input_nodes = [build_file_node(record, {}) for record in raw_records]
    output_nodes = [
        build_file_node(build_output_file_record(out_dir, {"path": str(path), "role": role}), {})
        for path, role in output_paths
    ]
    operation_id = stable_operation_id("geo_bulk_prepare", study_id, [str(node["id"]) for node in input_nodes])
    operation = build_analysis_node(
        analysis_id=operation_id,
        method="geo_bulk_prepare",
        name=f"prepare_{study_id}",
        description="Standardize a public GEO bulk RNA count matrix, MINiML sample metadata, and gene annotation.",
        parameters=parameters,
        command=None,
        entrypoint="geneset-extractors workflows geo_bulk_prepare",
        repo_url=REPO_URL,
        module="geneset_extractors.workflows.geo_bulk_prepare",
        script_url=REPO_URL,
        version=_resolve_git_commit(),
        dcc_url=parameters.get("landing_page_url") or REPO_URL,
        drc_url=REPO_URL,
    )
    primary_output = output_nodes[0]
    extra_outputs = output_nodes[1:]
    payload = {
        "geo_bulk_inputs": {
            "nodes": input_nodes + [operation] + output_nodes,
            "edges": build_edges(input_nodes, str(operation["id"]), str(primary_output["id"]), extra_outputs),
        }
    }
    path = out_dir / "geo_bulk_inputs.provenance_graph.json"
    write_canonical_json(path, payload)
    return path


def prepare_geo_bulk_inputs(
    *,
    counts_file: str,
    miniml_file: str,
    annotation_file: str,
    out_dir: str,
    study_id: str,
    sample_id_field: str = "title",
    group_characteristic: str,
    condition_a_values: str,
    condition_b_values: str,
    condition_a_label: str,
    condition_b_label: str,
    annotation_source_column: str,
    annotation_target_column: str,
    counts_source_url: str | None = None,
    miniml_source_url: str | None = None,
    annotation_source_url: str | None = None,
    landing_page_url: str | None = None,
) -> dict[str, object]:
    output_dir = Path(out_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    counts_path = Path(counts_file).resolve()
    miniml_path = Path(miniml_file).resolve()
    annotation_path = Path(annotation_file).resolve()
    for path in (counts_path, miniml_path, annotation_path):
        if not path.exists() or not path.is_file():
            raise ValueError(f"Missing GEO input file: {path}")

    rows = _parse_samples(
        miniml_path,
        sample_id_field=sample_id_field,
        group_characteristic=group_characteristic,
        condition_a_values=_split_csv(condition_a_values),
        condition_b_values=_split_csv(condition_b_values),
        condition_a_label=condition_a_label,
        condition_b_label=condition_b_label,
    )
    counts_out = output_dir / "counts.tsv"
    metadata_out = output_dir / "sample_metadata.tsv"
    mapping_out = output_dir / "feature_mapping.tsv"
    count_summary = _write_counts(counts_path, counts_out, {row["sample_id"] for row in rows})
    _write_metadata(rows, metadata_out)
    n_mappings = _write_feature_mapping(
        annotation_path,
        mapping_out,
        source_column=annotation_source_column,
        target_column=annotation_target_column,
    )
    summary = {
        "workflow": "geo_bulk_prepare",
        "study_id": study_id,
        "sample_id_field": sample_id_field,
        "group_characteristic": group_characteristic,
        "condition_a_values": sorted(_split_csv(condition_a_values)),
        "condition_b_values": sorted(_split_csv(condition_b_values)),
        "condition_a_label": condition_a_label,
        "condition_b_label": condition_b_label,
        "n_samples": len(rows),
        "n_condition_a": sum(row["condition"] == condition_a_label for row in rows),
        "n_condition_b": sum(row["condition"] == condition_b_label for row in rows),
        "n_feature_mappings": n_mappings,
        **count_summary,
    }
    summary_path = output_dir / "prepare_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    raw_records = [
        _source_record(counts_path, "geo_counts_archive", source_url=counts_source_url, landing_page_url=landing_page_url),
        _source_record(miniml_path, "geo_miniml_archive", source_url=miniml_source_url, landing_page_url=landing_page_url),
        _source_record(
            annotation_path,
            "geo_gene_annotation_archive",
            source_url=annotation_source_url,
            landing_page_url=landing_page_url,
        ),
    ]
    provenance_path = _write_provenance(
        out_dir=output_dir,
        study_id=study_id,
        raw_records=raw_records,
        output_paths=[
            (counts_out, "prepared_counts_tsv"),
            (metadata_out, "prepared_sample_metadata_tsv"),
            (mapping_out, "prepared_feature_mapping_tsv"),
            (summary_path, "prepare_summary_json"),
        ],
        parameters={**summary, "landing_page_url": landing_page_url or ""},
    )
    return {
        **summary,
        "out_dir": str(output_dir),
        "counts_tsv": str(counts_out),
        "sample_metadata_tsv": str(metadata_out),
        "feature_mapping_tsv": str(mapping_out),
        "provenance_graph": str(provenance_path),
    }


def run(args) -> dict[str, object]:
    return prepare_geo_bulk_inputs(
        counts_file=args.counts_file,
        miniml_file=args.miniml_file,
        annotation_file=args.annotation_file,
        out_dir=args.out_dir,
        study_id=args.study_id,
        sample_id_field=getattr(args, "sample_id_field", "title"),
        group_characteristic=args.group_characteristic,
        condition_a_values=args.condition_a_values,
        condition_b_values=args.condition_b_values,
        condition_a_label=args.condition_a_label,
        condition_b_label=args.condition_b_label,
        annotation_source_column=args.annotation_source_column,
        annotation_target_column=args.annotation_target_column,
        counts_source_url=getattr(args, "counts_source_url", None),
        miniml_source_url=getattr(args, "miniml_source_url", None),
        annotation_source_url=getattr(args, "annotation_source_url", None),
        landing_page_url=getattr(args, "landing_page_url", None),
    )
