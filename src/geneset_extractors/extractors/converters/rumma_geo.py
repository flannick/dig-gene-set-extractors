"""Reconstruct RummaGEO signed perturbation libraries from recorded source inputs.

The historical notebooks queried Harmonizome for a selection table before
processing GMT memberships.  That query result is a scientific input, rather
than something that can be recovered reliably from GEO sample titles, so this
converter requires it as a tab-delimited manifest.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context


DRUG_FALSE_POSITIVES = frozenset({"1B", "ATPA", "AVA", "C-1", "CDC", "FIT", "ITE", "PIT", "RITA", "TRIM", "compe", "iq", "niacin", "pen", "rutin"})
VALID_STATUSES = frozenset({"signature", "reversed"})
MODEL_FAMILIES = {"HZ1": "drug_perturbations", "HZ2": "gene_perturbations"}


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def _rows(path: Path) -> Iterable[dict[str, str]]:
    with _open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def _mappings(human_gene_info: Path, mouse_gene_info: Path, gene_orthologs: Path) -> tuple[set[str], dict[str, str]]:
    human_rows = [row for row in _rows(human_gene_info) if row.get("#tax_id") == "9606"]
    human_id_to_symbol = {row["GeneID"]: row["Symbol"] for row in human_rows if row.get("GeneID") and row.get("Symbol")}
    human_protein_coding = {row["Symbol"] for row in human_rows if row.get("type_of_gene") == "protein-coding" and row.get("Symbol")}
    mouse_symbol_to_id = {row["Symbol"]: row["GeneID"] for row in _rows(mouse_gene_info) if row.get("#tax_id") == "10090" and row.get("GeneID") and row.get("Symbol")}
    orthologs = {
        row["Other_GeneID"]: row["GeneID"]
        for row in _rows(gene_orthologs)
        if row.get("#tax_id") == "9606" and row.get("Other_tax_id") == "10090" and row.get("GeneID") and row.get("Other_GeneID")
    }
    return human_protein_coding, {
        symbol: human_id_to_symbol[orthologs[gene_id]]
        for symbol, gene_id in mouse_symbol_to_id.items()
        if gene_id in orthologs and orthologs[gene_id] in human_id_to_symbol
    }


def _source_memberships(path: Path) -> Iterable[tuple[str, str, str]]:
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"invalid GMT row {path}:{line_number}")
            try:
                source_term, direction = fields[0].rsplit(" ", 1)
            except ValueError as exc:
                raise ValueError(f"RummaGEO GMT name must end in a direction at {path}:{line_number}") from exc
            if direction not in {"up", "dn"}:
                raise ValueError(f"unexpected RummaGEO direction {direction!r} at {path}:{line_number}")
            for gene in fields[2:]:
                if gene:
                    yield source_term, direction, gene


def _normalized_term(row: dict[str, str], model_family: str) -> str:
    explicit = row.get("normalized_term", "").strip()
    if explicit:
        return explicit.replace("/", "-").replace(",", " ")
    required = ("gse", "search_term", "condition_1", "condition_2", "species")
    missing = [name for name in required if not row.get(name, "").strip()]
    if missing:
        raise ValueError("selection manifest row requires normalized_term or " + ", ".join(missing))
    context = row.get("context", "").strip()
    expression = row.get("expression", "").strip()
    pieces = [row["gse"].replace(",", "_"), row["search_term"]]
    if model_family == "gene_perturbations" and expression:
        pieces.append(expression)
    if context:
        pieces.append(context)
    pieces.extend([f"{row['condition_1']}_v_{row['condition_2']}", row["species"]])
    return "_".join(pieces).replace("/", "-").replace(",", " ")


def _selection(path: Path, model_id: str, model_family: str) -> dict[str, tuple[str, bool]]:
    selected: dict[str, tuple[str, bool]] = {}
    for row in _rows(path):
        source_term = row.get("source_term", "").strip()
        if not source_term:
            raise ValueError("selection manifest requires source_term")
        if row.get("model_id", "").strip() != model_id or row.get("status", "").strip() not in VALID_STATUSES:
            continue
        if model_family == "drug_perturbations" and row.get("search_term", "").strip() in DRUG_FALSE_POSITIVES:
            continue
        if source_term in selected:
            raise ValueError(f"selection manifest has duplicate selected source_term: {source_term}")
        selected[source_term] = (_normalized_term(row, model_family), row["status"].strip() == "reversed")
    if not selected:
        raise ValueError(f"selection manifest selected no {model_id} records")
    return selected


def _source_records(path: Path, roles: list[str]) -> dict[str, dict[str, object]]:
    """Validate pinned source locations/version labels before emitting metadata."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    records = payload.get("sources") if isinstance(payload, dict) else None
    if not isinstance(records, dict):
        raise ValueError("source manifest must be a JSON object with a sources object")
    result: dict[str, dict[str, object]] = {}
    for role in roles:
        record = records.get(role)
        if not isinstance(record, dict) or not str(record.get("url", "")).strip() or not str(record.get("version", "")).strip():
            raise ValueError(f"source manifest requires sources.{role}.url and sources.{role}.version")
        result[role] = record
    return result


def _validate_against_legacy(generated: dict[str, set[str]], legacy_gmt: Path) -> dict[str, object]:
    legacy: dict[str, set[str]] = {}
    with _open_text(legacy_gmt) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                legacy[fields[0]] = set(filter(None, fields[2:]))
    matched = sorted(set(generated) & set(legacy))
    intersection = sum(len(generated[name] & legacy[name]) for name in matched)
    generated_total = sum(len(generated[name]) for name in matched)
    legacy_total = sum(len(legacy[name]) for name in matched)
    return {
        "legacy_term_count": len(legacy), "generated_term_count": len(generated), "normalized_term_matches": len(matched),
        "membership_precision": intersection / generated_total if generated_total else 0.0,
        "membership_recall": intersection / legacy_total if legacy_total else 0.0,
        "membership_jaccard": intersection / (generated_total + legacy_total - intersection) if generated_total + legacy_total - intersection else 0.0,
        "per_set_jaccard": {name: len(generated[name] & legacy[name]) / len(generated[name] | legacy[name]) for name in matched},
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    activate_runtime_context("rumma_geo", getattr(args, "provenance_overlay_json", None))
    model_id = args.model_id
    if model_id not in MODEL_FAMILIES:
        raise ValueError(f"unsupported RummaGEO model_id: {model_id}")
    model_family = MODEL_FAMILIES[model_id]
    paths = [Path(getattr(args, name)).resolve() for name in ("human_gmt", "mouse_gmt", "selection_manifest", "human_gene_info", "mouse_gene_info", "gene_orthologs")]
    if any(not path.is_file() for path in paths):
        raise FileNotFoundError("all RummaGEO inputs must exist")
    roles = ["human_rummageo_gmt", "mouse_rummageo_gmt", "recorded_selection_manifest", "ncbi_human_gene_info", "ncbi_mouse_gene_info", "ncbi_gene_orthologs"]
    source_manifest = Path(args.source_manifest).resolve()
    if not source_manifest.is_file():
        raise FileNotFoundError(f"missing RummaGEO source manifest: {source_manifest}")
    source_records = _source_records(source_manifest, roles)
    human_pc, mouse_to_human = _mappings(paths[3], paths[4], paths[5])
    selected = _selection(paths[2], model_id, model_family)
    raw: list[tuple[str, str, str]] = []
    source_count = 0
    for source_path in paths[:2]:
        for source_term, direction, gene in _source_memberships(source_path):
            if source_term not in selected:
                continue
            source_count += 1
            term, reversed_signature = selected[source_term]
            # The published notebook's expression accidentally mapped both branches
            # to dn.  The standardized RummaGEO GMT reverses the sign correctly;
            # preserve that production behavior rather than perpetuating the bug.
            if reversed_signature:
                direction = "dn" if direction == "up" else "up"
            gene = mouse_to_human.get(gene, gene)
            if gene in human_pc:
                raw.append((term, direction, gene))
    duplicate_counts = Counter((term, gene) for term, _, gene in raw)
    filtered = [(term, direction, gene) for term, direction, gene in raw if duplicate_counts[(term, gene)] == 1]
    grouped: dict[tuple[str, str], set[str]] = defaultdict(set)
    for term, direction, gene in filtered:
        grouped[(term, direction)].add(gene)
    emitted = {f"{term}_{direction}": genes for (term, direction), genes in grouped.items() if len(genes) >= args.min_genes}
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "genesets.gmt").open("w", encoding="utf-8", newline="\n") as handle, (out_dir / "geneset.tsv").open("w", encoding="utf-8", newline="") as table, (out_dir / "signature_summary.tsv").open("w", encoding="utf-8", newline="") as summary:
        writer = csv.DictWriter(table, fieldnames=["term", "direction", "gene_id", "gene_symbol", "score", "signed_score", "sign", "rank"], delimiter="\t", lineterminator="\n")
        summary_writer = csv.DictWriter(summary, fieldnames=["set_name", "gene_count"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); summary_writer.writeheader()
        for name in sorted(emitted):
            genes = sorted(emitted[name])
            handle.write("\t".join([name, args.gmt_description, *genes]) + "\n")
            summary_writer.writerow({"set_name": name, "gene_count": len(genes)})
            term, direction = name.rsplit("_", 1)
            sign = 1 if direction == "up" else -1
            for rank, gene in enumerate(genes, start=1):
                writer.writerow({"term": term, "direction": direction, "gene_id": gene, "gene_symbol": gene, "score": 1, "signed_score": sign, "sign": sign, "rank": rank})
    diagnostics: dict[str, object] = {"model_id": model_id, "model_family": model_family, "n_selected_source_terms": len(selected), "n_source_memberships_selected": source_count, "n_memberships_after_human_filter": len(raw), "n_memberships_removed_as_duplicate_term_gene": len(raw) - len(filtered), "n_sets_before_min_genes": len(grouped), "n_sets_emitted": len(emitted)}
    if getattr(args, "legacy_gmt", None):
        diagnostics["legacy_validation"] = _validate_against_legacy(emitted, Path(args.legacy_gmt).resolve())
    (out_dir / "reconstruction_diagnostics.json").write_text(json.dumps(diagnostics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    metadata = make_metadata(converter_name="rumma_geo", parameters={"model_id": model_id, "model_family": model_family, "min_genes": args.min_genes, "selection_method": "recorded_harmonizome_query_manifest", "duplicate_membership_policy": "drop_all_duplicate_term_gene_pairs", "mouse_mapping": "NCBI_mouse_symbol_to_human_ortholog_symbol", "human_filter": "NCBI_type_of_gene_protein-coding_exact_symbol"}, data_type="transcriptomics", assay="rna_seq", organism="human", genome_build=args.genome_build, files=[input_file_record(path, role, resource_record=source_records[role]) for path, role in zip(paths, roles)] + [input_file_record(source_manifest, "input_provenance_manifest")], gene_annotation={"mode": "exact_symbol_filter", "source": "NCBI Gene", "gene_id_field": "gene_symbol", "synonym_rescue": False}, weights={"weight_type": "signed_unweighted", "normalization": {"method": "none"}, "aggregation": "RummaGEO notebook ternary membership"}, summary={"n_input_features": source_count, "n_genes": len(set().union(*emitted.values())) if emitted else 0, "n_features_assigned": len(filtered), "fraction_features_assigned": len(filtered) / source_count if source_count else 0.0, "n_gene_sets": len(emitted), **diagnostics}, output_files=[{"path": "genesets.gmt", "role": "gmt_library"}, {"path": "geneset.tsv", "role": "selected_program"}, {"path": "signature_summary.tsv", "role": "signature_summary"}, {"path": "reconstruction_diagnostics.json", "role": "reconstruction_diagnostics"}, {"path": "geneset.meta.json", "role": "metadata_json"}], gmt={"written": True, "path": "genesets.gmt", "prefer_symbol": True, "min_genes": args.min_genes, "max_genes": None, "plans": [{"name": "historical_rummageo_notebook", "method": "signed_ternary_membership", "parameters": {"deterministic_sort": "term then gene"}, "n_genes_emitted": sum(map(len, emitted.values())), "token_type": "gene_symbol", "n_duplicates_dropped": len(raw) - len(filtered)}]}, gene_set_description=f"RummaGEO {model_family.replace('_', ' ')} signatures reconstructed with model {model_id} from recorded source and NCBI mapping snapshots.", provenance_mirror_local_prefix=getattr(args, "provenance_mirror_local_prefix", None), provenance_mirror_remote_prefix=getattr(args, "provenance_mirror_remote_prefix", None))
    write_metadata(out_dir / "geneset.meta.json", metadata)
    return {"n_peaks": len(filtered), "n_genes": len(set().union(*emitted.values())) if emitted else 0, "n_sets": len(emitted), "out_dir": str(out_dir)}
