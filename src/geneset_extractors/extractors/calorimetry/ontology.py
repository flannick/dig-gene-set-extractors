from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from importlib.resources import files
from pathlib import Path
from typing import Iterable


def _open_text(path):
    if isinstance(path, (str, Path)):
        p = Path(path)
        if p.suffix.lower() == ".gz":
            return gzip.open(p, "rt", encoding="utf-8")
        return p.open("r", encoding="utf-8")
    name = str(path)
    if name.lower().endswith(".gz"):
        return gzip.open(path.open("rb"), "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_term_templates_tsv(path: str | Path) -> tuple[list[dict[str, str]], dict[str, object]]:
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Term template table missing header: {path}")
        rows = [{str(k): str(v or "") for k, v in row.items()} for row in reader]
    return rows, {"path": str(path), "n_rows": len(rows)}


def read_phenotype_gene_edges_tsv(path: str | Path) -> tuple[list[dict[str, str]], dict[str, object]]:
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Phenotype gene edge table missing header: {path}")
        rows = [{str(k): str(v or "") for k, v in row.items()} for row in reader]
    return rows, {"path": str(path), "n_rows": len(rows)}


def read_term_hierarchy_tsv(path: str | Path | None) -> tuple[list[dict[str, str]], dict[str, object] | None]:
    if not path:
        return [], None
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Term hierarchy table missing header: {path}")
        rows = [{str(k): str(v or "") for k, v in row.items()} for row in reader]
    return rows, {"path": str(path), "n_rows": len(rows)}


def build_term_template_index(rows: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    out: dict[str, dict[str, object]] = {}
    for row in rows:
        term_id = str(row.get("term_id", row.get("mp_id", ""))).strip()
        feature = str(row.get("feature", "")).strip()
        if not term_id or not feature:
            continue
        entry = out.setdefault(
            term_id,
            {
                "term_id": term_id,
                "term_label": str(row.get("term_label", row.get("name", term_id))).strip() or term_id,
                "program": str(row.get("program", "global")).strip() or "global",
                "weights": {},
            },
        )
        entry["weights"][feature] = float(row.get("weight", 0.0) or 0.0)
    return out


def build_gene_edge_index(rows: list[dict[str, str]]) -> dict[str, list[dict[str, object]]]:
    out: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        term_id = str(row.get("term_id", row.get("mp_id", ""))).strip()
        gene_id = str(row.get("gene_id", row.get("marker_symbol", ""))).strip()
        gene_symbol = str(row.get("gene_symbol", row.get("marker_symbol", gene_id))).strip() or gene_id
        if not term_id or not gene_id:
            continue
        out[term_id].append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "weight": float(row.get("weight", 1.0) or 1.0),
            }
        )
    return dict(out)


def build_term_neighbors(rows: list[dict[str, str]]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        parent = str(row.get("parent_term_id", row.get("parent", ""))).strip()
        child = str(row.get("child_term_id", row.get("child", ""))).strip()
        if not parent or not child:
            continue
        out[parent].add(child)
        out[child].add(parent)
    return {key: set(value) for key, value in out.items()}


def score_terms(
    *,
    feature_scores: dict[str, float],
    templates: dict[str, dict[str, object]],
    program_name: str,
) -> list[dict[str, object]]:
    scored: list[dict[str, object]] = []
    for term_id, entry in templates.items():
        if str(entry.get("program", "global")) not in {"global", program_name}:
            continue
        weights = entry.get("weights", {})
        if not isinstance(weights, dict):
            continue
        total_abs = sum(abs(float(weight)) for weight in weights.values())
        if total_abs <= 0.0:
            continue
        raw = 0.0
        supporting = 0
        for feature_name, weight in weights.items():
            if feature_name not in feature_scores:
                continue
            raw += float(weight) * float(feature_scores[feature_name])
            supporting += 1
        if supporting == 0:
            continue
        normalized = max(0.0, float(raw) / float(total_abs))
        scored.append(
            {
                "term_id": term_id,
                "term_label": str(entry.get("term_label", term_id)),
                "score": normalized,
                "raw_score": raw,
                "supporting_features": supporting,
            }
        )
    scored.sort(key=lambda row: (-float(row["score"]), str(row["term_id"])))
    return scored


def route_terms_to_genes(
    *,
    scored_terms: list[dict[str, object]],
    gene_edges: dict[str, list[dict[str, object]]],
    term_neighbors: dict[str, set[str]] | None,
    mode: str,
) -> tuple[dict[str, float], dict[str, str], dict[str, object]]:
    gene_scores: dict[str, float] = defaultdict(float)
    gene_symbols: dict[str, str] = {}
    matched_terms = 0
    expanded_terms = 0
    for term in scored_terms:
        base_score = float(term.get("score", 0.0) or 0.0)
        if base_score <= 0.0:
            continue
        term_id = str(term.get("term_id", ""))
        targets = list(gene_edges.get(term_id, []))
        if mode == "expanded" and term_neighbors:
            for neighbor in sorted(term_neighbors.get(term_id, set())):
                for edge in gene_edges.get(neighbor, []):
                    targets.append({**edge, "weight": float(edge.get("weight", 1.0) or 1.0) * 0.5})
                expanded_terms += 1
        if not targets:
            continue
        matched_terms += 1
        degree_penalty = 1.0 / max(1.0, float(len(gene_edges.get(term_id, [])) or 1))
        for edge in targets:
            gene_id = str(edge.get("gene_id", "")).strip()
            gene_symbol = str(edge.get("gene_symbol", gene_id)).strip() or gene_id
            if not gene_id:
                continue
            edge_weight = float(edge.get("weight", 1.0) or 1.0)
            gene_scores[gene_id] += base_score * degree_penalty * edge_weight
            gene_symbols[gene_id] = gene_symbol
    return dict(gene_scores), gene_symbols, {
        "n_scored_terms": len(scored_terms),
        "n_terms_with_gene_support": matched_terms,
        "n_expanded_neighbor_terms_used": expanded_terms,
    }


def top_terms(scored_terms: list[dict[str, object]], limit: int = 10) -> list[dict[str, object]]:
    return [dict(row) for row in scored_terms[: max(0, int(limit))]]


def packaged_resource_path(filename: str):
    return files("geneset_extractors.resources").joinpath(filename)


def read_packaged_term_templates_tsv(organism: str = "mouse") -> tuple[list[dict[str, str]], dict[str, object]]:
    if str(organism).strip().lower() != "mouse":
        raise ValueError("Packaged calorimetry term templates are currently mouse-first only; provide --term_templates_tsv for other organisms.")
    return read_term_templates_tsv(packaged_resource_path("calr_term_templates_mouse_v1.tsv"))


def read_packaged_phenotype_gene_edges_tsv(organism: str = "mouse") -> tuple[list[dict[str, str]], dict[str, object]]:
    if str(organism).strip().lower() != "mouse":
        raise ValueError(
            "Packaged calorimetry phenotype-gene edges are currently mouse-first only; provide --phenotype_gene_edges_tsv for other organisms."
        )
    return read_phenotype_gene_edges_tsv(packaged_resource_path("calr_phenotype_gene_edges_mouse_v1.tsv"))


def read_packaged_term_hierarchy_tsv(organism: str = "mouse") -> tuple[list[dict[str, str]], dict[str, object] | None]:
    if str(organism).strip().lower() != "mouse":
        return [], None
    return read_term_hierarchy_tsv(packaged_resource_path("calr_term_hierarchy_mouse_v1.tsv"))
