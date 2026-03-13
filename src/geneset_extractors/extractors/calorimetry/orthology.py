from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from importlib.resources import files
from pathlib import Path
from typing import Iterable

ORTHOLOG_RESOURCE_ID = "calorimetry_mouse_human_orthologs_v1"


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


def packaged_resource_path(filename: str):
    return files("geneset_extractors.resources").joinpath(filename)


def ortholog_resource_record(*, path: str, method: str) -> dict[str, str]:
    return {
        "id": ORTHOLOG_RESOURCE_ID,
        "method": method,
        "path": str(path),
        "provider": "packaged_default_or_bundle",
        "stable_id": f"package:{ORTHOLOG_RESOURCE_ID}",
        "version": "v1",
        "sha256": "",
        "license": "Derived from MGI HOM_MouseHumanSequence.rpt; provenance documented in calorimetry assay docs",
    }


def read_mouse_human_orthologs_tsv(path: str | Path):
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Mouse-human ortholog table missing header: {path}")
        rows = [{str(k): str(v or "") for k, v in row.items()} for row in reader]
    by_mouse_symbol: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_mouse_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    unique_mouse_symbols: set[str] = set()
    for row in rows:
        mouse_symbol = str(row.get("mouse_gene_symbol", "")).strip()
        mouse_id = str(row.get("mouse_entrez_id", row.get("mouse_mgi_id", ""))).strip()
        if mouse_symbol:
            by_mouse_symbol[mouse_symbol].append(row)
            unique_mouse_symbols.add(mouse_symbol)
        if mouse_id:
            by_mouse_id[mouse_id].append(row)
    return {
        "path": str(path),
        "rows": rows,
        "by_mouse_symbol": dict(by_mouse_symbol),
        "by_mouse_id": dict(by_mouse_id),
        "summary": {"path": str(path), "n_rows": len(rows), "n_mouse_symbols": len(unique_mouse_symbols)},
    }


def read_packaged_mouse_human_orthologs_tsv():
    return read_mouse_human_orthologs_tsv(packaged_resource_path("calr_mouse_human_orthologs_v1.tsv.gz"))


def _preferred_human_targets(entries: list[dict[str, str]], policy: str) -> list[dict[str, str]]:
    if not entries:
        return []
    unique_targets = {}
    for row in entries:
        symbol = str(row.get("human_gene_symbol", "")).strip()
        if not symbol:
            continue
        unique_targets[symbol] = row
    if policy == "expand_all":
        return [unique_targets[key] for key in sorted(unique_targets)]
    if len(unique_targets) == 1:
        return [next(iter(unique_targets.values()))]
    return []


def map_mouse_gene_to_human(
    *,
    mouse_gene_id: str,
    mouse_gene_symbol: str,
    orthologs: dict[str, object] | None,
    policy: str,
) -> tuple[list[dict[str, str]], str]:
    if orthologs is None:
        return [], "no_ortholog_resource"
    entries: list[dict[str, str]] = []
    if mouse_gene_symbol:
        entries = list(orthologs.get("by_mouse_symbol", {}).get(mouse_gene_symbol, []))
    if not entries and mouse_gene_id:
        entries = list(orthologs.get("by_mouse_id", {}).get(mouse_gene_id, []))
    if not entries:
        return [], "unmapped"
    preferred = _preferred_human_targets(entries, policy)
    if preferred:
        return preferred, "mapped"
    return [], "ambiguous"


def map_gene_scores_to_human(
    *,
    gene_scores: dict[str, float],
    gene_symbols: dict[str, str],
    orthologs: dict[str, object] | None,
    policy: str = "unique_only",
    fallback_to_source: bool = False,
) -> tuple[dict[str, float], dict[str, str], dict[str, object]]:
    mapped_scores: dict[str, float] = defaultdict(float)
    mapped_symbols: dict[str, str] = {}
    n_mapped = 0
    n_unmapped = 0
    n_ambiguous = 0
    for mouse_gene_id, score in gene_scores.items():
        mouse_symbol = str(gene_symbols.get(mouse_gene_id, mouse_gene_id)).strip()
        targets, status = map_mouse_gene_to_human(
            mouse_gene_id=str(mouse_gene_id),
            mouse_gene_symbol=mouse_symbol,
            orthologs=orthologs,
            policy=policy,
        )
        if targets:
            n_mapped += 1
            weight = float(score) / float(len(targets))
            for target in targets:
                human_symbol = str(target.get("human_gene_symbol", "")).strip()
                human_gene_id = human_symbol or str(target.get("human_hgnc_id", "")).strip()
                human_symbol = human_symbol or human_gene_id
                if not human_gene_id:
                    continue
                mapped_scores[human_gene_id] += weight
                mapped_symbols[human_gene_id] = human_symbol
            continue
        if status == "ambiguous":
            n_ambiguous += 1
        else:
            n_unmapped += 1
        if fallback_to_source:
            source_id = str(mouse_gene_id)
            mapped_scores[source_id] += float(score)
            mapped_symbols[source_id] = mouse_symbol or source_id
    summary = {
        "output_gene_species": "human",
        "ortholog_policy": policy,
        "ortholog_resource": None if orthologs is None else orthologs.get("summary", {}),
        "n_source_genes": len(gene_scores),
        "n_mapped_source_genes": n_mapped,
        "n_unmapped_source_genes": n_unmapped,
        "n_ambiguous_source_genes": n_ambiguous,
        "fallback_to_source": fallback_to_source,
    }
    return dict(mapped_scores), mapped_symbols, summary
