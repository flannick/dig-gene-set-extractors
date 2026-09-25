"""Write DAPPER-ID-1 provenance YAML from DIG's legacy provenance graph.

The structural conversion is adapted from DAPPER's
``schema/converter/geneset_to_dapper.py``. The identifier implementation is
the DAPPER-ID-1 algorithm for the DAPPER node classes emitted here, pinned to
DAPPER release ``0.2.0-a1`` at revision
``c0cfce549baded068aa9e82f81a1400025614b51``.
"""
from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import yaml
from rdflib import Graph, Literal, URIRef
from rdflib.compare import to_canonical_graph


_INPUT_LABELS = {"data input", "metadata input"}
_EDGE_ROLE = {"data input": "data_input", "metadata input": "metadata_input"}
_OUTPUT_LABEL = "data output"
_BUCKETS = (
    "files",
    "c2m2_files",
    "activities",
    "gene_sets",
    "gene_set_collections",
    "used_edges",
    "was_generated_by_edges",
)
DAPPER_RELEASE = "0.2.0-a1"
DAPPER_SCHEMA_REVISION = "c0cfce549baded068aa9e82f81a1400025614b51"
_CLASS_BY_BUCKET = {
    "files": "File",
    "c2m2_files": "C2M2File",
    "activities": "Activity",
    "gene_sets": "GeneSet",
    "gene_set_collections": "GeneSetCollection",
}
_HASHABLE_SLOTS = {
    "File": ("description", "filename", "md5", "mime_type", "name", "sha256", "size_in_bytes"),
    "C2M2File": ("c2m2_uuid", "dcc_url", "description", "drc_url", "filename", "local_id", "md5", "mime_type", "name", "persistent_id", "sha256", "size_in_bytes"),
    "Activity": ("activity_type", "code_version", "command", "container_image", "dcc_url", "description", "drc_url", "entrypoint", "has_agentic_workspace", "has_lineage_step", "name", "observed_command", "repo_url", "script_url", "software_name", "software_version"),
    "GeneSet": ("access_level", "alternate_identifier", "assay", "controlled_access", "data_type", "dcc_url", "description", "drc_url", "funded_by", "genome_build", "has_contributor", "has_creator", "has_gmt_file", "has_license", "has_recommended_citation", "is_described_by", "member_type", "members", "n_genes", "n_members", "n_sets", "name", "organism", "term", "term_prefix", "was_attributed_to", "was_derived_from", "was_generated_by"),
    "GeneSetCollection": ("access_level", "alternate_identifier", "assay", "controlled_access", "data_type", "dcc_url", "description", "drc_url", "funded_by", "genome_build", "has_contributor", "has_creator", "has_gmt_file", "has_license", "has_recommended_citation", "is_described_by", "member_type", "members", "n_genes", "n_members", "n_sets", "name", "organism", "term_prefix", "was_attributed_to", "was_derived_from", "was_generated_by"),
}
_RELATIONSHIP_SLOTS = {
    "Activity": {"has_agentic_workspace", "has_lineage_step"},
    "GeneSet": {
        "funded_by",
        "has_contributor",
        "has_creator",
        "has_gmt_file",
        "has_license",
        "has_recommended_citation",
        "is_described_by",
        "members",
        "in_gene_set_collection",
        "in_gmt_file",
        "was_attributed_to",
        "was_derived_from",
        "was_generated_by",
    },
    "GeneSetCollection": {
        "funded_by",
        "has_contributor",
        "has_creator",
        "has_gmt_file",
        "has_license",
        "has_recommended_citation",
        "is_described_by",
        "members",
        "was_attributed_to",
        "was_derived_from",
        "was_generated_by",
    },
}
# Some DAPPER references are deliberately unhashable representation or inverse
# links. They must be rewritten, but must not become minting dependencies: a
# row's inverse collection link would otherwise form a circular hash with the
# collection's identity-bearing ``members`` list.
_HASHABLE_RELATIONSHIP_SLOTS = {
    class_name: slots - {"in_gene_set_collection", "in_gmt_file"}
    for class_name, slots in _RELATIONSHIP_SLOTS.items()
}
_SELF = URIRef("urn:dapper:self")
_CLASS_PREDICATE = URIRef("urn:dapper:class")
_SLOT_NAMESPACE = "urn:dapper:slot:"
_SELF_REFERENCE = "urn:dapper:self-reference"


def _sha512t24u(payload: bytes) -> str:
    return base64.urlsafe_b64encode(hashlib.sha512(payload).digest()[:24]).decode("ascii")


def _dapper_digest(value: Any) -> str | None:
    if not isinstance(value, str) or not value.startswith("dapper:"):
        return None
    _, _, suffix = value.partition(":")
    _, dot, digest = suffix.partition(".")
    return digest if dot and digest else None


def _replace_self(value: Any, identifier: str | None, *, substring: bool) -> Any:
    """Apply DAPPER 0.2.0-a1's slot-aware self-reference rule."""
    if not identifier:
        return value
    if isinstance(value, str):
        return value.replace(identifier, _SELF_REFERENCE) if substring else (
            _SELF_REFERENCE if value == identifier else value
        )
    if isinstance(value, list):
        return [_replace_self(item, identifier, substring=substring) for item in value]
    if isinstance(value, dict):
        return {key: _replace_self(item, identifier, substring=substring) for key, item in value.items()}
    return value


def _literal(value: Any) -> Literal:
    digest = _dapper_digest(value)
    if digest is not None:
        value = digest
    if isinstance(value, (dict, list)):
        value = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return Literal(value)


def _compute_id(node: dict[str, Any], class_name: str, identifier: str | None) -> str:
    """Compute DAPPER-ID-1 exactly for DIG's emitted DAPPER node classes."""
    graph = Graph()
    graph.add((_SELF, _CLASS_PREDICATE, Literal(class_name)))
    for slot in _HASHABLE_SLOTS[class_name]:
        raw_value = node.get(slot)
        if slot in _RELATIONSHIP_SLOTS.get(class_name, set()) or isinstance(raw_value, (list, dict)):
            value = _replace_self(raw_value, identifier, substring=False)
        elif isinstance(raw_value, str) and (" " in raw_value or "\n" in raw_value):
            value = _replace_self(raw_value, identifier, substring=True)
        else:
            # DAPPER 0.2.0-a1 deliberately leaves literal scalars alone. A
            # scalar external identifier must not be mistaken for a reference
            # merely because it equals the node's pre-mint identifier.
            value = raw_value
        if value is None:
            continue
        if isinstance(value, list):
            for index, item in enumerate(value):
                if item is not None:
                    graph.add((_SELF, URIRef(f"{_SLOT_NAMESPACE}{slot}[{index}]"), _literal(item)))
        else:
            graph.add((_SELF, URIRef(f"{_SLOT_NAMESPACE}{slot}"), _literal(value)))
    canonical = to_canonical_graph(graph)
    lines = sorted(line for line in canonical.serialize(format="nt").split("\n") if line.strip())
    digest = _sha512t24u("\n".join(lines).encode("utf-8"))
    return f"dapper:{class_name}.{digest}"


def _rewrite_substrings(value: Any, identifiers: dict[str, str]) -> Any:
    if isinstance(value, str):
        for old, new in sorted(identifiers.items(), key=lambda item: len(item[0]), reverse=True):
            value = value.replace(old, new)
        return value
    if isinstance(value, list):
        return [_rewrite_substrings(item, identifiers) for item in value]
    if isinstance(value, dict):
        return {key: _rewrite_substrings(item, identifiers) for key, item in value.items()}
    return value


def _rewrite_node(
    node: dict[str, Any], class_name: str, identifiers: dict[str, str]
) -> dict[str, Any]:
    """Mirror DAPPER 0.2.0-a1's slot-aware rewrite rule."""
    rewritten: dict[str, Any] = {}
    for key, value in node.items():
        if key == "id":
            rewritten[key] = identifiers.get(value, value)
        elif key in _RELATIONSHIP_SLOTS.get(class_name, set()):
            rewritten[key] = _rewrite_exact(value, identifiers)
        elif isinstance(value, str) and (" " in value or "\n" in value):
            rewritten[key] = _rewrite_substrings(value, identifiers)
        elif isinstance(value, (list, dict)):
            # DIG's emitted DAPPER node slots are literals, not relationships;
            # list/dict members therefore use exact rather than substring rewrites.
            rewritten[key] = _rewrite_exact(value, identifiers)
        else:
            rewritten[key] = value
    return rewritten


def _rewrite_exact(value: Any, identifiers: dict[str, str]) -> Any:
    if isinstance(value, str):
        return identifiers.get(value, value)
    if isinstance(value, list):
        return [_rewrite_exact(item, identifiers) for item in value]
    if isinstance(value, dict):
        return {key: _rewrite_exact(item, identifiers) for key, item in value.items()}
    return value


def _local_relationship_ids(value: Any, known_ids: set[str]) -> set[str]:
    """Return local DAPPER-node references carried by a relationship slot."""
    if isinstance(value, str):
        return {value} if value in known_ids else set()
    if isinstance(value, list):
        return set().union(*(_local_relationship_ids(item, known_ids) for item in value))
    if isinstance(value, dict):
        return set().union(*(_local_relationship_ids(item, known_ids) for item in value.values()))
    return set()


def _mint_dapper_ids(document: dict[str, list[dict[str, Any]]]) -> None:
    """Replace DIG node IDs and every edge reference with DAPPER-ID-1 IDs.

    DAPPER 0.2.0-a1 makes a collection's ``has_gmt_file`` identity-bearing.
    Mint dependency nodes first, so the collection digests the minted file ID,
    exactly as DAPPER-ID-1 specifies.
    """
    nodes: dict[str, tuple[str, dict[str, Any]]] = {}
    for bucket, class_name in _CLASS_BY_BUCKET.items():
        for node in document.get(bucket, []):
            old_id = node.get("id")
            if isinstance(old_id, str):
                nodes[old_id] = (class_name, node)

    identifiers: dict[str, str] = {}
    state: dict[str, int] = {}

    def mint(old_id: str, trail: tuple[str, ...]) -> None:
        if state.get(old_id) == 2:
            return
        if state.get(old_id) == 1:
            cycle = " -> ".join(trail[trail.index(old_id):] + (old_id,))
            raise ValueError(f"cycle in DAPPER hashable references: {cycle}")
        state[old_id] = 1
        class_name, node = nodes[old_id]
        for slot in _HASHABLE_RELATIONSHIP_SLOTS.get(class_name, set()):
            for dependency in sorted(_local_relationship_ids(node.get(slot), set(nodes))):
                if dependency != old_id:
                    mint(dependency, trail + (old_id,))
        resolved = _rewrite_node(node, class_name, identifiers)
        identifiers[old_id] = _compute_id(resolved, class_name, old_id)
        state[old_id] = 2

    for old_id in sorted(nodes):
        mint(old_id, ())
    for bucket in _CLASS_BY_BUCKET:
        class_name = _CLASS_BY_BUCKET[bucket]
        document[bucket] = [
            _rewrite_node(node, class_name, identifiers) for node in document.get(bucket, [])
        ]
    for bucket in ("used_edges", "was_generated_by_edges"):
        # Edge fields carry whole identifiers. DAPPER 0.2.0-a1 prohibits
        # context-free substring replacement here.
        document[bucket] = [_rewrite_exact(edge, identifiers) for edge in document.get(bucket, [])]


def _clean(payload: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in payload.items() if value not in (None, "", [], {})}


def _sha256_by_local_id(metadata: dict[str, Any]) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for record in (metadata.get("input", {}) or {}).get("files", []) or []:
        if not isinstance(record, dict) or not record.get("sha256"):
            continue
        for key in ("path", "local_path"):
            value = record.get(key)
            if value:
                checksums[str(value)] = str(record["sha256"])
    return checksums


def _dapper_uri(value: object) -> str | None:
    """Keep only absolute URI values in DAPPER URI/CURIE slots.

    Legacy DIG graphs historically used local filesystem paths in C2M2's
    ``local_id``/``dcc_url``/``drc_url`` fields. DAPPER 0.2 treats those as
    identifiers, not locations.  Preserve the path in ``location`` instead
    of exporting an invalid URI-like value.
    """
    if not isinstance(value, str) or not value.strip():
        return None
    candidate = value.strip()
    return candidate if urlparse(candidate).scheme else None


def _c2m2_file(node: dict[str, Any], checksums: dict[str, str]) -> dict[str, Any]:
    c2m2 = node.get("c2m2_properties", {}) or {}
    raw_local_id = c2m2.get("local_id")
    location = node.get("location") or node.get("path") or node.get("local_path") or raw_local_id
    return _clean(
        {
            "id": node.get("id"),
            "name": node.get("name"),
            "description": node.get("description"),
            "filename": c2m2.get("filename"),
            "persistent_id": c2m2.get("persistent_id"),
            "local_id": _dapper_uri(raw_local_id),
            "c2m2_uuid": c2m2.get("_uuid"),
            "md5": c2m2.get("md5"),
            "sha256": checksums.get(str(raw_local_id)),
            "size_in_bytes": c2m2.get("size_in_bytes"),
            "dcc_url": _dapper_uri(node.get("dcc_url")),
            "drc_url": _dapper_uri(node.get("drc_url")),
            "location": location,
        }
    )


def _file(node: dict[str, Any], checksums: dict[str, str]) -> dict[str, Any]:
    """Map a DIG file node without C2M2 properties to DAPPER's generic File."""
    location = node.get("location") or node.get("path") or node.get("local_path")
    return _clean(
        {
            "id": node.get("id"),
            "name": node.get("name"),
            "description": node.get("description"),
            "filename": node.get("filename"),
            "md5": node.get("md5"),
            "sha256": node.get("sha256") or checksums.get(str(location)),
            "size_in_bytes": node.get("size_in_bytes") or node.get("size_bytes"),
            "mime_type": node.get("mime_type"),
            # Location and DRS registration are intentionally not identity
            # bearing in DAPPER 0.2.0-a1, but remain useful retrieval metadata.
            "location": location,
            "drs_representation": node.get("drs_representation"),
        }
    )


def _activity(node: dict[str, Any]) -> dict[str, Any]:
    analysis = node.get("analysis", {}) or {}
    environment = analysis.get("environment", {}) or {}
    c2m2 = node.get("c2m2_properties", {}) or {}
    result = _clean(
        {
            "id": node.get("id"),
            "name": node.get("name"),
            "description": node.get("description"),
            "command": analysis.get("command"),
            "observed_command": analysis.get("observed_command"),
            "script_url": analysis.get("script_url"),
            "repo_url": environment.get("repo_url"),
            "code_version": analysis.get("version"),
            "entrypoint": environment.get("entrypoint"),
            "container_image": environment.get("container_image"),
            "dcc_url": _dapper_uri(node.get("dcc_url")),
            "drc_url": _dapper_uri(node.get("drc_url")),
        }
    )
    if c2m2.get("synonyms"):
        result["aliases"] = list(c2m2["synonyms"])
    return result


def _is_gene_set_collection(node: dict[str, Any], metadata: dict[str, Any]) -> bool:
    summary = metadata.get("summary", {}) or {}
    return (
        node.get("type") == "GeneSetCollection"
        or summary.get("n_sets_emitted") is not None
        or node.get("n_sets") is not None
    )


def _gene_set(node: dict[str, Any], metadata: dict[str, Any]) -> dict[str, Any]:
    gene_set = metadata.get("gene_set", {}) or {}
    summary = metadata.get("summary", {}) or {}
    parameters = (metadata.get("converter", {}) or {}).get("parameters", {}) or {}
    is_collection = _is_gene_set_collection(node, metadata)
    return _clean(
        {
            "id": node.get("id"),
            "name": node.get("name"),
            "description": node.get("description") or gene_set.get("description"),
            "member_type": "gene_set" if is_collection else "gene",
            "assay": gene_set.get("assay"),
            "data_type": gene_set.get("data_type"),
            "organism": gene_set.get("organism"),
            "genome_build": gene_set.get("genome_build"),
            "n_genes": (
                gene_set.get("n_genes")
                if gene_set.get("n_genes") is not None
                else summary.get("n_genes")
            ),
            "n_sets": summary.get("n_sets_emitted", node.get("n_sets")),
            "term_prefix": parameters.get("term_prefix"),
            "dcc_url": node.get("dcc_url"),
            "drc_url": node.get("drc_url"),
        }
    )


def _link_collection_to_unique_gmt(document: dict[str, list[dict[str, Any]]]) -> None:
    """Link a library to a uniquely identified GMT from its own activity.

    A DAPPER collection may name a physical GMT via ``has_gmt_file``. We only
    infer that relationship when exactly one generated GMT shares the
    collection's producing activity; multiple candidates remain intentionally
    unmapped rather than guessing a filename.
    """
    generated = document["was_generated_by_edges"]
    gmt_ids = {
        node["id"]
        for bucket in ("files", "c2m2_files")
        for node in document[bucket]
        if str(node.get("filename", "")).lower().endswith(".gmt")
    }
    for collection in document["gene_set_collections"]:
        producers = {
            edge["object"]
            for edge in generated
            if edge.get("subject") == collection.get("id")
        }
        candidates = {
            edge["subject"]
            for edge in generated
            if edge.get("object") in producers and edge.get("subject") in gmt_ids
        }
        if len(candidates) == 1:
            collection["has_gmt_file"] = next(iter(candidates))


def _convert_graph(graph: dict[str, Any], metadata: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """Map one DIG ``{nodes, edges}`` graph to DAPPER document buckets."""
    document: dict[str, list[dict[str, Any]]] = {bucket: [] for bucket in _BUCKETS}
    checksums = _sha256_by_local_id(metadata)
    for node in graph.get("nodes", []):
        if not isinstance(node, dict):
            continue
        match node.get("type"):
            case "File":
                if node.get("c2m2_properties"):
                    document["c2m2_files"].append(_c2m2_file(node, checksums))
                else:
                    document["files"].append(_file(node, checksums))
            case "AnalysisType":
                document["activities"].append(_activity(node))
            case "GeneSet" | "GeneSetCollection":
                bucket = (
                    "gene_set_collections"
                    if _is_gene_set_collection(node, metadata)
                    else "gene_sets"
                )
                document[bucket].append(_gene_set(node, metadata))

    for edge in graph.get("edges", []):
        if not isinstance(edge, dict):
            continue
        label = edge.get("label")
        if label in _INPUT_LABELS:
            document["used_edges"].append(
                _clean(
                    {
                        "subject": edge.get("target"),
                        "predicate": "prov:used",
                        "object": edge.get("source"),
                        "edge_role": _EDGE_ROLE[str(label)],
                    }
                )
            )
        elif label == _OUTPUT_LABEL:
            document["was_generated_by_edges"].append(
                _clean(
                    {
                        "subject": edge.get("target"),
                        "predicate": "prov:wasGeneratedBy",
                        "object": edge.get("source"),
                    }
                )
            )
    _link_collection_to_unique_gmt(document)
    _mint_dapper_ids(document)
    return {key: value for key, value in document.items() if value}


def _dapper_export_config(metadata: dict[str, Any]) -> dict[str, Any] | None:
    """Return explicit row-export configuration, never guessing gene namespaces."""
    config = metadata.get("dapper")
    if not isinstance(config, dict):
        return None
    prefix = config.get("gene_member_prefix")
    prefix_uri = config.get("gene_member_prefix_uri")
    if not isinstance(prefix, str) or not prefix.strip():
        return None
    if not isinstance(prefix_uri, str) or not prefix_uri.strip():
        raise ValueError(
            "DAPPER row export requires dapper.gene_member_prefix_uri for "
            "the declared dapper.gene_member_prefix"
        )
    return config


def _declared_gmt_path(metadata: dict[str, Any], output_dir: Path) -> Path | None:
    """Resolve the single declared GMT used for an opt-in DAPPER row export."""
    candidates: list[Path] = []
    for record in ((metadata.get("output") or {}).get("files") or []):
        if not isinstance(record, dict):
            continue
        raw_path = record.get("path")
        if not isinstance(raw_path, str) or not raw_path.lower().endswith(".gmt"):
            continue
        path = Path(raw_path)
        path = path if path.is_absolute() else output_dir / path
        if path.exists() and path.is_file():
            candidates.append(path)
    unique = sorted({path.resolve() for path in candidates})
    if not unique:
        return None
    if len(unique) != 1:
        raise ValueError(
            "DAPPER row export requires exactly one declared existing GMT output; "
            f"found {len(unique)}"
        )
    return unique[0]


def _parse_gmt_rows(path: Path) -> list[tuple[str, list[str], bytes]]:
    """Parse a GMT while retaining every byte after the first tab-delimited field."""
    rows: list[tuple[str, list[str], bytes]] = []
    labels: set[str] = set()
    for line_number, raw_line in enumerate(path.read_bytes().splitlines(keepends=True), start=1):
        body = raw_line.rstrip(b"\r\n")
        parts = body.decode("utf-8").split("\t")
        if len(parts) < 3 or not parts[0] or any(not token for token in parts[2:]):
            raise ValueError(f"{path}: invalid GMT row {line_number}")
        label = parts[0]
        if label in labels:
            raise ValueError(f"{path}: duplicate GMT row label {label!r}")
        labels.add(label)
        rows.append((label, parts[2:], raw_line))
    if not rows:
        raise ValueError(f"{path}: no GMT rows found")
    return rows


def _base64_md5(data: bytes) -> str:
    return base64.b64encode(hashlib.md5(data).digest()).decode("ascii")


def _dapper_gmt_path(gmt_path: Path) -> Path:
    return gmt_path.with_name(f"{gmt_path.stem}.dapper-ids{gmt_path.suffix}")


def _find_gmt_file_id(document: dict[str, Any], gmt_path: Path) -> str | None:
    matches = [
        node.get("id")
        for bucket in ("files", "c2m2_files")
        for node in document.get(bucket, [])
        if node.get("filename") == gmt_path.name and isinstance(node.get("id"), str)
    ]
    return matches[0] if len(matches) == 1 else None


def _add_dapper_gmt_export(
    document: dict[str, Any], metadata: dict[str, Any], output_dir: Path
) -> None:
    """Add DAPPER's row-level GMT representation without altering the source GMT.

    The source GMT stays the extractor's compatibility artifact. The additive
    ``*.dapper-ids.gmt`` export replaces only its first column with each minted
    ``GeneSet`` ID, preserving every remaining byte on each row.
    """
    config = _dapper_export_config(metadata)
    if config is None:
        return
    gmt_path = _declared_gmt_path(metadata, output_dir)
    if gmt_path is None:
        return
    collections = document.get("gene_set_collections") or []
    if len(collections) != 1:
        raise ValueError(
            "DAPPER row export requires exactly one GeneSetCollection in the provenance graph"
        )
    source_gmt_id = _find_gmt_file_id(document, gmt_path)
    if source_gmt_id is None:
        raise ValueError(
            f"DAPPER row export could not identify one provenance File for {gmt_path.name}"
        )
    source_rows = _parse_gmt_rows(gmt_path)
    prefix = str(config["gene_member_prefix"]).strip()
    row_labels = config.get("row_display_names")
    if not isinstance(row_labels, dict):
        raise ValueError(
            "DAPPER row export requires dapper.row_display_names for every GMT row; "
            "do not reuse original GMT labels as human-readable GeneSet names"
        )
    producer_by_collection = {
        edge.get("object")
        for edge in document.get("was_generated_by_edges", [])
        if edge.get("subject") == collections[0].get("id")
    }
    if len(producer_by_collection) != 1:
        raise ValueError("DAPPER row export requires one producing Activity for the collection")
    producer_id = next(iter(producer_by_collection))

    row_nodes: list[dict[str, Any]] = []
    for index, (label, genes, _raw_line) in enumerate(source_rows, start=1):
        display_name = row_labels.get(label)
        if not isinstance(display_name, str) or not display_name.strip() or display_name == label:
            raise ValueError(
                "DAPPER row export requires a distinct non-empty readable name for "
                f"GMT label {label!r} in dapper.row_display_names"
            )
        members = [gene if ":" in gene else f"{prefix}:{gene}" for gene in genes]
        row = {
            "id": f"urn:dig:dapper-row:{index}",
            "name": display_name,
            "alternate_identifier": [label],
            "member_type": "gene",
            "members": members,
            "n_genes": len(set(members)),
            "was_generated_by": producer_id,
        }
        row["id"] = _compute_id(row, "GeneSet", row["id"])
        row_nodes.append(row)

    export_path = _dapper_gmt_path(gmt_path)
    rendered_lines = []
    for row, (_label, _genes, raw_line) in zip(row_nodes, source_rows, strict=True):
        first_tab = raw_line.find(b"\t")
        rendered_lines.append(row["id"].encode("utf-8") + raw_line[first_tab:])
    export_bytes = b"".join(rendered_lines)
    export_path.write_bytes(export_bytes)

    export_activity = {
        "id": "urn:dig:dapper-gmt-export",
        "name": "Export GMT with DAPPER GeneSet identifiers",
        "description": (
            "Replace original GMT first-column labels with minted DAPPER GeneSet identifiers; "
            "preserve descriptions, gene members, row order, and line endings."
        ),
        "entrypoint": "geneset_extractors.core.dapper_provenance",
    }
    export_activity["id"] = _compute_id(export_activity, "Activity", export_activity["id"])
    export_file = {
        "id": "urn:dig:dapper-gmt-file",
        "name": f"{gmt_path.stem} with DAPPER GeneSet identifiers",
        "description": "GMT export whose first-column names are DAPPER GeneSet identifiers.",
        "filename": export_path.name,
        "location": export_path.name,
        "md5": _base64_md5(export_bytes),
        "sha256": hashlib.sha256(export_bytes).hexdigest(),
        "size_in_bytes": len(export_bytes),
    }
    export_file["id"] = _compute_id(export_file, "File", export_file["id"])

    old_collection = collections[0]
    collection = dict(old_collection)
    collection["id"] = "urn:dig:dapper-collection"
    collection["members"] = [row["id"] for row in row_nodes]
    collection["n_members"] = len(row_nodes)
    collection["n_sets"] = len(row_nodes)
    collection["n_genes"] = len({gene for row in row_nodes for gene in row["members"]})
    collection["has_gmt_file"] = export_file["id"]
    collection["member_type"] = "gene_set"
    collection["id"] = _compute_id(collection, "GeneSetCollection", collection["id"])

    for row in row_nodes:
        row["in_gene_set_collection"] = [collection["id"]]
        row["in_gmt_file"] = export_file["id"]
        row["gmt_entry"] = row["id"]

    old_collection_id = old_collection["id"]
    document["gene_set_collections"] = [collection]
    document.setdefault("gene_sets", []).extend(row_nodes)
    document.setdefault("files", []).append(export_file)
    document.setdefault("activities", []).append(export_activity)
    for edge in document.get("was_generated_by_edges", []):
        if edge.get("subject") == old_collection_id:
            edge["subject"] = collection["id"]
    document.setdefault("used_edges", []).append(
        {"subject": export_activity["id"], "predicate": "prov:used", "object": source_gmt_id, "edge_role": "data_input"}
    )
    document.setdefault("was_generated_by_edges", []).append(
        {"subject": export_file["id"], "predicate": "prov:wasGeneratedBy", "object": export_activity["id"]}
    )
    document["prefixes"] = {
        str(config["gene_member_prefix"]): str(config["gene_member_prefix_uri"])
    }


def build_dapper_provenance(
    legacy_payload: dict[str, Any],
    metadata: dict[str, Any],
    *,
    output_dir: Path | None = None,
) -> dict[str, Any]:
    """Convert all graphs in a legacy DIG payload into one DAPPER YAML document."""
    document: dict[str, Any] = {bucket: [] for bucket in _BUCKETS}
    for graph in legacy_payload.values():
        if not isinstance(graph, dict) or "nodes" not in graph:
            continue
        for bucket, values in _convert_graph(graph, metadata).items():
            document[bucket].extend(values)
    if output_dir is not None:
        _add_dapper_gmt_export(document, metadata, output_dir)
    return {key: value for key, value in document.items() if value}


def write_dapper_provenance(
    path: str | Path, legacy_payload: dict[str, Any], metadata: dict[str, Any]
) -> Path:
    """Write a deterministic, readable DAPPER YAML sidecar."""
    output_path = Path(path)
    output_path.write_text(
        yaml.safe_dump(
            build_dapper_provenance(legacy_payload, metadata, output_dir=output_path.parent),
            sort_keys=False,
            default_flow_style=False,
            allow_unicode=True,
        ),
        encoding="utf-8",
    )
    return output_path
