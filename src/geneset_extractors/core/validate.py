from __future__ import annotations

import csv
import json
import re
from pathlib import Path


def _validate_required_fallback(payload: object, schema: dict[str, object], path: str = "$") -> None:
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: expected object")
    required = schema.get("required", [])
    if isinstance(required, list):
        for key in required:
            if key not in payload:
                raise ValueError(f"{path}: missing required key '{key}'")
    props = schema.get("properties", {})
    if not isinstance(props, dict):
        return
    for key, sub_schema in props.items():
        if key in payload and isinstance(sub_schema, dict) and sub_schema.get("type") == "object":
            _validate_required_fallback(payload[key], sub_schema, f"{path}.{key}")


def _validate_provenance_graph_fallback(payload: object, path: str = "$") -> None:
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: expected top-level provenance object")
    for graph_key, graph in payload.items():
        if not isinstance(graph, dict):
            raise ValueError(f"{path}.{graph_key}: expected provenance graph object")
        nodes = graph.get("nodes")
        edges = graph.get("edges")
        if not isinstance(nodes, list) or not nodes:
            raise ValueError(f"{path}.{graph_key}.nodes: expected non-empty array")
        if not isinstance(edges, list) or not edges:
            raise ValueError(f"{path}.{graph_key}.edges: expected non-empty array")
        for index, node in enumerate(nodes):
            if not isinstance(node, dict):
                raise ValueError(f"{path}.{graph_key}.nodes[{index}]: expected object")
            for required_key in ("id", "type", "name", "description", "dcc_url", "drc_url"):
                if required_key not in node:
                    raise ValueError(f"{path}.{graph_key}.nodes[{index}]: missing required key '{required_key}'")
        for index, edge in enumerate(edges):
            if not isinstance(edge, dict):
                raise ValueError(f"{path}.{graph_key}.edges[{index}]: expected object")
            for required_key in ("id", "source", "target", "label", "description"):
                if required_key not in edge:
                    raise ValueError(f"{path}.{graph_key}.edges[{index}]: missing required key '{required_key}'")


def validate_metadata_schema(meta_path: Path, schema_path: Path) -> None:
    payload = json.loads(meta_path.read_text(encoding="utf-8"))
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    try:
        import jsonschema  # type: ignore
    except ModuleNotFoundError:
        # Compatibility fallback when jsonschema is not yet installed.
        _validate_required_fallback(payload, schema)
    else:
        jsonschema.validate(payload, schema)


def validate_provenance_schema(provenance_path: Path, schema_path: Path) -> None:
    payload = json.loads(provenance_path.read_text(encoding="utf-8"))
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    try:
        import jsonschema  # type: ignore
    except ModuleNotFoundError:
        if "$defs" in schema and isinstance(schema.get("additionalProperties"), dict):
            _validate_provenance_graph_fallback(payload)
        else:
            _validate_required_fallback(payload, schema)
    else:
        jsonschema.validate(payload, schema)


def validate_geneset_tsv(geneset_path: Path) -> None:
    with geneset_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"gene_id", "score"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise ValueError("geneset.tsv missing required columns gene_id and score")
        seen: set[str] = set()
        for row in reader:
            gid = row["gene_id"]
            if gid in seen:
                raise ValueError(f"duplicate gene_id in geneset.tsv: {gid}")
            seen.add(gid)
            float(row["score"])
            if "weight" in row and str(row["weight"]).strip() != "":
                float(row["weight"])
            if "rank" in row and str(row["rank"]).strip() != "":
                int(float(row["rank"]))


def validate_gmt(gmt_path: Path) -> None:
    with gmt_path.open("r", encoding="utf-8") as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line:
                raise ValueError(f"{gmt_path}: line {line_no} is empty")
            parts = line.split("\t")
            if len(parts) == 2:
                name, genes_field = parts
                tokens = genes_field.split(" ")
            elif len(parts) >= 3:
                name = parts[0]
                genes_field = "\t".join(parts[2:])
                tokens = parts[2:]
            else:
                raise ValueError(f"{gmt_path}: line {line_no} must contain at least one set name and one gene token")
            if not name.strip():
                raise ValueError(f"{gmt_path}: line {line_no} has empty set name")
            if not genes_field:
                raise ValueError(f"{gmt_path}: line {line_no} has empty gene list")
            if any(tok == "" for tok in tokens):
                raise ValueError(f"{gmt_path}: line {line_no} has empty gene token")


def _resolve_manifest_path(root: Path, manifest_value: str) -> Path:
    candidate = Path(manifest_value)
    if candidate.is_absolute():
        return candidate
    joined = root / candidate
    if joined.exists():
        return joined
    return candidate


def _validate_grouped_output_dir(out_dir: Path, schema_path: Path) -> dict[str, object]:
    manifest = out_dir / "manifest.tsv"
    if not manifest.exists():
        raise FileNotFoundError("output dir must contain geneset.tsv and geneset.meta.json, or manifest.tsv for grouped outputs")
    with manifest.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "path" not in reader.fieldnames:
            raise ValueError("manifest.tsv must contain a path column")
        rows = list(reader)
    if not rows:
        raise ValueError("manifest.tsv contains no group rows")
    failures: list[str] = []
    for row in rows:
        group_path = _resolve_manifest_path(out_dir, str(row["path"]))
        try:
            _validate_single_output_dir(group_path, schema_path)
        except Exception as exc:
            failures.append(f"{group_path}: {exc}")
    if failures:
        raise ValueError("grouped validation failed: " + "; ".join(failures))
    root_gmt = out_dir / "genesets.gmt"
    if root_gmt.exists():
        validate_gmt(root_gmt)
    return {"mode": "grouped", "n_groups": len(rows)}


def _validate_single_output_dir(out: Path, schema_path: Path) -> None:
    geneset = out / "geneset.tsv"
    geneset_full = out / "geneset.full.tsv"
    gmt = out / "genesets.gmt"
    meta = out / "geneset.meta.json"
    provenance = out / "geneset.provenance.json"
    if not geneset.exists() or not meta.exists():
        raise FileNotFoundError("output dir must contain geneset.tsv and geneset.meta.json")
    validate_geneset_tsv(geneset)
    if geneset_full.exists():
        validate_geneset_tsv(geneset_full)
    if gmt.exists():
        validate_gmt(gmt)
    validate_metadata_schema(meta, schema_path)
    provenance_schema = schema_path.with_name("geneset_provenance.schema.json")
    meta_payload = json.loads(meta.read_text(encoding="utf-8"))
    if provenance.exists():
        validate_provenance_schema(provenance, provenance_schema)
    elif isinstance(meta_payload.get("provenance"), dict):
        raise FileNotFoundError("metadata references geneset.provenance.json but file is missing")


def validate_output_dir(out_dir: str | Path, schema_path: str | Path) -> dict[str, object]:
    out = Path(out_dir)
    schema = Path(schema_path)
    geneset = out / "geneset.tsv"
    meta = out / "geneset.meta.json"
    if geneset.exists() and meta.exists():
        _validate_single_output_dir(out, schema)
        return {"mode": "single", "n_groups": 1}
    return _validate_grouped_output_dir(out, schema)


# ---------------------------------------------------------------------------
# Submission-level gate: cross-file invariants over a <library>_all_models tree.
#
# Complements per-output validate_output_dir with the invariants that distinguish
# an accepted submission from a rejected one (calibrated against an accepted
# reference package): clean single-token source URIs, integrity fields on every
# file node, rerunnable commands carrying real version SHAs, meta<->gmt name
# parity, and uniform .orig snapshots. Library-specific coverage rules stay with
# the caller; these checks are framework-generic.
# ---------------------------------------------------------------------------

_SHA40 = re.compile(r"^[0-9a-f]{40}$")
_INTERNAL_NAME_TOKENS = ("deg_long__", "condition=")


def _provenance_graph(prov_payload: object) -> dict:
    if not isinstance(prov_payload, dict) or not prov_payload:
        return {"nodes": [], "edges": []}
    graph = prov_payload[next(iter(prov_payload))]
    return graph if isinstance(graph, dict) else {"nodes": [], "edges": []}


def _gmt_first_column(gmt_path: Path) -> list[str]:
    names: list[str] = []
    if not gmt_path.exists():
        return names
    with gmt_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                names.append(line.rstrip("\n").split("\t")[0])
    return names


def _is_unclean_identifier(value: object) -> bool:
    # A clean source identifier is a single whitespace-free token (URI / path /
    # filename). Prose descriptions contain whitespace; an accepted reference
    # submission has zero file nodes whose dcc_url/drc_url contains whitespace.
    return isinstance(value, str) and bool(value.strip()) and " " in value.strip()


def validate_submission_tree(
    root: str | Path,
    *,
    require_md5: bool = True,
    require_version_sha: bool = True,
    require_uniform_orig: bool = True,
) -> dict[str, object]:
    """Check submission invariants across a <library>_all_models output tree.

    Returns ``{"failures": [str, ...], "counts": {...}, "models": [...]}``.
    An empty ``failures`` list means the tree passes; callers decide whether a
    non-empty list should block bundling.
    """
    root = Path(root)
    prov_files = sorted(root.glob("genesets/**/geneset.provenance.json"))
    counts = {
        "file_nodes": 0,
        "prose_dcc_url": 0,
        "missing_md5": 0,
        "prose_commands": 0,
        "nonsha_versions": 0,
        "name_mismatch": 0,
        "dup_row_files": 0,
        "provenance_files": len(prov_files),
    }
    examples: dict[str, str] = {}
    orig_by_model: dict[str, int] = {}
    models: set[str] = set()

    for prov_path in prov_files:
        out_dir = prov_path.parent
        rel = prov_path.relative_to(root).as_posix()
        parts = rel.split("/")
        model = parts[3] if len(parts) > 3 else "?"
        models.add(model)
        graph = _provenance_graph(json.loads(prov_path.read_text(encoding="utf-8")))
        for node in graph.get("nodes", []):
            if not isinstance(node, dict):
                continue
            if node.get("type") == "File":
                counts["file_nodes"] += 1
                if _is_unclean_identifier(node.get("dcc_url")) or _is_unclean_identifier(node.get("drc_url")):
                    counts["prose_dcc_url"] += 1
                    examples.setdefault("inv1", f"{rel}: {str(node.get('dcc_url'))[:70]}")
                c2m2 = node.get("c2m2_properties")
                has_integrity = (
                    isinstance(c2m2, dict)
                    and c2m2.get("md5")
                    and c2m2.get("size_in_bytes") is not None
                )
                if require_md5 and not has_integrity:
                    counts["missing_md5"] += 1
            analysis = node.get("analysis") if isinstance(node.get("analysis"), dict) else None
            if analysis is not None:
                command = analysis.get("command")
                if isinstance(command, str) and "(inputs:" in command:
                    counts["prose_commands"] += 1
                    examples.setdefault("inv3", f"{rel}: {command[:70]}")
                version = analysis.get("version")
                if require_version_sha and isinstance(version, str) and not _SHA40.match(version):
                    counts["nonsha_versions"] += 1
                    examples.setdefault("inv3v", f"{rel}: version={version!r}")

        gmt_names = _gmt_first_column(out_dir / "genesets.gmt")
        if gmt_names and len(gmt_names) != len(set(gmt_names)):
            counts["dup_row_files"] += 1
        meta_path = out_dir / "geneset.meta.json"
        if meta_path.exists() and gmt_names:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            gmt_meta = meta.get("gmt", {}) if isinstance(meta, dict) else {}
            mismatch = False
            for field in ("emitted_outputs", "requested_outputs", "plans"):
                names = [o.get("name") for o in gmt_meta.get(field, []) if isinstance(o, dict)]
                names = [n for n in names if n]
                if names and sorted(names) != sorted(gmt_names):
                    mismatch = True
                    examples.setdefault("inv4", f"{rel}: {field}={names} vs gmt={gmt_names}")
                    break
            emitted = [o.get("name", "") for o in gmt_meta.get("emitted_outputs", []) if isinstance(o, dict)]
            if any(any(tok in (n or "") for tok in _INTERNAL_NAME_TOKENS) for n in emitted):
                mismatch = True
                examples.setdefault("inv4", f"{rel}: internal token in emitted_outputs {emitted}")
            if mismatch:
                counts["name_mismatch"] += 1

    for orig in root.glob("genesets/*/models/*/extractor/**/*.orig"):
        mdl = orig.relative_to(root).as_posix().split("/")[3]
        orig_by_model[mdl] = orig_by_model.get(mdl, 0) + 1
    models_with_orig = {m for m, n in orig_by_model.items() if n}
    orig_uniform = (not models_with_orig) or (models_with_orig == models)

    failures: list[str] = []
    if counts["prose_dcc_url"]:
        failures.append(
            f"INV-1 source identifiers: {counts['prose_dcc_url']}/{counts['file_nodes']} "
            f"file nodes carry prose/whitespace dcc_url (e.g. {examples.get('inv1', '')})"
        )
    if require_md5 and counts["missing_md5"]:
        failures.append(
            f"INV-2 integrity: {counts['missing_md5']}/{counts['file_nodes']} file nodes missing md5/size_in_bytes"
        )
    if counts["prose_commands"]:
        failures.append(
            f"INV-3 command fidelity: {counts['prose_commands']} analysis nodes with prose '(inputs:' "
            f"command (e.g. {examples.get('inv3', '')})"
        )
    if require_version_sha and counts["nonsha_versions"]:
        failures.append(
            f"INV-3 version: {counts['nonsha_versions']} analysis nodes with non-SHA version "
            f"(e.g. {examples.get('inv3v', '')})"
        )
    if counts["name_mismatch"]:
        failures.append(
            f"INV-4 name parity: {counts['name_mismatch']} outputs where meta names != gmt col1 "
            f"(e.g. {examples.get('inv4', '')})"
        )
    if counts["dup_row_files"]:
        failures.append(f"INV-5 gmt: {counts['dup_row_files']} gmt files with duplicate row names")
    if require_uniform_orig and not orig_uniform:
        failures.append(
            f"INV-7 .orig not uniform: present for {sorted(models_with_orig)} of models {sorted(models)}"
        )

    return {"failures": failures, "counts": counts, "models": sorted(models)}
