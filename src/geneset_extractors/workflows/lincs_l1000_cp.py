"""Export individual LINCS L1000 chemical-perturbation CD signatures."""
from __future__ import annotations

import csv
import math
import shutil
from pathlib import Path
from urllib.parse import unquote, urlparse
from urllib.request import urlopen

from geneset_extractors.workflows.gtex_runtime_common import write_tsv, write_workflow_provenance_graph


LIBRARY_ID = "54198d6e-fe17-5ef8-91ac-02b425761653"
FILENAME_PREFIX = "L1000_LINCS_DCIC_"


def persistent_id_to_term(persistent_id: str) -> str:
    """Convert a SigCom persistent ID (filename or URL) to a legacy GMT term."""
    value = str(persistent_id or "").strip()
    filename = Path(unquote(urlparse(value).path or value)).name
    if filename.startswith(FILENAME_PREFIX):
        filename = filename[len(FILENAME_PREFIX) :]
    if filename.endswith(".tsv"):
        filename = filename[: -len(".tsv")]
    if not filename:
        raise ValueError(f"Cannot derive a term from persistent_id={persistent_id!r}")
    return filename


def _source_url(row: dict[str, str], source_url_base: str) -> str:
    explicit = str(row.get("source_url", "")).strip()
    if explicit:
        return explicit
    persistent_id = str(row.get("persistent_id", "")).strip()
    if persistent_id.startswith(("https://", "http://")):
        return persistent_id
    if not persistent_id:
        raise ValueError("Signature manifest row has neither source_url nor persistent_id")
    return source_url_base.rstrip("/") + "/" + persistent_id


def _read_signature(path: Path, top_n: int) -> tuple[list[tuple[str, float]], list[tuple[str, float]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"symbol", "CD-coefficient"}.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required column(s): {', '.join(sorted(missing))}")
        rows: list[tuple[str, float]] = []
        symbols: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            symbol = str(row.get("symbol", "")).strip()
            if not symbol:
                raise ValueError(f"{path}:{line_number} has an empty symbol")
            if symbol in symbols:
                raise ValueError(f"{path}:{line_number} has duplicate symbol {symbol!r}")
            try:
                coefficient = float(str(row.get("CD-coefficient", "")).strip())
            except ValueError as error:
                raise ValueError(f"{path}:{line_number} has non-numeric CD-coefficient") from error
            if not math.isfinite(coefficient):
                raise ValueError(f"{path}:{line_number} has non-finite CD-coefficient")
            symbols.add(symbol)
            rows.append((symbol, coefficient))
    if len(rows) < top_n * 2:
        raise ValueError(f"{path} has {len(rows)} unique genes; need at least {top_n * 2}")
    ranked = sorted(rows, key=lambda item: (-item[1], item[0]))
    return ranked[:top_n], ranked[-top_n:]


def _materialize_source(row: dict[str, str], cache_dir: Path, source_url_base: str, timeout: int) -> Path:
    source_path = str(row.get("source_path", "")).strip()
    if source_path:
        path = Path(source_path).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Missing source_path from signature manifest: {path}")
        return path
    url = _source_url(row, source_url_base)
    filename = Path(unquote(urlparse(url).path)).name
    if not filename:
        raise ValueError(f"Cannot derive cache filename from source URL {url!r}")
    destination = cache_dir / filename
    if destination.is_file() and destination.stat().st_size > 0:
        return destination
    temporary = destination.with_suffix(destination.suffix + ".partial")
    try:
        with urlopen(url, timeout=timeout) as response, temporary.open("wb") as output:
            shutil.copyfileobj(response, output)
        temporary.replace(destination)
    finally:
        if temporary.exists():
            temporary.unlink()
    return destination


def _manifest_rows(path: Path, limit_signatures: int | None) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "persistent_id" not in rows[0]:
        raise ValueError(f"Signature manifest must be non-empty and have a persistent_id column: {path}")
    rows = sorted(rows, key=lambda row: persistent_id_to_term(row.get("persistent_id", "")))
    terms = [persistent_id_to_term(row.get("persistent_id", "")) for row in rows]
    if len(terms) != len(set(terms)):
        raise ValueError("Signature manifest resolves multiple persistent IDs to the same GMT term")
    return rows[:limit_signatures] if limit_signatures is not None else rows


def run(args) -> dict[str, object]:
    manifest = Path(args.signature_manifest_tsv).resolve()
    if not manifest.is_file():
        raise FileNotFoundError(f"Missing signature manifest: {manifest}")
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(args.cache_dir).resolve() if args.cache_dir else out_dir / "source_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    top_n = int(args.top_n)
    if top_n <= 0:
        raise ValueError("top_n must be positive")
    manifest_rows = _manifest_rows(manifest, args.limit_signatures)
    signed_rows: list[dict[str, object]] = []
    gmt_rows: list[tuple[str, list[str]]] = []
    source_rows: list[dict[str, object]] = []
    for row in manifest_rows:
        persistent_id = str(row["persistent_id"])
        term = persistent_id_to_term(persistent_id)
        source = _materialize_source(row, cache_dir, args.source_url_base, int(args.request_timeout))
        up, down = _read_signature(source, top_n)
        gmt_rows.extend([(f"{term} up", [symbol for symbol, _ in up]), (f"{term} down", [symbol for symbol, _ in down])])
        for sign, records in ((1, up), (-1, down)):
            for symbol, coefficient in records:
                signed_rows.append({"term": term, "gene_id": symbol, "gene_symbol": symbol, "score": abs(coefficient), "sign": sign})
        source_rows.append({"persistent_id": persistent_id, "term": term, "source_path": str(source), "source_url": _source_url(row, args.source_url_base)})
    signed_path = out_dir / "lincs_l1000_cp_signed_term_gene.tsv"
    write_tsv(signed_path, signed_rows, ["term", "gene_id", "gene_symbol", "score", "sign"])
    source_path = out_dir / "lincs_l1000_cp_sources.tsv"
    write_tsv(source_path, source_rows, ["persistent_id", "term", "source_path", "source_url"])
    workflow_gmt = out_dir / "l1000_cp.gmt"
    with workflow_gmt.open("w", encoding="utf-8", newline="\n") as handle:
        for term, genes in gmt_rows:
            handle.write("\t".join([term, "", *genes]) + "\n")
    write_workflow_provenance_graph(
        workflow_name="lincs_l1000_cp",
        module_name="geneset_extractors.workflows.lincs_l1000_cp",
        output_dir=out_dir,
        focus_output_path=signed_path,
        output_paths=[(signed_path, "signed_term_gene_tsv"), (source_path, "source_manifest_tsv"), (workflow_gmt, "per_signature_gmt")],
        input_paths=[(manifest, "sigcom_signature_manifest")],
        parameters={"sigcom_library_id": LIBRARY_ID, "top_n": top_n, "ranking": "CD-coefficient descending; symbol ascending", "n_signatures": len(manifest_rows), "n_sets": len(gmt_rows)},
    )
    return {"n_rows": len(signed_rows), "n_sets": len(gmt_rows), "out_dir": str(out_dir)}
