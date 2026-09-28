"""Deterministic scientific-style white-paper sidecars for GMT outputs.

The Markdown document is the canonical source.  The companion PDF renders the
same document with a deliberately small built-in PDF writer so normal DIG
generation does not acquire a document-rendering runtime dependency.
"""
from __future__ import annotations

import hashlib
import json
import statistics
import textwrap
import unicodedata
from pathlib import Path
from typing import Any


WHITE_PAPER_VERSION = "1.0.0"
WHITE_PAPER_MARKDOWN_FILENAME = "geneset.whitepaper.md"
WHITE_PAPER_PDF_FILENAME = "geneset.whitepaper.pdf"


def _resolve_declared_gmts(metadata: dict[str, Any], output_dir: Path) -> list[Path]:
    candidates: list[Path] = []
    for record in ((metadata.get("output") or {}).get("files") or []):
        if not isinstance(record, dict):
            continue
        value = record.get("path")
        if not isinstance(value, str) or not value.lower().endswith(".gmt"):
            continue
        candidate = Path(value)
        candidate = candidate if candidate.is_absolute() else output_dir / candidate
        if candidate.is_file():
            candidates.append(candidate.resolve())
    unique = sorted(set(candidates))
    return unique


def _gmt_statistics(gmt_path: Path) -> dict[str, int | float]:
    gene_sets = 0
    genes: set[str] = set()
    members_per_set: list[int] = []
    for line_number, raw in enumerate(gmt_path.read_text(encoding="utf-8").splitlines(), start=1):
        fields = raw.split("\t")
        if len(fields) < 3 or not fields[0] or any(not gene for gene in fields[2:]):
            raise ValueError(f"{gmt_path}: invalid GMT row {line_number}")
        row_genes = set(fields[2:])
        gene_sets += 1
        genes.update(row_genes)
        members_per_set.append(len(row_genes))
    if not members_per_set:
        raise ValueError(f"{gmt_path}: no GMT rows found")
    return {
        "n_gene_sets": gene_sets,
        "n_distinct_genes": len(genes),
        "n_gene_memberships": sum(members_per_set),
        "min_genes_per_set": min(members_per_set),
        "median_genes_per_set": statistics.median(members_per_set),
        "max_genes_per_set": max(members_per_set),
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _table(rows: list[tuple[str, object]]) -> str:
    lines = ["| Field | Value |", "| --- | --- |"]
    lines.extend(f"| {field} | {str(value).replace('|', '\\|')} |" for field, value in rows)
    return "\n".join(lines)


def _compact_json(value: object) -> str:
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False)


def render_white_paper(
    *,
    metadata: dict[str, Any],
    gmt_path: Path,
    legacy_provenance_path: Path,
    dapper_provenance_path: Path,
) -> str:
    """Render one canonical scientific-style Markdown report from output facts."""
    statistics_payload = _gmt_statistics(gmt_path)
    gene_set = metadata.get("gene_set") or {}
    converter = metadata.get("converter") or {}
    execution = converter.get("execution") or {}
    code = converter.get("code") or {}
    summary = metadata.get("summary") or {}
    inputs = (metadata.get("input") or {}).get("files") or []
    input_rows = []
    for index, record in enumerate(inputs, start=1):
        if not isinstance(record, dict):
            continue
        label = record.get("role") or f"input_{index}"
        identifier = record.get("canonical_uri") or record.get("local_path") or record.get("path") or "not declared"
        input_rows.append((str(label), identifier))
    provenance_rows = [
        ("Legacy provenance", f"`{legacy_provenance_path.name}` ({_sha256(legacy_provenance_path)})"),
        ("DAPPER provenance", f"`{dapper_provenance_path.name}` ({_sha256(dapper_provenance_path)})"),
        ("DIG repository", code.get("repo_url") or "not declared"),
        ("DIG commit", code.get("git_commit") or "not declared"),
    ]
    command = execution.get("observed_command") or execution.get("command") or "not declared"
    if isinstance(command, list):
        command = " ".join(str(part) for part in command)
    return "\n".join(
        [
            "---",
            f"white_paper_version: {WHITE_PAPER_VERSION}",
            f"gmt_sha256: {_sha256(gmt_path)}",
            "---",
            "",
            f"# {gene_set.get('name') or metadata.get('geneset_id') or gmt_path.stem}",
            "",
            "## Abstract",
            "",
            str(gene_set.get("description") or "This report describes a DIG-generated gene-set artifact and its recorded provenance."),
            "",
            "## Scope and artifact",
            "",
            _table(
                [
                    ("Gene-set identifier", metadata.get("geneset_id") or "not declared"),
                    ("Assay", gene_set.get("assay") or "not declared"),
                    ("Data type", gene_set.get("data_type") or "not declared"),
                    ("Organism", gene_set.get("organism") or "not declared"),
                    ("Genome build", gene_set.get("genome_build") or "not declared"),
                    ("GMT artifact", gmt_path.name),
                    ("GMT SHA-256", _sha256(gmt_path)),
                ]
            ),
            "",
            "## Source materials",
            "",
            _table(input_rows or [("Inputs", "No input records were declared.")]),
            "",
            "## Computational method",
            "",
            _table(
                [
                    ("Converter or workflow", converter.get("name") or "not declared"),
                    ("Software version", converter.get("version") or "not declared"),
                    ("Entrypoint", execution.get("entrypoint") or "not declared"),
                ]
            ),
            "",
            "### Recorded parameters",
            "",
            "```json",
            _compact_json(converter.get("parameters") or {}),
            "```",
            "",
            "### Recorded command",
            "",
            "```text",
            str(command),
            "```",
            "",
            "## Gene-set statistics",
            "",
            _table(
                [
                    ("Named gene sets", statistics_payload["n_gene_sets"]),
                    ("Distinct genes", statistics_payload["n_distinct_genes"]),
                    ("Gene memberships", statistics_payload["n_gene_memberships"]),
                    ("Genes per set (min / median / max)", f"{statistics_payload['min_genes_per_set']} / {statistics_payload['median_genes_per_set']} / {statistics_payload['max_genes_per_set']}"),
                    ("Metadata-reported genes", summary.get("n_genes") or "not declared"),
                ]
            ),
            "",
            "## Provenance and reproducibility",
            "",
            _table(provenance_rows),
            "",
            "The paired legacy JSON and DAPPER YAML sidecars are the authoritative machine-readable provenance records. This white paper is a human-readable summary derived from those records and the GMT; it does not replace them.",
            "",
        ]
    )


def _pdf_escape(value: str) -> str:
    # Base-14 PDF fonts are WinAnsi. Keep all source content in the Markdown;
    # normalize only the display fallback when a glyph cannot be represented.
    rendered = unicodedata.normalize("NFKD", value).encode("latin-1", "replace").decode("latin-1")
    return rendered.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def _pdf_stream(markdown: str) -> bytes:
    lines: list[str] = []
    for source_line in markdown.splitlines():
        lines.extend(textwrap.wrap(source_line, width=92, replace_whitespace=False) or [""])
    page_lines = 48
    pages = [lines[index : index + page_lines] for index in range(0, len(lines), page_lines)] or [[]]
    objects: list[bytes] = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"",  # Filled after all page object numbers are known.
    ]
    page_ids: list[int] = []
    for page in pages:
        content = "BT\n/F1 9 Tf\n50 760 Td\n12 TL\n" + "\n".join(
            f"({_pdf_escape(line)}) Tj\nT*" for line in page
        ) + "\nET\n"
        content_bytes = content.encode("latin-1")
        stream_id = len(objects) + 1
        page_id = stream_id + 1
        objects.append(f"<< /Length {len(content_bytes)} >>\nstream\n".encode("ascii") + content_bytes + b"endstream")
        objects.append(b"")  # Filled after the shared font object is known.
        page_ids.append(page_id)
    font_id = len(objects) + 1
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    for page_id, stream_id in zip(page_ids, range(3, 3 + len(page_ids) * 2, 2), strict=True):
        objects[page_id - 1] = (
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] "
            f"/Resources << /Font << /F1 {font_id} 0 R >> >> /Contents {stream_id} 0 R >>"
        ).encode("ascii")
    objects[1] = (
        f"<< /Type /Pages /Count {len(page_ids)} "
        f"/Kids [{' '.join(f'{page_id} 0 R' for page_id in page_ids)}] >>"
    ).encode("ascii")
    source_hash = hashlib.sha256(markdown.encode("utf-8")).hexdigest()
    objects.append(f"<< /Title (DIG gene-set white paper) /Subject (Canonical Markdown SHA-256: {source_hash}) >>".encode("ascii"))
    encoded = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offsets = [0]
    for index, body in enumerate(objects, start=1):
        offsets.append(len(encoded))
        encoded.extend(f"{index} 0 obj\n".encode("ascii"))
        encoded.extend(body)
        encoded.extend(b"\nendobj\n")
    startxref = len(encoded)
    encoded.extend(f"xref\n0 {len(objects) + 1}\n0000000000 65535 f \n".encode("ascii"))
    encoded.extend(b"".join(f"{offset:010d} 00000 n \n".encode("ascii") for offset in offsets[1:]))
    encoded.extend(f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R /Info {len(objects)} 0 R >>\nstartxref\n{startxref}\n%%EOF\n".encode("ascii"))
    return bytes(encoded)


def write_white_paper(
    *, metadata: dict[str, Any], output_dir: str | Path,
    legacy_provenance_path: str | Path, dapper_provenance_path: str | Path,
) -> list[tuple[Path, Path]]:
    """Write matching Markdown/PDF sidecars for every declared existing GMT."""
    directory = Path(output_dir)
    gmt_paths = _resolve_declared_gmts(metadata, directory)
    legacy_path = Path(legacy_provenance_path)
    dapper_path = Path(dapper_provenance_path)
    if gmt_paths and (not legacy_path.is_file() or not dapper_path.is_file()):
        raise FileNotFoundError("white-paper generation requires both provenance sidecars")
    paths: list[tuple[Path, Path]] = []
    for gmt_path in gmt_paths:
        markdown = render_white_paper(
            metadata=metadata,
            gmt_path=gmt_path,
            legacy_provenance_path=legacy_path,
            dapper_provenance_path=dapper_path,
        )
        if len(gmt_paths) == 1:
            markdown_path = directory / WHITE_PAPER_MARKDOWN_FILENAME
            pdf_path = directory / WHITE_PAPER_PDF_FILENAME
        else:
            markdown_path = directory / f"{gmt_path.stem}.whitepaper.md"
            pdf_path = directory / f"{gmt_path.stem}.whitepaper.pdf"
        markdown_path.write_text(markdown, encoding="utf-8", newline="\n")
        pdf_path.write_bytes(_pdf_stream(markdown))
        paths.append((markdown_path, pdf_path))
    return paths


def write_white_paper_from_metadata(metadata_path: str | Path) -> list[tuple[Path, Path]]:
    """Rebuild sidecars for an existing metadata/provenance/GMT combination."""
    path = Path(metadata_path)
    metadata = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(metadata, dict):
        raise ValueError("metadata payload must be a JSON object")
    provenance = metadata.get("provenance") or {}
    legacy_name = str(provenance.get("path") or "geneset.provenance.legacy.json")
    dapper_name = str(provenance.get("dapper_path") or "geneset.provenance.dapper.yaml")
    return write_white_paper(
        metadata=metadata,
        output_dir=path.parent,
        legacy_provenance_path=path.parent / legacy_name,
        dapper_provenance_path=path.parent / dapper_name,
    )
