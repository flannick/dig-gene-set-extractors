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
    for field, value in rows:
        escaped_value = str(value).replace("|", "\\|")
        lines.append(f"| {field} | {escaped_value} |")
    return "\n".join(lines)


def _compact_json(value: object) -> str:
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False)


def _human_list(values: list[str]) -> str:
    if not values:
        return "none"
    if len(values) == 1:
        return values[0]
    if len(values) == 2:
        return f"{values[0]} and {values[1]}"
    return f"{', '.join(values[:-1])}, and {values[-1]}"


def _input_description(record: dict[str, Any], index: int) -> str:
    role = str(record.get("role") or f"input {index}")
    source = _display_input_source(record)
    details = []
    for key, label in (
        ("provider", "provider"),
        ("version", "release/version"),
        ("access_level", "access level"),
        ("license", "license"),
        ("sha256", "SHA-256"),
    ):
        value = record.get(key)
        if value not in (None, ""):
            details.append(f"{label}={value}")
    suffix = f" ({'; '.join(details)})" if details else ""
    return f"The **{role}** input was obtained from `{source}`{suffix}."


def _display_input_source(record: dict[str, Any]) -> str:
    for key in ("canonical_uri", "download_url", "persistent_id", "stable_id"):
        value = record.get(key)
        if value not in (None, ""):
            return str(value)
    local_value = record.get("local_path") or record.get("path")
    if local_value not in (None, ""):
        # White papers are typically published. Keep the data artifact's
        # filename while avoiding a contributor-specific workspace path.
        return f"local file {Path(str(local_value)).name} (workspace path omitted)"
    return "an undeclared location"


def _parameter_sentence(parameters: object) -> str:
    if not isinstance(parameters, dict) or not parameters:
        return "No converter parameters were recorded."
    items = [f"`{key}={json.dumps(value, sort_keys=True, ensure_ascii=False)}`" for key, value in sorted(parameters.items())]
    return "The recorded parameterization was " + _human_list(items) + "."


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
        identifier = _display_input_source(record)
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
    input_paragraphs = [
        _input_description(record, index)
        for index, record in enumerate(inputs, start=1)
        if isinstance(record, dict)
    ]
    data_paragraph = " ".join(input_paragraphs) if input_paragraphs else (
        "No input records were declared in the metadata; the provenance records should be consulted before interpreting this artifact."
    )
    workflow_paragraph = (
        f"Gene sets were generated by `{converter.get('name') or 'an undeclared converter'}` "
        f"(software version `{converter.get('version') or 'not declared'}`) through "
        f"`{execution.get('entrypoint') or 'an undeclared entrypoint'}`. "
        f"{_parameter_sentence(converter.get('parameters'))}"
    )
    output_paragraph = (
        f"The workflow produced `{gmt_path.name}`, containing {statistics_payload['n_gene_sets']} named gene sets, "
        f"{statistics_payload['n_distinct_genes']} distinct genes, and {statistics_payload['n_gene_memberships']} total gene memberships. "
        f"Set sizes range from {statistics_payload['min_genes_per_set']} to {statistics_payload['max_genes_per_set']} genes "
        f"with a median of {statistics_payload['median_genes_per_set']}."
    )
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
            "## Data used",
            "",
            data_paragraph,
            "",
            "The library describes "
            f"{gene_set.get('assay') or 'an undeclared assay'} {gene_set.get('data_type') or 'data'} "
            f"for {gene_set.get('organism') or 'an undeclared organism'} "
            f"on {gene_set.get('genome_build') or 'an undeclared genome build'}.",
            "",
            "### Input inventory",
            "",
            _table(input_rows or [("Inputs", "No input records were declared.")]),
            "",
            "## Workflow and parameterization",
            "",
            workflow_paragraph,
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
            "## Output generated",
            "",
            output_paragraph,
            "",
            _table(
                [
                    ("Gene-set identifier", metadata.get("geneset_id") or "not declared"),
                    ("GMT artifact", gmt_path.name),
                    ("GMT SHA-256", _sha256(gmt_path)),
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


def _pdf_plain_text(value: str) -> str:
    return value.replace("`", "").replace("**", "").replace("\\|", "|")


def _pdf_lines(markdown: str) -> list[tuple[str, str, float, float]]:
    """Translate the report's small Markdown subset into styled PDF lines."""
    source = markdown.splitlines()
    lines: list[tuple[str, str, float, float]] = []
    in_front_matter = False
    in_code = False
    table_header = False
    for index, raw_line in enumerate(source):
        if index == 0 and raw_line == "---":
            in_front_matter = True
            continue
        if in_front_matter:
            if raw_line == "---":
                in_front_matter = False
            continue
        if raw_line.startswith("```"):
            in_code = not in_code
            continue
        if raw_line.startswith("# "):
            lines.append((_pdf_plain_text(raw_line[2:]), "F2", 18.0, 27.0))
            continue
        if raw_line.startswith("## "):
            lines.append((_pdf_plain_text(raw_line[3:]), "F2", 13.0, 22.0))
            continue
        if raw_line.startswith("### "):
            lines.append((_pdf_plain_text(raw_line[4:]), "F2", 10.5, 17.0))
            continue
        if raw_line.startswith("| ") and raw_line.endswith(" |"):
            cells = [cell.strip() for cell in raw_line.strip("|").split("|")]
            if cells == ["Field", "Value"]:
                table_header = True
                continue
            if table_header and all(cell.replace("-", "") == "" for cell in cells):
                continue
            table_header = False
            text = f"{_pdf_plain_text(cells[0])}: {_pdf_plain_text(' | '.join(cells[1:]))}"
            lines.extend((part, "F1", 8.5, 11.0) for part in textwrap.wrap(text, width=102) or [""])
            continue
        if not raw_line:
            lines.append(("", "F1", 9.0, 7.0))
            continue
        font, size, leading = ("F3", 8.0, 10.5) if in_code else ("F1", 9.5, 14.0)
        lines.extend(
            (_pdf_plain_text(part), font, size, leading)
            for part in textwrap.wrap(raw_line, width=96, replace_whitespace=False) or [""]
        )
    return lines


def _paginate_pdf_lines(lines: list[tuple[str, str, float, float]]) -> list[list[tuple[str, str, float, float]]]:
    pages: list[list[tuple[str, str, float, float]]] = [[]]
    height = 756.0
    for line in lines:
        if height - line[3] < 48.0 and pages[-1]:
            pages.append([])
            height = 756.0
        pages[-1].append(line)
        height -= line[3]
    return pages


def _pdf_page_content(page: list[tuple[str, str, float, float]], page_number: int) -> bytes:
    y = 756.0
    commands: list[str] = []
    for text, font, size, leading in page:
        if text:
            commands.append(f"BT /{font} {size:g} Tf 50 {y:g} Td ({_pdf_escape(text)}) Tj ET")
        y -= leading
    commands.append("0.4 w 50 38 m 562 38 l S")
    commands.append(f"BT /F1 8 Tf 50 25 Td (DIG gene-set white paper | Page {page_number}) Tj ET")
    return ("\n".join(commands) + "\n").encode("latin-1")


def _pdf_stream(markdown: str) -> bytes:
    pages = _paginate_pdf_lines(_pdf_lines(markdown))
    objects: list[bytes] = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"",  # Filled after all page object numbers are known.
    ]
    page_ids: list[int] = []
    for page_number, page in enumerate(pages, start=1):
        content_bytes = _pdf_page_content(page, page_number)
        stream_id = len(objects) + 1
        page_id = stream_id + 1
        objects.append(f"<< /Length {len(content_bytes)} >>\nstream\n".encode("ascii") + content_bytes + b"endstream")
        objects.append(b"")  # Filled after the shared font object is known.
        page_ids.append(page_id)
    regular_font_id = len(objects) + 1
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    bold_font_id = len(objects) + 1
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>")
    mono_font_id = len(objects) + 1
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Courier >>")
    for page_id, stream_id in zip(page_ids, range(3, 3 + len(page_ids) * 2, 2), strict=True):
        objects[page_id - 1] = (
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] "
            f"/Resources << /Font << /F1 {regular_font_id} 0 R /F2 {bold_font_id} 0 R /F3 {mono_font_id} 0 R >> >> "
            f"/Contents {stream_id} 0 R >>"
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
