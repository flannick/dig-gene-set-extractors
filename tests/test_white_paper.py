from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.white_paper import (
    WHITE_PAPER_MARKDOWN_FILENAME,
    WHITE_PAPER_PDF_FILENAME,
    write_white_paper_from_metadata,
)


def _metadata_with_gmt(tmp_path: Path) -> dict[str, object]:
    source = tmp_path / "source.tsv"
    source.write_text("gene_id\tscore\nGENE1\t1\n", encoding="utf-8")
    (tmp_path / "geneset.tsv").write_text("gene_id\tscore\nGENE1\t1\n", encoding="utf-8")
    (tmp_path / "genesets.gmt").write_text(
        "set_a\tSet A\tGENE1\tGENE2\nset_b\tSet B\tGENE2\tGENE3\tGENE4\n",
        encoding="utf-8",
    )
    return make_metadata(
        converter_name="toy_converter",
        parameters={"threshold": 0.05},
        data_type="expression",
        assay="bulk_rna",
        organism="human",
        genome_build="hg38",
        files=[input_file_record(source, "source_tsv")],
        gene_annotation={"mode": "none", "source": "toy", "gene_id_field": "symbol"},
        weights={"weight_type": "score", "normalization": {}, "aggregation": "none"},
        summary={
            "n_input_features": 4,
            "n_genes": 4,
            "n_features_assigned": 4,
            "fraction_features_assigned": 1.0,
            "n_sets_emitted": 2,
        },
        output_files=[{"path": "genesets.gmt", "role": "gmt_library"}],
        gene_set_description="A compact, synthetic white-paper test library.",
    )


def test_metadata_write_emits_matching_white_paper_sidecars(tmp_path: Path):
    write_metadata(tmp_path / "geneset.meta.json", _metadata_with_gmt(tmp_path))

    markdown_path = tmp_path / WHITE_PAPER_MARKDOWN_FILENAME
    pdf_path = tmp_path / WHITE_PAPER_PDF_FILENAME
    markdown = markdown_path.read_text(encoding="utf-8")
    pdf = pdf_path.read_bytes()
    assert "# toy_converter:" in markdown
    assert "## Computational method" in markdown
    assert "Named gene sets | 2" in markdown
    assert "Distinct genes | 4" in markdown
    assert "Genes per set (min / median / max) | 2 / 2.5 / 3" in markdown
    assert pdf.startswith(b"%PDF-1.4")
    assert pdf.rstrip().endswith(b"%%EOF")
    assert hashlib.sha256(markdown.encode("utf-8")).hexdigest().encode("ascii") in pdf


def test_white_paper_cli_rebuilds_an_existing_pair(tmp_path: Path):
    write_metadata(tmp_path / "geneset.meta.json", _metadata_with_gmt(tmp_path))
    (tmp_path / WHITE_PAPER_MARKDOWN_FILENAME).unlink()
    (tmp_path / WHITE_PAPER_PDF_FILENAME).unlink()
    result = subprocess.run(
        [sys.executable, "-m", "geneset_extractors.cli", "provenance", "white-paper", str(tmp_path / "geneset.meta.json")],
        capture_output=True,
        text=True,
        env={**__import__("os").environ, "PYTHONPATH": "src"},
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["status"] == "ok"
    assert len(payload["white_papers"]) == 1
    assert Path(payload["white_papers"][0]["markdown_path"]).is_file()
    assert Path(payload["white_papers"][0]["pdf_path"]).is_file()


def test_multiple_declared_gmts_receive_distinct_stemmed_sidecars(tmp_path: Path):
    metadata = _metadata_with_gmt(tmp_path)
    (tmp_path / "secondary.gmt").write_text("set_c\tSet C\tGENE5\n", encoding="utf-8")
    metadata["output"]["files"].append({"path": "secondary.gmt", "role": "secondary_gmt"})  # type: ignore[index]
    write_metadata(tmp_path / "geneset.meta.json", metadata)

    paths = write_white_paper_from_metadata(tmp_path / "geneset.meta.json")
    assert len(paths) == 2
    assert (tmp_path / "genesets.whitepaper.md").is_file()
    assert (tmp_path / "secondary.whitepaper.pdf").is_file()
