from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import pytest

from geneset_extractors.workflows.lincs_l1000_cp import persistent_id_to_term, run


def _source(path: Path, *, duplicate: bool = False) -> None:
    lines = ["symbol\tCD-coefficient"]
    for index in range(500):
        symbol = "GENE000" if duplicate and index == 1 else f"GENE{index:03d}"
        lines.append(f"{symbol}\t{500 - index}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _args(manifest: Path, out_dir: Path) -> Namespace:
    return Namespace(signature_manifest_tsv=str(manifest), out_dir=str(out_dir), cache_dir=None, top_n=250, source_url_base="https://example.invalid/", request_timeout=1, limit_signatures=None)


def test_persistent_id_term_handles_filename_and_url() -> None:
    value = "L1000_LINCS_DCIC_ABY001_A375_XH_A13_afatinib_10uM.tsv"
    assert persistent_id_to_term(value) == "ABY001_A375_XH_A13_afatinib_10uM"
    assert persistent_id_to_term("https://example.org/path/" + value) == "ABY001_A375_XH_A13_afatinib_10uM"


def test_exports_two_exact_250_gene_sets(tmp_path: Path) -> None:
    source = tmp_path / "signature.tsv"
    _source(source)
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(f"persistent_id\tsource_path\nL1000_LINCS_DCIC_test.tsv\t{source}\n", encoding="utf-8")
    out_dir = tmp_path / "out"
    result = run(_args(manifest, out_dir))
    assert result["n_sets"] == 2
    lines = (out_dir / "l1000_cp.gmt").read_text(encoding="utf-8").splitlines()
    assert [line.split("\t", 1)[0] for line in lines] == ["test up", "test down"]
    assert all(len(line.split("\t")) == 252 for line in lines)
    assert lines[0].split("\t")[2:5] == ["GENE000", "GENE001", "GENE002"]
    assert lines[1].split("\t")[-3:] == ["GENE497", "GENE498", "GENE499"]


def test_rejects_duplicate_symbols(tmp_path: Path) -> None:
    source = tmp_path / "signature.tsv"
    _source(source, duplicate=True)
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(f"persistent_id\tsource_path\nL1000_LINCS_DCIC_test.tsv\t{source}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate symbol"):
        run(_args(manifest, tmp_path / "out"))
