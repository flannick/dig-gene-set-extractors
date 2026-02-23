import gzip
from pathlib import Path

import pytest

from omics2geneset.io.reference_tables import read_atlas_gene_stats_tsv, read_ref_ubiquity_tsv


def _write_gz(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write(text)


def test_read_ref_ubiquity_accepts_bundle_ccre_schema(tmp_path: Path):
    bed_path = tmp_path / "encode_ccre_hg38.bed.gz"
    _write_gz(
        bed_path,
        "\n".join(
            [
                "chr1\t100\t200\tEH38D1\tEH38E1\tpELS",
                "chr2\t300\t450\tEH38D2\tEH38E2\tCA",
                "",
            ]
        ),
    )
    ubiq_path = tmp_path / "ccre_ubiquity_hg38.tsv.gz"
    _write_gz(
        ubiq_path,
        "\n".join(
            [
                "ccre_id\tdf\tN\tidf\tidf_norm",
                "EH38E1\t10\t100\t2.5\t0.5",
                "EH38E2\t20\t100\t1.5\t0.3",
                "EH38EX\t8\t100\t2.0\t0.4",
                "",
            ]
        ),
    )

    rows = read_ref_ubiquity_tsv(ubiq_path)
    assert len(rows) == 2
    assert rows[0] == {"chrom": "chr1", "start": 100, "end": 200, "idf_ref": 2.5}
    assert rows[1] == {"chrom": "chr2", "start": 300, "end": 450, "idf_ref": 1.5}


def test_read_ref_ubiquity_rejects_placeholder_na_table(tmp_path: Path):
    bed_path = tmp_path / "encode_ccre_mm10.bed.gz"
    _write_gz(
        bed_path,
        "\n".join(
            [
                "chr1\t10\t20\tEM10D1\tEM10E1\tCA",
                "",
            ]
        ),
    )
    ubiq_path = tmp_path / "ccre_ubiquity_mm10.tsv.gz"
    _write_gz(
        ubiq_path,
        "\n".join(
            [
                "ccre_id\tdf\tN\tidf\tidf_norm",
                "EM10E1\tNA\tNA\tNA\tNA",
                "",
            ]
        ),
    )
    with pytest.raises(ValueError, match="no usable rows"):
        read_ref_ubiquity_tsv(ubiq_path)


def test_read_atlas_gene_stats_accepts_long_profile_schema(tmp_path: Path):
    atlas_path = tmp_path / "atac_reference_profiles_hg38.tsv.gz"
    _write_gz(
        atlas_path,
        "\n".join(
            [
                "context_id\tcontext_label\tgene_id\tscore",
                "C1\tctx1\tGENE1\t1.0",
                "C2\tctx2\tGENE1\t3.0",
                "C3\tctx3\tGENE1\t5.0",
                "C1\tctx1\tGENE2\t2.0",
                "C2\tctx2\tGENE2\t2.0",
                "",
            ]
        ),
    )

    stats = read_atlas_gene_stats_tsv(atlas_path)
    assert stats["GENE1"] == (3.0, 2.0)
    assert stats["GENE2"] == (2.0, 0.0)
