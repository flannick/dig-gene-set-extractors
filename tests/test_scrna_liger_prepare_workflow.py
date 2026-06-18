import csv
import json
from pathlib import Path
from types import SimpleNamespace

from geneset_extractors.cli import main
from geneset_extractors.workflows.scrna_liger_prepare import run as run_scrna_liger_prepare


def _make_args(out_dir: Path) -> SimpleNamespace:
    return SimpleNamespace(
        matrix_tsv="tests/data/toy_scrna_matrix.tsv",
        h5ad=None,
        seurat_rds=None,
        mtx_dir=None,
        matrix_orientation="auto",
        matrix_cell_id_column="cell_id",
        matrix_gene_id_column=None,
        matrix_delim="\t",
        meta_tsv="tests/data/toy_scrna_meta.tsv",
        meta_cell_id_column="cell_id",
        dataset_column="donor_id",
        cell_type_column="cell_type",
        split_by_cell_type=True,
        cell_type_allowlist=None,
        min_cells_per_cell_type=1,
        bucket_columns=None,
        max_cells_per_bucket=200,
        max_cells_total=20000,
        seed=1,
        matrix_value_type="logcounts",
        min_total_per_cell=None,
        min_total_per_gene=None,
        out_dir=str(out_dir),
        keep_tmp=False,
        execute=False,
        organism="human",
        genome_build="hg38",
        liger_k_grid="10,12,14",
        liger_n_reps=3,
        liger_fixed_k=None,
        liger_top_n_genes=250,
        liger_min_cells_per_dataset=1,
        liger_min_features=0,
        liger_min_umi=0.0,
        liger_max_mito=100.0,
        extractor_top_k=250,
    )


def test_scrna_liger_prepare_creates_split_subsets_and_scripts(tmp_path: Path):
    out_dir = tmp_path / "prep"
    args = _make_args(out_dir)
    result = run_scrna_liger_prepare(args)

    assert int(result["n_subsets"]) == 2
    assert (out_dir / "prepare_summary.json").exists()
    assert (out_dir / "subsets_manifest.tsv").exists()

    with (out_dir / "subsets_manifest.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 2

    for row in rows:
        subset_dir = out_dir / "subsets" / row["subset_id"]
        assert (subset_dir / "counts_prefiltered.tsv").exists()
        assert (subset_dir / "meta.tsv").exists()
        run_liger = subset_dir / "run_liger.sh"
        run_convert = subset_dir / "run_geneset_extractors_from_liger.sh"
        assert run_liger.exists()
        assert run_convert.exists()
        assert run_liger.stat().st_mode & 0o111
        assert run_convert.stat().st_mode & 0o111
        text = run_convert.read_text(encoding="utf-8")
        assert "--liger_gene_loadings_tsv" in text

    summary = json.loads((out_dir / "prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["workflow"] == "scrna_liger_prepare"
    assert summary["n_subsets"] == 2


def test_scrna_liger_prepare_cli_entrypoint(tmp_path: Path):
    out_dir = tmp_path / "prep_cli"
    code = main(
        [
            "workflows",
            "scrna_liger_prepare",
            "--matrix_tsv",
            "tests/data/toy_scrna_matrix.tsv",
            "--meta_tsv",
            "tests/data/toy_scrna_meta.tsv",
            "--meta_cell_id_column",
            "cell_id",
            "--dataset_column",
            "donor_id",
            "--cell_type_column",
            "cell_type",
            "--min_cells_per_cell_type",
            "1",
            "--out_dir",
            str(out_dir),
            "--organism",
            "human",
            "--genome_build",
            "hg38",
        ]
    )
    assert code == 0
    assert (out_dir / "prepare_summary.json").exists()
    assert (out_dir / "subsets_manifest.tsv").exists()
