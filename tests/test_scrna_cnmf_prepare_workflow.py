import csv
import json
from pathlib import Path
from types import SimpleNamespace

from omics2geneset.cli import main
from omics2geneset.workflows.scrna_cnmf_prepare import run as run_scrna_cnmf_prepare


def _make_args(out_dir: Path) -> SimpleNamespace:
    return SimpleNamespace(
        matrix_tsv="tests/data/toy_scrna_matrix.tsv",
        matrix_cell_id_column="cell_id",
        matrix_delim="\t",
        meta_tsv="tests/data/toy_scrna_meta.tsv",
        meta_cell_id_column="cell_id",
        donor_column="donor_id",
        phenotype_column="phenotype",
        cell_type_column="cell_type",
        split_by_cell_type=None,
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
        cnmf_name=None,
        cnmf_k_list="5 8",
        cnmf_n_iter=20,
        cnmf_numgenes=50,
        cnmf_total_workers=1,
        cnmf_densify=True,
        write_postprocess_template=True,
        execute=False,
        organism="human",
        genome_build="hg38",
    )


def _read_matrix(path: Path) -> tuple[list[str], list[list[float]]]:
    with path.open("r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        rows: list[list[float]] = []
        for row in reader:
            values = [float(x) for x in row[1:]]
            rows.append(values)
    return header, rows


def test_scrna_cnmf_prepare_creates_split_subsets_and_filters_zero_totals(tmp_path: Path):
    out_dir = tmp_path / "prep"
    args = _make_args(out_dir)
    result = run_scrna_cnmf_prepare(args)

    assert int(result["n_subsets"]) == 2
    assert (out_dir / "prepare_summary.json").exists()
    assert (out_dir / "subsets_manifest.tsv").exists()

    with (out_dir / "subsets_manifest.tsv").open("r", encoding="utf-8") as fh:
        manifest_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(manifest_rows) == 2
    subset_ids = {row["subset_id"] for row in manifest_rows}
    assert subset_ids == {"cell_type=A", "cell_type=B"}

    for row in manifest_rows:
        subset_dir = out_dir / "subsets" / row["subset_id"]
        counts_path = subset_dir / "counts_prefiltered.tsv"
        meta_path = subset_dir / "meta.tsv"
        script_path = subset_dir / "run_cnmf.sh"
        assert counts_path.exists()
        assert meta_path.exists()
        assert script_path.exists()
        assert script_path.stat().st_mode & 0o111

        header, values = _read_matrix(counts_path)
        assert len(header) >= 2
        assert values
        # No zero-total cells
        assert all(sum(v) > 0.0 for v in values)
        # No zero-total genes
        n_genes = len(header) - 1
        gene_totals = [0.0 for _ in range(n_genes)]
        for row_vals in values:
            for i, val in enumerate(row_vals):
                gene_totals[i] += val
        assert all(total > 0.0 for total in gene_totals)

    summary = json.loads((out_dir / "prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["workflow"] == "scrna_cnmf_prepare"
    assert summary["split_by_cell_type"] is True


def test_scrna_cnmf_prepare_downsampling_is_deterministic(tmp_path: Path):
    out_dir1 = tmp_path / "prep1"
    out_dir2 = tmp_path / "prep2"

    args1 = _make_args(out_dir1)
    args1.max_cells_per_bucket = 2
    args1.max_cells_total = 3
    args1.seed = 11
    run_scrna_cnmf_prepare(args1)

    args2 = _make_args(out_dir2)
    args2.max_cells_per_bucket = 2
    args2.max_cells_total = 3
    args2.seed = 11
    run_scrna_cnmf_prepare(args2)

    manifest1 = (out_dir1 / "subsets_manifest.tsv").read_text(encoding="utf-8")
    manifest2 = (out_dir2 / "subsets_manifest.tsv").read_text(encoding="utf-8")
    assert manifest1 == manifest2

    for subset in ["cell_type=A", "cell_type=B"]:
        counts1 = (out_dir1 / "subsets" / subset / "counts_prefiltered.tsv").read_text(encoding="utf-8")
        counts2 = (out_dir2 / "subsets" / subset / "counts_prefiltered.tsv").read_text(encoding="utf-8")
        assert counts1 == counts2


def test_scrna_cnmf_prepare_cli_entrypoint(tmp_path: Path):
    out_dir = tmp_path / "prep_cli"
    code = main(
        [
            "workflows",
            "scrna_cnmf_prepare",
            "--matrix_tsv",
            "tests/data/toy_scrna_matrix.tsv",
            "--meta_tsv",
            "tests/data/toy_scrna_meta.tsv",
            "--meta_cell_id_column",
            "cell_id",
            "--cell_type_column",
            "cell_type",
            "--donor_column",
            "donor_id",
            "--min_cells_per_cell_type",
            "1",
            "--out_dir",
            str(out_dir),
        ]
    )
    assert code == 0
    assert (out_dir / "prepare_summary.json").exists()
    assert (out_dir / "subsets_manifest.tsv").exists()
