import csv
import json
from pathlib import Path
from types import SimpleNamespace

from geneset_extractors.cli import main
from geneset_extractors.workflows.scrna_liger_prepare import run as run_scrna_liger_prepare
from geneset_extractors.workflows.scrna_liger_prepare import write_runtime_provenance


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
    assert (out_dir / "prepare_summary.provenance_graph.json").exists()
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
        liger_text = run_liger.read_text(encoding="utf-8")
        text = run_convert.read_text(encoding="utf-8")
        assert "scrna_liger_runtime_provenance" in liger_text
        assert "liger_run.provenance_graph.json" in liger_text
        assert "--liger_gene_loadings_tsv" in text
        assert "--upstream_provenance_graph_json" in text
        assert "liger_run.provenance_graph.json" in text

    summary = json.loads((out_dir / "prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["workflow"] == "scrna_liger_prepare"
    assert summary["n_subsets"] == 2
    assert summary["prepare_provenance_graph_path"] == "prepare_summary.provenance_graph.json"


def test_scrna_liger_runtime_provenance_connects_liger_outputs(tmp_path: Path):
    out_dir = tmp_path / "prep_runtime"
    args = _make_args(out_dir)
    run_scrna_liger_prepare(args)

    subset_dir = out_dir / "subsets" / "cell_type=B"
    program_dir = subset_dir / "liger_out" / "B"
    program_dir.mkdir(parents=True, exist_ok=True)
    (program_dir / "gene_loadings.tsv").write_text("gene_id\tFactor_1\nGENE1\t0.5\n", encoding="utf-8")
    (program_dir / "gene_programs.txt").write_text("Factor_1\nGENE1\n", encoding="utf-8")
    (program_dir / "cell_scores.tsv").write_text("cell_id\tFactor_1\ncell1\t0.8\n", encoding="utf-8")
    (program_dir / "metadata.txt").write_text("cell_type\tk\tmethod\ttimestamp\nB\t10\tLIGER_iNMF\tnow\n", encoding="utf-8")
    (program_dir / "k_stability.tsv").write_text("k\tstability\n10\t0.9\n", encoding="utf-8")
    (program_dir / "factor_importance.txt").write_text("Factor_1\n1.2\n", encoding="utf-8")

    result = write_runtime_provenance(
        SimpleNamespace(
            subset_dir=str(subset_dir),
            input_mode="matrix_tsv",
            input_path="counts_prefiltered.tsv",
            liger_output_dir="liger_out",
            run_liger_script="run_liger.sh",
            r_script=str(Path("src/geneset_extractors/preprocessing/rnaseq/liger_inmf.R").resolve()),
            runtime_graph_out="liger_run.provenance_graph.json",
            prepare_provenance_graph_json=str(out_dir / "prepare_summary.provenance_graph.json"),
            meta_path="meta.tsv",
            dataset_column="donor_id",
            cell_type_column="cell_type",
            cell_type_label="B",
            max_cells_total=20000,
            min_cells_per_cell_type=1,
            seed=1,
            liger_top_n_genes=250,
            liger_k_grid="10,12,14",
            liger_n_reps=3,
            liger_fixed_k="",
            liger_min_cells_per_dataset=1,
            liger_min_features=0,
            liger_min_umi=0.0,
            liger_max_mito=100.0,
        )
    )
    graph_path = Path(result["runtime_provenance_graph"])
    assert graph_path.exists()
    payload = json.loads(graph_path.read_text(encoding="utf-8"))
    graph = next(iter(payload.values()))
    paths = {
        str(node.get("c2m2_properties", {}).get("local_id", ""))
        for node in graph["nodes"]
        if str(node.get("type", "")) == "File"
    }
    assert str((program_dir / "gene_loadings.tsv").resolve()) in paths
    assert str((subset_dir / "run_liger.sh").resolve()) in paths
    assert str((out_dir / "prepare_summary.json").resolve()) in paths


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
    assert (out_dir / "prepare_summary.provenance_graph.json").exists()
    assert (out_dir / "subsets_manifest.tsv").exists()
