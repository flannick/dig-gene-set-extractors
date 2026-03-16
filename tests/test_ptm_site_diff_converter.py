import csv
import gzip
import json
import shutil
from pathlib import Path

from geneset_extractors.converters import ptm_site_diff
from geneset_extractors.core.validate import validate_output_dir
from geneset_extractors.workflows.ptm_prepare_reference_bundle import run as run_ptm_prepare_reference_bundle


class Args:
    ptm_tsv = "tests/data/toy_ptm_site_diff.tsv"
    out_dir = "tests/tmp/ptm_site_diff"
    organism = "human"
    genome_build = "human"
    signature_name = "toy_ptm"
    dataset_label = "toy_ptm_dataset"
    ptm_type = "phospho"
    site_id_column = None
    site_group_column = None
    gene_id_column = None
    gene_symbol_column = None
    protein_accession_column = None
    residue_column = None
    position_column = None
    score_column = None
    stat_column = None
    logfc_column = None
    padj_column = None
    pvalue_column = None
    localization_prob_column = None
    peptide_count_column = None
    protein_logfc_column = None
    protein_stat_column = None
    score_mode = "auto"
    score_transform = "signed"
    protein_adjustment = "subtract"
    protein_adjustment_lambda = 1.0
    confidence_weight_mode = "combined"
    min_localization_prob = 0.75
    site_dup_policy = "highest_confidence"
    gene_aggregation = "signed_topk_mean"
    gene_topk_sites = 3
    ambiguous_gene_policy = "drop"
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    use_reference_bundle = True
    site_alias_resource_id = None
    site_ubiquity_resource_id = None
    select = "top_k"
    top_k = 200
    quantile = 0.01
    min_score = 0.0
    normalize = "within_set_l1"
    emit_full = True
    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = True
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 1
    gmt_max_genes = 10
    gmt_topk_list = "2"
    gmt_mass_list = ""
    gmt_split_signed = True
    gmt_format = "dig2col"
    emit_small_gene_sets = True
    neglog10p_cap = 50.0
    neglog10p_eps = 1e-300


class BundleArgs:
    sources_tsv = "tests/data/toy_ptm_sources.tsv"
    out_dir = "tests/tmp/ptm_bundle"
    organism = "human"
    ptm_type = "phospho"
    bundle_id = "toy_phospho_bundle_v1"



def _copy_resources(resources_dir: Path) -> None:
    resources_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile("tests/data/phosphosite_aliases_human_v1.tsv.gz", resources_dir / "phosphosite_aliases_human_v1.tsv.gz")
    shutil.copyfile("tests/data/phosphosite_ubiquity_human_v1.tsv.gz", resources_dir / "phosphosite_ubiquity_human_v1.tsv.gz")



def test_ptm_site_diff_converter_end_to_end(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args = Args()
    args.out_dir = str(tmp_path / "ptm_run")
    args.resources_dir = str(resources_dir)
    result = ptm_site_diff.run(args)
    assert result["resolved_score_mode"] == "logfc_times_neglog10p"

    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)

    with (Path(args.out_dir) / "geneset.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    assert {row["gene_symbol"] for row in rows[:2]} >= {"KCNN4", "MAPK1"}
    assert abs(sum(float(r["weight"]) for r in rows) - 1.0) < 1e-9

    gmt_lines = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8").strip().splitlines()
    assert any("__pos__" in line for line in gmt_lines)
    assert any("__neg__" in line for line in gmt_lines)

    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    provenance = json.loads((Path(args.out_dir) / "geneset.provenance.json").read_text(encoding="utf-8"))
    resources_info = meta["converter"]["parameters"]["resources"]
    assert len(resources_info["used"]) == 2
    assert meta["summary"]["n_sites_matched_to_ubiquity_prior"] >= 1
    alias_nodes = [node for node in provenance["nodes"] if node.get("role") == "site_alias_table"]
    assert alias_nodes
    assert alias_nodes[0]["identifiers"].get("stable_id")
    assert alias_nodes[0]["metadata"].get("provider")
    assert alias_nodes[0]["access"]["local_path"]



def test_ptm_site_diff_runs_without_bundle(tmp_path: Path, capsys):
    args = Args()
    args.out_dir = str(tmp_path / "ptm_no_bundle")
    args.resources_dir = None
    result = ptm_site_diff.run(args)
    assert result["n_genes"] >= 1
    captured = capsys.readouterr()
    assert "Missing" not in captured.err



def test_ptm_site_diff_protein_adjustment_changes_scores(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args_sub = Args()
    args_sub.out_dir = str(tmp_path / "ptm_subtract")
    args_sub.resources_dir = str(resources_dir)
    args_sub.protein_adjustment = "subtract"
    ptm_site_diff.run(args_sub)

    args_none = Args()
    args_none.out_dir = str(tmp_path / "ptm_none")
    args_none.resources_dir = str(resources_dir)
    args_none.protein_adjustment = "none"
    ptm_site_diff.run(args_none)

    def _score_map(path: Path) -> dict[str, float]:
        with path.open("r", encoding="utf-8") as fh:
            return {row["gene_id"]: float(row["score"]) for row in csv.DictReader(fh, delimiter="\t")}

    sub_scores = _score_map(Path(args_sub.out_dir) / "geneset.full.tsv")
    none_scores = _score_map(Path(args_none.out_dir) / "geneset.full.tsv")
    assert sub_scores["G_KCNN4"] != none_scores["G_KCNN4"]



def test_ptm_prepare_reference_bundle_workflow(tmp_path: Path):
    args = BundleArgs()
    args.out_dir = str(tmp_path / "ptm_bundle")
    result = run_ptm_prepare_reference_bundle(args)
    assert result["bundle_id"] == "toy_phospho_bundle_v1"

    out_dir = Path(args.out_dir)
    alias_path = out_dir / "phosphosite_aliases_human_v1.tsv.gz"
    ubiquity_path = out_dir / "phosphosite_ubiquity_human_v1.tsv.gz"
    assert alias_path.exists()
    assert ubiquity_path.exists()
    assert (out_dir / "bundle_provenance.json").exists()
    assert (out_dir / "local_resources_manifest.json").exists()

    with gzip.open(ubiquity_path, "rt", encoding="utf-8") as fh:
        rows = {row["canonical_site_key"]: row for row in csv.DictReader(fh, delimiter="\t")}
    assert rows["P12345|S|123|phospho"]["df_ref"] == "2"
    assert rows["P12345|S|123|phospho"]["n_samples_ref"] == "2"


def test_ptm_site_diff_emits_composition_qc_warning(tmp_path: Path, capsys):
    ptm_path = tmp_path / "composition_ptm.tsv"
    ptm_path.write_text(
        "\n".join(
            [
                "site_id\tgene_id\tgene_symbol\tprotein_accession\tresidue\tposition\tlog2fc\tpvalue",
                "HBB|S|1|phospho\tG_HBB\tHBB\tHBB\tS\t1\t2.0\t0.001",
                "HBA1|S|2|phospho\tG_HBA1\tHBA1\tHBA1\tS\t2\t1.8\t0.001",
                "VWF|S|3|phospho\tG_VWF\tVWF\tVWF\tS\t3\t1.5\t0.001",
                "PECAM1|S|4|phospho\tG_PECAM1\tPECAM1\tPECAM1\tS\t4\t1.4\t0.001",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    args = Args()
    args.ptm_tsv = str(ptm_path)
    args.out_dir = str(tmp_path / "ptm_composition")
    args.resources_dir = None
    args.use_reference_bundle = False
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "4"
    result = ptm_site_diff.run(args)
    assert result["n_genes"] == 4

    captured = capsys.readouterr()
    assert "composition QC suggests possible sample-composition dominance" in captured.err
    summary = json.loads((Path(args.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["tumor_intrinsic_confidence"] == "low"
    assert summary["composition_warning"]
