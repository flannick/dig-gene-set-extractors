import csv
import json
from pathlib import Path

from geneset_extractors.converters import calr_ontology_mapper
from geneset_extractors.core.validate import validate_output_dir


class Args:
    calr_data_csv = "tests/data/toy_calr_data.csv"
    session_csv = "tests/data/toy_calr_session.csv"
    exclusions_tsv = "tests/data/toy_calr_exclusions.tsv"
    out_dir = "tests/tmp/calr_ontology"
    organism = "mouse"
    genome_build = "mm39"
    dataset_label = "toy_calr"
    signature_name = "toy_calr_sig"
    analysis_start_hour = None
    analysis_end_hour = None
    photoperiod_lights_on_hour = None
    photoperiod_hours_light = 12.0
    exploratory_without_session = True
    mass_covariate = None
    min_group_size = 2
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    reference_bundle_id = None
    term_templates_tsv = None
    phenotype_gene_edges_tsv = None
    term_hierarchy_tsv = None
    select = "top_k"
    top_k = 10
    quantile = 0.01
    min_score = 0.0
    normalize = "within_set_l1"
    emit_full = True
    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = False
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 1
    gmt_max_genes = 20
    gmt_topk_list = "5"
    gmt_mass_list = ""
    gmt_split_signed = False
    gmt_format = "classic"
    emit_small_gene_sets = True


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_calr_ontology_mapper_end_to_end(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "calr_ontology")
    result = calr_ontology_mapper.run(args)
    assert result["n_groups"] > 0

    out_dir = Path(args.out_dir)
    manifest_rows = _read_rows(out_dir / "manifest.tsv")
    assert manifest_rows
    assert any(row["program"] == "thermogenesis" for row in manifest_rows)
    assert (out_dir / "genesets.gmt").exists()

    thermogenesis_path = out_dir / "program=thermogenesis__mode=core__contrast=KO" / "geneset.tsv"
    rows = _read_rows(thermogenesis_path)
    gene_symbols = {row["gene_symbol"] for row in rows}
    assert "Ucp1" in gene_symbols or "Ppargc1a" in gene_symbols

    meta = json.loads((out_dir / "program=thermogenesis__mode=core__contrast=KO" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["summary"]["session_mode"] == "explicit"
    assert meta["summary"]["routing_summary"]["n_terms_with_gene_support"] >= 1
    validate_output_dir(out_dir, Path("src/geneset_extractors/schemas/geneset_metadata.schema.json"))


def test_calr_ontology_mapper_exploratory_without_session_warns(tmp_path: Path, capsys):
    args = Args()
    args.out_dir = str(tmp_path / "calr_ontology_exploratory")
    args.session_csv = None
    result = calr_ontology_mapper.run(args)
    assert result["n_groups"] > 0
    captured = capsys.readouterr()
    assert "not fully CalR-equivalent" in captured.err
    root_summary = json.loads((Path(args.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert root_summary["session_mode"] == "exploratory"
