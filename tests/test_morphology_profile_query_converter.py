import csv
import json
from pathlib import Path

from geneset_extractors.converters import morphology_profile_query
from geneset_extractors.core.validate import validate_output_dir


class Args:
    query_profiles_tsv = "tests/data/toy_morph_query_profiles.tsv"
    out_dir = "tests/tmp/morphology"
    organism = "human"
    genome_build = "hg38"
    dataset_label = "toy_morph"
    signature_name = "toy_query"

    query_id_column = "sample_id"
    query_profiles_delimiter = "\t"
    query_metadata_tsv = "tests/data/toy_morph_query_metadata.tsv"
    query_metadata_id_column = "sample_id"
    query_metadata_delimiter = "\t"
    group_query_by = "query_group"
    query_aggregate = "median"
    exclude_query_ids_from_reference = False

    reference_profiles_tsv = "tests/data/toy_morph_reference_profiles.tsv"
    reference_profiles_parquet = None
    reference_profiles_delimiter = "\t"
    reference_id_column = "perturbation_id"
    reference_metadata_tsv = "tests/data/toy_morph_reference_metadata.tsv"
    reference_metadata_id_column = "perturbation_id"
    reference_metadata_delimiter = "\t"
    compound_targets_tsv = "tests/data/toy_morph_compound_targets.tsv"
    compound_targets_delimiter = "\t"
    compound_id_column = "compound_id"
    compound_target_gene_symbol_column = "gene_symbol"
    compound_target_weight_column = "weight"
    feature_stats_tsv = "tests/data/toy_morph_feature_stats.tsv"
    feature_schema_tsv = "tests/data/toy_morph_feature_schema.tsv"

    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    reference_bundle_id = None

    similarity_metric = "cosine"
    similarity_power = 1.0
    polarity = "both"
    compound_weight = 0.5
    genetic_weight = 0.5

    select = "top_k"
    top_k = 2
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
    gmt_max_genes = 10
    gmt_topk_list = "2"
    gmt_mass_list = ""
    gmt_split_signed = False
    gmt_format = "classic"
    emit_small_gene_sets = True


def _manifest_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _geneset_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_morphology_profile_query_explicit_files(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "morph_explicit")
    result = morphology_profile_query.run(args)
    assert result["n_groups"] == 4

    out_dir = Path(args.out_dir)
    manifest = _manifest_rows(out_dir / "manifest.tsv")
    assert len(manifest) == 4
    assert (out_dir / "genesets.gmt").exists()

    q1_similar = out_dir / "program=Q1__polarity=similar" / "geneset.tsv"
    rows = _geneset_rows(q1_similar)
    assert [row["gene_id"] for row in rows] == ["GENE1", "GENE3"]
    weight_sum = sum(float(row["weight"]) for row in rows)
    assert abs(weight_sum - 1.0) < 1e-9

    meta = json.loads((out_dir / "program=Q1__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["input"]["assay"] == "morphology"
    assert meta["summary"]["mapping_summary"]["n_compound_refs_with_targets"] == 2
    assert meta["gmt"]["written"] is True

    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(out_dir, schema)


def test_morphology_profile_query_bundle_mode(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "morph_bundle")
    args.reference_profiles_tsv = None
    args.reference_metadata_tsv = None
    args.compound_targets_tsv = None
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.resources_dir = "tests/data/morphology_bundle"
    args.reference_bundle_id = "morphology_jump_target_pilot_u2os_48h_v1"
    result = morphology_profile_query.run(args)
    assert result["n_groups"] == 4

    out_dir = Path(args.out_dir)
    meta = json.loads((out_dir / "program=Q1__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    resources = meta["summary"]["resources"]
    assert resources is not None
    used_ids = {row["id"] for row in resources["used"]}
    assert "morphology_jump_target_pilot_u2os_48h_v1" in used_ids


def test_morphology_profile_query_missing_bundle_warns_or_fails(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "morph_missing_bundle")
    args.reference_profiles_tsv = None
    args.reference_metadata_tsv = None
    args.compound_targets_tsv = None
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.resources_dir = str(tmp_path / "empty_resources")
    args.reference_bundle_id = "morphology_jump_target_pilot_u2os_48h_v1"
    Path(args.resources_dir).mkdir(parents=True, exist_ok=True)
    try:
        morphology_profile_query.run(args)
    except ValueError as exc:
        assert "Could not resolve morphology reference bundle" in str(exc)
    else:
        raise AssertionError("Expected missing bundle resolution to fail")
