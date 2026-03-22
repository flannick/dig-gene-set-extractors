import csv
import json
from pathlib import Path

from geneset_extractors.converters import sc_rna_marker
from geneset_extractors.core.validate import validate_output_dir
from tests.provenance_helpers import assert_manifest_has_enriched_columns


class Args:
    counts_tsv = "tests/data/toy_sc_rna_counts.tsv"
    out_dir = "tests/tmp/sc_rna"
    organism = "human"
    genome_build = "hg38"
    groups_tsv = None
    gene_id_column = "gene_id"
    barcode_column = "barcode"
    value_column = "count"
    group_barcode_column = "barcode"
    group_column = "group"
    peak_summary = "sum_counts"
    normalize = "l1"
    provenance_overlay_json = None


def test_sc_rna_marker_without_groups(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_rna")
    args.groups_tsv = None
    sc_rna_marker.run(args)
    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)


def test_sc_rna_marker_with_groups(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "sc_rna_groups")
    args.groups_tsv = "tests/data/toy_sc_rna_groups.tsv"
    sc_rna_marker.run(args)
    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir) / "group=g1", schema)
    validate_output_dir(Path(args.out_dir) / "group=g2", schema)
    assert (Path(args.out_dir) / "manifest.tsv").exists()
    rows = list(csv.DictReader((Path(args.out_dir) / "manifest.tsv").open("r", encoding="utf-8"), delimiter="\t"))
    assert_manifest_has_enriched_columns(rows)


def test_sc_rna_marker_overlay_propagates(tmp_path: Path):
    overlay_path = tmp_path / "overlay.json"
    overlay_path.write_text(
        '{"inputs":{"role:counts_tsv":{"landing_page_url":"https://example.org/counts"}}}',
        encoding="utf-8",
    )
    args = Args()
    args.out_dir = str(tmp_path / "sc_rna_overlay")
    args.groups_tsv = "tests/data/toy_sc_rna_groups.tsv"
    args.provenance_overlay_json = str(overlay_path)
    sc_rna_marker.run(args)
    rows = list(csv.DictReader((Path(args.out_dir) / "manifest.tsv").open("r", encoding="utf-8"), delimiter="\t"))
    provenance = json.loads((Path(args.out_dir) / rows[0]["provenance_path"]).read_text(encoding="utf-8"))
    node = next(node for node in provenance["nodes"] if node.get("role") == "counts_tsv")
    assert node["access"]["landing_page_url"] == "https://example.org/counts"
