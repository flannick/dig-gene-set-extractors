import json
from pathlib import Path

from geneset_extractors.external_import import import_external_gmt


def test_external_import_copies_gmt_and_emits_paired_provenance(tmp_path: Path):
    source_gmt = tmp_path / "source.gmt"
    source_gmt.write_text("Example\tdescription\tGENE1\tGENE2\n", encoding="utf-8")
    source = tmp_path / "source.json"
    source.write_text(json.dumps({
        "name": "External release", "uri_or_identifier": "doi:test", "release": "v1",
        "license": "CC-BY-4.0", "access_restrictions": "public", "organism": "human",
        "genome_build": "hg38", "assay": "rna_seq", "data_type": "precomputed_gene_sets",
        "documentation": "https://example.org/methods",
    }), encoding="utf-8")
    out = tmp_path / "out"
    result = import_external_gmt(gmt=source_gmt, out_dir=out, source_record=source, library_id="External", model_id="M1", display_name="Model 1", description="Imported unchanged", expected_sha256="sha256:" + __import__("hashlib").sha256(source_gmt.read_bytes()).hexdigest())
    assert result["n_gene_sets"] == 1
    assert (out / "genesets.gmt").read_bytes() == source_gmt.read_bytes()
    assert (out / "geneset.provenance.legacy.json").exists()
    assert (out / "geneset.provenance.dapper.yaml").exists()
    metadata = json.loads((out / "geneset.meta.json").read_text(encoding="utf-8"))
    assert metadata["external_import"]["regeneration_status"] == "incomplete_code"
