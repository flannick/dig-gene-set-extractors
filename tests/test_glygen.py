import csv
import json
from argparse import Namespace

from geneset_extractors.extractors.converters.glygen import run_glycan_synthesizing_enzymes, run_glycosylated_proteins


def _csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["uniprotkb_canonical_ac", "glytoucan_ac"])
        writer.writeheader(); writer.writerows(rows)


def test_glycosylated_proteins_normalizes_filters_and_sorts(tmp_path):
    rows = [{"uniprotkb_canonical_ac": f"P{i}", "glytoucan_ac": "GKEEP"} for i in range(1, 6)]
    rows += [{"uniprotkb_canonical_ac": "P1", "glytoucan_ac": "GKEEP"}, {"uniprotkb_canonical_ac": "P404", "glytoucan_ac": "GKEEP"}]
    _csv(tmp_path / "unicarbkb.csv", rows)
    _csv(tmp_path / "harvard.csv", [{"uniprotkb_canonical_ac": f"P{i}", "glytoucan_ac": "GDROP"} for i in range(1, 5)])
    _csv(tmp_path / "glyconnect.csv", [{"uniprotkb_canonical_ac": "P1", "glytoucan_ac": ""}])
    with (tmp_path / "master.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["uniprotkb_canonical_ac", "gene_name"]); writer.writeheader()
        writer.writerows([{ "uniprotkb_canonical_ac": f"P{i}", "gene_name": f" gene{i} "} for i in range(1, 6)])
    out = tmp_path / "out"
    result = run_glycosylated_proteins(Namespace(unicarbkb=tmp_path / "unicarbkb.csv", harvard=tmp_path / "harvard.csv", glyconnect=tmp_path / "glyconnect.csv", masterlist=tmp_path / "master.csv", out_dir=out, model_id="HZ1", min_genes=5, genome_build="hg38", gmt_description="test", provenance_overlay_json=None))
    assert result["n_gene_sets"] == 1
    assert (out / "genesets.gmt").read_text() == "GKEEP\ttest\tGENE1\tGENE2\tGENE3\tGENE4\tGENE5\n"


def test_glycan_enzymes_uses_cached_human_records_only(tmp_path):
    cache = tmp_path / "G00024MO.json"
    cache.write_text(json.dumps({"enzyme": [{"tax_id": 9606, "gene": "poglut2"}, {"tax_id": "9606", "gene": "POGLUT1"}, {"tax_id": 9606, "gene": "POGLUT3"}, {"tax_id": 9606, "gene": "POGLUT1"}, {"tax_id": 10090, "gene": "Poglut1"}, {"tax_id": 9606, "gene": ""}]}))
    missing = tmp_path / "empty.json"; missing.write_text("{}")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("accession\tstatus\tcache_file\nG00024MO\tcached\tG00024MO.json\nGEMPTY\tcached\tempty.json\n")
    out = tmp_path / "out"
    result = run_glycan_synthesizing_enzymes(Namespace(cache_manifest=manifest, out_dir=out, model_id="HZ2", genome_build="hg38", gmt_description="", provenance_overlay_json=None))
    assert result["n_gene_sets"] == 1
    assert (out / "genesets.gmt").read_text() == "glytoucan:G00024MO\t\tPOGLUT1\tPOGLUT2\tPOGLUT3\n"
