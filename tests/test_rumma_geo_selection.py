from __future__ import annotations

import csv
import json
from argparse import Namespace
from pathlib import Path

from geneset_extractors.extractors.converters.rumma_geo_selection import run


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_selects_notebook_shaped_gse193739_and_gse216043_records(tmp_path: Path) -> None:
    """Use cached query records, never a legacy GMT, for the two regression GSEs."""
    records = [
        {
            "id": "gene-record",
            "search_term": "knockdown",
            "term": "GSE193739-2-vs-1-human dn",
            "title": "GSE193739 cga siRNA",
            "sampleGroups": {"titles": {
                "1": "vcr strains sgc7901/vcr transfection cga sirna human gastric cancer cell line",
                "2": "adr nc strains sgc7901/adr transfection scrambled control sirna human gastric cancer cell line",
            }},
        },
        {
            "id": "drug-record",
            "search_term": "cisplatin",
            "term": "GSE216043-2-vs-3-human up",
            "title": "GSE216043 cisplatin",
            "sampleGroups": {"titles": {
                "2": "jar mln4924 cisplatin cell line choriocarcinoma 1um 4um 48h",
                "3": "2102ep dmso dmf cell line embryonal carcinoma 0.0025% 0.006% 48h",
            }},
        },
    ]
    query = tmp_path / "query_records.json"
    drug_terms = tmp_path / "drug_terms.json"
    query.write_text(json.dumps(records), encoding="utf-8")
    drug_terms.write_text(json.dumps(["cisplatin"]), encoding="utf-8")

    gene_dir, drug_dir = tmp_path / "HZ2", tmp_path / "HZ1"
    run(Namespace(query_records_json=str(query), drug_terms_json=None, model_id="HZ2", out_dir=str(gene_dir), provenance_overlay_json=None))
    run(Namespace(query_records_json=str(query), drug_terms_json=str(drug_terms), model_id="HZ1", out_dir=str(drug_dir), provenance_overlay_json=None))

    assert _rows(gene_dir / "selection_manifest.tsv") == [{
        "uuid": "gene-record", "source_term": "GSE193739-2-vs-1-human", "model_id": "HZ2",
        "status": "reversed", "search_term": "knockdown", "gse": "GSE193739",
        "condition_1": "2", "condition_2": "1", "species": "human", "source_direction": "dn",
        "condition_1_title": records[0]["sampleGroups"]["titles"]["2"],
        "condition_2_title": records[0]["sampleGroups"]["titles"]["1"], "title": "GSE193739 cga siRNA",
    }]
    assert _rows(drug_dir / "selection_manifest.tsv") == [{
        "uuid": "drug-record", "source_term": "GSE216043-2-vs-3-human", "model_id": "HZ1",
        "status": "signature", "search_term": "cisplatin", "gse": "GSE216043",
        "condition_1": "2", "condition_2": "3", "species": "human", "source_direction": "up",
        "condition_1_title": records[1]["sampleGroups"]["titles"]["2"],
        "condition_2_title": records[1]["sampleGroups"]["titles"]["3"], "title": "GSE216043 cisplatin",
    }]
