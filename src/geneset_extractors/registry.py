from __future__ import annotations

import json
from pathlib import Path

from geneset_extractors.extractors.converters import (
    atac_bulk,
    atac_bulk_matrix,
    atac_sc_10x,
    chipseq_peak,
    cnv_gene_extractor,
    drug_response_screen,
    methylation_cpg_diff,
    methylation_dmr,
    methylation_dmr_regions,
    morphology_profile_query,
    proteomics_diff,
    ptm_site_diff,
    ptm_site_matrix,
    rna_deg,
    rna_deg_multi,
    rna_sc_programs,
    sc_rna_marker,
)


CONVERTERS = {
    "atac_bulk": atac_bulk,
    "atac_bulk_matrix": atac_bulk_matrix,
    "atac_sc_10x": atac_sc_10x,
    "chipseq_peak": chipseq_peak,
    "cnv_gene_extractor": cnv_gene_extractor,
    "drug_response_screen": drug_response_screen,
    "methylation_cpg_diff": methylation_cpg_diff,
    "methylation_dmr": methylation_dmr,
    "methylation_dmr_regions": methylation_dmr_regions,
    "morphology_profile_query": morphology_profile_query,
    "proteomics_diff": proteomics_diff,
    "ptm_site_diff": ptm_site_diff,
    "ptm_site_matrix": ptm_site_matrix,
    "rna_deg": rna_deg,
    "rna_deg_multi": rna_deg_multi,
    "rna_sc_programs": rna_sc_programs,
    "sc_rna_marker": sc_rna_marker,
}


def _schemas_dir() -> Path:
    return Path(__file__).parent / "schemas"


def _specs_dir() -> Path:
    return Path(__file__).parent / "extractors" / "converters" / "specs"


def _validate_converter_spec(payload: dict[str, object], schema: dict[str, object]) -> None:
    try:
        import jsonschema  # type: ignore
    except ModuleNotFoundError:
        required = ["name", "description", "inputs", "parameters", "outputs"]
        missing = [k for k in required if k not in payload]
        if missing:
            raise ValueError(f"converter spec missing required fields: {', '.join(missing)}")
    else:
        jsonschema.validate(payload, schema)


def list_converters() -> list[tuple[str, str]]:
    out = []
    for name in sorted(CONVERTERS):
        spec = get_converter_spec(name)
        out.append((name, str(spec.get("description", ""))))
    return out


def get_converter(name: str):
    if name not in CONVERTERS:
        raise KeyError(f"Unknown converter: {name}")
    return CONVERTERS[name]


def get_converter_spec(name: str) -> dict[str, object]:
    p = _specs_dir() / f"{name}.json"
    if not p.exists():
        raise FileNotFoundError(f"Missing converter spec: {p}")
    payload = json.loads(p.read_text(encoding="utf-8"))
    schema = json.loads((_schemas_dir() / "converter_spec.schema.json").read_text(encoding="utf-8"))
    _validate_converter_spec(payload, schema)
    return payload
