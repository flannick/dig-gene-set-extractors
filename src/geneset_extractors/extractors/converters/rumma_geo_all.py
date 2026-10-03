"""Convenience orchestration for cached RummaGEO acquisition and reconstruction."""
from __future__ import annotations
import hashlib, json, sys
from argparse import Namespace
from pathlib import Path
from urllib.request import urlopen
from . import rumma_geo_acquire, rumma_geo_selection, rumma_geo

SIGCOM_URL = "https://s3.dev.maayanlab.cloud/sigcom-lincs/ranker/signatures_meta.json"

def _sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def _sources(paths): return {role:{"url":Path(path).resolve().as_uri(),"version":"sha256:"+_sha(path)} for role,path in paths.items()}

def run(args):
    out=Path(args.out_dir).resolve(); provenance=out/"provenance"; provenance.mkdir(parents=True,exist_ok=True)
    meta=Path(args.signatures_meta_json).resolve() if args.signatures_meta_json else provenance/"signatures_meta.json"
    if args.refresh_sources or not meta.is_file():
        print("[rumma_geo_all] acquiring SigCom signatures metadata", file=sys.stderr, flush=True)
        with urlopen(SIGCOM_URL, timeout=120) as response: meta.write_bytes(response.read())
    else: print("[rumma_geo_all] reusing cached SigCom signatures metadata", file=sys.stderr, flush=True)
    raw=json.loads(meta.read_text(encoding="utf-8")); terms=sorted({str(x.get("pert_name","")).strip() for x in raw.values() if str(x.get("pert_name","")).strip()})
    drug_terms=provenance/"sigcom_lincs_drug_terms.json"; drug_terms.write_text(json.dumps(terms,indent=2)+"\n",encoding="utf-8")
    common={"human_gmt":args.human_gmt,"mouse_gmt":args.mouse_gmt,"human_gene_info":args.human_gene_info,"mouse_gene_info":args.mouse_gene_info,"gene_orthologs":args.gene_orthologs}
    requested = ("HZ2", "HZ1") if args.models == "all" or args.prepare_sources else tuple(args.models.split(","))
    for model, label in (("HZ2","gene_perturbations"),("HZ1","drug_perturbations")):
        if model not in requested: continue
        cache=provenance/label/"acquisition"; query=cache/"rummageo_query_records.json"
        if args.refresh_sources or not query.is_file():
            print(f"[rumma_geo_all] acquiring {label} query records", file=sys.stderr, flush=True)
            rumma_geo_acquire.run(Namespace(model_id=model,drug_terms_json=str(drug_terms) if model=="HZ1" else None,out_dir=str(cache),endpoint=args.endpoint,page_size=args.page_size))
        else: print(f"[rumma_geo_all] reusing cached {label} query records", file=sys.stderr, flush=True)
        if args.prepare_sources: continue
        print(f"[rumma_geo_all] selecting {label} signatures", file=sys.stderr, flush=True)
        selection=provenance/label/"selection"; rumma_geo_selection.run(Namespace(query_records_json=str(query),drug_terms_json=str(drug_terms) if model=="HZ1" else None,model_id=model,out_dir=str(selection),provenance_overlay_json=None))
        roles={"human_rummageo_gmt":args.human_gmt,"mouse_rummageo_gmt":args.mouse_gmt,"recorded_selection_manifest":str(selection/"selection_manifest.tsv"),"ncbi_human_gene_info":args.human_gene_info,"ncbi_mouse_gene_info":args.mouse_gene_info,"ncbi_gene_orthologs":args.gene_orthologs}
        manifest=provenance/label/"source_manifest.json"; manifest.write_text(json.dumps({"sources":_sources(roles)},indent=2)+"\n")
        print(f"[rumma_geo_all] reconstructing {label} GMT", file=sys.stderr, flush=True)
        rumma_geo.run(Namespace(**common,selection_manifest=str(selection/"selection_manifest.tsv"),source_manifest=str(manifest),model_id=model,out_dir=str(out/label),genome_build="hg38",min_genes=5,gmt_description="RummaGEO perturbation signature",legacy_gmt=None,provenance_overlay_json=None,provenance_mirror_local_prefix=None,provenance_mirror_remote_prefix=None))
        print(f"[rumma_geo_all] completed {label}", file=sys.stderr, flush=True)
    if args.prepare_sources: (provenance / "sources_prepared.ok").write_text("ok\n", encoding="utf-8")
    return {"n_peaks":0,"n_genes":0,"out_dir":str(out)}
