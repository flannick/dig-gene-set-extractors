#!/usr/bin/env python3
# Corrected (de-biased) TF/eCLIP specific_Up regulons ∩ MitoCarta3.0 -> per-factor SPECIFICALLY-BOUND
# mitochondrial genes (promiscuous/HOT-region binding removed). Union per assay = mito genes specifically
# bound by any factor. Complements the earlier raw mito_regulation with the corrected regulons.
import os, json, glob, csv
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"interop_intersections","output","corrected_regulon_mito"); os.makedirs(OUT,exist_ok=True)
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")} if os.path.exists(fp) else set()
mito={l.split("\t")[0] for i,l in enumerate(open(HOME+"/mito_genesets/output/mitocarta3_reference.tsv")) if i>0 and l.strip()}
print("MitoCarta genes:",len(mito))
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(name,desc,genes,srcs,extra):
    if not genes: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="Intersection of "+" & ".join(srcs)+" (NIH-funded, public); binding background-corrected."
    m={"standard_name":name,"library":"NIH_corrected_regulon_mito","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"intersection_of":srcs,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"corrected_regulon_x_mitocarta","inputs":srcs,"public":True,"source_citation":cite,"funding":"NIH/NHGRI (ENCODE) + NIH/NIGMS (MitoCarta)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1
tot=0
for assay,base in [("TF","ENCODE_TF_regulon_bgcorrected"),("eCLIP","ENCODE_eCLIP_regulon_bgcorrected")]:
    n=0; cov=set()
    for d in glob.glob(os.path.join(HOME,"encode_regulons_bgcontrast","output",base,"*_specific_Up")):
        try: tgt=json.load(open(os.path.join(d,"geneset.meta.json"))).get("target",os.path.basename(d))
        except Exception: tgt=os.path.basename(d)
        inter=load(os.path.join(d,"geneset.tsv"))&mito
        if not inter: continue
        cov|=inter
        n+=emit(f"ENCODE_{assay}_{safe(tgt)}_specific_mito_targets",f"{tgt}: specifically-bound ({assay}) mitochondrial (MitoCarta) genes",inter,[f"ENCODE {assay} corrected regulon: {tgt} (NHGRI)","MitoCarta3.0 (NIGMS)"],{"assay":assay,"target":tgt})
    n+=emit(f"ENCODE_{assay}_any_factor_specific_mito_targets",f"MitoCarta genes specifically bound by >=1 {assay} factor",cov,[f"ENCODE {assay} corrected regulons (NHGRI)","MitoCarta3.0 (NIGMS)"],{"assay":assay})
    print(f"{assay}: {n} sets | union mito targets: {len(cov)}")
    tot+=n
print("TOTAL corrected-regulon x MitoCarta sets:",tot)
