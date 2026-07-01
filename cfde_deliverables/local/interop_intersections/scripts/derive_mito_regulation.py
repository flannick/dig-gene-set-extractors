#!/usr/bin/env python3
# TF regulon ∩ MitoCarta -> mitochondrial-regulation gene sets: for each nuclear master-regulator of mito
# biogenesis, the MitoCarta genes it binds (ENCODE TF ChIP). NOTE: inherits the pending ChIP-contrast
# ruling (TF regulons are absolute peak->gene); not for submission until that's resolved.
import os, json, glob, csv
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"interop_intersections","output","mito_regulation")
TFDIR=os.path.join(HOME,"encode_regulons","output_tf","regulon")
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")} if os.path.exists(fp) else set()
mito={l.split("\t")[0] for i,l in enumerate(open(HOME+"/mito_genesets/output/mitocarta3_reference.tsv")) if i>0 and l.strip()}
# established nuclear regulators of mitochondrial biogenesis/function
MASTERS=["NRF1","GABPA","GABPB1","ESRRA","ESRRG","YY1","TFAM","TFB2M","PPARGC1A","PPARGC1B",
         "NFE2L2","MYC","MAX","CREB1","PPRC1","SP1","CEBPB","FOXO1","FOXO3","MEF2A"]
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(name,desc,genes,srcs,extra):
    if not genes: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="Intersection of "+" & ".join(srcs)+" (all NIH-funded, public). Mito-regulation view; TF regulons are ENCODE TF ChIP peak->gene (pending ChIP-contrast ruling)."
    m={"standard_name":name,"library":"NIH_mito_regulation","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"intersection_of":srcs,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"tf_regulon_x_mitocarta","inputs":srcs,"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE) + NIH/NIGMS (MitoCarta)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1
n=0; covered=set()
for tf in MASTERS:
    reg=load(os.path.join(TFDIR,f"ENCODE_TF_regulon_{tf}","geneset.tsv"))
    if not reg: continue
    inter=reg&mito
    if not inter: continue
    covered|=inter
    n+=emit(f"MitoCarta_x_{safe(tf)}_TFtargets",f"Mitochondrial (MitoCarta3.0) genes bound by {tf} (ENCODE TF ChIP)",inter,
            [f"ENCODE TF ChIP regulon: {tf} (NHGRI)","MitoCarta3.0 (NIGMS)"],{"tf":tf})
emit("MitoCarta_regulated_by_any_master_TF","Mitochondrial genes bound by >=1 mito master-regulator TF (ENCODE)",covered,
     ["ENCODE TF ChIP regulons (mito masters; NHGRI)","MitoCarta3.0 (NIGMS)"],{})
print(f"mito-regulation sets: {n} per-TF + 1 union | mito genes covered by >=1 master TF: {len(covered)}/{len(mito)}")
