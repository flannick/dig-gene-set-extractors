#!/usr/bin/env python3
# TRIPLE multi-evidence convergence per biosample: specifically ACCESSIBLE (corrected) ∩ EXPRESSED (RNA-seq,
# same biosample) ∩ GTEx eGene (mapped tissue) = genes with three independent cis-regulatory signals agreeing.
import os, re, json, glob, csv
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"interop_intersections","output","triple_convergence"); os.makedirs(OUT,exist_ok=True)
def norm(s): return re.sub(r'[^a-z0-9]','',s.lower())
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")}
expr={}
for d in glob.glob(HOME+"/encode_rnaseq/output/ENCODE_RNAseq_expressed_*"):
    mp=os.path.join(d,"geneset.meta.json"); gp=os.path.join(d,"geneset.tsv")
    if os.path.exists(mp) and os.path.exists(gp):
        bs=json.load(open(mp)).get("biosample")
        if bs: expr[norm(bs)]=load(gp)
egene={}
for d in glob.glob(HOME+"/gtex_qtl/output_eqtl/GTEx_eQTL_eGenes_*"):
    egene[norm(os.path.basename(d).replace("GTEx_eQTL_eGenes_",""))]=load(os.path.join(d,"geneset.tsv"))
def mtis(bs):
    nb=norm(bs)
    for k in egene:
        if k and (k in nb or nb in k): return k
    return None
tot=0
for assay,base in [("ATAC","accessibility_bgcontrast/output/ENCODE_ATAC_accessible_bgcontrast"),("DNase","accessibility_bgcontrast/output/ENCODE_DNase_accessible_bgcontrast")]:
    n=0
    for d in glob.glob(os.path.join(HOME,base,"*_accessible_Up")):
        try: bs=json.load(open(os.path.join(d,"geneset.meta.json"))).get("biosample")
        except Exception: bs=None
        if not bs: continue
        nb=norm(bs); k=mtis(bs)
        if nb not in expr or not k: continue
        tri=load(os.path.join(d,"geneset.tsv")) & expr[nb] & egene[k]
        if not tri: continue
        name=f"triple_{assay}_accessible_expressed_eGene_{safe(bs)}"; dd=os.path.join(OUT,name); os.makedirs(dd,exist_ok=True); g=sorted(tri)
        open(dd+"/geneset.tsv","w").write("gene\n"+"\n".join(g)+"\n")
        open(dd+"/genesets.gmt","w").write(f"{name}\t{bs}: specifically accessible AND expressed AND GTEx eGene ({assay})\t"+"\t".join(g)+"\n")
        cite=f"Triple convergence: ENCODE {assay} corrected accessibility ∩ ENCODE RNA-seq expressed ∩ GTEx cis-eQTL eGene, {bs}. NIH-public."
        json.dump({"standard_name":name,"library":"NIH_triple_convergence","description":f"{bs}: accessible+expressed+eGene ({assay})","version":"0.1","file_type":"geneset","n_genes":len(g),"organism":"human","biosample":bs,"assay":assay,"derived_in_this_work":True,"source":cite},open(dd+"/geneset.meta.json","w"),indent=1)
        json.dump({"focus":name,"operation":"triple_convergence","inputs":[f"ENCODE {assay} corrected accessibility","ENCODE RNA-seq expressed","GTEx cis-eQTL eGenes"],"public":True,"source_citation":cite,"funding":"NIH/NHGRI + NIH Common Fund"},open(dd+"/geneset.provenance.json","w"),indent=1)
        n+=1
    print(f"{assay}: {n} triple-convergence sets")
    tot+=n
print("TOTAL triple-convergence sets:",tot)
