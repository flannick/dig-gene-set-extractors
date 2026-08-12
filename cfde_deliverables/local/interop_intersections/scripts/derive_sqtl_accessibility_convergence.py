#!/usr/bin/env python3
# NEW cross: GTEx cis-sQTL genetic regulation  x  ENCODE chromatin accessibility (background-corrected).
# Maps each ENCODE accessibility biosample -> GTEx tissue, then per matched biosample:
#   convergent_regulated_open = GTEx sGene(tissue) AND specifically-accessible(biosample)  (two independent
#                               regulatory signals agree — strong candidate cis-regulated gene)
#   eqtl_only  = sGene, NOT specifically accessible   (genetic regulation without local promoter opening / distal)
#   access_only = specifically accessible, NOT an sGene (open, no detected common cis-sQTL)
# Uses background-corrected accessibility "Up" sets. Local only (GTEx sQTL sGenes + bgcontrast Up).
import os, re, json, glob, csv
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"interop_intersections","output","sqtl_accessibility_convergence"); os.makedirs(OUT,exist_ok=True)
def norm(s): return re.sub(r'[^a-z0-9]','',s.lower())
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")}
# GTEx sQTL sGenes per tissue
egene={}; elab={}
for d in glob.glob(HOME+"/gtex_qtl/output_sqtl/GTEx_sQTL_sGenes_*"):
    t=os.path.basename(d).replace("GTEx_sQTL_sGenes_","")
    gp=os.path.join(d,"geneset.tsv")
    if os.path.exists(gp): egene[norm(t)]=load(gp); elab[norm(t)]=t
def match(bs):
    nb=norm(bs)
    for k in egene:
        if k and (k in nb or nb in k): return k
    return None
def emit(name,desc,genes,srcs,extra):
    if not genes: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="Intersection of "+" & ".join(srcs)+" (NIH-funded, public); accessibility background-corrected."
    m={"standard_name":name,"library":"NIH_sqtl_accessibility_convergence","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"intersection_of":srcs,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"sqtl_accessibility_convergence","inputs":srcs,"public":True,"source_citation":cite,"funding":"NIH Common Fund (GTEx) + NIH/NHGRI (ENCODE)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1
tot=0
for assay,base in [("ATAC","accessibility_bgcontrast/output/ENCODE_ATAC_accessible_bgcontrast"),
                   ("DNase","accessibility_bgcontrast/output/ENCODE_DNase_accessible_bgcontrast")]:
    n=0; matched=0
    for d in glob.glob(os.path.join(HOME,base,"*_accessible_Up")):
        try: bs=json.load(open(os.path.join(d,"geneset.meta.json"))).get("biosample")
        except Exception: bs=None
        if not bs: continue
        k=match(bs)
        if not k: continue
        matched+=1; A=load(os.path.join(d,"geneset.tsv")); E=egene[k]; tis=elab[k]
        src=[f"ENCODE {assay} specifically-accessible ({bs}; NHGRI)",f"GTEx cis-sQTL sGenes ({tis}; Common Fund)"]
        n+=emit(f"sQTL_{assay}_convergent_regulated_open_{safe(bs)}",f"{bs}->{tis}: GTEx sGene AND specifically accessible ({assay})",A&E,src,{"assay":assay,"biosample":bs,"tissue":tis,"class":"convergent"})
        n+=emit(f"sQTL_{assay}_eqtl_only_{safe(bs)}",f"{bs}->{tis}: GTEx sGene but NOT specifically accessible ({assay})",E-A,src,{"assay":assay,"biosample":bs,"tissue":tis,"class":"eqtl_only"})
        n+=emit(f"sQTL_{assay}_access_only_{safe(bs)}",f"{bs}->{tis}: specifically accessible but NOT a GTEx sGene ({assay})",A-E,src,{"assay":assay,"biosample":bs,"tissue":tis,"class":"access_only"})
    print(f"{assay}: {matched} biosamples matched to GTEx tissue -> {n} sets")
    tot+=n
print("TOTAL sQTL x accessibility convergence sets:",tot)
