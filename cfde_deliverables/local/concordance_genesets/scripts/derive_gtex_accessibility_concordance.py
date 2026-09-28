#!/usr/bin/env python3
# Concordance/NONconcordance of GTEx tissue expression vs ENCODE chromatin accessibility, per tissue.
# Maps each ENCODE accessibility biosample to a GTEx tissue (normalized name containment), then per matched
# biosample (in BOTH raw and control accessibility modes):
#   concordant_active        = accessible AND GTEx-tissue-enriched
#   open_but_not_enriched     = accessible, NOT GTEx-enriched     (poised / accessibility ahead of expression)
#   enriched_but_not_accessible = GTEx-enriched, NOT accessible   (distal regulation / expression-accessibility latency)
# Local only: GTEx tissue-enriched sets + ENCODE accessibility (raw zips + bgcontrast Up dirs).
import os, re, json, glob, zipfile
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"concordance_genesets","output_gtex"); os.makedirs(OUT,exist_ok=True)
def norm(s): return re.sub(r'[^a-z0-9]','',s.lower())
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
# GTEx tissue-enriched
gt={}; gtlab={}
for d in glob.glob(HOME+"/gtex_tissue_enriched/output/GTEx_tissue_enriched_*"):
    t=os.path.basename(d).replace("GTEx_tissue_enriched_","")
    gp=os.path.join(d,"geneset.tsv")
    if os.path.exists(gp): gt[norm(t)]={ln.strip() for ln in open(gp).read().splitlines()[1:] if ln.strip()}; gtlab[norm(t)]=t
def match_tissue(bs):
    nb=norm(bs)
    for k in gt:
        if k and (k in nb or nb in k): return k
    return None
def load_zip(zp):
    acc={}; lab={}; z=zipfile.ZipFile(zp)
    for n in z.namelist():
        if not n.endswith("geneset.tsv"): continue
        try: bs=json.loads(z.read(n.rsplit("/",1)[0]+"/geneset.meta.json")).get("biosample")
        except: bs=None
        if bs: acc[bs]={ln.strip() for ln in z.read(n).decode('utf-8','replace').splitlines()[1:] if ln.strip()}
    return acc
def load_updirs(base):
    acc={}
    for d in glob.glob(base+"/*_accessible_Up"):
        gp=os.path.join(d,"geneset.tsv"); mp=os.path.join(d,"geneset.meta.json")
        if os.path.exists(gp) and os.path.exists(mp):
            bs=json.load(open(mp)).get("biosample")
            if bs: acc[bs]={ln.strip() for ln in open(gp).read().splitlines()[1:] if ln.strip()}
    return acc
def emit(mode,name,desc,genes,assay,bs,tis,cls):
    if not genes: return 0
    d=os.path.join(OUT,mode,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"GTEx expression x ENCODE {assay} accessibility concordance ({mode}); biosample {bs} -> GTEx {tis}. ENCODE/NHGRI + GTEx/Common Fund; public."
    json.dump({"standard_name":name,"library":f"GTEx_x_accessibility_concordance_{mode}","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","assay":assay,"biosample":bs,"gtex_tissue":tis,"control_mode":mode,"class":cls,"derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"gtex_accessibility_concordance","control_mode":mode,"inputs":[f"ENCODE {assay} accessibility ({mode})","GTEx tissue-enriched"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE) + NIH Common Fund (GTEx)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1
SRC={"ATAC":("ENCODE_ATAC_accessible_genes_20260630.zip","accessibility_bgcontrast/output/ENCODE_ATAC_accessible_bgcontrast"),
     "DNase":("ENCODE_DNase_accessible_genes_20260630.zip","accessibility_bgcontrast/output/ENCODE_DNase_accessible_bgcontrast")}
summary={}
for assay,(zp,updir) in SRC.items():
    for mode,acc in [("raw",load_zip(os.path.join(HOME,zp))),("control",load_updirs(os.path.join(HOME,updir)))]:
        n=0; matched=0
        for bs,A in acc.items():
            k=match_tissue(bs)
            if not k: continue
            matched+=1; E=gt[k]; tis=gtlab[k]
            n+=emit(mode,f"GTEx_{assay}_{mode}_concordant_{safe(bs)}",f"{bs}->{tis}: accessible AND GTEx-enriched ({assay} {mode})",A&E,assay,bs,tis,"concordant")
            n+=emit(mode,f"GTEx_{assay}_{mode}_open_not_enriched_{safe(bs)}",f"{bs}->{tis}: accessible but NOT GTEx-enriched — poised ({assay} {mode})",A-E,assay,bs,tis,"open_not_enriched")
            n+=emit(mode,f"GTEx_{assay}_{mode}_enriched_not_accessible_{safe(bs)}",f"{bs}->{tis}: GTEx-enriched but NOT accessible — distal/latency ({assay} {mode})",E-A,assay,bs,tis,"enriched_not_accessible")
        summary[f"{assay}/{mode}"]=(matched,n); print(f"{assay} {mode}: {matched} biosamples matched to GTEx -> {n} sets")
print("summary:",summary)
