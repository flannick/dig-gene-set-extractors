#!/usr/bin/env python3
# Catch "background sneaking": per biosample, the genes whose accessibility is present in RAW but removed by
# the background CONTROL (A_raw - A_control), split by expression to disambiguate:
#   background_sneak_silent    = (raw-only accessible) AND NOT expressed  -> true ubiquitous/background open (control correctly removes)
#   background_absorbed_expressed = (raw-only accessible) AND expressed   -> real accessible+expressed genes the control over-removed
# Comparing these tells you how much discordance is background artifact vs the control over-correcting real signal.
import os, re, json, glob, zipfile
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"concordance_genesets","output","background_delta"); os.makedirs(OUT,exist_ok=True)
def norm(s): return re.sub(r'[^a-z0-9]','',s.lower())
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
expr={}
for d in glob.glob(HOME+"/encode_rnaseq/output/ENCODE_RNAseq_expressed_*"):
    mp=os.path.join(d,"geneset.meta.json"); gp=os.path.join(d,"geneset.tsv")
    if os.path.exists(mp) and os.path.exists(gp):
        bs=json.load(open(mp)).get("biosample")
        if bs: expr[norm(bs)]={ln.strip() for ln in open(gp).read().splitlines()[1:] if ln.strip()}
def load_zip(zp):
    acc={}; lab={}; z=zipfile.ZipFile(zp)
    for n in z.namelist():
        if not n.endswith("geneset.tsv"): continue
        try: bs=json.loads(z.read(n.rsplit("/",1)[0]+"/geneset.meta.json")).get("biosample")
        except: bs=None
        if bs: acc[norm(bs)]={ln.strip() for ln in z.read(n).decode('utf-8','replace').splitlines()[1:] if ln.strip()}; lab[norm(bs)]=bs
    return acc,lab
def load_updirs(base):
    acc={}; lab={}
    for d in glob.glob(base+"/*_accessible_Up"):
        gp=os.path.join(d,"geneset.tsv"); mp=os.path.join(d,"geneset.meta.json")
        if os.path.exists(gp) and os.path.exists(mp):
            bs=json.load(open(mp)).get("biosample")
            if bs: acc[norm(bs)]={ln.strip() for ln in open(gp).read().splitlines()[1:] if ln.strip()}; lab[norm(bs)]=bs
    return acc,lab
def emit(name,desc,genes,assay,bs,cls):
    if not genes: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"Raw-minus-control accessibility delta (background sneak) x expression, {assay} {bs}; ENCODE/NHGRI public."
    json.dump({"standard_name":name,"library":"ENCODE_background_sneak_delta","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","assay":assay,"biosample":bs,"class":cls,"derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"background_sneak_delta","inputs":[f"ENCODE {assay} raw accessible","ENCODE {assay} bg-corrected accessible","ENCODE RNA-seq expressed"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1
SRC={"ATAC":("ENCODE_ATAC_accessible_genes_20260630.zip","accessibility_bgcontrast/output/ENCODE_ATAC_accessible_bgcontrast"),
     "DNase":("ENCODE_DNase_accessible_genes_20260630.zip","accessibility_bgcontrast/output/ENCODE_DNase_accessible_bgcontrast")}
tot=0
for assay,(zp,updir) in SRC.items():
    raw,lab=load_zip(os.path.join(HOME,zp)); ctrl,_=load_updirs(os.path.join(HOME,updir))
    n=0
    for k in set(raw)&set(ctrl)&set(expr):
        delta=raw[k]-ctrl[k]; bs=lab[k]; E=expr[k]
        n+=emit(f"ENCODE_{assay}_bgsneak_silent_{safe(bs)}",f"{bs}: accessible in RAW only (removed by control) and NOT expressed — background sneak ({assay})",delta-E,assay,bs,"background_sneak_silent")
        n+=emit(f"ENCODE_{assay}_bgabsorbed_expressed_{safe(bs)}",f"{bs}: accessible in RAW only (removed by control) but EXPRESSED — real signal absorbed by background ({assay})",delta&E,assay,bs,"background_absorbed_expressed")
    print(f"{assay}: {len(set(raw)&set(ctrl)&set(expr))} biosamples -> {n} delta sets")
    tot+=n
print("TOTAL background-sneak delta sets:",tot)
