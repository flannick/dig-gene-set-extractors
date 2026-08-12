#!/usr/bin/env python3
# Per cell-type/tissue: chromatin accessibility vs expression CONCORDANCE, matched on the SAME ENCODE
# biosample. Run in TWO modes for comparison:
#   nocontrol = raw accessibility (promoter peak present)              vs expression
#   control   = background-corrected "specifically accessible" (Up)   vs expression
# Per biosample+mode:
#   concordant_active     = accessible AND expressed
#   open_but_silent       = accessible, NOT expressed        (poised / accessibility precedes expression)
#   expressed_but_closed  = expressed, NOT accessible        (distal regulation / expression-accessibility latency)
# Local only: ENCODE accessibility zips (raw) + bgcontrast Up dirs (control) + ENCODE RNA-seq expressed.
import os, re, json, glob, zipfile
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"concordance_genesets","output"); os.makedirs(OUT,exist_ok=True)
def norm(s): return re.sub(r'[^a-z0-9]','',s.lower())
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]

# expression: ENCODE RNA-seq expressed per biosample
expr={}
for d in glob.glob(HOME+"/encode_rnaseq/output/ENCODE_RNAseq_expressed_*"):
    mp=os.path.join(d,"geneset.meta.json"); gp=os.path.join(d,"geneset.tsv")
    if not (os.path.exists(mp) and os.path.exists(gp)): continue
    bs=json.load(open(mp)).get("biosample")
    if bs: expr[norm(bs)]={ln.strip() for ln in open(gp).read().splitlines()[1:] if ln.strip()}
print("RNA-seq biosamples:",len(expr))

def load_zip(zp):   # raw accessible per biosample
    out={},{}; acc={}; lab={}; z=zipfile.ZipFile(zp)
    for n in z.namelist():
        if not n.endswith("geneset.tsv"): continue
        try: bs=json.loads(z.read(n.rsplit("/",1)[0]+"/geneset.meta.json")).get("biosample")
        except: bs=None
        if bs: acc[norm(bs)]={ln.strip() for ln in z.read(n).decode("utf-8","replace").splitlines()[1:] if ln.strip()}; lab[norm(bs)]=bs
    return acc,lab
def load_updirs(base):   # control: specifically-accessible Up per biosample
    acc={}; lab={}
    for d in glob.glob(base+"/*_accessible_Up"):
        mp=os.path.join(d,"geneset.meta.json"); gp=os.path.join(d,"geneset.tsv")
        if not (os.path.exists(mp) and os.path.exists(gp)): continue
        bs=json.load(open(mp)).get("biosample")
        if bs: acc[norm(bs)]={ln.strip() for ln in open(gp).read().splitlines()[1:] if ln.strip()}; lab[norm(bs)]=bs
    return acc,lab

def emit(mode,name,desc,genes,srcs,extra):
    if not genes: return 0
    d=os.path.join(OUT,mode,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"Accessibility×expression concordance ({mode}; same ENCODE biosample): "+" & ".join(srcs)+" (NIH/NHGRI, public)."
    m={"standard_name":name,"library":f"ENCODE_accessibility_expression_concordance_{mode}","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","control_mode":mode,"derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"accessibility_expression_concordance","control_mode":mode,"inputs":srcs,"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1

SRC={"ATAC":("ENCODE_ATAC_accessible_genes_20260630.zip","accessibility_bgcontrast/output/ENCODE_ATAC_accessible_bgcontrast"),
     "DNase":("ENCODE_DNase_accessible_genes_20260630.zip","accessibility_bgcontrast/output/ENCODE_DNase_accessible_bgcontrast")}
summary={}
for assay,(zp,updir) in SRC.items():
    for mode,loader in [("raw",     lambda: load_zip(os.path.join(HOME,zp))),
                        ("control", lambda: load_updirs(os.path.join(HOME,updir)))]:
        acc,lab=loader()
        matched=set(acc)&set(expr); n=0
        for k in matched:
            A=acc[k]; E=expr[k]; bs=lab[k]
            src=[f"ENCODE {assay} {'raw-accessible (uncontrolled)' if mode=='raw' else 'bg-corrected specifically-accessible'} ({bs}; NHGRI)",f"ENCODE RNA-seq expressed ({bs}; NHGRI)"]
            n+=emit(mode,f"ENCODE_{assay}_{mode}_concordant_active_{safe(bs)}",f"{bs}: accessible AND expressed ({assay} {mode})",A&E,src,{"assay":assay,"biosample":bs,"class":"concordant_active"})
            n+=emit(mode,f"ENCODE_{assay}_{mode}_open_but_silent_{safe(bs)}",f"{bs}: accessible but NOT expressed — poised ({assay} {mode})",A-E,src,{"assay":assay,"biosample":bs,"class":"open_silent"})
            n+=emit(mode,f"ENCODE_{assay}_{mode}_expressed_but_closed_{safe(bs)}",f"{bs}: expressed but NOT accessible — distal/latency ({assay} {mode})",E-A,src,{"assay":assay,"biosample":bs,"class":"expressed_closed"})
        summary[f"{assay}/{mode}"]=(len(matched),n)
        print(f"{assay} {mode}: {len(matched)} matched -> {n} sets")
print("summary:",summary)
