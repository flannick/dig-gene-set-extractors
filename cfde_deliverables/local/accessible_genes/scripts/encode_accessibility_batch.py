#!/usr/bin/env python3
# ENCODE accessibility -> per-biosample accessible-gene sets (promoter TSS+/-1kb peak overlap).
# PORTABLE: pure Python stdlib, runs on Python 3.9+ (no toolchain). Lean: download->derive->DELETE peak.
# Runs LOCALLY or ON THE CLUSTER (only needs internet to encodeproject.org + UCSC).
# Config via env:
#   ENC_ASSAY   (default "ATAC-seq"; also "DNase-seq")
#   ENC_OUTPUT  (default "IDR thresholded peaks"; for DNase use "peaks")
#   ENC_LIB     (default "ENCODE_ATAC_accessible")
#   OUTDIR      (default ./output)   TMPDIR (default ./tmp)
#   REFGENE     (default $TMPDIR/refGene_hg38.txt.gz; auto-downloaded if missing)
import gzip, json, os, collections, urllib.request, urllib.parse
ASSAY=os.environ.get("ENC_ASSAY","ATAC-seq"); OUTPUT=os.environ.get("ENC_OUTPUT","IDR thresholded peaks")
LIB=os.environ.get("ENC_LIB","ENCODE_ATAC_accessible")
OUT=os.environ.get("OUTDIR","./output"); TMP=os.environ.get("TMPDIR","./tmp")
os.makedirs(OUT,exist_ok=True); os.makedirs(TMP,exist_ok=True)
REFGENE=os.environ.get("REFGENE",os.path.join(TMP,"refGene_hg38.txt.gz"))
WIN=BIN=1000
LOG=open(os.path.join(OUT,"batch_log.txt"),"a")
def log(m): LOG.write(m+"\n"); LOG.flush()

if not os.path.exists(REFGENE):
    log("downloading UCSC refGene hg38 ..."); urllib.request.urlretrieve(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz", REFGENE)
MODE=os.environ.get("ENC_MODE","promoter")   # promoter (TSS+/-WIN) | body (full gene span)
DESC=os.environ.get("ENC_DESC","with promoter-proximal")
gene_reg=collections.defaultdict(set)
for line in gzip.open(REFGENE,'rt'):
    f=line.rstrip('\n').split('\t')
    if len(f)<13: continue
    ch,st,a,b,g=f[2],f[3],f[4],f[5],f[12]
    if '_' in ch or not a.isdigit() or not b.isdigit(): continue
    a=int(a); b=int(b)
    if MODE=="body": gene_reg[g].add((ch,a,b))
    else: tss=a if st=='+' else b; gene_reg[g].add((ch,tss-WIN,tss+WIN))
log(f"refGene genes: {len(gene_reg)} (mode={MODE})")

_params={"type":"File","assay_title":ASSAY,"file_format":"bed","output_type":OUTPUT,
   "assembly":"GRCh38","status":"released","limit":"all","format":"json"}
if os.environ.get("ENC_TARGET"): _params["target.label"]=os.environ["ENC_TARGET"]  # e.g. H3K4me3
q=urllib.parse.urlencode(_params)
d=json.load(urllib.request.urlopen(urllib.request.Request("https://www.encodeproject.org/search/?"+q,
   headers={"Accept":"application/json","User-Agent":"x"}),timeout=300))
seen={}
for f in d.get("@graph",[]):
    bo=f.get("biosample_ontology",{}); bs=bo.get("term_name") if isinstance(bo,dict) else None
    if bs and bs not in seen: seen[bs]=(f.get("accession"),f.get("dataset","").strip("/").split("/")[-1],f.get("href"))
log(f"{ASSAY} distinct biosamples: {len(seen)}")

def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def accessible(p):
    cov=collections.defaultdict(set)
    for line in gzip.open(p,'rt'):
        f=line.rstrip('\n').split('\t')
        if len(f)<3 or not f[1].isdigit(): continue
        for bb in range(int(f[1])//BIN,int(f[2])//BIN+1): cov[f[0]].add(bb)
    out=[]
    for g,regs in gene_reg.items():
        for c,s,e in regs:
            cb=cov.get(c)
            if cb and any(bb in cb for bb in range(s//BIN,e//BIN+1)): out.append(g); break
    return sorted(set(out))
done=fail=0
for bs,(acc,exp,href) in seen.items():
    if not href: continue
    name=f"{LIB}_{safe(bs)}"; dd=os.path.join(OUT,name)
    if os.path.exists(os.path.join(dd,"geneset.tsv")): done+=1; continue
    pk=os.path.join(TMP,acc+".bed.gz")
    try:
        urllib.request.urlretrieve("https://www.encodeproject.org"+href, pk)
        genes=accessible(pk); os.remove(pk); os.makedirs(dd,exist_ok=True)
        cite=f"Derived from ENCODE {ASSAY} peak {acc} (exp {exp}, biosample {bs}), GRCh38, ENCODE/NHGRI, public; UCSC refGene hg38."
        open(os.path.join(dd,"geneset.tsv"),"w").write("gene\n"+"\n".join(genes)+"\n")
        open(os.path.join(dd,"genesets.gmt"),"w").write(f"{name}\tGenes {DESC} ENCODE {ASSAY} peaks in {bs}\t"+"\t".join(genes)+"\n")
        json.dump({"standard_name":name,"library":LIB,"description":f"Genes {DESC} ENCODE {ASSAY} peaks in {bs} (DERIVED from ENCODE public peaks; mode={MODE}).","version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","assembly":"GRCh38","biosample":bs,"assay":ASSAY,"mark":os.environ.get("ENC_TARGET",""),"overlap_mode":MODE,"derived_in_this_work":True,"source":cite},open(os.path.join(dd,"geneset.meta.json"),"w"),indent=1)
        json.dump({"focus":name,"operation":"derive_accessible_genes","inputs":[f"ENCODE {acc} {ASSAY} peaks (NHGRI; public)","UCSC refGene hg38"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"},open(os.path.join(dd,"geneset.provenance.json"),"w"),indent=1)
        done+=1; log(f"[{done}] {bs}: {len(genes)} genes")
    except Exception as e:
        fail+=1; log(f"FAIL {bs} {acc}: {str(e)[:80]}")
        try: os.remove(pk)
        except OSError: pass
log(f"=== DONE {ASSAY}: {done} sets, {fail} failures ===")
LOG.close()
