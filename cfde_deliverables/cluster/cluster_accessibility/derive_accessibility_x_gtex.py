#!/usr/bin/env python3
# CLUSTER: ENCODE accessibility (ATAC/DNase) consensus x GTEx tissue-enriched. Runs against the
# already-uploaded accessible-gene set dirs (no re-download). PORTABLE pure stdlib (Python 3.9+).
# Env:
#   ACCESS_DIR  = dir containing per-biosample accessible-gene set subdirs (each has geneset.tsv)
#   ENRICHED_DIR= dir containing GTEx tissue-enriched set subdirs (each has geneset.tsv) [preferred],
#                 OR set GTEX_TSTAT to a GTEx t-stat matrix (gene + per-tissue columns) instead.
#   OUTDIR      = output dir            FRAC = consensus fraction (default 0.25)   TSTAT_THR (default 4)
#   LIBNAME     = label prefix (e.g. ENCODE_ATAC_accessible)
import os, csv, json, glob, collections
ACCESS=os.environ["ACCESS_DIR"]; OUT=os.environ.get("OUTDIR","./accessibility_x_gtex")
ENRICHED=os.environ.get("ENRICHED_DIR"); TSTAT=os.environ.get("GTEX_TSTAT")
FRAC=float(os.environ.get("FRAC","0.25")); THR=float(os.environ.get("TSTAT_THR","4"))
LIB=os.environ.get("LIBNAME","ENCODE_accessible"); os.makedirs(OUT,exist_ok=True)
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")}
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]

# consensus over accessibility per-biosample sets
files=glob.glob(os.path.join(ACCESS,"*","geneset.tsv"))
freq=collections.Counter()
for fp in files:
    for g in load(fp): freq[g]+=1
n=len(files); need=max(2,int(round(FRAC*n)))
consensus={g for g,c in freq.items() if c>=need}
print(f"accessibility biosamples: {n} | consensus genes (>= {FRAC}): {len(consensus)}")

# tissue-enriched gene sets, either from enriched dirs or from a t-stat matrix
tissues={}
if ENRICHED:
    for d in glob.glob(os.path.join(ENRICHED,"*")):
        gp=os.path.join(d,"geneset.tsv")
        if os.path.exists(gp):
            t=os.path.basename(d)
            for pre in ("GTEx_tissue_enriched_",): t=t.replace(pre,"")
            tissues[t]=load(gp)
elif TSTAT:
    rows=list(csv.reader(open(TSTAT),delimiter='\t')); cols=rows[0][1:]
    for ti,t in enumerate(cols):
        tissues[t]={r[0] for r in rows[1:] if len(r)>ti+1 and r[ti+1] not in("","NA") and float(r[ti+1])>=THR}
else:
    raise SystemExit("set ENRICHED_DIR or GTEX_TSTAT")
print(f"GTEx tissues: {len(tissues)}")

def emit(name,desc,genes,t):
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"ENCODE accessibility consensus (>= {FRAC} of biosamples; NHGRI) intersected with GTEx tissue-enrichment ({t}; NIH Common Fund); public."
    json.dump({"standard_name":name,"library":LIB+"_x_GTEx","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","tissue":t,"derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"accessibility_consensus_x_gtex","inputs":["ENCODE accessibility sets (NHGRI; public)","GTEx tissue-enrichment (NIH Common Fund; public)"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE) + NIH Common Fund (GTEx)"},open(d+"/geneset.provenance.json","w"),indent=1)
nx=0
for t,enr in tissues.items():
    inter=consensus&enr
    if not inter: continue
    emit(f"{LIB}_consensus_x_GTEx_enriched_{safe(t)}",f"{LIB} consensus genes GTEx-enriched in {t}",inter,t); nx+=1
print(f"accessibility x_GTEx sets: {nx}")
