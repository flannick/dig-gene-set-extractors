#!/usr/bin/env python3
# D: ENCODE RNA-seq consensus expressed genes (expressed in >= FRAC of biosamples) x GTEx tissue-enriched.
import csv, json, os, glob, collections
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
GTEX="/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv"
OUT=os.path.join(HOME,"encode_rnaseq","x_GTEx"); os.makedirs(OUT,exist_ok=True)
FRAC=float(os.environ.get("FRAC","0.25")); THR=float(os.environ.get("TSTAT_THR","4"))
files=glob.glob(os.path.join(HOME,"encode_rnaseq","output","ENCODE_RNAseq_expressed_*","geneset.tsv"))
freq=collections.Counter()
for fp in files:
    for r in csv.DictReader(open(fp),delimiter='\t'):
        if r.get("gene"): freq[r["gene"]]+=1
n=len(files); need=max(2,int(round(FRAC*n)))
consensus={g for g,c in freq.items() if c>=need}
rows=list(csv.reader(open(GTEX),delimiter='\t')); tissues=rows[0][1:]
tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(d,name,desc,genes,extra):
    os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    gx="GTEx" in name
    cite=f"Derived from ENCODE RNA-seq consensus expressed (>= {FRAC} of {n} biosamples)"+(f"; intersected with GTEx tissue-enrichment (t>={THR}; NIH Common Fund)" if gx else "")+"; ENCODE/NHGRI public."
    m={"standard_name":name,"library":"ENCODE_RNAseq_x_GTEx","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"consensus_fraction":FRAC,"n_biosamples":n,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"rnaseq_consensus_x_gtex","inputs":["ENCODE RNA-seq expressed-gene sets (NHGRI; public)"]+(["GTEx.tstat.hgnc.tsv (NIH Common Fund)"] if gx else []),"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE) + NIH Common Fund (GTEx)"},open(d+"/geneset.provenance.json","w"),indent=1)
emit(os.path.join(OUT,"consensus","ENCODE_RNAseq_consensus_expressed"),"ENCODE_RNAseq_consensus_expressed",
     f"Genes expressed in >= {FRAC} of ENCODE RNA-seq biosamples (n={n})",consensus,{})
nx=0
for ti,t in enumerate(tissues):
    e=sorted(g for g in consensus if g in tmap and tmap[g][ti]>=THR)
    if not e: continue
    sn=f"ENCODE_RNAseq_consensus_x_GTEx_enriched_{safe(t)}"
    emit(os.path.join(OUT,"x_GTEx",sn),sn,f"Consensus-expressed genes GTEx-enriched (t>={THR}) in {t}",e,{"tissue":t}); nx+=1
print(f"RNA-seq biosamples: {n} | consensus expressed: {len(consensus)} | x_GTEx tissue sets: {nx}")
