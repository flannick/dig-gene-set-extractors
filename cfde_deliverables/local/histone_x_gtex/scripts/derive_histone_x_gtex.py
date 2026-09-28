#!/usr/bin/env python3
# Histone marks x GTEx: per-mark CONSENSUS gene set (genes marked in >= FRAC of that mark's ENCODE
# biosamples) intersected with GTEx tissue-enriched genes per tissue -> tissue-resolved histone sets.
import csv, json, os, glob, collections
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
GTEX="/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv"
OUT=os.path.join(HOME,"histone_x_gtex","output"); os.makedirs(OUT,exist_ok=True)
FRAC=float(os.environ.get("FRAC","0.25")); THR=float(os.environ.get("TSTAT_THR","4"))
MARKS={
 "H3K4me3": HOME+"/h3k4me3_genesets/output",
 "H3K27me3": HOME+"/histone_genesets/H3K27me3",
 "H3K36me3": HOME+"/histone_genesets/H3K36me3",
 "H3K27ac": HOME+"/histone_genesets/H3K27ac",
 "H3K4me1": HOME+"/histone_genesets/H3K4me1",
 "H3K9me3": HOME+"/histone_genesets/H3K9me3",
 "H3K9ac": HOME+"/histone_genesets/H3K9ac",
}
rows=list(csv.reader(open(GTEX),delimiter='\t')); tissues=rows[0][1:]
tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(d,name,desc,genes,extra):
    os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"Derived from ENCODE {extra.get('mark','')} histone ChIP-seq consensus (>= {FRAC} of biosamples), GRCh38, ENCODE/NHGRI, public" + (f"; intersected with GTEx tissue-enrichment (t>={THR}; NIH Common Fund)" if "GTEx" in name else "") + "."
    m={"standard_name":name,"library":"ENCODE_histone_x_GTEx","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"consensus_fraction":FRAC,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"histone_consensus_x_gtex","inputs":["ENCODE histone ChIP-seq peak->gene sets (NHGRI; public)"]+(["GTEx.tstat.hgnc.tsv (NIH Common Fund)"] if "GTEx" in name else []),"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE) + NIH Common Fund (GTEx)"},open(d+"/geneset.provenance.json","w"),indent=1)

summary=[]; nx=0
for mark,base in MARKS.items():
    files=glob.glob(os.path.join(base,"*","geneset.tsv"))
    if not files: continue
    freq=collections.Counter()
    for fp in files:
        gs={r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")}
        for g in gs: freq[g]+=1
    n=len(files); need=max(2,int(round(FRAC*n)))
    consensus={g for g,c in freq.items() if c>=need}
    emit(os.path.join(OUT,"consensus",f"ENCODE_{mark}_consensus"),f"ENCODE_{mark}_consensus",
         f"Genes marked by {mark} in >= {FRAC} of ENCODE biosamples (n={n})",consensus,{"mark":mark,"n_biosamples":n})
    for ti,t in enumerate(tissues):
        e=sorted(g for g in consensus if g in tmap and tmap[g][ti]>=THR)
        if not e: continue
        sn=f"ENCODE_{mark}_consensus_x_GTEx_enriched_{safe(t)}"
        emit(os.path.join(OUT,"x_GTEx",sn),sn,f"{mark} consensus genes GTEx-enriched (t>={THR}) in {t}",e,{"mark":mark,"tissue":t}); nx+=1
    summary.append((mark,n,len(consensus)))
print("mark | #biosamples | #consensus genes (>= %.0f%%):"%(FRAC*100))
for m,n,c in summary: print(f"  {m:9} {n:4} {c}")
print("histone x_GTEx tissue sets:",nx)
