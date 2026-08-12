#!/usr/bin/env python3
# #4: per-tissue 'tissue-enriched' gene sets DERIVED from local GTEx specificity t-stats.
# NOTE: GTEx t-stat = RELATIVE tissue specificity, NOT absolute expression (documented in meta).
import csv, json, os
G=os.environ.get("GTEX","/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv")
OUT=os.environ.get("OUTDIR",os.path.expanduser("~/Claude/proj-valiation-challenge/gtex_tissue_enriched/output"))
os.makedirs(OUT,exist_ok=True); THR=float(os.environ.get("TSTAT_THR","4"))
rows=list(csv.reader(open(G),delimiter='\t')); tissues=rows[0][1:]
data=[(r[0],[float(x) for x in r[1:]]) for r in rows[1:]]
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
for ti,t in enumerate(tissues):
    enr=sorted([(g,v[ti]) for g,v in data if v[ti]>=THR],key=lambda kv:-kv[1])
    genes=[g for g,_ in enr]; name=f"GTEx_tissue_enriched_{safe(t)}"; d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True)
    cite=f"Derived from GTEx tissue-specificity t-stats (GTEx.tstat.hgnc.tsv) for {t}; GTEx (NIH Common Fund), public aggregate."
    open(d+"/geneset.tsv","w").write("gene\tgtex_tstat\n"+"\n".join(f"{g}\t{s:.2f}" for g,s in enr)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\tGenes tissue-enriched (GTEx t-stat>={THR}) in {t}\t"+"\t".join(genes)+"\n")
    json.dump({"standard_name":name,"library":"GTEx_tissue_enriched","description":f"Genes with tissue-SPECIFIC enrichment (GTEx t-stat>={THR}; relative specificity, NOT absolute expression) in {t} (DERIVED from GTEx public t-stats).","version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","tissue":t,"tstat_threshold":THR,"derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"derive_tissue_enriched_genes","method":f"GTEx t-stat>={THR} (relative tissue specificity)","inputs":["GTEx.tstat.hgnc.tsv (NIH Common Fund; public aggregate)"],"source_citation":cite,"public":True,"funding":"NIH Common Fund (GTEx)"},open(d+"/geneset.provenance.json","w"),indent=1)
counts=sorted(((t,sum(1 for g,v in data if v[ti]>=THR)) for ti,t in enumerate(tissues)),key=lambda x:-x[1])
print(f"GTEx tissue-enriched: {len(tissues)} tissues, t-stat>={THR}")
for t,n in counts[:5]+counts[-2:]: print(f"  {t}: {n}")
