#!/usr/bin/env python3
# I: GTEx tissue-DEPLETED genes (t-stat <= -THR) per tissue — inverse of tissue-enriched (specifically LOW).
import csv, json, os
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
GTEX="/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv"
OUT=os.path.join(HOME,"gtex_tissue_depleted","output"); os.makedirs(OUT,exist_ok=True); THR=4.0
rows=list(csv.reader(open(GTEX),delimiter='\t')); tissues=rows[0][1:]
tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(d,name,desc,genes,t):
    os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"Derived from GTEx tissue t-statistics (t<=-{THR}; relative under-expression), NIH Common Fund, public."
    json.dump({"standard_name":name,"library":"GTEx_tissue_depleted","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","tissue":t,"derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"gtex_tissue_depleted","inputs":["GTEx.tstat.hgnc.tsv (NIH Common Fund; public)"],"source_citation":cite,"public":True,"funding":"NIH Common Fund (GTEx)"},open(d+"/geneset.provenance.json","w"),indent=1)
n=0
for ti,t in enumerate(tissues):
    g=[gene for gene,v in tmap.items() if v[ti]<=-THR]
    if not g: continue
    sn=f"GTEx_tissue_depleted_{safe(t)}"
    emit(os.path.join(OUT,sn),sn,f"Genes relatively DEPLETED (t<=-{THR}) in {t} (GTEx)",g,t); n+=1
print("GTEx tissue-depleted sets:",n)
