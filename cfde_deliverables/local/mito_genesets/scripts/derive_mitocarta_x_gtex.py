#!/usr/bin/env python3
# #3a: MitoCarta3.0 (Broad/NIGMS) -> DERIVED mito gene sets. We do NOT ship MitoCarta verbatim (published
# product); we ship derivations: MitoCarta x GTEx tissue-enrichment (tissue-resolved mito biology).
# The MitoCarta x TF-regulon intersection (mito *regulation*) is built separately once TF regulons exist.
import os, json, csv
import pandas as pd
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
TMP="/Users/gage/.claude/jobs/32851e29/tmp"
XLS=os.path.join(TMP,"MitoCarta3.xls")
GTEX="/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv"; THR=4.0
OUT=os.path.join(HOME,"mito_genesets","output"); os.makedirs(OUT,exist_ok=True)
df=pd.read_excel(XLS,sheet_name="A Human MitoCarta3.0",engine="xlrd")
mito=sorted({str(s).strip() for s in df["Symbol"].dropna() if str(s).strip()})
# reference table (NOT a geneset deliverable): symbol + sublocalization + pathways
with open(os.path.join(OUT,"mitocarta3_reference.tsv"),"w") as fh:
    fh.write("symbol\tsub_mito_localization\tmito_pathways\n")
    for _,r in df.iterrows():
        fh.write(f"{r.get('Symbol','')}\t{r.get('MitoCarta3.0_SubMitoLocalization','')}\t{str(r.get('MitoCarta3.0_MitoPathways','')).replace(chr(10),' ')}\n")
print("MitoCarta3.0 mito genes:",len(mito))
rows=list(csv.reader(open(GTEX),delimiter='\t')); tissues=rows[0][1:]
tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(d,name,desc,genes,extra):
    os.makedirs(d,exist_ok=True); genes=sorted(g for g in genes if g)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"Derived from MitoCarta3.0 (Broad Institute; NIH/NIGMS) intersected with GTEx tissue-enrichment (t>={THR}; NIH Common Fund); both public. Inventory cited, not redistributed verbatim."
    m={"standard_name":name,"library":"MitoCarta3_x_GTEx","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"mitocarta_x_gtex","inputs":["MitoCarta3.0 (Broad; NIH/NIGMS; public, cited)","GTEx.tstat.hgnc.tsv (NIH Common Fund)"],"source_citation":cite,"public":True,"funding":"NIH/NIGMS (MitoCarta) + NIH Common Fund (GTEx)"},open(d+"/geneset.provenance.json","w"),indent=1)
nx=0
for ti,t in enumerate(tissues):
    e=sorted(g for g in mito if g in tmap and tmap[g][ti]>=THR)
    if not e: continue
    sn=f"MitoCarta3_x_GTEx_enriched_{safe(t)}"
    emit(os.path.join(OUT,"x_GTEx",sn),sn,f"MitoCarta3.0 mitochondrial genes GTEx-enriched (t>={THR}) in {t}",e,{"tissue":t}); nx+=1
print("MitoCarta x_GTEx tissue sets:",nx)
