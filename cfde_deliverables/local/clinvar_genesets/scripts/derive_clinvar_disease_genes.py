#!/usr/bin/env python3
# F: ClinVar (NLM) per-condition gene sets. Genes with PATHOGENIC / LIKELY-PATHOGENIC germline variants
# for each condition (PhenotypeList). GRCh38 rows, human. Lean: stream variant_summary.txt.gz.
import os, json, gzip, urllib.request, collections
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
TMP=os.environ.get("TMPDIR","/Users/gage/.claude/jobs/32851e29/tmp")
OUT=os.path.join(HOME,"clinvar_genesets","output"); os.makedirs(OUT,exist_ok=True)
MIN_GENES=int(os.environ.get("MIN_GENES","5"))
VS=os.path.join(TMP,"variant_summary.txt.gz")
if not os.path.exists(VS):
    urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",VS)
cond=collections.defaultdict(set); allpath=set()
with gzip.open(VS,'rt') as fh:
    hdr=fh.readline().rstrip("\n").lstrip("#").split("\t")
    ix={c:i for i,c in enumerate(hdr)}
    iG=ix.get("GeneSymbol"); iS=ix.get("ClinicalSignificance"); iP=ix.get("PhenotypeList")
    iA=ix.get("Assembly"); iT=ix.get("Type")
    for line in fh:
        f=line.rstrip("\n").split("\t")
        if len(f)<=max(iG,iS,iP,iA): continue
        if f[iA]!="GRCh38": continue
        sig=f[iS].lower()
        if "pathogenic" not in sig or "conflict" in sig: continue   # keeps Pathogenic / Likely pathogenic
        gene=f[iG]
        if not gene or gene in("-","") or ";" in gene: continue     # skip multi-gene/intergenic
        allpath.add(gene)
        for ph in f[iP].split("|"):
            ph=ph.strip()
            if ph and ph.lower() not in("not provided","not specified","see cases",""):
                cond[ph].add(gene)
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(sub,name,desc,genes,extra):
    d=os.path.join(OUT,sub,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="Derived from ClinVar variant_summary (GRCh38; Pathogenic/Likely-pathogenic; gene with reported variant for condition), NCBI/NLM/NIH, public domain."
    m={"standard_name":name,"library":"ClinVar_condition_genes","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"clinvar_condition_genes","inputs":["ClinVar variant_summary (NCBI/NLM/NIH; public domain)"],"source_citation":cite,"public":True,"funding":"NIH/NLM (ClinVar)"},open(d+"/geneset.provenance.json","w"),indent=1)
emit("all","ClinVar_pathogenic_genes_all","Genes with any Pathogenic/Likely-pathogenic ClinVar variant (GRCh38)",allpath,{})
n=0
for ph,genes in cond.items():
    if len(genes)<MIN_GENES: continue
    sn=f"ClinVar_{safe(ph)}"
    emit("by_condition",sn,f"Genes with Pathogenic/Likely-pathogenic ClinVar variants for: {ph}",genes,{"condition":ph}); n+=1
print(f"ClinVar pathogenic genes: {len(allpath)} | conditions (>= {MIN_GENES} genes): {n}")
try: os.remove(VS)
except OSError: pass
