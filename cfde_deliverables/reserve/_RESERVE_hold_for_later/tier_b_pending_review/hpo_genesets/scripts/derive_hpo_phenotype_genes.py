#!/usr/bin/env python3
# HPO genes_to_phenotype -> per-phenotype gene sets (genes linked to each HPO term). NHGRI-supported; public.
import os, json, urllib.request, collections
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
TMP=os.environ.get("TMPDIR","/Users/gage/.claude/jobs/32851e29/tmp")
OUT=os.path.join(HOME,"hpo_genesets","output"); os.makedirs(OUT,exist_ok=True)
MIN=int(os.environ.get("MIN_GENES","5"))
F=os.path.join(TMP,"genes_to_phenotype.txt")
if not os.path.exists(F):
    urllib.request.urlretrieve("https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt",F)
pheno=collections.defaultdict(set)
with open(F) as fh:
    hdr=fh.readline().rstrip("\n").split("\t"); ix={c:i for i,c in enumerate(hdr)}
    gi=ix.get("gene_symbol", ix.get("entrez-gene-symbol",1))
    ni=ix.get("hpo_name", ix.get("HPO-Term-Name"))
    for line in fh:
        f=line.rstrip("\n").split("\t")
        if len(f)<=max(gi,ni or 0): continue
        g=f[gi]; ph=f[ni] if ni is not None else None
        if g and ph: pheno[ph].add(g)
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(name,desc,genes,ph):
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="Derived from HPO genes_to_phenotype (Human Phenotype Ontology; NIH/NHGRI-supported; public)."
    json.dump({"standard_name":name,"library":"HPO_phenotype_genes","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","phenotype":ph,"derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"hpo_phenotype_genes","inputs":["HPO genes_to_phenotype (NHGRI-supported; public)"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (HPO/Monarch)"},open(d+"/geneset.provenance.json","w"),indent=1)
n=0
for ph,genes in pheno.items():
    if len(genes)<MIN: continue
    emit(f"HPO_{safe(ph)}",f"Genes associated with HPO phenotype: {ph}",genes,ph); n+=1
print(f"HPO phenotypes (>= {MIN} genes): {n} | total phenotypes: {len(pheno)}")
try: os.remove(F)
except OSError: pass
