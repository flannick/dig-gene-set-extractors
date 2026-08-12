#!/usr/bin/env python3
# GlyGen disease-glyco + cancer-glyco-mutation gene sets. Decodes UniProt AC->gene via UniProt map.
import csv, json, os, io, gzip, urllib.request, collections
OUT=os.path.expanduser("~/Claude/proj-valiation-challenge/glygen_genesets/output")
MAP=os.environ.get("AC2GENE","/Users/gage/.claude/jobs/32851e29/tmp/uniprot_ac2gene.tsv.gz")
TMP="/Users/gage/.claude/jobs/32851e29/tmp"
B="https://data.glygen.org/ln2data/releases/data/current/reviewed/"
ac2g={}
for line in gzip.open(MAP,'rt'):
    p=line.rstrip('\n').split('\t')
    if len(p)>=2 and p[0]!="Entry" and p[1]: ac2g[p[0]]=p[1]
def G(ac): return ac2g.get((ac or "").split("-")[0])
def text(fname):
    return urllib.request.urlopen(B+fname,timeout=300).read().decode("utf-8","replace")
def emit(sub,name,desc,genes,extra):
    d=os.path.join(OUT,sub,name); os.makedirs(d,exist_ok=True); genes=sorted(g for g in genes if g)
    cite=f"Derived from GlyGen reviewed data (GlyGen, NIH/NIGMS Common Fund, public); UniProt AC->gene where needed."
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    m={"standard_name":name,"library":"GlyGen_disease_cancer_glyco","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"derive_glygen_disease_cancer","inputs":["GlyGen reviewed disease/mutation files (NIH/NIGMS; public)","UniProt AC->gene map (public)"],"source_citation":cite,"public":True,"funding":"NIH Common Fund / NIGMS (GlyGen)"},open(d+"/geneset.provenance.json","w"),indent=1)

# ---- disease-glyco ----
dis_all=set(); by_dis=collections.defaultdict(set)
r=csv.DictReader(io.StringIO(text("human_protein_disease_glyco.csv")))
for row in r:
    g=row.get("gene_symbol") or G(row.get("uniprotkb_canonical_ac"))
    if not g: continue
    dis_all.add(g);
    dn=(row.get("mondo_disease_name") or row.get("mondo_label") or "").strip()
    if dn: by_dis[dn].add(g)
r=csv.DictReader(io.StringIO(text("human_protein_disease_glycosmos.csv")))
for row in r:
    g=G(row.get("uniprotkb_canonical_ac"))
    if g:
        dis_all.add(g); dn=(row.get("ggd_name") or "").strip()
        if dn: by_dis[dn].add(g)
emit("disease_glyco","GlyGen_disease_associated_glycoproteins","Human glycoproteins with disease-associated glycosylation (GlyGen)",dis_all,{})
ndis=0
for dn,genes in by_dis.items():
    if len(genes)>=3:
        sn="GlyGen_disease_glyco_"+"".join(c if c.isalnum() else "_" for c in dn)[:45]
        emit("disease_glyco/by_disease",sn,f"Glycoproteins in disease '{dn}' (GlyGen)",genes,{"disease":dn}); ndis+=1
print(f"disease-glyco: {len(dis_all)} genes overall | {ndis} per-disease sets (>=3)")

# ---- cancer-glyco-mutation ----
for fname,name,desc in [
  ("human_protein_mutation_cancer_glycoeffect.csv","GlyGen_cancer_glyco_mutation_genes","Genes with cancer mutations affecting glycosylation (GlyGen)"),
  ("human_protein_mutation_cancer_glycosylation_loss.csv","GlyGen_cancer_glycosylation_loss_genes","Genes where cancer mutations cause loss of a glycosylation site (GlyGen)")]:
    fp=os.path.join(TMP,fname); urllib.request.urlretrieve(B+fname,fp)
    genes=set()
    with open(fp,encoding="utf-8",errors="replace") as fh:
        for row in csv.DictReader(fh):
            g=G(row.get("uniprotkb_canonical_ac"))
            if g: genes.add(g)
    os.remove(fp); emit("cancer_glyco",name,desc,genes,{}); print(f"{name}: {len(genes)} genes")
