#!/usr/bin/env python3
# #3: GlyGen-derived glyco gene sets, STANDALONE (#3a) and x GTEx tissue-enriched (#3b).
# Imports gene-centric GlyGen reviewed files (public, NIH/NIGMS Common Fund Glycoscience), cited.
# x GTEx = GlyGen glyco genes that are GTEx tissue-enriched (t-stat>=THR) -> granular tissue-resolved.
import csv, json, os, io, urllib.request
BASE=os.path.expanduser("~/Claude/proj-valiation-challenge/glygen_genesets"); OUT=os.path.join(BASE,"output")
os.makedirs(OUT,exist_ok=True)
GTEX=os.environ.get("GTEX","/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv")
BURL="https://data.glygen.org/ln2data/releases/data/current/reviewed/"; THR=float(os.environ.get("TSTAT_THR","4"))
FILES={"GlyGen_glycosyltransferases":"human_protein_glycosyltransferase.csv",
       "GlyGen_glycohydrolases":"human_protein_glycohydrolase.csv",
       "GlyGen_glycogenes":"human_protein_glycogenes.csv",
       "GlyGen_glycosylation_motif_proteins":"human_protein_glycosylation_motifs.csv"}
def gcol(fl):
    for c in ("gene_symbol","gene_name","gene","hgnc_symbol"):
        if c in (fl or []): return c
    return None
def fetch(fname):
    data=urllib.request.urlopen(BURL+fname,timeout=180).read().decode("utf-8","replace")
    r=csv.DictReader(io.StringIO(data)); gc=gcol(r.fieldnames)
    if not gc: return None,(r.fieldnames or [])[:8]
    return sorted({row[gc].strip() for row in r if row.get(gc) and row[gc].strip() not in ("","NA")}), gc

standalone={}
for name,fname in FILES.items():
    try: genes,info=fetch(fname)
    except Exception as e: print(f"FAIL {name}: {e}"); continue
    if genes is None: print(f"SKIP {name}: no gene column (cols {info})"); continue
    standalone[name]=genes; d=os.path.join(OUT,"standalone",name); os.makedirs(d,exist_ok=True)
    cite=f"Derived from GlyGen reviewed dataset {fname} (GlyGen, NIH/NIGMS Common Fund Glycoscience, public)."
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{name} (GlyGen-derived)\t"+"\t".join(genes)+"\n")
    json.dump({"standard_name":name,"library":"GlyGen_glyco","description":f"{name} (DERIVED from GlyGen {fname}).","version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"import_glygen_glyco_geneset","inputs":[f"GlyGen {fname} (NIH/NIGMS; public)"],"source_citation":cite,"public":True,"funding":"NIH Common Fund / NIGMS (GlyGen)"},open(d+"/geneset.provenance.json","w"),indent=1)
    print(f"{name}: {len(genes)} genes (col {info})")

rows=list(csv.reader(open(GTEX),delimiter='\t')); tissues=rows[0][1:]
tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
nx=0
for name,genes in standalone.items():
    for ti,t in enumerate(tissues):
        enr=sorted(g for g in genes if g in tmap and tmap[g][ti]>=THR)
        if not enr: continue
        sn=f"{name}_x_GTEx_enriched_{safe(t)}"; d=os.path.join(OUT,"x_GTEx",sn); os.makedirs(d,exist_ok=True)
        cite=f"GlyGen {name} genes tissue-enriched (GTEx t>={THR}) in {t}. GlyGen (NIH/NIGMS) x GTEx (NIH Common Fund), public."
        open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(enr)+"\n")
        open(d+"/genesets.gmt","w").write(f"{sn}\t{name} genes GTEx-enriched (t>={THR}) in {t}\t"+"\t".join(enr)+"\n")
        json.dump({"standard_name":sn,"library":"GlyGen_glyco_x_GTEx","description":f"{name} genes with GTEx tissue-enrichment (t>={THR}) in {t} (DERIVED; GlyGen x GTEx).","version":"0.1","file_type":"geneset","n_genes":len(enr),"organism":"human","tissue":t,"derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
        json.dump({"focus":sn,"operation":"glygen_glyco_x_gtex_tissue","inputs":[f"GlyGen {name}","GTEx.tstat.hgnc.tsv (NIH Common Fund)"],"source_citation":cite,"public":True},open(d+"/geneset.provenance.json","w"),indent=1)
        nx+=1
print(f"x_GTEx tissue-enriched sets written: {nx}")
