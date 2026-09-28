#!/usr/bin/env python3
# ClinGen dosage sensitivity -> haploinsufficient / triplosensitive gene sets. NHGRI; public.
import os, json, urllib.request, csv
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
TMP=os.environ.get("TMPDIR","/Users/gage/.claude/jobs/32851e29/tmp")
OUT=os.path.join(HOME,"clingen_dosage","output"); os.makedirs(OUT,exist_ok=True)
F=os.path.join(TMP,"ClinGen_gene_curation_list_GRCh38.tsv")
if not os.path.exists(F):
    urllib.request.urlretrieve("https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv",F)
# file has comment lines starting with '#'; find the header row (contains 'Gene Symbol')
lines=[l.rstrip("\n") for l in open(F)]
hi=set(); hi_some=set(); ts=set(); ts_some=set()
hdr=None
for i,l in enumerate(lines):
    if l.lstrip("#").startswith("Gene Symbol") and "\t" in l:   # real header (not the 'create link' comment)
        hdr=l.lstrip("#").split("\t"); start=i+1; break
ix={c.strip():j for j,c in enumerate(hdr)} if hdr else {}
gi=next((ix[c] for c in ix if c.lower()=="gene symbol"),0)
hii=next((ix[c] for c in ix if "haploinsufficiency score" in c.lower()),None)
tsi=next((ix[c] for c in ix if "triplosensitivity score" in c.lower()),None)
for l in lines[start:]:
    if not l.strip() or l.startswith("#"): continue
    f=l.split("\t")
    if len(f)<=gi: continue
    g=f[gi].strip()
    if not g: continue
    hv=f[hii].strip() if (hii is not None and len(f)>hii) else ""
    tv=f[tsi].strip() if (tsi is not None and len(f)>tsi) else ""
    if hv=="3": hi.add(g)
    if hv in ("1","2","3"): hi_some.add(g)
    if tv=="3": ts.add(g)
    if tv in ("1","2","3"): ts_some.add(g)
def emit(name,desc,genes,extra):
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="Derived from ClinGen dosage sensitivity curation (GRCh38; ClinGen/NHGRI; public)."
    m={"standard_name":name,"library":"ClinGen_dosage","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"clingen_dosage","inputs":["ClinGen dosage sensitivity (NHGRI; public)"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ClinGen)"},open(d+"/geneset.provenance.json","w"),indent=1)
emit("ClinGen_haploinsufficient_score3","ClinGen haploinsufficient genes (HI score 3, sufficient evidence)",hi,{"score":"HI=3"})
emit("ClinGen_haploinsufficient_some_evidence","ClinGen haploinsufficient genes (HI score 1-3)",hi_some,{"score":"HI in 1-3"})
emit("ClinGen_triplosensitive_score3","ClinGen triplosensitive genes (TS score 3)",ts,{"score":"TS=3"})
emit("ClinGen_triplosensitive_some_evidence","ClinGen triplosensitive genes (TS score 1-3)",ts_some,{"score":"TS in 1-3"})
print(f"ClinGen: HI=3 {len(hi)} | HI 1-3 {len(hi_some)} | TS=3 {len(ts)} | TS 1-3 {len(ts_some)}")
try: os.remove(F)
except OSError: pass
