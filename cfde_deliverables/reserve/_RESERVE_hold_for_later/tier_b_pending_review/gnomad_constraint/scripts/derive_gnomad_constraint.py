#!/usr/bin/env python3
# gnomAD v2.1.1 LoF constraint -> constrained-gene sets (LoF intolerance). Broad/NIH; public.
import os, json, gzip, urllib.request
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
TMP=os.environ.get("TMPDIR","/Users/gage/.claude/jobs/32851e29/tmp")
OUT=os.path.join(HOME,"gnomad_constraint","output"); os.makedirs(OUT,exist_ok=True)
F=os.path.join(TMP,"gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
if not os.path.exists(F):
    urllib.request.urlretrieve("https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",F)
rows=[]
with gzip.open(F,'rt') as fh:
    hdr=fh.readline().rstrip("\n").split("\t"); ix={c:i for i,c in enumerate(hdr)}
    for line in fh:
        f=line.rstrip("\n").split("\t"); rows.append(f)
def num(f,c):
    try: return float(f[ix[c]])
    except: return None
def emit(name,desc,genes,extra):
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(set(g for g in genes if g))
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="Derived from gnomAD v2.1.1 LoF constraint metrics by gene (Broad/NIH; public)."
    m={"standard_name":name,"library":"gnomAD_constraint","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"gnomad_constraint","inputs":["gnomAD v2.1.1 constraint (Broad/NIH; public)"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (gnomAD)"},open(d+"/geneset.provenance.json","w"),indent=1)
def col(f,c): return f[ix[c]] if c in ix and ix[c]<len(f) else ""
hi_loeuf=[col(f,"gene") for f in rows if (num(f,"oe_lof_upper") is not None and num(f,"oe_lof_upper")<0.35)]
loeuf06=[col(f,"gene") for f in rows if (num(f,"oe_lof_upper") is not None and num(f,"oe_lof_upper")<0.6)]
pli=[col(f,"gene") for f in rows if (num(f,"pLI") is not None and num(f,"pLI")>=0.9)]
mis_z=[col(f,"gene") for f in rows if (num(f,"mis_z") is not None and num(f,"mis_z")>=3.09)]
emit("gnomAD_LoF_constrained_LOEUF_lt0.35","Highly LoF-constrained genes (gnomAD LOEUF<0.35)",hi_loeuf,{"metric":"LOEUF<0.35"})
emit("gnomAD_LoF_constrained_LOEUF_lt0.6","LoF-constrained genes (gnomAD LOEUF<0.6)",loeuf06,{"metric":"LOEUF<0.6"})
emit("gnomAD_LoF_intolerant_pLI_ge0.9","LoF-intolerant genes (gnomAD pLI>=0.9)",pli,{"metric":"pLI>=0.9"})
emit("gnomAD_missense_constrained_misZ_ge3.09","Missense-constrained genes (gnomAD mis_z>=3.09)",mis_z,{"metric":"mis_z>=3.09"})
print(f"gnomAD: LOEUF<0.35={len(set(hi_loeuf))} LOEUF<0.6={len(set(loeuf06))} pLI>=0.9={len(set(pli))} misZ>=3.09={len(set(mis_z))}")
try: os.remove(F)
except OSError: pass
