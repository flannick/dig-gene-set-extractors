#!/usr/bin/env python3
# PTM enzyme + effector(substrate) gene lists per PTM type, standalone AND x GTEx tissue-enriched.
# Generalizes glyco A/B (enzymes) + C/D (effectors) to ALL PTM types (iPTMnet). Human only.
# iPTMnet resolves PTM TYPE (phospho/methyl/ubiq/acetyl/...), NOT ubiquitin linkage subtypes.
# Source iPTMnet ptm.txt (NIH/NIGMS). LICENSE CC BY-NC-SA 4.0 (non-commercial; attribute; cite source DBs).
import os, csv, json, collections
PTM=os.environ.get("IPTM","/Users/gage/.claude/jobs/32851e29/tmp/iptm_ptm.txt")
GTEX=os.environ.get("GTEX","/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv")
OUT=os.path.expanduser("~/Claude/proj-valiation-challenge/iptmnet_genesets/output/enzyme_effector_x_gtex")
THR=float(os.environ.get("TSTAT_THR","4")); LIC="CC BY-NC-SA 4.0 (non-commercial; attribution; cite iPTMnet + source DBs)"
os.makedirs(OUT,exist_ok=True)
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:50]

enz=collections.defaultdict(set); sub=collections.defaultdict(set)   # ptm_type -> genes
for line in open(PTM,encoding="utf-8",errors="replace"):
    c=line.rstrip("\n").split("\t")
    if len(c)<8 or "Human" not in c[4]: continue
    t=c[0]
    if c[3] and c[3] not in ("","-"): sub[t].add(c[3])
    if c[7] and c[7] not in ("","-"): enz[t].add(c[7])

rows=list(csv.reader(open(GTEX),delimiter='\t')); tissues=rows[0][1:]
tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}

def emit(d,name,desc,genes,extra):
    os.makedirs(d,exist_ok=True); genes=sorted(genes)
    cite="Derived from iPTMnet ptm.txt (PIR/U.Delaware; NIH/NIGMS), release 6.2, public."
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    m={"standard_name":name,"library":"iPTMnet_PTM_enzyme_effector","description":desc,"version":"6.2",
       "file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,
       "license":LIC,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"iptmnet_ptm_enzyme_effector_x_gtex","inputs":["iPTMnet ptm.txt (NIH/NIGMS; CC BY-NC-SA 4.0)"]+(["GTEx.tstat.hgnc.tsv (NIH Common Fund)"] if "GTEx" in name else []),"source_citation":cite,"license":LIC,"public":True,"funding":"NIH/NIGMS (iPTMnet)"},open(d+"/geneset.provenance.json","w"),indent=1)

nstand=nx=0
for role,dd in (("enzymes",enz),("substrates",sub)):
    for t,genes in dd.items():
        T=safe(t)
        emit(os.path.join(OUT,"standalone",f"iPTMnet_{T}_{role}"),f"iPTMnet_{T}_{role}",
             f"Human {role} of {t} (iPTMnet)",genes,{"ptm_type":t,"role":role}); nstand+=1
        for ti,tis in enumerate(tissues):
            e=sorted(g for g in genes if g in tmap and tmap[g][ti]>=THR)
            if not e: continue
            sn=f"iPTMnet_{T}_{role}_x_GTEx_enriched_{safe(tis)}"
            emit(os.path.join(OUT,"x_GTEx",sn),sn,f"{t} {role} GTEx-enriched (t>={THR}) in {tis}",e,{"ptm_type":t,"role":role,"tissue":tis}); nx+=1
print("PTM types:",sorted(enz)|set() if False else sorted(set(enz)|set(sub)))
print("enzyme set sizes:",{t:len(g) for t,g in sorted(enz.items(),key=lambda kv:-len(kv[1]))})
print(f"standalone enzyme+substrate sets: {nstand} | x_GTEx tissue sets: {nx}")
