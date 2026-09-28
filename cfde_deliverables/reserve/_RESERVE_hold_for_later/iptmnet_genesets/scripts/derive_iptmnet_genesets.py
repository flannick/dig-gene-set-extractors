#!/usr/bin/env python3
# #2: iPTMnet PTM gene sets (human). PTM-type substrate sets + per-enzyme substrate sets.
# Source: iPTMnet ptm.txt (PIR/U.Delaware; NIH/NIGMS). LICENSE: CC BY-NC-SA 4.0 (non-commercial,
# attribution, share-alike; cite iPTMnet + the per-row source databases). ptm.txt cols:
# 0 ptm_type 1 source 2 substrate_AC 3 substrate_genename 4 organism 5 site 6 enzyme_AC 7 enzyme_genename 8 note 9 pmid
import os, collections, json, urllib.request
OUT=os.path.expanduser("~/Claude/proj-valiation-challenge/iptmnet_genesets/output"); os.makedirs(OUT,exist_ok=True)
PTM=os.environ.get("IPTM","/Users/gage/.claude/jobs/32851e29/tmp/iptm_ptm.txt")
if not os.path.exists(PTM):
    urllib.request.urlretrieve("https://research.bioinformatics.udel.edu/iptmnet_data/files/current/ptm.txt", PTM)
LIC="CC BY-NC-SA 4.0 (non-commercial; attribution; share-alike; cite iPTMnet + per-row source databases)"
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:50]

by_type=collections.defaultdict(set); by_enz=collections.defaultdict(set); enz_roster=set()
for line in open(PTM, encoding="utf-8", errors="replace"):
    c=line.rstrip("\n").split("\t")
    if len(c)<8: continue
    ptm,sub_g,org,enz_g=c[0],c[3],c[4],c[7]
    if "Human" not in org: continue
    if sub_g and sub_g not in ("","-"): by_type[ptm].add(sub_g)
    if enz_g and enz_g not in ("","-"):
        enz_roster.add(enz_g)
        if sub_g and sub_g not in ("","-"): by_enz[(ptm,enz_g)].add(sub_g)

def emit(d,name,desc,genes,extra):
    os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    meta={"standard_name":name,"library":"iPTMnet_PTM","description":desc,"version":"6.2","file_type":"geneset",
          "n_genes":len(genes),"organism":"human","derived_in_this_work":True,"license":LIC,
          "source":"Derived from iPTMnet ptm.txt (PIR/U.Delaware; NIH/NIGMS), release 6.2, public."}
    meta.update(extra); json.dump(meta,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"derive_iptmnet_ptm_geneset","inputs":["iPTMnet ptm.txt (NIH/NIGMS; CC BY-NC-SA 4.0)"],
               "source_citation":"iPTMnet (Huang et al., NAR 2018); cite per-row source DBs.","license":LIC,
               "public":True,"funding":"NIH/NIGMS (iPTMnet, PIR)"},open(d+"/geneset.provenance.json","w"),indent=1)

# PTM-type substrate sets
for ptm,genes in by_type.items():
    emit(os.path.join(OUT,"by_ptm_type",f"iPTMnet_{safe(ptm)}_substrates"),
         f"iPTMnet_{safe(ptm)}_substrates", f"Human proteins undergoing {ptm} (iPTMnet)", genes, {"ptm_type":ptm})
print("PTM-type substrate sets:", {p:len(g) for p,g in sorted(by_type.items(),key=lambda kv:-len(kv[1]))})
# per-enzyme substrate sets (enzymes with >=10 substrates)
ne=0
for (ptm,enz),subs in by_enz.items():
    if len(subs)>=10:
        emit(os.path.join(OUT,"by_enzyme",f"iPTMnet_{safe(ptm)}_substrates_of_{safe(enz)}"),
             f"iPTMnet_{safe(ptm)}_substrates_of_{safe(enz)}", f"Human {ptm} substrates of {enz} (iPTMnet)", subs,
             {"ptm_type":ptm,"enzyme":enz}); ne+=1
print(f"per-enzyme substrate sets (>=10 subs): {ne} | enzyme roster size: {len(enz_roster)}")
