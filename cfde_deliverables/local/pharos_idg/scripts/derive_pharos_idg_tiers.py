#!/usr/bin/env python3
# H: Pharos / IDG (NIH Common Fund) Target Development Level tiers -> gene-symbol sets.
# Tclin (drugs w/ known MoA), Tchem (potent ligands), Tbio (biology known), Tdark (understudied).
import os, json, urllib.request, socket
socket.setdefaulttimeout(45)   # global socket timeout — the urlopen(timeout=) arg alone let it hang 19h at 0% CPU
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"pharos_idg","output"); os.makedirs(OUT,exist_ok=True)
API="https://pharos-api.ncats.io/graphql"
import time
def gql(q):
    for attempt in range(4):
        try:
            req=urllib.request.Request(API,data=json.dumps({"query":q}).encode(),headers={"Content-Type":"application/json","User-Agent":"x"})
            return json.load(urllib.request.urlopen(req,timeout=120))
        except Exception:
            if attempt==3: raise
            time.sleep(2)
def fetch(tdl):
    syms=set(); skip=0
    while True:
        q='{targets(filter:{facets:[{facet:"Target Development Level",values:["%s"]}]},top:1000,skip:%d){targets{sym}}}'%(tdl,skip)
        d=gql(q); t=d.get("data",{}).get("targets",{}).get("targets",[])
        if not t: break
        for x in t:
            if x.get("sym"): syms.add(x["sym"])
        skip+=len(t)   # page by ACTUAL returned size (API caps page length)
        if skip % 500 == 0: print(f"  {tdl}: {skip} fetched ({len(syms)} syms)",flush=True)
    return syms
def emit(name,desc,genes,tdl):
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"Derived from Pharos/IDG Target Development Level = {tdl} (NIH Common Fund IDG; public)."
    json.dump({"standard_name":name,"library":"Pharos_IDG_TDL","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","tdl":tdl,"derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"pharos_idg_tdl","inputs":["Pharos/IDG GraphQL API (NIH Common Fund; public)"],"source_citation":cite,"public":True,"funding":"NIH Common Fund (IDG)"},open(d+"/geneset.provenance.json","w"),indent=1)
total=0
for tdl in ("Tclin","Tchem","Tbio","Tdark"):
    g=fetch(tdl); total+=len(g)
    emit(f"Pharos_IDG_{tdl}",f"IDG {tdl} targets (Pharos Target Development Level)",g,tdl)
    print(f"{tdl}: {len(g)} genes")
print("total tier genes:",total)
