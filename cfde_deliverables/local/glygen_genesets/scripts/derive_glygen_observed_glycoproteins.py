#!/usr/bin/env python3
# #1: GlyGen OBSERVED glycoproteins (experimental site data), split by glycosylation_type (N/O-linked)
# and structure_glycan_type (high-mannose=IMMATURE vs complex/hybrid=MATURE). Decodes UniProt AC->gene.
# Aggregates the human_proteoform_glycosylation_sites_*.csv files (skips any that return non-CSV).
import csv, json, os, io, gzip, urllib.request, collections
OUT=os.path.expanduser("~/Claude/proj-valiation-challenge/glygen_genesets/output")
MAP="/Users/gage/.claude/jobs/32851e29/tmp/uniprot_ac2gene.tsv.gz"
GTEX="/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv"
B="https://data.glygen.org/ln2data/releases/data/current/reviewed/"; THR=4.0
FILES=["glyconnect","gptwiki","literature_mining_manually_verified","literature","embl","harvard","carbbank","glycomeatlas","c_man"]
ac2g={}
for line in gzip.open(MAP,'rt'):
    p=line.rstrip('\n').split('\t')
    if len(p)>=2 and p[0]!="Entry" and p[1]: ac2g[p[0]]=p[1]
def G(ac): return ac2g.get((ac or "").split("-")[0])

allg=set(); bylink=collections.defaultdict(set); bymat=collections.defaultdict(set)
def maturity(gt):
    s=(gt or "").lower()
    if "high" in s and "mannose" in s or "oligomannose" in s or "paucimannose" in s: return "immature_high_mannose"
    if any(k in s for k in("complex","hybrid","antenna","sial","fucos")): return "mature_complex_hybrid"
    return None
for short in FILES:
    fn=f"human_proteoform_glycosylation_sites_{short}.csv"
    try: data=urllib.request.urlopen(B+fn,timeout=120).read().decode("utf-8","replace")
    except Exception: continue
    if data[:1]=="<": continue   # SPA/HTML = file not present
    r=csv.DictReader(io.StringIO(data))
    for row in r:
        g=row.get("gene_symbol") or G(row.get("uniprotkb_canonical_ac"))
        if not g: continue
        allg.add(g)
        lt=(row.get("glycosylation_type") or "").strip().lower()
        if "n-link" in lt or lt=="n-linked": bylink["N_linked"].add(g)
        elif "o-link" in lt or lt=="o-linked": bylink["O_linked"].add(g)
        m=maturity(row.get("structure_glycan_type"))
        if m: bymat[m].add(g)
print("observed glycoproteins:",len(allg),"| by link:",{k:len(v) for k,v in bylink.items()},"| by maturity:",{k:len(v) for k,v in bymat.items()})

def emit(sub,name,desc,genes,extra):
    d=os.path.join(OUT,sub,name); os.makedirs(d,exist_ok=True); genes=sorted(g for g in genes if g)
    cite="Derived from GlyGen human_proteoform_glycosylation_sites_* (experimental; GlyGen NIH/NIGMS, public); UniProt AC->gene."
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    m={"standard_name":name,"library":"GlyGen_observed_glycoproteins","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"glygen_observed_glycoproteins","inputs":["GlyGen proteoform glycosylation_sites (NIH/NIGMS; public)","UniProt AC->gene map"]+(["GTEx.tstat.hgnc.tsv (NIH Common Fund)"] if "GTEx" in name else []),"source_citation":cite,"public":True,"funding":"NIH Common Fund / NIGMS (GlyGen)"},open(d+"/geneset.provenance.json","w"),indent=1)

emit("observed","GlyGen_observed_glycoproteins","Human proteins with experimentally-observed glycosylation (GlyGen)",allg,{})
for k,genes in bylink.items(): emit("observed",f"GlyGen_observed_{k}_glycoproteins",f"Human proteins with observed {k.replace('_','-')} glycosylation (GlyGen)",genes,{"link":k})
for k,genes in bymat.items(): emit("observed",f"GlyGen_observed_{k}_glycoproteins",f"Human glycoproteins observed carrying {k.replace('_',' ')} glycans (GlyGen; {'IMMATURE' if 'immature' in k else 'MATURE'})",genes,{"maturity":k})
# x GTEx for overall + N/O
rows=list(csv.reader(open(GTEX),delimiter='\t')); tissues=rows[0][1:]; tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
nx=0
for setname,genes in [("observed_glycoproteins",allg)]+[("observed_"+k,v) for k,v in bylink.items()]:
    for ti,t in enumerate(tissues):
        e=sorted(g for g in genes if g in tmap and tmap[g][ti]>=THR)
        if not e: continue
        sn=f"GlyGen_{setname}_x_GTEx_enriched_{safe(t)}"
        emit("observed_x_GTEx",sn,f"GlyGen {setname} GTEx-enriched (t>={THR}) in {t}",e,{"tissue":t}); nx+=1
print("observed x_GTEx sets:",nx)
