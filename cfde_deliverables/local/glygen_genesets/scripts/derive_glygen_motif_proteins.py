#!/usr/bin/env python3
# Decode GlyGen human_protein_glycosylation_motifs.csv (UniProt AC only) -> gene symbols via UniProt map.
# = proteins bearing glycosylation sequons/motifs (sequence-predicted glyco-capable). Standalone + xGTEx.
import csv, json, os, io, gzip, urllib.request, collections
OUT=os.path.expanduser("~/Claude/proj-valiation-challenge/glygen_genesets/output")
MAP=os.environ.get("AC2GENE","/Users/gage/.claude/jobs/32851e29/tmp/uniprot_ac2gene.tsv.gz")
GTEX=os.environ.get("GTEX","/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv")
BURL="https://data.glygen.org/ln2data/releases/data/current/reviewed/"; THR=float(os.environ.get("TSTAT_THR","4"))

ac2g={}
for line in gzip.open(MAP,'rt'):
    p=line.rstrip('\n').split('\t')
    if len(p)>=2 and p[0]!="Entry" and p[1]: ac2g[p[0]]=p[1]
print("AC->gene map entries:", len(ac2g))

data=urllib.request.urlopen(BURL+"human_protein_glycosylation_motifs.csv",timeout=120).read().decode("utf-8","replace")
r=csv.DictReader(io.StringIO(data))
allg=set(); bymotif=collections.defaultdict(set); n=0; mapped=0
for row in r:
    n+=1; ac=(row.get("uniprotkb_canonical_ac") or "").split("-")[0]
    g=ac2g.get(ac)
    if not g: continue
    mapped+=1; allg.add(g); bymotif[(row.get("motif") or "NA").strip()].add(g)
print(f"motif rows: {n} | mapped to gene: {mapped} | distinct proteins: {len(allg)} | motif types: {dict((k,len(v)) for k,v in bymotif.items())}")

def emit(d,name,desc,genes,extra):
    os.makedirs(d,exist_ok=True); genes=sorted(genes)
    cite=("Derived from GlyGen human_protein_glycosylation_motifs.csv (GlyGen, NIH/NIGMS Common Fund, public); "
          "UniProt AC->gene_symbol mapping (UniProt, public).")
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    m={"standard_name":name,"library":"GlyGen_glyco","description":desc,"version":"0.1","file_type":"geneset",
       "n_genes":len(genes),"organism":"human","derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"decode_glygen_motif_proteins","inputs":["GlyGen human_protein_glycosylation_motifs.csv (NIH/NIGMS; public)","UniProt AC->gene map (public)"],"source_citation":cite,"public":True,"funding":"NIH Common Fund / NIGMS (GlyGen)"},open(d+"/geneset.provenance.json","w"),indent=1)

emit(os.path.join(OUT,"standalone","GlyGen_glycosylation_motif_proteins"),
     "GlyGen_glycosylation_motif_proteins","Human proteins bearing glycosylation sequons/motifs (GlyGen; sequence-predicted glyco-capable)",allg,{})
for motif,genes in bymotif.items():
    if len(genes)>=20:
        sn="GlyGen_motif_"+"".join(c if c.isalnum() else "_" for c in motif)[:40]+"_proteins"
        emit(os.path.join(OUT,"standalone",sn),sn,f"Human proteins with glycosylation motif '{motif}' (GlyGen)",genes,{"motif":motif})

# x GTEx tissue-enriched (overall motif set)
rows=list(csv.reader(open(GTEX),delimiter='\t')); tissues=rows[0][1:]
tmap={r0[0]:[float(x) for x in r0[1:]] for r0 in rows[1:]}
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
nx=0
for ti,t in enumerate(tissues):
    enr=sorted(g for g in allg if g in tmap and tmap[g][ti]>=THR)
    if not enr: continue
    sn=f"GlyGen_glycosylation_motif_proteins_x_GTEx_enriched_{safe(t)}"
    emit(os.path.join(OUT,"x_GTEx",sn),sn,f"Glyco-motif proteins GTEx-enriched (t>={THR}) in {t}",enr,{"tissue":t}); nx+=1
print("motif-protein x_GTEx sets:", nx)
