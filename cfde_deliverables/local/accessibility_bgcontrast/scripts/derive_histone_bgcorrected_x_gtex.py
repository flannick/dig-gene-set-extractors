#!/usr/bin/env python3
# Rebuild histone×GTEx from BACKGROUND-CORRECTED sets. Per mark: union of genes that are specifically-marked
# vs background (the "_accessible_Up" sets) across biosamples, intersected with GTEx tissue-enriched (t>=4).
# Replaces the old absolute-consensus histone×GTEx with a de-biased version.
import os, json, glob, csv, collections
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
UPDIRS=[HOME+"/accessibility_bgcontrast/output/ENCODE_histone_bgcontrast",
        HOME+"/accessibility_bgcontrast/output/ENCODE_H3K9_bgcontrast"]
GTEX="/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv"; THR=4.0
OUT=os.path.join(HOME,"accessibility_bgcontrast","output","histone_bgcorrected_x_GTEx"); os.makedirs(OUT,exist_ok=True)
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")}
# gather Up sets grouped by mark
bymark=collections.defaultdict(set)
for base in UPDIRS:
    for d in glob.glob(base+"/*_accessible_Up"):
        mp=os.path.join(d,"geneset.meta.json")
        try: mark=json.load(open(mp)).get("group","?")
        except Exception: mark="?"
        bymark[mark]|=load(os.path.join(d,"geneset.tsv"))
rows=list(csv.reader(open(GTEX),delimiter='\t')); tissues=rows[0][1:]
tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(name,desc,genes,extra):
    if not genes: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=("Background-corrected histone×GTEx: genes specifically marked vs cross-biosample background "
          "(bgcontrast Up, union across biosamples) INTERSECTED with GTEx tissue-enrichment (t>=4). "
          "ENCODE Histone ChIP-seq (NHGRI) + GTEx (Common Fund); public.")
    m={"standard_name":name,"library":"ENCODE_histone_bgcorrected_x_GTEx","description":desc,"version":"0.2",
       "file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"histone_bgcorrected_x_gtex","inputs":["ENCODE histone bgcontrast Up sets (NHGRI; public)","GTEx tissue-enrichment (Common Fund)"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE) + NIH Common Fund (GTEx)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1
nx=0
for mark,genes in bymark.items():
    mk=mark.replace("ENCODE_","").split("_")[0] if mark!="?" else "histone"
    for ti,t in enumerate(tissues):
        e=sorted(g for g in genes if g in tmap and tmap[g][ti]>=THR)
        if not e: continue
        nx+=emit(f"ENCODE_{safe(mk)}_bgcorrected_x_GTEx_enriched_{safe(t)}",f"{mk} specifically-marked (vs background) genes GTEx-enriched (t>={THR}) in {t}",e,{"mark":mk,"tissue":t})
print("marks:",{m:len(g) for m,g in bymark.items()})
print("bgcorrected histone×GTEx sets:",nx)
