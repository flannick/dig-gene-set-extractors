#!/usr/bin/env python3
# Combinatorial chromatin-STATE gene sets per biosample, from the background-corrected histone marks.
# Uses specifically-marked ("Up") genes per mark, matched within a biosample:
#   bivalent      = H3K4me3 & H3K27me3            (poised / developmental)
#   active        = H3K4me3 & H3K27ac             (active promoter+enhancer)
#   pure_polycomb = H3K27me3 - H3K4me3            (stably repressed)
#   transcribed   = H3K36me3                       (gene-body / actively transcribed)
#   heterochromatin = H3K9me3                      (constitutive silencing)
import os, re, json, glob, collections, csv
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"interop_intersections","output","chromatin_states"); os.makedirs(OUT,exist_ok=True)
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")}
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
# gather corrected histone "Up" sets keyed by (biosample -> mark -> genes)
bio=collections.defaultdict(dict)
for base in ["accessibility_bgcontrast/output/ENCODE_histone_bgcontrast","accessibility_bgcontrast/output/ENCODE_H3K9_bgcontrast"]:
    for d in glob.glob(os.path.join(HOME,base,"*_accessible_Up")):
        mp=os.path.join(d,"geneset.meta.json")
        try: m=json.load(open(mp))
        except Exception: continue
        grp=m.get("group",""); bs=m.get("biosample")
        mk=re.sub(r'^ENCODE_','',grp).split("_")[0] if grp else None
        if bs and mk: bio[bs][mk]=load(os.path.join(d,"geneset.tsv"))
print("biosamples with corrected histone marks:",len(bio))
def emit(name,desc,genes,bs,cls,marks):
    if not genes: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"Combinatorial chromatin state from background-corrected ENCODE histone marks ({'+'.join(marks)}) in {bs}. ENCODE/NHGRI public."
    json.dump({"standard_name":name,"library":"ENCODE_chromatin_state","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","biosample":bs,"state":cls,"marks":marks,"derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"chromatin_state","inputs":[f"ENCODE corrected histone {m} Up ({bs})" for m in marks],"public":True,"source_citation":cite,"funding":"NIH/NHGRI (ENCODE)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1
nb=nc=collections.Counter()
tot=0; states=collections.Counter()
for bs,marks in bio.items():
    def g(m): return marks.get(m,set())
    combos=[]
    if g("H3K4me3") and g("H3K27me3"): combos.append(("bivalent",g("H3K4me3")&g("H3K27me3"),["H3K4me3","H3K27me3"]))
    if g("H3K4me3") and g("H3K27ac"):  combos.append(("active",g("H3K4me3")&g("H3K27ac"),["H3K4me3","H3K27ac"]))
    if g("H3K27me3"): combos.append(("pure_polycomb",g("H3K27me3")-g("H3K4me3"),["H3K27me3"]))
    if g("H3K36me3"): combos.append(("transcribed",g("H3K36me3"),["H3K36me3"]))
    if g("H3K9me3"):  combos.append(("heterochromatin",g("H3K9me3"),["H3K9me3"]))
    for cls,genes,mk in combos:
        n=emit(f"ENCODE_chromstate_{cls}_{safe(bs)}",f"{bs}: {cls} chromatin state ({'+'.join(mk)})",genes,bs,cls,mk)
        tot+=n; 
        if n: states[cls]+=1
print("chromatin-state sets:",tot,"| by state:",dict(states))
