#!/usr/bin/env python3
# Master-regulator-per-tissue: for each GTEx tissue, rank TFs by how many of the tissue's ENRICHED genes
# their CORRECTED (specifically-bound) targets cover; emit the tissue-enriched genes bound by the top-K TFs
# ("regulated tissue-specific core"), recording the top master-regulator TFs + counts in meta.
import os, json, glob, csv, collections
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"interop_intersections","output","master_regulator_per_tissue"); os.makedirs(OUT,exist_ok=True)
TOPK=int(os.environ.get("TOPK","10"))
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")} if os.path.exists(fp) else set()
# corrected TF specific targets
tf={}
for d in glob.glob(HOME+"/encode_regulons_bgcontrast/output/ENCODE_TF_regulon_bgcorrected/*_specific_Up"):
    try: t=json.load(open(os.path.join(d,"geneset.meta.json"))).get("target",os.path.basename(d))
    except Exception: t=os.path.basename(d)
    g=load(os.path.join(d,"geneset.tsv"))
    if g: tf[t]=g
print("corrected TF factors:",len(tf))
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
n=0
for d in glob.glob(HOME+"/gtex_tissue_enriched/output/GTEx_tissue_enriched_*"):
    tis=os.path.basename(d).replace("GTEx_tissue_enriched_","")
    enr=load(os.path.join(d,"geneset.tsv"))
    if not enr: continue
    ov=[(t,len(g&enr)) for t,g in tf.items()]
    ov=[x for x in ov if x[1]>0]; ov.sort(key=lambda x:-x[1])
    top=ov[:TOPK]
    core=set()
    for t,_ in top: core|=(tf[t]&enr)
    if not core: continue
    name=f"MasterReg_{safe(tis)}_TFbound_enriched"
    dd=os.path.join(OUT,name); os.makedirs(dd,exist_ok=True); core=sorted(core)
    desc=f"{tis}: tissue-enriched genes specifically bound by top-{TOPK} TFs (candidate master regulators: "+", ".join(f"{t}({c})" for t,c in top)+")"
    open(dd+"/geneset.tsv","w").write("gene\n"+"\n".join(core)+"\n")
    open(dd+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(core)+"\n")
    cite="GTEx tissue-enriched genes intersected with corrected ENCODE TF ChIP specific targets; top TFs = candidate master regulators. NIH-public."
    json.dump({"standard_name":name,"library":"NIH_master_regulator_per_tissue","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(core),"organism":"human","tissue":tis,"top_master_regulators":[{"tf":t,"n_enriched_targets":c} for t,c in top],"derived_in_this_work":True,"source":cite},open(dd+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"master_regulator_per_tissue","inputs":["GTEx tissue-enriched (Common Fund)","ENCODE TF corrected regulons (NHGRI)"],"public":True,"source_citation":cite,"funding":"NIH Common Fund (GTEx) + NIH/NHGRI (ENCODE)"},open(dd+"/geneset.provenance.json","w"),indent=1)
    n+=1
print("master-regulator-per-tissue sets:",n)
