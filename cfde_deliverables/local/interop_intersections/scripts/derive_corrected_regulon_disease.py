#!/usr/bin/env python3
# Corrected (background-de-biased) TF/eCLIP "specific_Up" regulons  ∩  ClinVar pathogenic genes
# -> per-factor "specifically-bound disease targets": disease genes each TF/RBP targets ABOVE background
# (promiscuous/HOT-region binding already removed). Also a union "disease genes specifically bound by any factor".
import os, json, glob, csv, collections
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"interop_intersections","output","corrected_regulon_disease"); os.makedirs(OUT,exist_ok=True)
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")} if os.path.exists(fp) else set()
clinvar=load(HOME+"/clinvar_genesets/output/all/ClinVar_pathogenic_genes_all/geneset.tsv")
print("ClinVar pathogenic genes:",len(clinvar))
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(name,desc,genes,srcs,extra):
    if not genes: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="Intersection of "+" & ".join(srcs)+" (NIH-funded, public); binding is background-corrected (specificity)."
    m={"standard_name":name,"library":"NIH_corrected_regulon_disease","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"intersection_of":srcs,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"corrected_regulon_x_clinvar","inputs":srcs,"public":True,"source_citation":cite,"funding":"NIH/NHGRI (ENCODE) + NIH/NLM (ClinVar)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1
tot=0
for assay,base in [("TF","ENCODE_TF_regulon_bgcorrected"),("eCLIP","ENCODE_eCLIP_regulon_bgcorrected")]:
    n=0; unioncov=set()
    for d in glob.glob(os.path.join(HOME,"encode_regulons_bgcontrast","output",base,"*_specific_Up")):
        try: tgt=json.load(open(os.path.join(d,"geneset.meta.json"))).get("target",os.path.basename(d))
        except Exception: tgt=os.path.basename(d)
        inter=load(os.path.join(d,"geneset.tsv"))&clinvar
        if not inter: continue
        unioncov|=inter
        n+=emit(f"ENCODE_{assay}_{safe(tgt)}_specific_disease_targets",f"{tgt}: specifically-bound ({assay}) ClinVar disease genes",inter,[f"ENCODE {assay} corrected regulon: {tgt} (NHGRI)","ClinVar (NIH/NLM)"],{"assay":assay,"target":tgt})
    n+=emit(f"ENCODE_{assay}_any_factor_specific_disease_targets",f"ClinVar disease genes specifically bound by >=1 {assay} factor",unioncov,[f"ENCODE {assay} corrected regulons (NHGRI)","ClinVar (NIH/NLM)"],{"assay":assay})
    print(f"{assay}: {n} sets | union disease targets: {len(unioncov)}")
    tot+=n
print("TOTAL corrected-regulon x disease sets:",tot)
