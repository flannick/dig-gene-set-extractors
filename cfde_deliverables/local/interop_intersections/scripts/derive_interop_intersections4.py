#!/usr/bin/env python3
# Batch 4 (Tier-A only; no gnomAD/HPO): ClinGen dosage x GTEx QTL, and GTEx eGene tissue-specificity.
import os, json, glob, csv, collections
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"interop_intersections","output")
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")} if os.path.exists(fp) else set()
def norm(t): return "".join(c for c in t.lower() if c.isalnum())
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(sub,name,desc,genes,srcs,extra=None):
    if not genes: return 0
    d=os.path.join(OUT,sub,name); os.makedirs(d,exist_ok=True); genes=sorted(g for g in genes if g)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="Intersection/derivation of "+" & ".join(srcs)+" (all NIH-funded, public)."
    m={"standard_name":name,"library":"NIH_interop_intersection","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"intersection_of":srcs,"source":cite}
    if extra: m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"cross_resource_intersection","inputs":srcs,"source_citation":cite,"public":True,"funding":"NIH (multiple resources)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1

hi=load(HOME+"/clingen_dosage/output/ClinGen_haploinsufficient_score3/geneset.tsv")
hi_some=load(HOME+"/clingen_dosage/output/ClinGen_haploinsufficient_some_evidence/geneset.tsv")
egenes={norm(os.path.basename(d).replace("GTEx_eQTL_eGenes_","")):load(d+"/geneset.tsv") for d in glob.glob(HOME+"/gtex_qtl/output_eqtl/GTEx_eQTL_eGenes_*")}
sgenes={norm(os.path.basename(d).replace("GTEx_sQTL_sGenes_","")):load(d+"/geneset.tsv") for d in glob.glob(HOME+"/gtex_qtl/output_sqtl/GTEx_sQTL_sGenes_*")}
labels={norm(os.path.basename(d).replace("GTEx_eQTL_eGenes_","")):os.path.basename(d).replace("GTEx_eQTL_eGenes_","") for d in glob.glob(HOME+"/gtex_qtl/output_eqtl/GTEx_eQTL_eGenes_*")}
eg_u=set().union(*egenes.values()) if egenes else set(); sg_u=set().union(*sgenes.values()) if sgenes else set()

c=0
c+=emit("dosage_qtl","ClinGen_HI_x_GTEx_eGenes","Haploinsufficient genes (ClinGen HI=3) with a GTEx cis-eQTL",hi&eg_u,["ClinGen dosage (NHGRI)","GTEx eQTL (Common Fund)"])
c+=emit("dosage_qtl","ClinGen_HI_x_GTEx_sGenes","Haploinsufficient genes (ClinGen HI=3) with a GTEx splicing QTL",hi&sg_u,["ClinGen dosage (NHGRI)","GTEx sQTL (Common Fund)"])
c+=emit("dosage_qtl","ClinGen_HI_someEvidence_x_GTEx_eGenes","ClinGen haploinsufficient (HI 1-3) genes with a GTEx cis-eQTL",hi_some&eg_u,["ClinGen dosage (NHGRI)","GTEx eQTL (Common Fund)"])
print("ClinGen x QTL crosses:",c)

# eGene tissue-specificity: how many tissues is each gene an eGene in?
cnt=collections.Counter()
for s in egenes.values():
    for g in s: cnt[g]+=1
ntis=len(egenes)
constitutive={g for g,k in cnt.items() if k>=max(2,int(0.8*ntis))}
emit("eqtl_specificity","GTEx_constitutive_eGenes",f"Genes that are a cis-eQTL eGene in >=80% of GTEx tissues ({ntis} tissues; broadly genetically regulated)",constitutive,["GTEx eQTL (Common Fund)"])
nspec=0
for k,genes in egenes.items():
    spec={g for g in genes if cnt[g]<=2}   # eGene in <=2 tissues total => tissue-restricted regulation
    if not spec: continue
    nspec+=emit("eqtl_specificity",f"GTEx_tissue_specific_eGenes_{safe(labels[k])}",f"Genes that are a cis-eQTL eGene in {labels[k]} and in <=2 tissues total (tissue-restricted genetic regulation)",spec,["GTEx eQTL (Common Fund)"],{"tissue":labels[k]})
print(f"eGene specificity: constitutive {len(constitutive)} | tissue-specific sets {nspec}")
print("TOTAL shippable interop sets now:",len(glob.glob(os.path.join(OUT,"*","*","geneset.tsv"))))
