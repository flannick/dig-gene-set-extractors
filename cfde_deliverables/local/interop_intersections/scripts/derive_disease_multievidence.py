#!/usr/bin/env python3
# Multi-evidence disease-gene integration: ClinVar pathogenic genes annotated by regulatory evidence axes:
#   eGene (GTEx cis-eQTL, any tissue) | specifically-BOUND (corrected TF/eCLIP, any factor) |
#   specifically-ACCESSIBLE (corrected ATAC/DNase/histone, any biosample).
# Emits per-axis disease sets + the multi-evidence intersection (highest-confidence regulated disease genes)
# + tiers by number of supporting axes.
import os, glob, json, csv, collections
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"interop_intersections","output","disease_multievidence"); os.makedirs(OUT,exist_ok=True)
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")} if os.path.exists(fp) else set()
def union(glob_pat):
    u=set()
    for f in glob.glob(glob_pat): u|=load(f)
    return u
clinvar=load(HOME+"/clinvar_genesets/output/all/ClinVar_pathogenic_genes_all/geneset.tsv")
egene=union(HOME+"/gtex_qtl/output_eqtl/GTEx_eQTL_eGenes_*/geneset.tsv")
bound=union(HOME+"/encode_regulons_bgcontrast/output/ENCODE_TF_regulon_bgcorrected/*_specific_Up/geneset.tsv") | \
      union(HOME+"/encode_regulons_bgcontrast/output/ENCODE_eCLIP_regulon_bgcorrected/*_specific_Up/geneset.tsv")
access=union(HOME+"/accessibility_bgcontrast/output/ENCODE_ATAC_accessible_bgcontrast/*_accessible_Up/geneset.tsv") | \
       union(HOME+"/accessibility_bgcontrast/output/ENCODE_DNase_accessible_bgcontrast/*_accessible_Up/geneset.tsv") | \
       union(HOME+"/accessibility_bgcontrast/output/ENCODE_histone_bgcontrast/*_accessible_Up/geneset.tsv")
print(f"ClinVar={len(clinvar)} eGene-union={len(egene)} bound-union={len(bound)} access-union={len(access)}")
axes={"eQTL_regulated":egene,"specifically_bound":bound,"specifically_accessible":access}
def emit(name,desc,genes,srcs,extra):
    if not genes: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="ClinVar disease genes x regulatory evidence: "+" & ".join(srcs)+" (all NIH-public)."
    m={"standard_name":name,"library":"NIH_disease_multievidence","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"disease_multievidence","inputs":srcs,"public":True,"source_citation":cite,"funding":"NIH (multiple)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1
# per-axis disease sets
for ax,s in axes.items():
    emit(f"ClinVar_disease_x_{ax}",f"ClinVar disease genes that are {ax.replace('_',' ')}",clinvar&s,["ClinVar (NIH/NLM)",ax],{"axis":ax})
# multi-evidence intersection + tiers
tier=collections.Counter()
allthree=set(); tier_sets=collections.defaultdict(set)
for g in clinvar:
    k=sum(g in s for s in axes.values())
    tier[k]+=1
    if k>=1: tier_sets[k].add(g)
    if k==3: allthree.add(g)
emit("ClinVar_disease_multievidence_all3",f"Disease genes supported by ALL 3 regulatory axes (eQTL + specific binding + specific accessibility)",allthree,["ClinVar","GTEx eQTL","ENCODE corrected binding","ENCODE corrected accessibility"],{"n_axes":3})
for k in (2,3):
    u=set().union(*[v for kk,v in tier_sets.items() if kk>=k]) if tier_sets else set()
    emit(f"ClinVar_disease_ge{k}_axes",f"Disease genes supported by >= {k} regulatory axes",u,["ClinVar","GTEx","ENCODE"],{"min_axes":k})
print("disease genes by #evidence axes:",dict(tier))
print("all-3-axes disease genes:",len(allthree))
