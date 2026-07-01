#!/usr/bin/env python3
# Cross-NIH-resource intersection gene sets (interoperability deliverable). Ready-now crosses that don't
# depend on still-running jobs: MitoCarta/eGenes/sGenes x ClinVar, and eGenes x GTEx tissue-enriched.
import os, json, glob, csv
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"interop_intersections","output"); os.makedirs(OUT,exist_ok=True)
def load_set(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")}
def norm(t): return "".join(c for c in t.lower() if c.isalnum())
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(sub,name,desc,genes,srcs,extra):
    d=os.path.join(OUT,sub,name); os.makedirs(d,exist_ok=True); genes=sorted(g for g in genes if g)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="Intersection of "+ " & ".join(srcs)+" (all NIH-funded, public)."
    m={"standard_name":name,"library":"NIH_interop_intersection","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"intersection_of":srcs,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"cross_resource_intersection","inputs":srcs,"source_citation":cite,"public":True,"funding":"NIH (multiple resources)"},open(d+"/geneset.provenance.json","w"),indent=1)

clinvar_all=load_set(HOME+"/clinvar_genesets/output/all/ClinVar_pathogenic_genes_all/geneset.tsv")
mito={l.split("\t")[0] for i,l in enumerate(open(HOME+"/mito_genesets/output/mitocarta3_reference.tsv")) if i>0 and l.strip()}
egenes={norm(os.path.basename(d).replace("GTEx_eQTL_eGenes_","")):load_set(d+"/geneset.tsv") for d in glob.glob(HOME+"/gtex_qtl/output_eqtl/GTEx_eQTL_eGenes_*")}
sgenes={norm(os.path.basename(d).replace("GTEx_sQTL_sGenes_","")):load_set(d+"/geneset.tsv") for d in glob.glob(HOME+"/gtex_qtl/output_sqtl/GTEx_sQTL_sGenes_*")}
tenr={norm(os.path.basename(d).replace("GTEx_tissue_enriched_","")):(os.path.basename(d).replace("GTEx_tissue_enriched_",""),load_set(d+"/geneset.tsv")) for d in glob.glob(HOME+"/gtex_tissue_enriched/output/GTEx_tissue_enriched_*")}
egenes_union=set().union(*egenes.values()) if egenes else set()
sgenes_union=set().union(*sgenes.values()) if sgenes else set()

# 1. MitoCarta ∩ ClinVar
emit("mito_disease","MitoCarta_x_ClinVar_disease_genes","Mitochondrial genes (MitoCarta3.0) with pathogenic ClinVar variants",
     mito&clinvar_all,["MitoCarta3.0 (NIH/NIGMS)","ClinVar (NIH/NLM)"],{})
# 2. eGenes ∩ ClinVar ; 3. sGenes ∩ ClinVar
emit("qtl_disease","GTEx_eQTL_eGenes_x_ClinVar_disease_genes","Disease genes (ClinVar pathogenic) that have a GTEx cis-eQTL in any tissue",
     egenes_union&clinvar_all,["GTEx cis-eQTL eGenes (NIH Common Fund)","ClinVar (NIH/NLM)"],{})
emit("qtl_disease","GTEx_sQTL_sGenes_x_ClinVar_disease_genes","Disease genes (ClinVar pathogenic) that have a GTEx splicing QTL in any tissue",
     sgenes_union&clinvar_all,["GTEx cis-sQTL sGenes (NIH Common Fund)","ClinVar (NIH/NLM)"],{})
print(f"MitoCarta∩ClinVar: {len(mito&clinvar_all)} | eGenes∩ClinVar: {len(egenes_union&clinvar_all)} | sGenes∩ClinVar: {len(sgenes_union&clinvar_all)}")
# 4. eGenes ∩ tissue-enriched (matched tissue)
n=0
for k,(label,enr) in tenr.items():
    if k in egenes:
        inter=egenes[k]&enr
        if not inter: continue
        sn=f"GTEx_eGenes_x_tissue_enriched_{safe(label)}"
        emit("eqtl_tissuespecific",sn,f"Genes both GTEx tissue-enriched AND cis-eQTL eGene in {label}",inter,
             ["GTEx cis-eQTL eGenes (NIH Common Fund)","GTEx tissue-enriched (NIH Common Fund)"],{"tissue":label}); n+=1
print(f"eGenes∩tissue-enriched matched-tissue sets: {n}")
