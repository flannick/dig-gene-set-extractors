#!/usr/bin/env python3
# Batch 3 intersections: constraint/dosage/glyco x disease/QTL/tissue, plus high-value TRIPLE intersections
# (disease + constrained + tissue-specific). All inputs already on disk; no dependence on running jobs.
import os, json, glob, csv
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
    cite="Intersection of "+" & ".join(srcs)+" (all NIH-funded, public)."
    m={"standard_name":name,"library":"NIH_interop_intersection","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"intersection_of":srcs,"source":cite}
    if extra: m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"cross_resource_intersection","inputs":srcs,"source_citation":cite,"public":True,"funding":"NIH (multiple resources)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1

clinvar=load(HOME+"/clinvar_genesets/output/all/ClinVar_pathogenic_genes_all/geneset.tsv")
mito={l.split("\t")[0] for i,l in enumerate(open(HOME+"/mito_genesets/output/mitocarta3_reference.tsv")) if i>0 and l.strip()}
constr=load(HOME+"/gnomad_constraint/output/gnomAD_LoF_constrained_LOEUF_lt0.35/geneset.tsv")
mis=load(HOME+"/gnomad_constraint/output/gnomAD_missense_constrained_misZ_ge3.09/geneset.tsv")
hi=load(HOME+"/clingen_dosage/output/ClinGen_haploinsufficient_score3/geneset.tsv")
glygen=load(HOME+"/glygen_genesets/output/observed/GlyGen_observed_glycoproteins/geneset.tsv")
egenes={norm(os.path.basename(d).replace("GTEx_eQTL_eGenes_","")):load(d+"/geneset.tsv") for d in glob.glob(HOME+"/gtex_qtl/output_eqtl/GTEx_eQTL_eGenes_*")}
sgenes={norm(os.path.basename(d).replace("GTEx_sQTL_sGenes_","")):load(d+"/geneset.tsv") for d in glob.glob(HOME+"/gtex_qtl/output_sqtl/GTEx_sQTL_sGenes_*")}
tenr={norm(os.path.basename(d).replace("GTEx_tissue_enriched_","")):(os.path.basename(d).replace("GTEx_tissue_enriched_",""),load(d+"/geneset.tsv")) for d in glob.glob(HOME+"/gtex_tissue_enriched/output/GTEx_tissue_enriched_*")}
eg_u=set().union(*egenes.values()) if egenes else set(); sg_u=set().union(*sgenes.values()) if sgenes else set()

c=0
c+=emit("dosage_cross","ClinGen_HI_x_ClinVar","ClinGen haploinsufficient genes with pathogenic ClinVar variants",hi&clinvar,["ClinGen dosage (NHGRI)","ClinVar (NIH/NLM)"])
c+=emit("dosage_cross","ClinGen_HI_x_MitoCarta","Haploinsufficient mitochondrial genes",hi&mito,["ClinGen dosage (NHGRI)","MitoCarta3.0 (NIGMS)"])
c+=emit("constraint_cross","gnomAD_constrained_x_eGenes","LoF-constrained genes with a GTEx cis-eQTL",constr&eg_u,["gnomAD constraint (Broad/NIH)","GTEx eQTL (Common Fund)"])
c+=emit("constraint_cross","gnomAD_constrained_x_sGenes","LoF-constrained genes with a GTEx splicing QTL",constr&sg_u,["gnomAD constraint (Broad/NIH)","GTEx sQTL (Common Fund)"])
c+=emit("constraint_cross","gnomAD_missense_constrained_x_ClinVar","Missense-constrained genes with pathogenic ClinVar variants",mis&clinvar,["gnomAD missense constraint (Broad/NIH)","ClinVar (NIH/NLM)"])
c+=emit("glyco_cross","GlyGen_observed_x_eGenes","Observed glycoproteins with a GTEx cis-eQTL",glygen&eg_u,["GlyGen observed (NIGMS)","GTEx eQTL (Common Fund)"])
c+=emit("glyco_cross","GlyGen_observed_x_MitoCarta","Observed glycoproteins that are mitochondrial",glygen&mito,["GlyGen observed (NIGMS)","MitoCarta3.0 (NIGMS)"])
print("global crosses:",c)

# per-tissue single crosses + triples
n=collections.Counter() if False else {"constr_t":0,"hi_t":0,"gly_t":0,"triple_cvt":0,"triple_egt":0}
for k,(label,enr) in tenr.items():
    n["constr_t"]+=emit("constraint_tissue",f"gnomAD_constrained_x_tissue_enriched_{safe(label)}",f"LoF-constrained genes tissue-enriched in {label}",constr&enr,["gnomAD constraint (Broad/NIH)","GTEx tissue-enriched (Common Fund)"],{"tissue":label})
    n["hi_t"]+=emit("dosage_tissue",f"ClinGen_HI_x_tissue_enriched_{safe(label)}",f"Haploinsufficient genes tissue-enriched in {label}",hi&enr,["ClinGen dosage (NHGRI)","GTEx tissue-enriched (Common Fund)"],{"tissue":label})
    n["gly_t"]+=emit("glyco_tissue",f"GlyGen_observed_x_tissue_enriched_{safe(label)}",f"Observed glycoproteins tissue-enriched in {label}",glygen&enr,["GlyGen observed (NIGMS)","GTEx tissue-enriched (Common Fund)"],{"tissue":label})
    # TRIPLE: disease + constrained + tissue-specific
    n["triple_cvt"]+=emit("triple_disease_constrained_tissue",f"ClinVar_AND_constrained_AND_enriched_{safe(label)}",f"Disease genes that are LoF-constrained AND tissue-enriched in {label}",clinvar&constr&enr,["ClinVar (NIH/NLM)","gnomAD constraint (Broad/NIH)","GTEx tissue-enriched (Common Fund)"],{"tissue":label})
    # TRIPLE: disease gene with eQTL in the very tissue it is specific to
    if k in egenes:
        n["triple_egt"]+=emit("triple_disease_eqtl_tissue",f"ClinVar_AND_eGene_AND_enriched_{safe(label)}",f"Disease genes that are GTEx eGenes AND tissue-enriched in {label}",clinvar&egenes[k]&enr,["ClinVar (NIH/NLM)","GTEx eQTL (Common Fund)","GTEx tissue-enriched (Common Fund)"],{"tissue":label})
print("per-tissue/triple:",dict(n))
import collections
print("TOTAL interop sets now:",len(glob.glob(os.path.join(OUT,"*","*","geneset.tsv"))))
