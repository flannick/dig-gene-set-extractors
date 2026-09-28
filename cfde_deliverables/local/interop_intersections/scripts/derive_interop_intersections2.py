#!/usr/bin/env python3
# Batch 2 of cross-NIH-resource intersections (all inputs ready). Adds constraint/glyco/QTL/disease crosses.
import os, json, glob, csv
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"interop_intersections","output")
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")} if os.path.exists(fp) else set()
def norm(t): return "".join(c for c in t.lower() if c.isalnum())
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(sub,name,desc,genes,srcs,extra=None):
    d=os.path.join(OUT,sub,name); os.makedirs(d,exist_ok=True); genes=sorted(g for g in genes if g)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="Intersection of "+" & ".join(srcs)+" (all NIH-funded, public)."
    m={"standard_name":name,"library":"NIH_interop_intersection","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"intersection_of":srcs,"source":cite}
    if extra: m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"cross_resource_intersection","inputs":srcs,"source_citation":cite,"public":True,"funding":"NIH (multiple resources)"},open(d+"/geneset.provenance.json","w"),indent=1)

clinvar=load(HOME+"/clinvar_genesets/output/all/ClinVar_pathogenic_genes_all/geneset.tsv")
mito={l.split("\t")[0] for i,l in enumerate(open(HOME+"/mito_genesets/output/mitocarta3_reference.tsv")) if i>0 and l.strip()}
constr=load(HOME+"/gnomad_constraint/output/gnomAD_LoF_constrained_LOEUF_lt0.35/geneset.tsv")
clingen_hi=load(HOME+"/clingen_dosage/output/ClinGen_haploinsufficient_score3/geneset.tsv")
glygen=load(HOME+"/glygen_genesets/output/observed/GlyGen_observed_glycoproteins/geneset.tsv")
egenes={norm(os.path.basename(d).replace("GTEx_eQTL_eGenes_","")):load(d+"/geneset.tsv") for d in glob.glob(HOME+"/gtex_qtl/output_eqtl/GTEx_eQTL_eGenes_*")}
sgenes={norm(os.path.basename(d).replace("GTEx_sQTL_sGenes_","")):load(d+"/geneset.tsv") for d in glob.glob(HOME+"/gtex_qtl/output_sqtl/GTEx_sQTL_sGenes_*")}
tenr={norm(os.path.basename(d).replace("GTEx_tissue_enriched_","")):(os.path.basename(d).replace("GTEx_tissue_enriched_",""),load(d+"/geneset.tsv")) for d in glob.glob(HOME+"/gtex_tissue_enriched/output/GTEx_tissue_enriched_*")}
eg_u=set().union(*egenes.values()) if egenes else set(); sg_u=set().union(*sgenes.values()) if sgenes else set()

C=lambda a,b: a&b
emit("constraint_disease","gnomAD_constrained_x_ClinVar","Highly LoF-constrained (LOEUF<0.35) genes with pathogenic ClinVar variants",C(constr,clinvar),["gnomAD constraint (Broad/NIH)","ClinVar (NIH/NLM)"])
emit("constraint_disease","gnomAD_constrained_x_ClinGen_HI","LoF-constrained genes that are ClinGen haploinsufficient (HI=3)",C(constr,clingen_hi),["gnomAD constraint (Broad/NIH)","ClinGen dosage (NHGRI)"])
emit("constraint_disease","gnomAD_constrained_x_MitoCarta","LoF-constrained mitochondrial genes",C(constr,mito),["gnomAD constraint (Broad/NIH)","MitoCarta3.0 (NIH/NIGMS)"])
emit("mito_qtl","MitoCarta_x_GTEx_eGenes","Mitochondrial genes with a GTEx cis-eQTL (any tissue)",C(mito,eg_u),["MitoCarta3.0 (NIH/NIGMS)","GTEx eQTL (Common Fund)"])
emit("mito_qtl","MitoCarta_x_GTEx_sGenes","Mitochondrial genes with a GTEx splicing QTL (any tissue)",C(mito,sg_u),["MitoCarta3.0 (NIH/NIGMS)","GTEx sQTL (Common Fund)"])
emit("glyco_cross","GlyGen_observed_x_ClinVar","Observed glycoproteins with pathogenic ClinVar variants",C(glygen,clinvar),["GlyGen observed glycoproteins (NIH/NIGMS)","ClinVar (NIH/NLM)"])
emit("glyco_cross","GlyGen_observed_x_gnomAD_constrained","Observed glycoproteins that are LoF-constrained",C(glygen,constr),["GlyGen observed glycoproteins (NIH/NIGMS)","gnomAD constraint (Broad/NIH)"])
print("constraint/mito/glyco crosses done")

# per-tissue: sQTL x tissue-enriched ; eGenes ∩ sGenes
n1=n2=n3=0
for k,(label,enr) in tenr.items():
    if k in sgenes and sgenes[k]&enr:
        emit("sqtl_tissuespecific",f"GTEx_sGenes_x_tissue_enriched_{safe(label)}",f"Genes both GTEx tissue-enriched AND sQTL sGene in {label}",sgenes[k]&enr,["GTEx sQTL (Common Fund)","GTEx tissue-enriched (Common Fund)"],{"tissue":label}); n1+=1
    if clinvar and enr&clinvar:
        emit("disease_tissuespecific",f"ClinVar_x_tissue_enriched_{safe(label)}",f"Disease genes (ClinVar) tissue-enriched in {label}",enr&clinvar,["ClinVar (NIH/NLM)","GTEx tissue-enriched (Common Fund)"],{"tissue":label}); n3+=1
for k in egenes:
    if k in sgenes and egenes[k]&sgenes[k]:
        lbl=k
        emit("eqtl_sqtl",f"GTEx_eGenes_AND_sGenes_{safe(k)}",f"Genes with BOTH a cis-eQTL and a splicing QTL in tissue {k}",egenes[k]&sgenes[k],["GTEx eQTL (Common Fund)","GTEx sQTL (Common Fund)"],{"tissue":k}); n2+=1
print(f"sQTL∩enriched: {n1} | eGenes∩sGenes: {n2} | ClinVar∩enriched: {n3}")

# per-RBP eCLIP regulon ∩ ClinVar (disease-relevant RBP targets)
nr=0
for d in glob.glob(HOME+"/encode_regulons/output_eclip/regulon/ENCODE_eCLIP_RBP_target_*"):
    rbp=os.path.basename(d).replace("ENCODE_eCLIP_RBP_target_","")
    inter=load(d+"/geneset.tsv")&clinvar
    if not inter: continue
    emit("rbp_disease",f"eCLIP_{safe(rbp)}_targets_x_ClinVar",f"{rbp} eCLIP target genes that are ClinVar disease genes",inter,["ENCODE eCLIP (NHGRI)","ClinVar (NIH/NLM)"],{"rbp":rbp}); nr+=1
print(f"eCLIP-RBP ∩ ClinVar sets: {nr}")
