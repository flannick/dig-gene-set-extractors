import os, json, glob, csv
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")} if os.path.exists(fp) else set()
targets={
 "ClinGenHI":(load(HOME+"/clingen_dosage/output/ClinGen_haploinsufficient_score3/geneset.tsv"),"ClinGen haploinsufficient (HI=3)"),
 "GlyGenObserved":(load(HOME+"/glygen_genesets/output/observed/GlyGen_observed_glycoproteins/geneset.tsv"),"GlyGen observed glycoproteins"),
}
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(sub,name,desc,genes,srcs,extra):
    if not genes: return 0
    d=os.path.join(HOME,"interop_intersections","output",sub,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite="Intersection of "+" & ".join(srcs)+" (NIH-public); binding background-corrected."
    m={"standard_name":name,"library":"NIH_corrected_regulon_"+sub,"description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"corrected_regulon_cross","inputs":srcs,"public":True,"source_citation":cite,"funding":"NIH"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1
for tkey,(tset,tdesc) in targets.items():
    for assay,base in [("TF","ENCODE_TF_regulon_bgcorrected"),("eCLIP","ENCODE_eCLIP_regulon_bgcorrected")]:
        n=0; cov=set()
        for d in glob.glob(os.path.join(HOME,"encode_regulons_bgcontrast","output",base,"*_specific_Up")):
            try: tgt=json.load(open(os.path.join(d,"geneset.meta.json"))).get("target",os.path.basename(d))
            except Exception: tgt=os.path.basename(d)
            inter=load(os.path.join(d,"geneset.tsv"))&tset
            if not inter: continue
            cov|=inter
            n+=emit(f"corrected_regulon_{tkey}",f"ENCODE_{assay}_{safe(tgt)}_specific_x_{tkey}",f"{tgt}: specifically-bound ({assay}) {tdesc}",inter,[f"ENCODE {assay} corrected regulon: {tgt}",tdesc],{"assay":assay,"target":tgt})
        n+=emit(f"corrected_regulon_{tkey}",f"ENCODE_{assay}_any_specific_x_{tkey}",f"{tdesc} specifically bound by >=1 {assay} factor",cov,[f"ENCODE {assay} corrected regulons",tdesc],{"assay":assay})
        print(f"{tkey} {assay}: {n} sets | union {len(cov)}")
