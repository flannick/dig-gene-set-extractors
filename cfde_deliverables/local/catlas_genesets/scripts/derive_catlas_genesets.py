#!/usr/bin/env python3
# CATLAS (BICCN human scATAC) cell-type accessible-gene sets. Inputs = ALREADY background-controlled
# specifically-accessible regions (*_Up.bed.gz; prevalence-bg / LOO / <0.25; 1kb bins) from the tissue-
# validation pipeline. This does the region->gene step (promoter TSS+/-1kb overlap) -> contract sets.
import os, gzip, glob, json, collections
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
CB="/Users/gage/Codex/PIGEAN_EAGGL/repos/factor-pgs-tissue-sc-validation/data/processed/catlas_bgcontrol"
REF="/Users/gage/.claude/jobs/32851e29/tmp/refGene_hg38_clean.txt.gz"
OUT=os.path.join(HOME,"catlas_genesets","output"); os.makedirs(OUT,exist_ok=True); BIN=1000; WIN=1000
LOG=open(os.path.join(OUT,"batch_log.txt"),"a")
def log(m): LOG.write(m+"\n"); LOG.flush(); print(m)
promo=collections.defaultdict(set)
for line in gzip.open(REF,'rt'):
    f=line.rstrip('\n').split('\t')
    if len(f)<13 or '_' in f[2] or not f[4].isdigit() or not f[5].isdigit(): continue
    tss=int(f[4]) if f[3]=='+' else int(f[5])
    for b in range((tss-WIN)//BIN,(tss+WIN)//BIN+1): promo[f[12]].add((f[2],b))
log(f"refGene genes: {len(promo)}")
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def covered(fp):
    cov=set()
    for line in gzip.open(fp,'rt'):
        f=line.rstrip('\n').split('\t')
        if len(f)<3 or not f[1].isdigit(): continue
        for b in range(int(f[1])//BIN,int(f[2])//BIN+1): cov.add((f[0],b))
    return cov
n=0
for fp in sorted(glob.glob(os.path.join(CB,"*_Up.bed.gz"))):
    ct=os.path.basename(fp).replace("_Up.bed.gz","").replace("___"," / ")
    cov=covered(fp)
    genes=sorted(g for g,pb in promo.items() if pb & cov)
    if not genes: continue
    name=f"CATLAS_{safe(os.path.basename(fp).replace('_Up.bed.gz',''))}_accessible_Up"
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\tGenes with promoters in specifically-accessible regions ({ct}); CATLAS scATAC, background-controlled\t"+"\t".join(genes)+"\n")
    cite=("Genes mapped (promoter TSS+/-1kb) from background-prevalence-controlled cell-type specifically-accessible "
          "regions (LOO, prevalence<0.25, 1kb bins) derived from CATLAS human scATAC (Zhang et al. 2021 Cell; BICCN/NIH; "
          "hg38) [+ ENCODE tissue comparators; see chromatin_manifest.py]. Same control method as our ENCODE accessibility. Public.")
    json.dump({"standard_name":name,"library":"CATLAS_scATAC_accessible_bgcontrolled","description":f"Background-controlled specifically-accessible genes in {ct} (CATLAS scATAC)","version":"1.0","file_type":"geneset","n_genes":len(genes),"organism":"human","assembly":"GRCh38","cell_type":ct,"method":"prevalence_background_LOO_region_then_promoter_map","derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"catlas_accessible_genes","inputs":["CATLAS scATAC cell-type cCREs (Zhang 2021 Cell; BICCN/NIH; public)","background-prevalence control (LOO)","UCSC refGene hg38 promoters"],"source_citation":cite,"public":True,"funding":"NIH BRAIN Initiative / BICCN (CATLAS)"},open(d+"/geneset.provenance.json","w"),indent=1)
    n+=1; log(f"[{n}] {ct}: {len(genes)} genes")
log(f"=== DONE CATLAS: {n} accessible-gene sets ===")
LOG.close()
