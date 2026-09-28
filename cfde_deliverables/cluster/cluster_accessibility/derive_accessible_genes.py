#!/usr/bin/env python
"""Derive an 'accessible genes' gene set FROM an ENCODE/4DN ATAC-seq peak BED.
A gene is 'accessible' if an ATAC peak overlaps its promoter (TSS +/- WIN bp).
PROVENANCE: derived from the provider's PUBLIC peaks (cited); not their data re-generated.
Generalized: call once per peak file/biosample; loop over many to build the full deliverable on-cluster.
Usage: derive_accessible_genes.py <peak.bed[.gz]> <refGene_hg38.txt.gz> <biosample> <file_acc> <exp_acc> <outdir>"""
import gzip, csv, json, os, sys, collections
PEAK,REFGENE,BIOSAMPLE,FILEACC,EXPACC,OUTDIR = sys.argv[1:7]
WIN=1000; BIN=1000   # promoter half-window; genome binning for fast overlap

def opn(p): return gzip.open(p,'rt') if p.endswith('.gz') else open(p)

# refGene hg38 -> gene -> {(chrom, TSS)}
gene_tss=collections.defaultdict(set)
for line in opn(REFGENE):
    f=line.rstrip('\n').split('\t')
    if len(f)<13: continue
    chrom,strand,txS,txE,gene=f[2],f[3],f[4],f[5],f[12]
    if '_' in chrom or not txS.isdigit(): continue
    tss=int(txS) if strand=='+' else int(txE)
    gene_tss[gene].add((chrom,tss))

# peaks -> covered 1kb bins per chrom
covered=collections.defaultdict(set); npk=0
for line in opn(PEAK):
    f=line.rstrip('\n').split('\t')
    if len(f)<3 or not f[1].isdigit(): continue
    chrom=f[0]; s=int(f[1]); e=int(f[2])
    for b in range(s//BIN, e//BIN+1): covered[chrom].add(b)
    npk+=1

# accessible = promoter bin overlaps a peak bin
acc=[]
for gene,tsss in gene_tss.items():
    for chrom,tss in tsss:
        cb=covered.get(chrom)
        if cb and any(b in cb for b in range((tss-WIN)//BIN,(tss+WIN)//BIN+1)):
            acc.append(gene); break
acc=sorted(set(acc))

name=f"ENCODE_ATAC_accessible_{BIOSAMPLE}"
cite=(f"Derived from ENCODE ATAC-seq peak file {FILEACC} (experiment {EXPACC}, biosample {BIOSAMPLE}), "
      f"GRCh38, ENCODE/NHGRI, public; gene TSS from UCSC refGene hg38, public.")
d=os.path.join(OUTDIR,name); os.makedirs(d,exist_ok=True)
open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(acc)+"\n")
open(d+"/genesets.gmt","w").write(name+"\t"+f"Genes with promoter-proximal (TSS+/-{WIN}bp) ENCODE ATAC accessibility in {BIOSAMPLE}"+"\t"+"\t".join(acc)+"\n")
json.dump({"standard_name":name,"library":"ENCODE_ATAC_accessible","description":
   f"Genes with promoter-proximal ENCODE ATAC-seq accessibility in {BIOSAMPLE} (DERIVED from ENCODE public peaks).",
   "version":"0.1","file_type":"geneset","n_genes":len(acc),"organism":"human","assembly":"GRCh38",
   "biosample":BIOSAMPLE,"derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
json.dump({"focus":name,"operation":"derive_accessible_genes_from_ATAC",
   "method":f"ATAC peak overlap with promoter (TSS +/-{WIN}bp), {BIN}bp binning",
   "inputs":[f"ENCODE {FILEACC} peaks (GRCh38; ENCODE/NHGRI; public)","UCSC refGene hg38 TSS (public)"],
   "source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"},open(d+"/geneset.provenance.json","w"),indent=1)
print(f"peaks={npk} | genes-with-TSS={len(gene_tss)} | ACCESSIBLE genes={len(acc)}")
print("sample:", acc[:15])
print("wrote", d)
