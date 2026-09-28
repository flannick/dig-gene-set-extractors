#!/usr/bin/env python3
# Build hg38 promoter->gene bin map for CATLAS fragment aggregation.
# refGene: bin(0) name(1) chrom(2) strand(3) txStart(4) txEnd(5) ... name2(12)=gene symbol.
# TSS = txStart if '+' else txEnd. Promoter window = TSS +/- WIN. Genome binned into BIN bp.
# Output: promoter_bins.tsv with "chrom:binidx \t gene1;gene2;..." (one line per occupied bin).
import os, gzip, urllib.request, collections
HERE=os.path.dirname(os.path.abspath(__file__))
WORK=os.environ.get("WORK", HERE)
WIN=int(os.environ.get("PROMOTER_WIN","1000"))     # +/- around TSS
BIN=int(os.environ.get("BIN_BP","1000"))
REF=os.path.join(WORK,"refGene.hg38.txt.gz")
OUT=os.path.join(WORK,"promoter_bins.tsv")
URL="http://hgdownload.soc.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
if not os.path.exists(REF):
    print("downloading refGene hg38 ...", flush=True)
    tmp=REF+".part"; urllib.request.urlretrieve(URL,tmp); os.replace(tmp,REF)
bin2genes=collections.defaultdict(set); ngene=set()
op=gzip.open(REF,"rt")
for line in op:
    f=line.rstrip("\n").split("\t")
    if len(f)<13: continue
    chrom=f[2]; st=f[3]
    if "_" in chrom or chrom=="chrM": continue     # keep primary chromosomes
    try: tss=int(f[4]) if st=="+" else int(f[5])
    except ValueError: continue
    gene=f[12]
    if not gene: continue
    ngene.add(gene)
    lo=tss-WIN; hi=tss+WIN
    for b in range(lo//BIN, hi//BIN + 1):
        bin2genes[(chrom,b)].add(gene)
with open(OUT,"w") as o:
    for (chrom,b),genes in bin2genes.items():
        o.write(f"{chrom}:{b}\t{';'.join(sorted(genes))}\n")
print(f"refGene genes: {len(ngene)} | occupied promoter bins: {len(bin2genes)} -> {OUT}", flush=True)
