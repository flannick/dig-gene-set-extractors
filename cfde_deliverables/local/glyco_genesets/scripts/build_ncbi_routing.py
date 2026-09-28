#!/usr/bin/env python
"""NIH-native (NCBI Gene/RefSeq) secretory-routing layer to corroborate UniProt for C/D.
Streams NCBI gene2go (human, Cellular Component) — nothing large is stored; keeps a small
table restricted to our gene universe (T2D effectors + A/B enzymes)."""
import gzip, csv, io, os, re, urllib.request
BASE=os.path.expanduser("~/Claude/proj-valiation-challenge/glyco_genesets")
EFF=BASE+"/output/T2D_effectors.tsv"
A="DPAGT1 ALG1 ALG2 ALG3 ALG5 ALG6 ALG8 ALG9 ALG10 ALG10B ALG11 ALG12 ALG13 ALG14 DPM1 DPM2 DPM3 MPDU1 DOLK SRD5A3 RFT1 STT3A STT3B RPN1 RPN2 DDOST DAD1 OSTC TUSC3 MAGT1 OST4 KRTCAP2 MOGS GANAB PRKCSH MAN1B1".split()
B="MAN1A1 MAN1A2 MAN1C1 MAN2A1 MAN2A2 MGAT1 MGAT2 MGAT3 MGAT4A MGAT4B MGAT4C MGAT5 MGAT5B B4GALT1 B4GALT2 B4GALT3 B4GALT4 B4GALT5 B4GALT6 B4GALT7 ST3GAL1 ST3GAL2 ST3GAL3 ST3GAL4 ST3GAL5 ST3GAL6 ST6GAL1 ST6GAL2 ST6GALNAC1 ST6GALNAC2 ST6GALNAC3 ST6GALNAC4 ST6GALNAC6 ST8SIA1 ST8SIA2 ST8SIA4 FUT1 FUT2 FUT3 FUT4 FUT5 FUT6 FUT7 FUT8 FUT9 FUT10 FUT11 POFUT1 POFUT2 GALNT1 GALNT2 GALNT3 GALNT4 GALNT6 GALNT7 GALNT10 GALNT12 C1GALT1 C1GALT1C1 GCNT1 GCNT3".split()
universe=set(A)|set(B)|{r["gene"] for r in csv.DictReader(open(EFF),delimiter='\t')}
print("gene universe:",len(universe))

def stream(url):
    req=urllib.request.Request(url,headers={"User-Agent":"x"})
    return gzip.GzipFile(fileobj=urllib.request.urlopen(req,timeout=300))

# 1) symbol/synonym -> GeneID (human gene_info)
print("downloading Homo_sapiens.gene_info.gz ...")
sym2id={}
for line in io.TextIOWrapper(stream("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"),encoding="utf-8"):
    if line.startswith("#"): continue
    f=line.rstrip("\n").split("\t")
    gid=f[1]; sym=f[2]; syns=f[4]
    if sym in universe: sym2id[sym]=gid
    for s in syns.split("|"):
        if s in universe and s not in sym2id: sym2id[s]=gid
id2sym={v:k for k,v in sym2id.items()}
print("mapped to GeneID:",len(sym2id),"/",len(universe))

# 2) gene2go (human, Component) -> secretory-pathway CC flag
SECRETORY=re.compile(r"extracellular|plasma membrane|cell surface|endoplasmic reticulum|golgi|lysosom|endosom|secretory|external side|cell-surface",re.I)
print("streaming gene2go.gz (filtering human Component) ...")
ncbi={}  # symbol -> set(GO terms matching secretory CC)
n=0
for line in io.TextIOWrapper(stream("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"),encoding="utf-8"):
    if line.startswith("#"): continue
    f=line.rstrip("\n").split("\t")
    if f[0]!="9606": continue
    gid=f[1]
    if gid not in id2sym: continue
    if len(f)<8 or f[7]!="Component": continue
    term=f[5]
    if SECRETORY.search(term):
        ncbi.setdefault(id2sym[gid],set()).add(term)
    n+=1
print("human Component rows for our genes:",n,"| genes with secretory CC:",len(ncbi))

out=BASE+"/output/ncbi_routing.tsv"
with open(out,'w',newline='') as fo:
    w=csv.writer(fo,delimiter='\t',lineterminator='\n'); w.writerow(["gene","ncbi_secretory_CC","ncbi_CC_terms"])
    for g in sorted(universe):
        terms=ncbi.get(g)
        w.writerow([g, "yes" if terms else "no", ";".join(sorted(terms))[:200] if terms else ""])
print("wrote",out)
