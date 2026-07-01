#!/usr/bin/env python
"""Batch: one 'accessible genes' set per ENCODE human ATAC biosample (IDR peaks).
Lean: download peak -> derive -> DELETE peak. Resume-safe. Cites each ENCODE accession."""
import gzip, json, os, collections, urllib.request
BASE=os.path.expanduser("~/Claude/proj-valiation-challenge/accessible_genes")
OUT=BASE+"/output"; os.makedirs(OUT,exist_ok=True)
REFGENE=os.path.expanduser("~/.claude/jobs/32851e29/tmp/refGene_hg38.txt.gz")
TMP=os.path.expanduser("~/.claude/jobs/32851e29/tmp")
WIN=BIN=1000
LOG=open(OUT+"/batch_log.txt","a")
def log(m): LOG.write(m+"\n"); LOG.flush()

gene_tss=collections.defaultdict(set)
for line in gzip.open(REFGENE,'rt'):
    f=line.rstrip('\n').split('\t')
    if len(f)<13: continue
    chrom,strand,txS,txE,gene=f[2],f[3],f[4],f[5],f[12]
    if '_' in chrom or not txS.isdigit(): continue
    gene_tss[gene].add((chrom,int(txS) if strand=='+' else int(txE)))
log(f"refGene genes: {len(gene_tss)}")

url=("https://www.encodeproject.org/search/?type=File&assay_title=ATAC-seq&file_format=bed"
     "&output_type=IDR+thresholded+peaks&assembly=GRCh38&status=released&limit=all&format=json")
d=json.load(urllib.request.urlopen(urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"x"}),timeout=180))
seen={}
for f in d.get("@graph",[]):
    bo=f.get("biosample_ontology",{}); bs=bo.get("term_name") if isinstance(bo,dict) else None
    if bs and bs not in seen:
        seen[bs]=(f.get("accession"), f.get("dataset","").strip("/").split("/")[-1], f.get("href"))
log(f"distinct biosamples: {len(seen)}")

def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def accessible(p):
    cov=collections.defaultdict(set); n=0
    for line in gzip.open(p,'rt'):
        f=line.rstrip('\n').split('\t')
        if len(f)<3 or not f[1].isdigit(): continue
        for b in range(int(f[1])//BIN,int(f[2])//BIN+1): cov[f[0]].add(b)
        n+=1
    acc=[]
    for g,tsss in gene_tss.items():
        for c,t in tsss:
            cb=cov.get(c)
            if cb and any(b in cb for b in range((t-WIN)//BIN,(t+WIN)//BIN+1)): acc.append(g); break
    return sorted(set(acc)), n

done=fail=0
for bs,(acc,exp,href) in seen.items():
    name=f"ENCODE_ATAC_accessible_{safe(bs)}"; dd=os.path.join(OUT,name)
    if os.path.exists(dd+"/geneset.tsv"): done+=1; continue
    if not href: continue
    pk=os.path.join(TMP,acc+".bed.gz")
    try:
        urllib.request.urlretrieve("https://www.encodeproject.org"+href, pk)
        genes,npk=accessible(pk); os.remove(pk)
        os.makedirs(dd,exist_ok=True)
        cite=f"Derived from ENCODE ATAC peak {acc} (exp {exp}, biosample {bs}), GRCh38, ENCODE/NHGRI, public; UCSC refGene hg38."
        open(dd+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
        open(dd+"/genesets.gmt","w").write(f"{name}\tGenes promoter-proximal to ENCODE ATAC peaks in {bs}\t"+"\t".join(genes)+"\n")
        json.dump({"standard_name":name,"library":"ENCODE_ATAC_accessible","description":f"Genes with promoter-proximal ENCODE ATAC accessibility in {bs} (DERIVED from ENCODE public peaks).","version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","assembly":"GRCh38","biosample":bs,"derived_in_this_work":True,"source":cite},open(dd+"/geneset.meta.json","w"),indent=1)
        json.dump({"focus":name,"operation":"derive_accessible_genes_from_ATAC","inputs":[f"ENCODE {acc} peaks (NHGRI; public)","UCSC refGene hg38"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"},open(dd+"/geneset.provenance.json","w"),indent=1)
        done+=1; log(f"[{done}] {bs}: {len(genes)} genes / {npk} peaks")
    except Exception as e:
        fail+=1; log(f"FAIL {bs} {acc}: {str(e)[:80]}")
        try: os.remove(pk)
        except: pass
log(f"=== DONE: {done} sets, {fail} failures ===")
LOG.close()
