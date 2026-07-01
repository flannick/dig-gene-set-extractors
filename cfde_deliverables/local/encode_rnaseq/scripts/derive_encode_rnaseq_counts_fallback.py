#!/usr/bin/env python3
# Fallback for biosamples whose ENCODE "gene quantifications" file is a STAR ReadsPerGene.out.tab
# (raw counts, no TPM/FPKM, no length). Use counts -> CPM>=THR as the "expressed" call.
# STAR cols: gene_id, unstranded, strand1, strand2 (after N_* header rows). Pick correct strand
# = column with max total assigned reads. Reprocesses ONLY the failed biosamples from batch_log.txt.
import gzip, json, os, csv, re, urllib.request
OUT=os.path.expanduser("~/Claude/proj-valiation-challenge/encode_rnaseq/output")
TMP="/Users/gage/.claude/jobs/32851e29/tmp"
GI=os.path.join(TMP,"Homo_sapiens.gene_info.gz"); CPM_THR=float(os.environ.get("CPM_THR","1"))

ensg2sym={}
for line in gzip.open(GI,'rt'):
    if line.startswith("#"): continue
    f=line.rstrip("\n").split("\t")
    for x in f[5].split("|"):
        if x.startswith("Ensembl:"): ensg2sym[x.split(":")[1]]=f[2]

# parse failed (biosample, accession) pairs from the run log
fails=[]
for line in open(os.path.join(OUT,"batch_log.txt")):
    m=re.match(r"FAIL (.+) (ENCFF\w+): no TPM column",line.strip())
    if m: fails.append((m.group(1),m.group(2)))
print("failed biosamples to retry:",len(fails))
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]

done=fail=0
for bs,acc in fails:
    name=f"ENCODE_RNAseq_expressed_{safe(bs)}"; dd=os.path.join(OUT,name)
    if os.path.exists(os.path.join(dd,"geneset.tsv")): done+=1; continue
    fp=os.path.join(TMP,acc+".tsv")
    try:
        urllib.request.urlretrieve(f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.tsv",fp)
        rows=[]
        with open(fp,encoding="utf-8",errors="replace") as fh:
            for r in csv.reader(fh,delimiter='\t'):
                if not r or r[0].startswith("N_") or r[0]=="gene_id": continue
                if len(r)>=4 and r[0].startswith("ENSG"): rows.append(r)
        os.remove(fp)
        if not rows: raise ValueError("no ENSG count rows")
        # choose strand column (1,2,3) by max total assigned reads
        tot=[0,0,0]
        for r in rows:
            for j in (1,2,3):
                try: tot[j-1]+=int(float(r[j]))
                except: pass
        col=[1,2,3][tot.index(max(tot))]; libsum=max(tot)
        if libsum<=0: raise ValueError("zero library counts")
        genes=set()
        for r in rows:
            try: c=float(r[col])
            except: continue
            if c/libsum*1e6>=CPM_THR:
                g=ensg2sym.get(r[0].split(".")[0])
                if g: genes.add(g)
        genes=sorted(genes); os.makedirs(dd,exist_ok=True)
        cite=(f"Derived from ENCODE RNA-seq STAR gene read counts {acc} (biosample {bs}; CPM>={CPM_THR} "
              f"on strand-col {col}; no TPM/length in source so CPM used in lieu of TPM), GRCh38, ENCODE/NHGRI, "
              f"public; ENSG->symbol via NCBI gene_info (NLM/NIH).")
        open(dd+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
        open(dd+"/genesets.gmt","w").write(f"{name}\tGenes expressed (CPM>={CPM_THR}) in {bs} by ENCODE RNA-seq (STAR counts)\t"+"\t".join(genes)+"\n")
        json.dump({"standard_name":name,"library":"ENCODE_RNAseq_expressed","description":f"Genes expressed (CPM>={CPM_THR}) in {bs} (DERIVED from ENCODE RNA-seq STAR read counts; CPM used in lieu of TPM).","version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","assembly":"GRCh38","biosample":bs,"quant_method":"STAR_counts_CPM","cpm_threshold":CPM_THR,"derived_in_this_work":True,"source":cite},open(dd+"/geneset.meta.json","w"),indent=1)
        json.dump({"focus":name,"operation":"derive_expressed_genes_from_rnaseq_counts","inputs":[f"ENCODE {acc} STAR gene read counts (NHGRI; public)","NCBI gene_info ENSG->symbol (NLM/NIH; public)"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"},open(dd+"/geneset.provenance.json","w"),indent=1)
        done+=1; print(f"[{done}] {bs}: {len(genes)} expressed (CPM>={CPM_THR}, strand-col {col})")
    except Exception as e:
        fail+=1; print(f"FAIL {bs} {acc}: {str(e)[:90]}")
        try: os.remove(fp)
        except OSError: pass
print(f"=== fallback DONE: {done} recovered, {fail} failures ===")
