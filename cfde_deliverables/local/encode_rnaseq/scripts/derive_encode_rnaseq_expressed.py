#!/usr/bin/env python3
# #3: per-biosample "expressed genes" from ENCODE RNA-seq gene quantifications (TPM>=THR).
# ENSG->symbol via NCBI gene_info (NLM/NIH). Lean: download tsv -> parse -> DELETE. Resume-safe.
import gzip, json, os, csv, io, urllib.request, urllib.parse, collections
OUT=os.environ.get("OUTDIR",os.path.expanduser("~/Claude/proj-valiation-challenge/encode_rnaseq/output"))
TMP=os.environ.get("TMPDIR","/Users/gage/.claude/jobs/32851e29/tmp"); os.makedirs(OUT,exist_ok=True)
THR=float(os.environ.get("TPM_THR","1")); LOG=open(os.path.join(OUT,"batch_log.txt"),"a")
def log(m): LOG.write(m+"\n"); LOG.flush()

# ENSG -> symbol from NCBI gene_info
GI=os.path.join(TMP,"Homo_sapiens.gene_info.gz")
if not os.path.exists(GI):
    urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",GI)
ensg2sym={}
for line in gzip.open(GI,'rt'):
    if line.startswith("#"): continue
    f=line.rstrip("\n").split("\t"); sym=f[2]
    for x in f[5].split("|"):
        if x.startswith("Ensembl:"): ensg2sym[x.split(":")[1]]=sym
log(f"ENSG->symbol map: {len(ensg2sym)}")

params=[("type","File"),("assay_title","total RNA-seq"),("assay_title","polyA plus RNA-seq"),
        ("output_type","gene quantifications"),("file_format","tsv"),("assembly","GRCh38"),
        ("status","released"),("limit","all"),("format","json")]
d=json.load(urllib.request.urlopen(urllib.request.Request("https://www.encodeproject.org/search/?"+urllib.parse.urlencode(params),
   headers={"Accept":"application/json","User-Agent":"x"}),timeout=300))
seen={}
for f in d.get("@graph",[]):
    bo=f.get("biosample_ontology",{}); bs=bo.get("term_name") if isinstance(bo,dict) else None
    if bs and bs not in seen: seen[bs]=(f.get("accession"),f.get("dataset","").strip("/").split("/")[-1],f.get("href"))
log(f"distinct biosamples: {len(seen)}")
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]

done=fail=0
for bs,(acc,exp,href) in seen.items():
    name=f"ENCODE_RNAseq_expressed_{safe(bs)}"; dd=os.path.join(OUT,name)
    if os.path.exists(os.path.join(dd,"geneset.tsv")): done+=1; continue
    if not href: continue
    fp=os.path.join(TMP,acc+".tsv")
    try:
        urllib.request.urlretrieve("https://www.encodeproject.org"+href,fp)
        genes=set()
        with open(fp,encoding="utf-8",errors="replace") as fh:
            rd=csv.reader(fh,delimiter='\t'); hdr=next(rd)
            gi=hdr.index("gene_id") if "gene_id" in hdr else 0
            ti=hdr.index("TPM") if "TPM" in hdr else None
            if ti is None: raise ValueError("no TPM column")
            for row in rd:
                if len(row)<=ti: continue
                try: tpm=float(row[ti])
                except: continue
                if tpm>=THR:
                    g=ensg2sym.get(row[gi].split(".")[0])
                    if g: genes.add(g)
        os.remove(fp); genes=sorted(genes); os.makedirs(dd,exist_ok=True)
        cite=f"Derived from ENCODE RNA-seq gene quantifications {acc} (exp {exp}, biosample {bs}; TPM>={THR}), GRCh38, ENCODE/NHGRI, public; ENSG->symbol via NCBI gene_info (NLM/NIH)."
        open(dd+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
        open(dd+"/genesets.gmt","w").write(f"{name}\tGenes expressed (TPM>={THR}) in {bs} by ENCODE RNA-seq\t"+"\t".join(genes)+"\n")
        json.dump({"standard_name":name,"library":"ENCODE_RNAseq_expressed","description":f"Genes expressed (TPM>={THR}) in {bs} (DERIVED from ENCODE RNA-seq).","version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","assembly":"GRCh38","biosample":bs,"tpm_threshold":THR,"derived_in_this_work":True,"source":cite},open(dd+"/geneset.meta.json","w"),indent=1)
        json.dump({"focus":name,"operation":"derive_expressed_genes_from_rnaseq","inputs":[f"ENCODE {acc} gene quantifications (NHGRI; public)","NCBI gene_info ENSG->symbol (NLM/NIH; public)"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"},open(dd+"/geneset.provenance.json","w"),indent=1)
        done+=1; log(f"[{done}] {bs}: {len(genes)} expressed genes")
    except Exception as e:
        fail+=1; log(f"FAIL {bs} {acc}: {str(e)[:80]}")
        try: os.remove(fp)
        except OSError: pass
log(f"=== DONE: {done} sets, {fail} failures ===")
LOG.close()
