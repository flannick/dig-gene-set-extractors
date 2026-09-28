#!/usr/bin/env python3
# E/G: GTEx per-tissue QTL genes. eQTL eGenes (cis-eQTL) or sQTL sGenes (splicing) with qval < QVAL.
# Lean: download the GTEx QTL tar -> stream-extract only the per-gene .txt.gz members -> DELETE tar.
# Env: QTL_TAR_URL, QTL_SUFFIX (".v8.egenes.txt.gz" | ".v8.sgenes.txt.gz"), QTL_LIB, QTL_KIND, OUTDIR, QVAL.
import os, json, gzip, io, csv, tarfile, urllib.request, collections
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
TMP=os.environ.get("TMPDIR","/Users/gage/.claude/jobs/32851e29/tmp")
URL=os.environ["QTL_TAR_URL"]; SUF=os.environ.get("QTL_SUFFIX",".v8.egenes.txt.gz")
LIB=os.environ.get("QTL_LIB","GTEx_eQTL_eGenes"); KIND=os.environ.get("QTL_KIND","cis-eQTL eGene")
OUT=os.environ.get("OUTDIR",os.path.join(HOME,"gtex_qtl","output_eqtl")); QVAL=float(os.environ.get("QVAL","0.05"))
os.makedirs(OUT,exist_ok=True); LOG=open(os.path.join(OUT,"batch_log.txt"),"a")
def log(m): LOG.write(m+"\n"); LOG.flush()
GI=os.path.join(TMP,"Homo_sapiens.gene_info.gz")
if not os.path.exists(GI):
    urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",GI)
ensg2sym={}
for line in gzip.open(GI,'rt'):
    if line.startswith("#"): continue
    f=line.rstrip("\n").split("\t")
    for x in f[5].split("|"):
        if x.startswith("Ensembl:"): ensg2sym[x.split(":")[1]]=f[2]
tar=os.path.join(TMP,os.path.basename(URL))
if not os.path.exists(tar):
    log(f"downloading {URL} ..."); urllib.request.urlretrieve(URL,tar)
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(d,name,desc,genes,t):
    os.makedirs(d,exist_ok=True); genes=sorted(g for g in genes if g)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"Derived from GTEx v8 {KIND} per-gene results ({t}; qval<{QVAL}), NIH Common Fund, public; gene_name from GTEx (ENSG->symbol via NCBI gene_info fallback)."
    json.dump({"standard_name":name,"library":LIB,"description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","tissue":t,"qtl":KIND,"qval_threshold":QVAL,"derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"gtex_qtl_genes","inputs":[f"GTEx v8 {KIND} (NIH Common Fund; public)","NCBI gene_info (NLM/NIH)"],"source_citation":cite,"public":True,"funding":"NIH Common Fund (GTEx)"},open(d+"/geneset.provenance.json","w"),indent=1)
n=0
with tarfile.open(tar) as tf:
    for m in tf.getmembers():
        if not m.name.endswith(SUF): continue
        tissue=os.path.basename(m.name)[:-len(SUF)]
        raw=tf.extractfile(m).read()
        rd=csv.reader(io.StringIO(gzip.decompress(raw).decode("utf-8","replace")),delimiter='\t')
        hdr=next(rd)
        qi=hdr.index("qval") if "qval" in hdr else None
        gni=hdr.index("gene_name") if "gene_name" in hdr else None
        gi=hdr.index("gene_id") if "gene_id" in hdr else (hdr.index("group_id") if "group_id" in hdr else 0)
        if qi is None: log(f"SKIP {tissue}: no qval"); continue
        genes=set()
        for r in rd:
            if len(r)<=qi: continue
            try: q=float(r[qi])
            except: continue
            if q<QVAL:
                g=r[gni] if (gni is not None and len(r)>gni and r[gni] not in("","NA")) else ensg2sym.get(r[gi].split(".")[0])
                if g: genes.add(g)
        if not genes: continue
        sn=f"{LIB}_{safe(tissue)}"
        emit(os.path.join(OUT,sn),sn,f"Genes with a significant GTEx {KIND} (qval<{QVAL}) in {tissue}",genes,tissue); n+=1
        log(f"[{n}] {tissue}: {len(genes)} genes")
try: os.remove(tar); log("deleted tar")
except OSError: pass
log(f"=== DONE {KIND}: {n} tissue sets ===")
print(f"{KIND}: {n} tissue sets")
