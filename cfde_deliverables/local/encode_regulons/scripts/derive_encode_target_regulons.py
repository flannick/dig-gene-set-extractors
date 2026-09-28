#!/usr/bin/env python3
# A/B: ENCODE per-TARGET regulons. For each target factor (TF for TF ChIP-seq; RBP for eCLIP),
# union the genes whose promoter (TSS+/-1kb) overlaps that factor's peaks across all its experiments
# -> one "regulon"/target-gene set per factor, + x GTEx tissue-enriched. PORTABLE pure stdlib (py3.9+).
# Lean: download each peak bed -> derive -> DELETE. Resume-safe (skips a factor already written).
# Env: ENC_ASSAY ("TF ChIP-seq" | "eCLIP"), ENC_LIB, OUTDIR, TMPDIR, REFGENE, GTEX_TSTAT, FRAC, TSTAT_THR.
import gzip, json, os, collections, urllib.request, urllib.parse
ASSAY=os.environ.get("ENC_ASSAY","TF ChIP-seq")
LIB=os.environ.get("ENC_LIB","ENCODE_TF_regulon")
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.environ.get("OUTDIR",os.path.join(HOME,"encode_regulons","output_tf"))
TMP=os.environ.get("TMPDIR","/Users/gage/.claude/jobs/32851e29/tmp")
GTEX=os.environ.get("GTEX_TSTAT","/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv")
THR=float(os.environ.get("TSTAT_THR","4")); os.makedirs(OUT,exist_ok=True); os.makedirs(TMP,exist_ok=True)
REFGENE=os.environ.get("REFGENE",os.path.join(TMP,"refGene_hg38.txt.gz"))
WIN=BIN=1000; LOG=open(os.path.join(OUT,"batch_log.txt"),"a")
def log(m): LOG.write(m+"\n"); LOG.flush()
PEAK_TYPES={"IDR thresholded peaks","conservative IDR thresholded peaks","optimal IDR thresholded peaks",
            "pseudoreplicated IDR thresholded peaks","peaks"}

if not os.path.exists(REFGENE):
    log("downloading UCSC refGene hg38 ..."); urllib.request.urlretrieve(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz", REFGENE)
gene_reg=collections.defaultdict(set)
for line in gzip.open(REFGENE,'rt'):
    f=line.rstrip('\n').split('\t')
    if len(f)<13: continue
    ch,st,a,b,g=f[2],f[3],f[4],f[5],f[12]
    if '_' in ch or not a.isdigit() or not b.isdigit(): continue
    tss=int(a) if st=='+' else int(b); gene_reg[g].add((ch,tss-WIN,tss+WIN))
log(f"refGene genes: {len(gene_reg)}")

params=[("type","File"),("assay_title",ASSAY),("file_format","bed"),("assembly","GRCh38"),
        ("status","released"),("limit","all"),("format","json"),
        ("field","accession"),("field","href"),("field","output_type"),("field","target")]
d=json.load(urllib.request.urlopen(urllib.request.Request("https://www.encodeproject.org/search/?"+urllib.parse.urlencode(params),
   headers={"Accept":"application/json","User-Agent":"x"}),timeout=600))
by_target=collections.defaultdict(list)
for f in d.get("@graph",[]):
    if f.get("output_type") not in PEAK_TYPES: continue
    tl=f.get("target")
    if isinstance(tl,dict): tl=tl.get("label") or (tl.get("genes") or [{}])[0].get("symbol")
    elif isinstance(tl,str): tl=tl.strip("/").split("/")[-1].replace("-human","")
    href=f.get("href")
    if tl and href: by_target[tl].append((f.get("accession"),href))
log(f"{ASSAY} targets: {len(by_target)}")

def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def peak_genes(p):
    cov=collections.defaultdict(set)
    for line in gzip.open(p,'rt'):
        f=line.rstrip('\n').split('\t')
        if len(f)<3 or not f[1].isdigit(): continue
        for bb in range(int(f[1])//BIN,int(f[2])//BIN+1): cov[f[0]].add(bb)
    out=set()
    for g,regs in gene_reg.items():
        for c,s,e in regs:
            cb=cov.get(c)
            if cb and any(bb in cb for bb in range(s//BIN,e//BIN+1)): out.add(g); break
    return out

def emit(d,name,desc,genes,extra,gx=False):
    os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"Derived from ENCODE {ASSAY} peaks (promoter TSS+/-{WIN}bp overlap, unioned across experiments), GRCh38, ENCODE/NHGRI, public; UCSC refGene hg38"+(f"; intersected with GTEx tissue-enrichment (t>={THR}; NIH Common Fund)" if gx else "")+"."
    m={"standard_name":name,"library":LIB+("_x_GTEx" if gx else ""),"description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","assembly":"GRCh38","assay":ASSAY,"derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"encode_target_regulon"+("_x_gtex" if gx else ""),"inputs":[f"ENCODE {ASSAY} peaks (NHGRI; public)","UCSC refGene hg38"]+(["GTEx.tstat.hgnc.tsv (NIH Common Fund)"] if gx else []),"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"+(" + NIH Common Fund (GTEx)" if gx else "")},open(d+"/geneset.provenance.json","w"),indent=1)

# GTEx map
rows=[l.rstrip('\n').split('\t') for l in open(GTEX)]; tissues=rows[0][1:]
tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}

done=fail=nx=0
for tgt,files in sorted(by_target.items()):
    name=f"{LIB}_{safe(tgt)}"; dd=os.path.join(OUT,"regulon",name)
    if os.path.exists(os.path.join(dd,"geneset.tsv")):
        done+=1
    else:
        genes=set(); ok=False
        for acc,href in files:
            pk=os.path.join(TMP,acc+".bed.gz")
            try:
                urllib.request.urlretrieve("https://www.encodeproject.org"+href,pk)
                genes|=peak_genes(pk); os.remove(pk); ok=True
            except Exception as e:
                log(f"  peak FAIL {tgt} {acc}: {str(e)[:60]}")
                try: os.remove(pk)
                except OSError: pass
        if not ok: fail+=1; log(f"FAIL {tgt}: no peaks"); continue
        emit(dd,name,f"Genes with promoter-proximal ENCODE {ASSAY} binding for {tgt} ({len(files)} expt)",genes,{"target":tgt,"n_experiments":len(files)})
        done+=1; log(f"[{done}] {tgt}: {len(genes)} target genes ({len(files)} expt)")
    # x GTEx
    gset={r["gene"] for r in __import__("csv").DictReader(open(os.path.join(dd,"geneset.tsv")),delimiter='\t') if r.get("gene")}
    for ti,t in enumerate(tissues):
        e=sorted(g for g in gset if g in tmap and tmap[g][ti]>=THR)
        if not e: continue
        sn=f"{LIB}_{safe(tgt)}_x_GTEx_enriched_{safe(t)}"
        gd=os.path.join(OUT,"x_GTEx",sn)
        if os.path.exists(os.path.join(gd,"geneset.tsv")): continue
        emit(gd,sn,f"{tgt} ENCODE {ASSAY} target genes GTEx-enriched (t>={THR}) in {t}",e,{"target":tgt,"tissue":t},gx=True); nx+=1
log(f"=== DONE {ASSAY}: {done} regulons, {fail} failed, {nx} x_GTEx sets ===")
LOG.close()
print(f"{ASSAY}: {done} regulons, {fail} failed, {nx} x_GTEx sets")
