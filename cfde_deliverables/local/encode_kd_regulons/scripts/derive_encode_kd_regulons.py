#!/usr/bin/env python3
# #2: ENCODE shRNA/CRISPR knockdown + RNA-seq -> FUNCTIONAL regulons. For each knocked-down factor,
# average expression over KD replicate files vs its matched control experiments, compute log2FC, and call
# regulated genes (|log2FC|>=LFC and max-expr>=EXPR). Splits up/down. Complements binding-based eCLIP/TF.
# PORTABLE stdlib (py3.9+). Lean: download each quant tsv -> parse -> DELETE. Resume-safe per factor.
import gzip, json, os, csv, math, collections, urllib.request, urllib.parse, socket
socket.setdefaulttimeout(120)   # so urlretrieve/urlopen can't hang forever on a stalled socket (was hanging at 0% CPU)
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.environ.get("OUTDIR",os.path.join(HOME,"encode_kd_regulons","output")); os.makedirs(OUT,exist_ok=True)
TMP=os.environ.get("TMPDIR","/Users/gage/.claude/jobs/32851e29/tmp/kd"); os.makedirs(TMP,exist_ok=True)
LFC=float(os.environ.get("LFC","1.0")); EXPR=float(os.environ.get("EXPR","1.0"))
LOG=open(os.path.join(OUT,"batch_log.txt"),"a")
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

def J(url):
    return json.load(urllib.request.urlopen(urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"x"}),timeout=600))
def exp_quant_files(exp_id):
    # return list of hrefs for "gene quantifications" tsv GRCh38 of an experiment
    d=J("https://www.encodeproject.org"+exp_id+"?format=json")
    out=[]
    for f in d.get("files",[]):
        if not isinstance(f,dict): continue
        if f.get("output_type","").startswith("gene quant") and f.get("file_format")=="tsv" and f.get("assembly")=="GRCh38" and f.get("status")=="released":
            out.append((f.get("accession"),f.get("href")))
    return out,d.get("possible_controls",[])

def parse_expr(fp):
    # return dict ENSG(noversion)->expression. Prefer TPM; else FPKM; else STAR counts->CPM.
    rows=[]; hdr=None
    with open(fp,encoding="utf-8",errors="replace") as fh:
        rd=csv.reader(fh,delimiter='\t')
        first=next(rd,None)
        if first and ("gene_id" in first or "TPM" in first or "FPKM" in first):
            hdr=first
            for r in rd: rows.append(r)
        else:
            if first: rows.append(first)
            for r in rd: rows.append(r)
    vals={}
    if hdr:
        gi=hdr.index("gene_id") if "gene_id" in hdr else 0
        col=hdr.index("TPM") if "TPM" in hdr else (hdr.index("FPKM") if "FPKM" in hdr else None)
        if col is None: return {}
        for r in rows:
            if len(r)<=col or len(r)<=gi or not r[gi].startswith("ENSG"): continue
            try: vals[r[gi].split(".")[0]]=float(r[col])
            except: pass
    else:  # STAR ReadsPerGene
        data=[r for r in rows if r and r[0].startswith("ENSG") and len(r)>=4]
        if not data: return {}
        tot=[0,0,0]
        for r in data:
            for j in (1,2,3):
                try: tot[j-1]+=int(float(r[j]))
                except: pass
        c=[1,2,3][tot.index(max(tot))]; lib=max(tot) or 1
        for r in data:
            try: vals[r[0].split(".")[0]]=float(r[c])/lib*1e6
            except: pass
    return vals

def mean_expr(hrefs):
    acc=collections.defaultdict(float); n=0
    for a,h in hrefs:
        fp=os.path.join(TMP,a+".tsv")
        try:
            urllib.request.urlretrieve("https://www.encodeproject.org"+h,fp)
            v=parse_expr(fp); os.remove(fp)
            if v:
                n+=1
                for g,x in v.items(): acc[g]+=x
        except Exception as e:
            log(f"  file FAIL {a}: {str(e)[:50]}")
            try: os.remove(fp)
            except OSError: pass
    if n==0: return {}
    return {g:s/n for g,s in acc.items()}

# collect KD experiments by target
by_target=collections.defaultdict(lambda:{"kd":[],"ctrl":set()})
ASSAYS=os.environ.get("ENC_KD_ASSAYS","shRNA RNA-seq").split("|")
PERT=os.environ.get("ENC_KD_PERT","shRNA knockdown")
LIB=os.environ.get("ENC_KD_LIB","ENCODE_shRNA_KD_regulon")
PREFIX=os.environ.get("ENC_KD_PREFIX","ENCODE_shRNA_KD")
for assay in ASSAYS:
    params=[("type","Experiment"),("assay_title",assay),("assembly","GRCh38"),("status","released"),
            ("limit","all"),("format","json"),("field","accession"),("field","target.label"),
            ("field","@id"),("field","possible_controls")]
    try: d=J("https://www.encodeproject.org/search/?"+urllib.parse.urlencode(params))
    except Exception as ex: log(f"assay query skip {assay}: {str(ex)[:50]}"); continue
    for e in d.get("@graph",[]):
        tgt=(e.get("target") or {}).get("label") if isinstance(e.get("target"),dict) else None
        if not tgt: continue
        by_target[tgt]["kd"].append(e["@id"])
        for c in e.get("possible_controls",[]):
            cid=c.get("@id") if isinstance(c,dict) else c
            if cid: by_target[tgt]["ctrl"].add(cid)
log(f"KD targets: {len(by_target)}")

def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(name,desc,genes,extra):
    d=os.path.join(OUT,"regulon",name); os.makedirs(d,exist_ok=True); genes=sorted(g for g in genes if g)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"Derived from ENCODE {PERT} RNA-seq in K562/HepG2 (KD vs matched control, |log2FC|>={LFC}, expr>={EXPR}), GRCh38, ENCODE/NHGRI, public; ENSG->symbol NCBI gene_info (NLM/NIH)."
    m={"standard_name":name,"library":LIB,"description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","derived_in_this_work":True,
       "evidence_type":"experimental_perturbation","biosample":"K562/HepG2 (immortalized cancer cell lines)","perturbation":PERT,
       "caveat":f"Cell-line {PERT} differential expression (direct + indirect/secondary effects). NOT tissue-level human biology. For direct targets, intersect with binding evidence (eCLIP/TF ChIP).",
       "note":"FC-based signature (no formal replicate statistics)","source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"encode_kd_regulon","evidence_type":"experimental_perturbation","biosample":"K562/HepG2","inputs":[f"ENCODE {PERT} RNA-seq gene quantifications (NHGRI; public)","NCBI gene_info (NLM/NIH)"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"},open(d+"/geneset.provenance.json","w"),indent=1)

done=skip=0
for tgt,info in sorted(by_target.items()):
    base=f"{PREFIX}_{safe(tgt)}"
    if os.path.exists(os.path.join(OUT,"regulon",base+"_regulated","geneset.tsv")): done+=1; continue
    if not info["ctrl"]: skip+=1; log(f"SKIP {tgt}: no control"); continue
    kd_h=[]; [kd_h.extend(exp_quant_files(e)[0]) for e in info["kd"]]
    ct_h=[]; [ct_h.extend(exp_quant_files(c)[0]) for c in info["ctrl"]]
    if not kd_h or not ct_h: skip+=1; log(f"SKIP {tgt}: no quant files"); continue
    kd=mean_expr(kd_h); ct=mean_expr(ct_h)
    if not kd or not ct: skip+=1; log(f"SKIP {tgt}: empty expr"); continue
    up=set(); dn=set()
    for g in set(kd)|set(ct):
        a=kd.get(g,0.0); b=ct.get(g,0.0)
        if max(a,b)<EXPR: continue
        lfc=math.log2((a+1)/(b+1))
        if lfc>=LFC: up.add(ensg2sym.get(g))
        elif lfc<=-LFC: dn.add(ensg2sym.get(g))
    up.discard(None); dn.discard(None); reg=up|dn
    emit(base+"_regulated",f"Genes differentially expressed on {tgt} knockdown (ENCODE)",reg,{"target":tgt,"direction":"both"})
    emit(base+"_up",f"Genes UP on {tgt} knockdown (de-repressed; ENCODE)",up,{"target":tgt,"direction":"up"})
    emit(base+"_down",f"Genes DOWN on {tgt} knockdown ({tgt}-dependent; ENCODE)",dn,{"target":tgt,"direction":"down"})
    done+=1; log(f"[{done}] {tgt}: regulated {len(reg)} (up {len(up)} / down {len(dn)})")
log(f"=== DONE KD: {done} factors, {skip} skipped ===")
print(f"ENCODE KD regulons: {done} factors, {skip} skipped")
LOG.close()
