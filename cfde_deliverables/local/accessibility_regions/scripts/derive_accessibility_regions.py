#!/usr/bin/env python3
# REGION-LEVEL background-corrected accessibility (advisor design: region ± vs background -> map to genes).
# Re-downloads ENCODE peak BEDs and KEEPS them (regions preserved for reuse in other projects).
# Pass1: download+keep each experiment's peaks, bin to fixed-width, accumulate per-bin prevalence.
# Pass2: per experiment, call region UP (peak here & prevalence<LOW) / DOWN (no peak & prevalence>HIGH),
#        write the ± region BEDs, and map ± regions -> genes (promoter TSS+/-1kb) -> Up/Down gene sets.
# PORTABLE stdlib. Env: ENC_ASSAY, ENC_OUTPUT, ENC_TARGET, LIB, OUTDIR, REGIONDIR, REFGENE, TMPDIR, BIN, LOW, HIGH
import os, gzip, json, collections, urllib.request, urllib.parse
ASSAY=os.environ.get("ENC_ASSAY","ATAC-seq"); OUTPUT=os.environ.get("ENC_OUTPUT","IDR thresholded peaks")
LIB=os.environ.get("LIB","ENCODE_ATAC_region"); TMP=os.environ.get("TMPDIR","/Users/gage/.claude/jobs/32851e29/tmp/regions")
OUT=os.environ.get("OUTDIR",os.path.expanduser("~/Claude/proj-valiation-challenge/accessibility_regions/output")); OUT=os.path.join(OUT,LIB)
REGIONDIR=os.environ.get("REGIONDIR",os.path.expanduser("~/Claude/proj-valiation-challenge/accessibility_regions/regions"))+"/"+LIB
BIN=int(os.environ.get("BIN","1000")); LOW=float(os.environ.get("LOW","0.25")); HIGH=float(os.environ.get("HIGH","0.75")); WIN=1000
for d in (OUT,REGIONDIR,TMP): os.makedirs(d,exist_ok=True)
REFGENE=os.environ.get("REFGENE",os.path.join(TMP,"refGene_hg38.txt.gz"))
LOG=open(os.path.join(OUT,"batch_log.txt"),"a")
def log(m): LOG.write(m+"\n"); LOG.flush(); print(m,flush=True)
if not os.path.exists(REFGENE):
    urllib.request.urlretrieve("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz",REFGENE)
# promoter bins per gene
promo=collections.defaultdict(set)
for line in gzip.open(REFGENE,'rt'):
    f=line.rstrip('\n').split('\t')
    if len(f)<13 or '_' in f[2] or not f[4].isdigit() or not f[5].isdigit(): continue
    st=f[3]; tss=int(f[4]) if st=='+' else int(f[5])   # refGene: f[2]=chrom f[3]=strand f[4]=txStart f[5]=txEnd f[12]=symbol
    for b in range((tss-WIN)//BIN,(tss+WIN)//BIN+1): promo[f[12]].add((f[2],b))
log(f"refGene genes: {len(promo)}")
# ENCODE peak file list
params=[("type","File"),("assay_title",ASSAY),("file_format","bed"),("output_type",OUTPUT),
        ("assembly","GRCh38"),("status","released"),("limit","all"),("format","json"),
        ("field","accession"),("field","href"),("field","biosample_ontology")]
if os.environ.get("ENC_TARGET"): params.append(("target.label",os.environ["ENC_TARGET"]))
d=json.load(urllib.request.urlopen(urllib.request.Request("https://www.encodeproject.org/search/?"+urllib.parse.urlencode(params),headers={"Accept":"application/json","User-Agent":"x"}),timeout=300))
seen={}
for f in d.get("@graph",[]):
    bo=f.get("biosample_ontology",{}); bs=bo.get("term_name") if isinstance(bo,dict) else None
    if bs and bs not in seen and f.get("href"): seen[bs]=(f["accession"],f["href"])
log(f"{ASSAY}: {len(seen)} experiments")
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def bins_of(path):
    covered=set()
    for line in gzip.open(path,'rt'):
        f=line.rstrip('\n').split('\t')
        if len(f)<3 or not f[1].isdigit(): continue
        for b in range(int(f[1])//BIN,int(f[2])//BIN+1): covered.add((f[0],b))
    return covered
# PASS 1: download+keep, bin, accumulate prevalence
expbins={}; prev=collections.Counter()
import socket; socket.setdefaulttimeout(120)
for bs,(acc,href) in seen.items():
    bed=os.path.join(REGIONDIR,acc+".bed.gz")
    if not os.path.exists(bed):
        try:                                                # ATOMIC: .part then replace -> no partial files
            urllib.request.urlretrieve("https://www.encodeproject.org"+href,bed+".part")
            os.replace(bed+".part",bed)
        except Exception as e:
            try: os.remove(bed+".part")
            except OSError: pass
            log(f"  dl FAIL {bs} {acc}: {str(e)[:60]}"); continue
    try: cb=bins_of(bed)
    except Exception as e: log(f"  parse FAIL {acc}: {str(e)[:50]}"); continue
    expbins[bs]=(acc,cb)
    for b in cb: prev[b]+=1
N=len(expbins); log(f"kept {N} peak BEDs in {REGIONDIR}; binned")
# PASS 2: region ± -> genes
def emit(name,desc,genes,extra):
    if not genes: return 0
    dd=os.path.join(OUT,name); os.makedirs(dd,exist_ok=True); genes=sorted(genes)
    open(dd+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(dd+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=(f"REGION-LEVEL background-corrected ENCODE {ASSAY}: {BIN}bp regions called +/- vs cross-experiment "
          f"prevalence (N={N}, leave-one-out; Up<{LOW}, Down>{HIGH}), regions mapped to genes by promoter "
          f"overlap (TSS+/-{WIN}bp). Region BEDs kept. GRCh38, ENCODE/NHGRI public; UCSC refGene.")
    m={"standard_name":name,"library":LIB,"description":desc,"version":"1.0","file_type":"geneset","n_genes":len(genes),
       "organism":"human","assembly":"GRCh38","assay":ASSAY,"method":"region_level_prevalence_background",
       "bin_bp":BIN,"low":LOW,"high":HIGH,"background_n":N,"derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(dd+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"accessibility_region_contrast","inputs":[f"ENCODE {ASSAY} peaks (kept BEDs; NHGRI public)","cross-experiment region prevalence background","UCSC refGene promoters"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"},open(dd+"/geneset.provenance.json","w"),indent=1)
    return 1
allbins=set(prev); nu=nd=0
# precompute high-prevalence bins ONCE (DOWN candidates); per-exp DOWN is then a fast set-difference.
# (a DOWN bin is not covered by the exp, so leave-one-out prevalence = prev[b]/(N-1) for all such exps.)
highbins={b for b in allbins if N>1 and prev[b]/(N-1) > HIGH}
log(f"high-prevalence (DOWN-candidate) bins: {len(highbins)}")
for bs,(acc,cb) in expbins.items():
    up_bins={b for b in cb if (prev[b]-1)/(N-1) < LOW} if N>1 else set()
    dn_bins=highbins - cb
    # write region BEDs (kept)
    for tag,bset in (("Up",up_bins),("Down",dn_bins)):
        with open(os.path.join(REGIONDIR,f"{acc}_{tag}.regions.bed"),"w") as fh:
            for c,b in sorted(bset): fh.write(f"{c}\t{b*BIN}\t{(b+1)*BIN}\n")
    up_g={g for g,pb in promo.items() if pb & up_bins}
    dn_g={g for g,pb in promo.items() if pb & dn_bins}
    nu+=emit(f"{LIB}_{safe(bs)}_accessible_Up",f"Genes at regions specifically accessible ({ASSAY}) in {bs} vs background",up_g,{"biosample":bs,"direction":"Up"})
    nd+=emit(f"{LIB}_{safe(bs)}_inaccessible_Down",f"Genes at regions specifically INaccessible ({ASSAY}) in {bs} vs background",dn_g,{"biosample":bs,"direction":"Down"})
log(f"=== DONE {LIB}: {nu} Up + {nd} Down gene sets; region BEDs kept in {REGIONDIR} ===")
