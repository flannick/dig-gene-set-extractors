#!/usr/bin/env python3
# Same background-prevalence correction as ATAC/DNase, applied to TF ChIP / eCLIP regulons.
# Background = per-gene BINDING prevalence across all factors of the assay (how many TFs/RBPs target it);
# a gene bound by many factors is a promiscuous/HOT-region artifact (binding analog of a ubiquitously-open
# promoter). Leave-one-out. Per factor:
#   specific_Up   = bound by THIS factor AND prevalence(excl self) < LOW   (specifically-bound targets)
#   promiscuous_Down = NOT bound here AND prevalence > HIGH                (commonly-bound elsewhere, not by this factor)
import os, json, glob, csv, collections
REGDIR=os.environ["REGULON_DIR"]                 # dir of <name>/geneset.tsv (per-factor regulons)
LIB=os.environ.get("LIB","ENCODE_TF_regulon_bgcorrected"); ASSAY=os.environ.get("ASSAY","TF ChIP-seq")
OUT=os.environ.get("OUTDIR",os.path.expanduser("~/Claude/proj-valiation-challenge/encode_regulons_bgcontrast/output")); OUT=os.path.join(OUT,LIB)
LOW=float(os.environ.get("LOW","0.25")); HIGH=float(os.environ.get("HIGH","0.75")); os.makedirs(OUT,exist_ok=True)
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def load(fp): return {r["gene"] for r in csv.DictReader(open(fp),delimiter='\t') if r.get("gene")}
sets={}; lab={}
for d in glob.glob(os.path.join(REGDIR,"*")):
    gp=os.path.join(d,"geneset.tsv")
    if not os.path.exists(gp): continue
    g=load(gp)
    if not g: continue
    tgt=os.path.basename(d)
    try: tgt=json.load(open(os.path.join(d,"geneset.meta.json"))).get("target",tgt)
    except Exception: pass
    sets[os.path.basename(d)]=g; lab[os.path.basename(d)]=tgt
N=len(sets); count=collections.Counter()
for g in sets.values():
    for x in g: count[x]+=1
highprev={x for x in count if N>1 and count[x]/(N-1) > HIGH}   # promiscuous / HOT genes
print(f"{ASSAY}: {N} factors | genes={len(count)} | promiscuous(>{HIGH}) genes={len(highprev)}")
def emit(name,desc,genes,tgt,cls):
    if not genes: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=(f"Background-corrected ENCODE {ASSAY} regulon: per-factor targets vs cross-factor BINDING prevalence "
          f"(N={N} factors, leave-one-out; Up<{LOW}, Down>{HIGH}); removes promiscuous/HOT-region binding. "
          f"Relative binding specificity — same correction as ATAC/DNase. ENCODE/NHGRI public.")
    m={"standard_name":name,"library":LIB,"description":desc,"version":"1.0","file_type":"geneset","n_genes":len(genes),
       "organism":"human","assay":ASSAY,"target":tgt,"method":"binding_prevalence_background_leave_one_out",
       "background_n_factors":N,"low":LOW,"high":HIGH,"class":cls,
       "caveat":"Cross-factor binding-prevalence background (specificity); removes promiscuous/HOT-region genes.",
       "derived_in_this_work":True,"source":cite}
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"regulon_background_contrast","inputs":[f"ENCODE {ASSAY} per-factor regulons (NHGRI; public)","cross-factor binding prevalence background"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1
nu=nd=0
for k,g in sets.items():
    tgt=lab[k]
    up={x for x in g if (count[x]-1)/(N-1) < LOW} if N>1 else set()
    dn=highprev - g
    nu+=emit(f"{LIB}_{safe(tgt)}_specific_Up",f"{tgt}: specifically-bound targets vs background ({ASSAY})",up,tgt,"specific_Up")
    nd+=emit(f"{LIB}_{safe(tgt)}_promiscuous_Down",f"{tgt}: commonly-bound-elsewhere genes NOT bound by {tgt} ({ASSAY})",dn,tgt,"promiscuous_Down")
print(f"{LIB}: {nu} specific_Up + {nd} promiscuous_Down sets")
