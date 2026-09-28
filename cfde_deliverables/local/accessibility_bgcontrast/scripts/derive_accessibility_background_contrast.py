#!/usr/bin/env python3
# Background-corrected chromatin-accessibility gene sets. For a collection of per-biosample accessible-gene
# lists (an existing submission zip), build an AGGREGATED BACKGROUND = per-gene accessibility PREVALENCE
# across all experiments, then per experiment emit matched Up/Down sets vs that background (LEAVE-ONE-OUT):
#   Up   = genes accessible in THIS experiment AND prevalence(excl self) < LOW  (specifically accessible)
#   Down = genes NOT accessible here AND prevalence(excl self) > HIGH           (specifically closed)
# This removes constitutively-open (housekeeping) promoters (high prevalence) -> de-biases the sets.
# HONEST SCOPE: this is a PEAK-CALL prevalence background (relative accessibility specificity), NOT a
# normalized read-count differential (no GC/library/DESeq2 — that needs raw reads). Thresholds documented.
import os, json, zipfile, collections
IN=os.environ["IN_ZIP"]                                   # existing per-biosample submission zip
LIB=os.environ.get("LIB","ENCODE_ATAC_accessible_bgcontrast")
ASSAY=os.environ.get("ASSAY","ATAC-seq")
LOW=float(os.environ.get("LOW","0.25")); HIGH=float(os.environ.get("HIGH","0.75"))
OUT=os.environ.get("OUTDIR",os.path.expanduser("~/Claude/proj-valiation-challenge/accessibility_bgcontrast/output"))
OUT=os.path.join(OUT,LIB); os.makedirs(OUT,exist_ok=True)

# GROUP_KEY: if set (e.g. "library"), the prevalence background is computed WITHIN each group (needed for
# histone marks — H3K4me3 must not be pooled with H3K27me3). LIB_FILTER: only process sets whose library
# contains this substring (so a mixed zip like batch3 only gets its histone sets contrasted).
GROUP_KEY=os.environ.get("GROUP_KEY"); LIB_FILTER=os.environ.get("LIB_FILTER")
# SYM_UNIVERSE: path to a valid-symbol reference (e.g. GTEx tstat); if set, genes are filtered to real
# HGNC symbols — needed for rE2G whose target column is polluted with ENSG/enhancer-element IDs.
VALID=None
if os.environ.get("SYM_UNIVERSE"):
    rows=[l.split("\t")[0].strip() for l in open(os.environ["SYM_UNIVERSE"])]
    VALID=set(rows[1:])  # skip header
    print(f"symbol universe: {len(VALID)} valid symbols")
z=zipfile.ZipFile(IN)
recs=[]  # (setdir, genes, label, group)
for n in z.namelist():
    if not n.endswith("geneset.tsv"): continue
    setdir=n.split("/")[-2]
    genes={ln.strip() for ln in z.read(n).decode("utf-8","replace").splitlines()[1:] if ln.strip()}
    if VALID is not None: genes={g for g in genes if g in VALID}
    if not genes: continue
    label=setdir; group="ALL"
    try:
        m=json.loads(z.read(n.rsplit("/",1)[0]+"/geneset.meta.json"))
        label=m.get("biosample",setdir)
        if LIB_FILTER and LIB_FILTER not in str(m.get("library","")): continue
        if GROUP_KEY: group=str(m.get(GROUP_KEY,"ALL"))
    except Exception:
        if LIB_FILTER: continue
    recs.append((setdir,genes,label,group))
# per-group background
groups=collections.defaultdict(list)
for r in recs: groups[r[3]].append(r)
print(f"{ASSAY}: {len(recs)} sets in {len(groups)} group(s): "+", ".join(f"{g}={len(v)}" for g,v in groups.items()))

def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
def emit(name,desc,genes,extra):
    if not genes: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=(f"Background-corrected ENCODE {ASSAY} accessibility: gene promoter accessible-vs-aggregated-background "
          f"(prevalence across {N} ENCODE experiments, leave-one-out; Up<{LOW}, Down>{HIGH}). "
          f"Relative accessibility SPECIFICITY at ENCODE peak-call resolution — NOT a normalized read-count "
          f"differential (no GC/library/DESeq2). GRCh38, ENCODE/NHGRI, public; UCSC refGene promoters.")
    m={"standard_name":name,"library":LIB,"description":desc,"version":"0.2","file_type":"geneset",
       "n_genes":len(genes),"organism":"human","assembly":"GRCh38","assay":ASSAY,
       "method":"prevalence_background_contrast_leave_one_out","background_n_experiments":N,
       "low_prevalence_threshold":LOW,"high_prevalence_threshold":HIGH,
       "caveat":"Peak-call prevalence background (relative specificity); not a GC/library-normalized read-count differential.",
       "derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"accessibility_background_contrast","inputs":[f"ENCODE {ASSAY} per-biosample accessible-gene sets (NHGRI; public)","aggregated cross-experiment prevalence background (leave-one-out)"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"},open(d+"/geneset.provenance.json","w"),indent=1)
    return 1

nu=nd=0
for grp,members in groups.items():
    N=len(members)                                          # per-group background (e.g. per histone mark)
    count=collections.Counter()
    for _,gs,_,_ in members:
        for g in gs: count[g]+=1
    universe=set(count)
    pre=LIB if grp=="ALL" else f"{LIB}_{safe(grp)}"
    for setdir,gs,lab,_ in members:
        up=set(); dn=set()
        for g in gs:
            prev=(count[g]-1)/(N-1) if N>1 else 0.0         # leave-one-out (g is in this sample)
            if prev<LOW: up.add(g)
        for g in universe:
            if g in gs: continue
            prev=count[g]/(N-1) if N>1 else 0.0             # g absent here; prevalence elsewhere
            if prev>HIGH: dn.add(g)
        nu+=emit(f"{pre}_{safe(lab)}_accessible_Up",f"Genes specifically accessible ({ASSAY}{'' if grp=='ALL' else ' '+grp}) in {lab} vs background (prevalence<{LOW})",up,{"biosample":lab,"group":grp,"direction":"Up"})
        nd+=emit(f"{pre}_{safe(lab)}_inaccessible_Down",f"Genes specifically INaccessible ({ASSAY}{'' if grp=='ALL' else ' '+grp}) in {lab} vs background (elsewhere-open>{HIGH})",dn,{"biosample":lab,"group":grp,"direction":"Down"})
    print(f"  {grp}: N={N} universe={len(universe)}")
print(f"{LIB}: {nu} Up + {nd} Down sets")
