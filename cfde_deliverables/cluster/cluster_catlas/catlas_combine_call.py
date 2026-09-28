#!/usr/bin/env python3
# Combine per-sample CATLAS partials -> per-cell-type promoter accessibility -> gene sets.
# Signal(ct,gene) = sum over samples of promoter fragment count, normalized per cell:
#     sig = total_count / ncells(ct)   (mean promoter fragments/cell)
# accessible(ct,gene) = sig >= MINSIG AND total_count >= MINCOUNT
# Background-prevalence control across all cell types (same method as ENCODE ATAC/DNase):
#     prevalence(gene) = fraction of cell types where accessible
#     specific_Up(ct)  = accessible in ct AND prevalence(excl ct) < LOW   (cell-type-specific open)
#     closed_Down(ct)  = NOT accessible in ct AND prevalence > HIGH       (broadly-open but closed here)
# Also emits RAW (uncontrolled) accessible set per cell type, labeled 'raw'.
# Honest scope: relative cell-type accessibility SPECIFICITY (promoter fragment CPM/cell), NOT a
# GC/library-normalized read-count differential.
import os, glob, gzip, json, collections
HERE=os.path.dirname(os.path.abspath(__file__))
WORK=os.environ.get("WORK",HERE)
PARTDIR=os.environ.get("PARTDIR",os.path.join(WORK,"partials"))
OUT=os.environ.get("OUTDIR_SETS",os.path.join(WORK,"catlas_genesets_full")); os.makedirs(OUT,exist_ok=True)
META=os.environ.get("META",os.path.join(WORK,"GSE184462_metadata.tsv.gz"))
MINSIG=float(os.environ.get("MINSIG","0.15"))     # mean promoter fragments/cell to call accessible
MINCOUNT=int(os.environ.get("MINCOUNT","10"))     # min raw fragments (guards tiny cell types)
LOW=float(os.environ.get("LOW","0.25")); HIGH=float(os.environ.get("HIGH","0.75"))
MIN_CELLS=int(os.environ.get("MIN_CELLS","25"))   # drop cell types with too few cells
CIT=("Zhang K, Hocker JD, Miller M, et al. A single-cell atlas of chromatin accessibility in the human "
     "genome. Cell 2021;184(24):5985-6001 (GSE184462; CATlas). NIH-funded, public.")

def safe(s): return "".join(c if (c.isalnum() or c in "._-") else "_" for c in s)[:70]

# cells per cell type (global, from metadata)
ncells=collections.Counter()
with gzip.open(META,"rt") as m:
    next(m)
    for line in m:
        f=line.rstrip("\n").split("\t")
        if len(f)>=7: ncells[f[6]]+=1
print(f"cell types in metadata: {len(ncells)} | total cells: {sum(ncells.values())}",flush=True)

# sum partials
count=collections.defaultdict(int)   # (ct,gene)->count
parts=sorted(glob.glob(os.path.join(PARTDIR,"*.partial.tsv")))
print(f"partial files: {len(parts)}",flush=True)
for p in parts:
    with open(p) as fh:
        for line in fh:
            a=line.rstrip("\n").split("\t")
            if len(a)!=3: continue
            ct,g,c=a[0],a[1],a[2]
            try: count[(ct,g)]+=int(c)
            except ValueError: pass

# accessible calls per cell type
celltypes=[c for c in ncells if ncells[c]>=MIN_CELLS]
accessible={ct:set() for ct in celltypes}
for (ct,g),c in count.items():
    if ct not in accessible: continue
    n=ncells[ct]
    if n<=0: continue
    if c>=MINCOUNT and (c/n)>=MINSIG:
        accessible[ct].add(g)
for ct in celltypes:
    print(f"  {ct}: {len(accessible[ct])} accessible genes ({ncells[ct]} cells)",flush=True)

# prevalence across cell types
prev=collections.Counter()
for ct in celltypes:
    for g in accessible[ct]: prev[g]+=1
N=len(celltypes)
allgenes=set(prev)
highprev={g for g in allgenes if N>1 and prev[g]/(N-1)>HIGH}
print(f"N cell types used: {N} | genes accessible somewhere: {len(allgenes)} | broadly-open(>{HIGH}): {len(highprev)}",flush=True)

def emit(name,desc,genes,cls,ct,extra):
    if not genes: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    m={"standard_name":name,"library":"CATLAS_scATAC_accessibility","description":desc,"version":"1.0",
       "file_type":"geneset","n_genes":len(genes),"organism":"human","assay":"scATAC-seq (CATlas)",
       "cell_type":ct,"n_cells":ncells.get(ct,0),"class":cls,
       "method":"promoter fragment/cell; background = cross-cell-type accessibility prevalence (leave-one-out)",
       "min_sig_per_cell":MINSIG,"min_count":MINCOUNT,"low":LOW,"high":HIGH,
       "caveat":"Relative cell-type accessibility specificity (promoter fragment CPM/cell); NOT a GC/library-normalized read-count differential.",
       "derived_in_this_work":True,"source":CIT}
    m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"catlas_promoter_accessibility_bgcontrol","inputs":
               ["CATlas scATAC fragments (GSE184462; NIH-funded, public)","hg38 refGene promoters (UCSC)",
                "cross-cell-type accessibility prevalence background"],
               "public":True,"source_citation":CIT,"funding":"NIH (CATlas / Zhang et al. 2021 Cell)"},
              open(d+"/geneset.provenance.json","w"),indent=1)
    return 1

tot=0; raw=0; up=0; dn=0
for ct in celltypes:
    A=accessible[ct]; s=safe(ct)
    raw+=emit(f"CATLAS_{s}_accessible_raw",f"{ct}: promoter-accessible genes (RAW, no background control)",A,"raw",ct,{})
    U={g for g in A if N>1 and (prev[g]-1)/(N-1)<LOW}
    up+=emit(f"CATLAS_{s}_specific_Up",f"{ct}: cell-type-specifically accessible genes (background-corrected)",U,"specific_Up",ct,{})
    D=highprev-A
    dn+=emit(f"CATLAS_{s}_closed_Down",f"{ct}: broadly-open genes that are CLOSED in {ct} (background-corrected)",D,"closed_Down",ct,{})
tot=raw+up+dn
print(f"\nCELL TYPES: {N}\nSETS: raw={raw} specific_Up={up} closed_Down={dn} TOTAL={tot}\n-> {OUT}",flush=True)
