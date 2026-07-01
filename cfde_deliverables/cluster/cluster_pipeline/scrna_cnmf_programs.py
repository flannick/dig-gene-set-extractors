#!/usr/bin/env python3
# SELF-CONTAINED scRNA -> cNMF gene programs -> contract gene sets. No external toolchain.
# Uses the installed `cnmf` + `anndata`. Steps: balanced cell subsample -> AnnData(counts) ->
# cNMF prepare/factorize/combine/consensus -> top-N genes per program -> geneset.tsv/gmt/meta/provenance.
# Env: MATRIX_TSV META_TSV OUTDIR NAME CELL_ID_COL CELL_TYPE_COL DONOR_COL MAX_CELLS K NHVG NITER TOPN SEED ORGANISM CITATION
import os, json, glob, gc
import numpy as np, pandas as pd, anndata as ad
from cnmf import cNMF

M=os.environ["MATRIX_TSV"]; E=os.environ["META_TSV"]; OUT=os.environ["OUTDIR"]; NAME=os.environ.get("NAME","scrna")
CID=os.environ.get("CELL_ID_COL","cell_id"); CT=os.environ.get("CELL_TYPE_COL","cell_type"); DN=os.environ.get("DONOR_COL","donor_id")
MAXC=int(os.environ.get("MAX_CELLS","20000")); K=int(os.environ.get("K","10")); NHVG=int(os.environ.get("NHVG","2000"))
NITER=int(os.environ.get("NITER","20")); TOPN=int(os.environ.get("TOPN","100")); SEED=int(os.environ.get("SEED","14"))
ORG=os.environ.get("ORGANISM","human"); CITE=os.environ.get("CITATION",NAME)
os.makedirs(OUT,exist_ok=True); GS=os.path.join(OUT,"genesets"); os.makedirs(GS,exist_ok=True)

# 1) metadata + balanced subsample across cell types
meta=pd.read_csv(E,sep="\t",dtype=str).dropna(subset=[CID]).set_index(CID)
rng=np.random.RandomState(SEED)
if CT in meta.columns and len(meta)>MAXC:
    keep=[]
    for ct,idx in meta.groupby(CT).groups.items():
        idx=list(map(str,idx)); rng.shuffle(idx); keep+=idx[:max(1,MAXC//max(1,meta[CT].nunique()))]
    keep=keep[:MAXC]
else:
    keep=list(meta.index.astype(str))[:MAXC]
keep=set(map(str,keep))
print(f"[prep] metadata cells={len(meta)} selected={len(keep)}",flush=True)

# 2) stream the (cells x genes) matrix, keep only selected cells
parts=[]
for ch in pd.read_csv(M,sep="\t",index_col=0,chunksize=5000):
    ch.index=ch.index.astype(str)
    sub=ch[ch.index.isin(keep)]
    if len(sub): parts.append(sub.astype("float32"))
mat=pd.concat(parts); del parts; gc.collect()
mat=mat[~mat.index.duplicated()]
print(f"[prep] matrix subset: {mat.shape[0]} cells x {mat.shape[1]} genes",flush=True)

# 3) AnnData(counts)
common=[c for c in (CT,DN) if c in meta.columns]
obs=meta.loc[mat.index, common].copy() if common else pd.DataFrame(index=mat.index)
adata=ad.AnnData(X=mat.values, obs=obs, var=pd.DataFrame(index=mat.columns.astype(str)))
del mat; gc.collect()
# QC: cNMF divides by per-gene std, so all-zero cells / zero-variance genes -> NaN. Drop them
# (needed for sparse 10x data; harmless for dense SMART-seq).
import numpy as np
ct=np.asarray(adata.X.sum(1)).ravel(); adata=adata[ct>0].copy()
gt=np.asarray(adata.X.sum(0)).ravel(); adata=adata[:,gt>0].copy()
gv=np.asarray(adata.X.var(axis=0)).ravel(); adata=adata[:,gv>0].copy()
print(f"[prep] after QC: {adata.n_obs} cells x {adata.n_vars} genes (dropped zero cells/genes)",flush=True)
if adata.n_obs<50 or adata.n_vars<200:
    raise SystemExit(f"too few cells/genes after QC ({adata.n_obs}x{adata.n_vars}) — input matrix likely truncated/corrupt; clear work/ and re-run")
h5=os.path.join(OUT,NAME+".counts.h5ad"); adata.write(h5)
print(f"[prep] wrote {h5} ({adata.n_obs} x {adata.n_vars})",flush=True)

# 4) cNMF
c=cNMF(output_dir=OUT,name=NAME)
c.prepare(counts_fn=h5, components=[K], n_iter=NITER, seed=SEED, num_highvar_genes=NHVG)
c.factorize(worker_i=0, total_workers=1)
c.combine()
c.consensus(k=K, density_threshold=2.0, show_clustering=False)
print("[cnmf] consensus done",flush=True)

# 5) top-N genes per program from the gene_spectra_score matrix
hits=glob.glob(os.path.join(OUT,NAME,f"*gene_spectra_score.k_{K}.*txt")) or \
     glob.glob(os.path.join(OUT,NAME,f"*spectra_score*k_{K}*txt"))
if not hits: raise SystemExit(f"no gene_spectra_score file for k={K} under {os.path.join(OUT,NAME)}")
scores=pd.read_csv(hits[0],sep="\t",index_col=0)
if scores.shape[0]!=K and scores.shape[1]==K: scores=scores.T   # ensure programs x genes
def emit(prog,genes):
    nm=f"{NAME}_cNMF_program_{prog}"; d=os.path.join(GS,nm); os.makedirs(d,exist_ok=True); genes=[str(g) for g in genes]
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{nm}\tcNMF gene program {prog} ({ORG}); {CITE}\t"+"\t".join(genes)+"\n")
    cite=f"Consensus NMF (cNMF) gene program (K={K}; top {TOPN} genes by spectra score) from {CITE} ({ORG} scRNA-seq); NIH-funded, public."
    json.dump({"standard_name":nm,"library":f"scRNA_cNMF_programs_{NAME}","description":f"cNMF gene program {prog} ({ORG} scRNA-seq; {CITE})","version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":ORG,"method":"cNMF","K":K,"program":int(prog),"derived_in_this_work":True,"source":cite},open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":nm,"operation":"scrna_cnmf_program","inputs":[f"{CITE} (NIH; public scRNA-seq)"],"source_citation":cite,"public":True,"funding":"NIH"},open(d+"/geneset.provenance.json","w"),indent=1)
for prog in scores.index:
    emit(prog, scores.loc[prog].sort_values(ascending=False).head(TOPN).index.tolist())

# decision/inputs MANIFEST (zipped with the gene sets)
import datetime
manifest=f"""# Submission MANIFEST — {NAME}

Generated: {datetime.date.today().isoformat()}
Deliverable: scRNA-seq cNMF gene programs ({scores.shape[0]} gene sets)

## Data / inputs
- Source: {CITE}
- Organism: {ORG}
- Matrix: {os.path.basename(M)}  | Metadata: {os.path.basename(E)}
- Cells available: {len(meta)}  | Cells used (balanced subsample across {CT}): {adata.n_obs}  | Genes: {adata.n_vars}

## Compliance
- Access: publicly accessible (anonymous download).
- Funding: NIH-funded (BICCN / NIH BRAIN Initiative for the wired Allen datasets).
- Each gene set carries source citation + funding in its geneset.provenance.json.

## Method & decisions
- Consensus NMF (cNMF): prepare -> factorize ({NITER} iters) -> combine -> consensus.
- K (programs): {K}   | density_threshold: 2.0 (no replicate filtering)
- Highly-variable genes for factorization: {NHVG}
- Gene set per program = top {TOPN} genes by cNMF gene-spectra score.
- Cell subsample cap: {MAXC} (balanced across cell types)  | seed: {SEED}
- NOTE: programs are global (across the dataset), not split per cell type.

## Outputs
- {scores.shape[0]} programs, each a contract gene set (geneset.tsv / genesets.gmt / geneset.meta.json / geneset.provenance.json).
"""
open(os.path.join(GS,"MANIFEST.md"),"w").write(manifest)
print(f"[done] wrote {scores.shape[0]} cNMF program gene sets + MANIFEST.md to {GS}",flush=True)
