#!/usr/bin/env python3
# catlas_fetal_genesets.py — fetal cell-type accessible-gene sets from GSE149683 Cicero gene activity scores.
# Source: "A human cell atlas of fetal chromatin accessibility" (GSE149683; BICCN/NIH; sci-ATAC-seq3; public).
# Same background-prevalence control as adult CATLAS (catlas_combine_call.py).
#
# File_S6 format: long-format sparse triplet CSV
#   row = gene symbol
#   col = {barcode}_{cell_type}_{tissue}   (barcode has no underscores; tissue is last element)
#   x   = Cicero gene activity score (or NA = present but unscored)
# Presence of any entry for (gene, cell) = accessible in that cell (NA treated as present).
#
# Pipeline:
#   1. Parse File_S2 metadata -> total cells per cell_type (denominator)
#   2. Stream File_S6 -> count (gene, cell_type) co-occurrences
#   3. accessible(ct, gene) if count/total_cells >= MIN_FRAC
#   4. LOO prevalence control (same thresholds as adult CATLAS)
#
# Outputs (same format as catlas_combine_call.py):
#   CATLAS_FETAL_<ct>_accessible_raw  — fraction >= MIN_FRAC (uncontrolled)
#   CATLAS_FETAL_<ct>_specific_Up     — raw AND LOO prevalence < LOW (background-controlled)
import os, gzip, json, csv, re, collections
WORK=os.environ.get("WORK",os.path.expanduser("~/CFDE/catlas_work"))
S6=os.environ.get("FILE_S6",os.path.join(WORK,"fetal","GSE149683_File_S6.csv.gz"))
S2=os.environ.get("FILE_S2",os.path.join(WORK,"fetal","GSE149683_File_S2.metadata.txt.gz"))
OUT=os.environ.get("OUTDIR_FETAL",os.path.join(WORK,"fetal","catlas_fetal_genesets_full"))
os.makedirs(OUT,exist_ok=True)
LOG=open(os.path.join(OUT,"batch_log.txt"),"a")
def log(m): LOG.write(m+"\n"); LOG.flush(); print(m)
MIN_FRAC=float(os.environ.get("MIN_FRAC","0.10"))   # fraction of cells in CT that must have gene (raw)
LOW=float(os.environ.get("LOW","0.25"))             # LOO prevalence threshold (controlled)
MIN_GENES=int(os.environ.get("MIN_GENES","10"))
CIT=("GSE149683: A human cell atlas of fetal chromatin accessibility (BICCN/NIH; sci-ATAC-seq3; public). "
     "Cicero gene activity scores used as accessibility proxy; presence in sparse matrix = accessible in that cell. "
     "Same LOO prevalence control as adult CATLAS pipeline (Zhang et al. 2021 Cell).")

def safe(s): return "".join(c if (c.isalnum() or c in "._-") else "_" for c in s)[:70]

# --- Step 1: total cells per cell_type from File_S2 metadata ---
# File_S2 expected: TSV with header; look for a column named 'cell_type' or similar
log(f"Loading cell counts from File_S2: {S2}")
ct_total = collections.Counter()
with gzip.open(S2,"rt") as f:
    header=f.readline().rstrip("\n").split("\t")
    # find cell_type column (case-insensitive scan)
    ct_col=next((i for i,h in enumerate(header) if "cell_type" in h.lower() or "celltype" in h.lower()),None)
    if ct_col is None:
        # fallback: print header and try last column
        log(f"WARNING: no cell_type column found in File_S2 header: {header[:10]}")
        ct_col=len(header)-1
    log(f"Using column {ct_col} ('{header[ct_col]}') as cell_type")
    for line in f:
        f2=line.rstrip("\n").split("\t")
        if len(f2)>ct_col: ct_total[f2[ct_col]]+=1
log(f"Cell types in File_S2: {len(ct_total)} | total cells: {sum(ct_total.values())}")

# --- Step 2: stream File_S6, count (gene, cell_type) occurrences ---
# col format: {barcode}_{cell_type}_{tissue}
# barcode contains no underscores; tissue is single word (last element)
# cell_type = col.split('_',1)[1].rsplit('_',1)[0]
log(f"Streaming File_S6: {S6}")
gene_ct = collections.Counter()   # (gene, cell_type) -> cell count
rows_read=0; rows_skipped=0
with gzip.open(S6,"rt") as f:
    reader=csv.reader(f)
    next(reader)  # skip header: row,col,x
    for row in reader:
        if len(row)<2: rows_skipped+=1; continue
        gene=row[0].strip(); col=row[1].strip()
        if not gene or not col: rows_skipped+=1; continue
        # extract cell_type from col
        try:
            after_barcode=col.split('_',1)[1]           # drop barcode prefix
            cell_type=after_barcode.rsplit('_',1)[0]    # drop tissue suffix
        except IndexError: rows_skipped+=1; continue
        gene_ct[(gene,cell_type)]+=1
        rows_read+=1
        if rows_read % 5_000_000 == 0: log(f"  ...{rows_read//1_000_000}M rows read")
log(f"Rows read: {rows_read} | skipped: {rows_skipped}")
log(f"Unique (gene, cell_type) pairs: {len(gene_ct)}")

# --- Step 3: threshold -> accessible sets ---
# Fall back to counting unique cells from File_S6 if File_S2 didn't give us a cell_type column
if not ct_total:
    log("WARNING: using gene_ct row-counts as cell totals (File_S2 parse failed)")
    ct_total_fallback=collections.Counter()
    for (g,ct),c in gene_ct.items(): ct_total_fallback[ct]+=c
    ct_total=ct_total_fallback

accessible=collections.defaultdict(set)
for (gene,ct),count in gene_ct.items():
    total=ct_total.get(ct,0)
    if total>0 and count/total>=MIN_FRAC:
        accessible[ct].add(gene)

celltypes=sorted(accessible)
for ct in celltypes:
    log(f"  {ct}: {len(accessible[ct])} accessible genes ({ct_total.get(ct,0)} cells)")

# --- Step 4: LOO prevalence control ---
N=len(celltypes)
prev=collections.Counter()
for ct in celltypes:
    for g in accessible[ct]: prev[g]+=1
allgenes=set(prev)
log(f"Cell types: {N} | genes accessible somewhere: {len(allgenes)}")
log(f"Broadly active (>0.75 prevalence): {sum(1 for c in prev.values() if N>1 and c/(N-1)>0.75)}")

def emit(name,desc,genes,cls,ct,extra={}):
    if len(genes)<MIN_GENES: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    m={"standard_name":name,"library":"CATLAS_fetal_scATAC_accessibility","description":desc,
       "version":"1.0","file_type":"geneset","n_genes":len(genes),"organism":"human",
       "assembly":"hg19/hg38","assay":"scATAC-seq (sci-ATAC-seq3; fetal CATlas GSE149683)",
       "cell_type":ct,"class":cls,"n_cells":ct_total.get(ct,0),
       "method":f"Cicero gene activity: fraction of cells >= {MIN_FRAC}; LOO prevalence background (same as adult CATLAS)",
       "min_frac":MIN_FRAC,"low":LOW,
       "caveat":("Cicero gene activity used as accessibility proxy; NOT raw fragment-based. "
                 "Cell-type matching to adult CATLAS is name-based inference."),
       "derived_in_this_work":True,"source":CIT}
    m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"catlas_fetal_accessible_genes",
               "inputs":["GSE149683 File_S6 Cicero gene activity (NIH/BICCN; public)",
                         f"File_S2 cell metadata (cell type counts)",
                         "cross-cell-type LOO prevalence background"],
               "public":True,"source_citation":CIT,
               "funding":"NIH BRAIN Initiative / BICCN (fetal CATlas)"},
              open(d+"/geneset.provenance.json","w"),indent=1)
    return 1

raw_n=0; ctrl_n=0
for ct in celltypes:
    s=safe(ct); A=accessible[ct]
    raw_n+=emit(f"CATLAS_FETAL_{s}_accessible_raw",
                f"Fetal {ct}: Cicero-active genes (RAW, >= {MIN_FRAC} of cells; GSE149683)",
                A,"raw",ct)
    U={g for g in A if N>1 and (prev[g]-1)/(N-1)<LOW}
    ctrl_n+=emit(f"CATLAS_FETAL_{s}_specific_Up",
                 f"Fetal {ct}: cell-type-specifically active genes (LOO prevalence<{LOW}; GSE149683)",
                 U,"specific_Up",ct)

log(f"=== DONE fetal gene sets: raw={raw_n} specific_Up={ctrl_n} TOTAL={raw_n+ctrl_n} ===")
LOG.close()
