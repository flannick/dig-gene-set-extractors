#!/usr/bin/env python3
# catlas_maturation_contrast.py — adult vs fetal CATLAS maturation contrast gene sets.
# Matches adult cell types (GSE184462 Batch 13) to fetal cell types (GSE149683) by normalized name.
# For each matched pair produces 6 gene set types (raw + controlled x 3 directions):
#   CATLAS_<ct>_raw_matureOpen_Up     — adult_raw - fetal_raw  (opens during maturation)
#   CATLAS_<ct>_raw_fetalOpen_Up      — fetal_raw - adult_raw  (closes during maturation)
#   CATLAS_<ct>_raw_developmentalStable — adult_raw & fetal_raw (stable through maturation)
#   CATLAS_<ct>_ctrl_matureOpen_Up    — adult_ctrl - fetal_ctrl (controlled)
#   CATLAS_<ct>_ctrl_fetalOpen_Up     — fetal_ctrl - adult_ctrl (controlled)
#   CATLAS_<ct>_ctrl_developmentalStable — adult_ctrl & fetal_ctrl (controlled)
import os, json, re, collections
WORK=os.environ.get("WORK",os.path.expanduser("~/CFDE/catlas_work"))
ADULT_SETS=os.environ.get("ADULT_SETS",os.path.join(WORK,"catlas_genesets_full"))
FETAL_SETS=os.environ.get("FETAL_SETS",os.path.join(WORK,"fetal","catlas_fetal_genesets_full"))
OUT=os.environ.get("OUTDIR_CONTRAST",os.path.join(WORK,"fetal","catlas_maturation_contrast")); os.makedirs(OUT,exist_ok=True)
HERE=os.path.dirname(os.path.abspath(__file__))
MAP_FILE=os.environ.get("CT_MAP",os.path.join(HERE,"catlas_celltype_map.tsv"))
LOG=open(os.path.join(OUT,"batch_log.txt"),"a")
def log(m): LOG.write(m+"\n"); LOG.flush(); print(m)
MIN_GENES=int(os.environ.get("MIN_GENES","10"))
CIT_A=("Zhang K et al. A single-cell atlas of chromatin accessibility in the human genome. "
       "Cell 2021;184(24):5985-6001 (GSE184462; adult CATlas). NIH-funded, public.")
CIT_F=("GSE149683: A human cell atlas of fetal chromatin accessibility (BICCN/NIH; sci-ATAC-seq3; public).")

def safe(s): return "".join(c if (c.isalnum() or c in "._-") else "_" for c in s)[:70]

def load_map(path):
    """Load fetal->adult mapping from TSV; adult may be comma-separated for multi-type union."""
    m={}
    if not os.path.isfile(path): return m
    for line in open(path):
        line=line.strip()
        if not line or line.startswith('#'): continue
        parts=line.split('\t')
        if len(parts)<2: continue
        fetal=parts[0].strip(); adult=[a.strip() for a in parts[1].split(',') if a.strip()]
        m[fetal]=adult
    return m

def load_sets(setdir,suffix):
    """Load gene sets with given suffix from a directory; returns {cell_type: set(genes)}."""
    result={}
    if not os.path.isdir(setdir): return result
    for name in os.listdir(setdir):
        if not name.endswith(suffix): continue
        meta_f=os.path.join(setdir,name,"geneset.meta.json")
        tsv_f=os.path.join(setdir,name,"geneset.tsv")
        if not os.path.isfile(tsv_f): continue
        ct=None
        if os.path.isfile(meta_f):
            try: ct=json.load(open(meta_f)).get("cell_type")
            except Exception: pass
        if not ct: ct=name  # fallback
        genes=set(l.strip() for l in open(tsv_f) if l.strip() and l.strip()!="gene")
        result[ct]=genes
    return result

log(f"Loading cell type map from {MAP_FILE}")
ct_map=load_map(MAP_FILE)
log(f"Mapped fetal cell types: {len(ct_map)}")

log(f"Loading adult sets from {ADULT_SETS}")
adult_raw =load_sets(ADULT_SETS,"_accessible_raw")
adult_ctrl=load_sets(ADULT_SETS,"_specific_Up")
log(f"Adult: raw={len(adult_raw)} controlled={len(adult_ctrl)}")

log(f"Loading fetal sets from {FETAL_SETS}")
fetal_raw =load_sets(FETAL_SETS,"_accessible_raw")
fetal_ctrl=load_sets(FETAL_SETS,"_specific_Up")
log(f"Fetal: raw={len(fetal_raw)} controlled={len(fetal_ctrl)}")

# Build direct lookup for adult (by cell_type field from meta.json)
adult_raw_by_ct ={k:v for k,v in adult_raw.items()}
adult_ctrl_by_ct={k:v for k,v in adult_ctrl.items()}

def get_union(ct_list, lookup):
    """Return union of gene sets for a list of adult cell types."""
    result=set()
    found=[]
    for ct in ct_list:
        if ct in lookup: result|=lookup[ct]; found.append(ct)
    return result, found

def emit(name,desc,genes,cls,adult_ct,fetal_ct,extra={}):
    if len(genes)<MIN_GENES: return 0
    d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True); genes=sorted(genes)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    m={"standard_name":name,"library":"CATLAS_maturation_contrast","description":desc,
       "version":"1.0","file_type":"geneset","n_genes":len(genes),"organism":"human",
       "class":cls,"adult_cell_type":adult_ct,"fetal_cell_type":fetal_ct,
       "method":"set contrast: adult CATlas (GSE184462) vs fetal CATlas (GSE149683); name-matched cell types",
       "caveat":("Cell type matching is name-based inference; adult and fetal cell types may not be exact "
                 "developmental counterparts. Maturation direction inferred, not measured longitudinally."),
       "derived_in_this_work":True,"source_adult":CIT_A,"source_fetal":CIT_F}
    m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"catlas_maturation_contrast",
               "inputs":[f"Adult: {adult_ct} (GSE184462)",f"Fetal: {fetal_ct} (GSE149683)"],
               "public":True,"funding":"NIH (adult CATlas + BICCN fetal CATlas)"},
              open(d+"/geneset.provenance.json","w"),indent=1)
    return 1

matched=0; skipped=0; n_sets=0
for fetal_ct,f_raw in sorted(fetal_raw.items()):
    if fetal_ct not in ct_map:
        log(f"  NO MAP: {fetal_ct}")
        skipped+=1; continue
    adult_cts=ct_map[fetal_ct]
    a_raw,found_raw=get_union(adult_cts,adult_raw_by_ct)
    a_ctrl,found_ctrl=get_union(adult_cts,adult_ctrl_by_ct)
    if not found_raw:
        log(f"  NO ADULT DATA: {fetal_ct} -> {adult_cts}")
        skipped+=1; continue
    f_ctrl=fetal_ctrl.get(fetal_ct,set())
    adult_ct=",".join(found_raw)
    s=safe(fetal_ct); matched+=1

    contrasts=[
        # (adult_genes, fetal_genes, suffix, class, description_template)
        (a_raw, f_raw, "raw_matureOpen_Up",      "raw_matureOpen",
         f"{adult_ct}: opens during maturation — adult-accessible, fetal-closed (raw)"),
        (f_raw, a_raw, "raw_fetalOpen_Up",        "raw_fetalOpen",
         f"{adult_ct}: closes during maturation — fetal-accessible, adult-closed (raw)"),
        (None,  None,  "raw_developmentalStable", "raw_stable",
         f"{adult_ct}: chromatin stable through maturation — accessible in both fetal and adult (raw)"),
        (a_ctrl,f_ctrl,"ctrl_matureOpen_Up",      "ctrl_matureOpen",
         f"{adult_ct}: opens during maturation — adult-specific, fetal-not-specific (controlled)"),
        (f_ctrl,a_ctrl,"ctrl_fetalOpen_Up",        "ctrl_fetalOpen",
         f"{adult_ct}: closes during maturation — fetal-specific, adult-not-specific (controlled)"),
        (None,  None,  "ctrl_developmentalStable","ctrl_stable",
         f"{adult_ct}: stable specific chromatin — specific in both fetal and adult (controlled)"),
    ]

    for a,f,suffix,cls,desc in contrasts:
        if "Stable" in suffix:
            # Stable = intersection of the relevant pair
            base="raw" if suffix.startswith("raw") else "ctrl"
            if base=="raw": genes=a_raw & f_raw
            else:           genes=a_ctrl & f_ctrl
        else:
            genes=a - f   # first set minus second
        n_sets+=emit(f"CATLAS_{s}_{suffix}",desc,genes,cls,adult_ct,fetal_ct)

    log(f"  matched: {adult_ct} <-> {fetal_ct}")

log(f"\nMatched pairs: {matched} | unmatched fetal: {skipped}")
log(f"=== DONE maturation contrast: {n_sets} sets ===")
LOG.close()
