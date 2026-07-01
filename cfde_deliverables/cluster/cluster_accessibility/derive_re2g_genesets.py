#!/usr/bin/env python
"""One 'regulatory-target genes' set per ENCODE biosample, DERIVED from ENCODE-rE2G
thresholded element-gene links (cCRE/enhancer -> gene; distal+proximal). Public, no-auth.
Lean: download bed -> parse unique TargetGene (col 6) -> DELETE. Resume-safe. Cites each accession."""
import gzip, json, os, urllib.request
OUT=os.environ.get("OUTDIR", os.path.expanduser("~/Claude/proj-valiation-challenge/re2g_genesets/output"))
TMP=os.environ.get("TMPDIR","./tmp")
os.makedirs(OUT,exist_ok=True); os.makedirs(TMP,exist_ok=True)
LOG=open(OUT+"/batch_log.txt","a")
def log(m): LOG.write(m+"\n"); LOG.flush()

url=("https://www.encodeproject.org/search/?type=File&output_type=thresholded+element+gene+links"
     "&file_format=bed&assembly=GRCh38&status=released&limit=all&format=json")
d=json.load(urllib.request.urlopen(urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"x"}),timeout=180))
seen={}
for f in d.get("@graph",[]):
    bo=f.get("biosample_ontology",{}); bs=bo.get("term_name") if isinstance(bo,dict) else None
    if bs and bs not in seen: seen[bs]=(f.get("accession"), f.get("dataset","").strip("/").split("/")[-1], f.get("href"))
log(f"distinct biosamples (thresholded element-gene links): {len(seen)}")

def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
done=fail=0
for bs,(acc,ds,href) in seen.items():
    name=f"ENCODE_rE2G_regulatory_targets_{safe(bs)}"; dd=os.path.join(OUT,name)
    if os.path.exists(dd+"/geneset.tsv"): done+=1; continue
    if not href: continue
    pk=os.path.join(TMP,acc+".bed.gz")
    try:
        urllib.request.urlretrieve("https://www.encodeproject.org"+href, pk)
        genes=set()
        for line in gzip.open(pk,'rt'):
            if line.startswith("#"): continue
            c=line.rstrip("\n").split("\t")
            if len(c)>5 and c[5] and c[5]!="NA": genes.add(c[5])
        os.remove(pk); genes=sorted(genes)
        os.makedirs(dd,exist_ok=True)
        cite=f"Derived from ENCODE-rE2G thresholded element-gene links {acc} (annotation {ds}, biosample {bs}), GRCh38, ENCODE/NHGRI, public."
        open(dd+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
        open(dd+"/genesets.gmt","w").write(f"{name}\tConfident cCRE/enhancer regulatory-target genes (ENCODE-rE2G, distal+proximal) in {bs}\t"+"\t".join(genes)+"\n")
        json.dump({"standard_name":name,"library":"ENCODE_rE2G_regulatory_targets","description":f"Genes that are confident regulatory targets in {bs} via ENCODE-rE2G cCRE->gene links (DERIVED from ENCODE public predictions; distal+proximal).","version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","assembly":"GRCh38","biosample":bs,"derived_in_this_work":True,"source":cite},open(dd+"/geneset.meta.json","w"),indent=1)
        json.dump({"focus":name,"operation":"derive_regulatory_target_genes_from_rE2G","inputs":[f"ENCODE-rE2G {acc} thresholded element-gene links (NHGRI; public)"],"source_citation":cite,"public":True,"funding":"NIH/NHGRI (ENCODE)"},open(dd+"/geneset.provenance.json","w"),indent=1)
        done+=1; log(f"[{done}] {bs}: {len(genes)} regulatory-target genes")
    except Exception as e:
        fail+=1; log(f"FAIL {bs} {acc}: {str(e)[:80]}")
        try: os.remove(pk)
        except: pass
log(f"=== DONE: {done} sets, {fail} failures ===")
LOG.close()
