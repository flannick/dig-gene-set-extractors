#!/usr/bin/env python3
# #1: GENCODE cis-antisense -> protein-coding target map. A lncRNA/antisense gene whose locus overlaps a
# protein-coding gene on the OPPOSITE strand is treated as a cis-antisense regulator of that coding gene.
# Emits: mapping table + gene set of coding genes under cis-antisense regulation + x GTEx. PORTABLE stdlib.
import gzip, json, os, collections, urllib.request, bisect, csv
HOME=os.path.expanduser("~/Claude/proj-valiation-challenge")
OUT=os.path.join(HOME,"gencode_antisense","output"); os.makedirs(OUT,exist_ok=True)
TMP=os.environ.get("TMPDIR","/Users/gage/.claude/jobs/32851e29/tmp")
GTEX="/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv"; THR=4.0
REL=os.environ.get("GENCODE_REL","46")
GTF=os.path.join(TMP,f"gencode.v{REL}.annotation.gtf.gz")
if not os.path.exists(GTF):
    urllib.request.urlretrieve(f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{REL}/gencode.v{REL}.annotation.gtf.gz",GTF)
def attr(s,k):
    i=s.find(k+' "')
    if i<0: return None
    i+=len(k)+2; return s[i:s.find('"',i)]
coding=collections.defaultdict(list); ncrna=[]
for line in gzip.open(GTF,'rt'):
    if line.startswith("#"): continue
    f=line.split("\t")
    if len(f)<9 or f[2]!="gene": continue
    ch,st,en,strand,a=f[0],int(f[3]),int(f[4]),f[6],f[8]
    gt=attr(a,"gene_type"); nm=attr(a,"gene_name")
    if gt=="protein_coding": coding[ch].append((st,en,strand,nm))
    elif gt in ("lncRNA","antisense","antisense_RNA"): ncrna.append((ch,st,en,strand,nm))
for ch in coding: coding[ch].sort()
starts={ch:[x[0] for x in v] for ch,v in coding.items()}
pairs=[]; targets=set()
for ch,st,en,strand,nm in ncrna:
    arr=coding.get(ch);
    if not arr: continue
    j=bisect.bisect_left(starts[ch],en)
    k=j-1
    while k>=0 and arr[k][1]>=st:   # coding.start<=en and coding.end>=st => overlap
        cs,ce,cstrand,cnm=arr[k]
        if cstrand!=strand and cnm and nm:
            pairs.append((nm,cnm,ch,strand)); targets.add(cnm)
        k-=1
print(f"lncRNA/antisense genes: {len(ncrna)} | cis-antisense pairs: {len(pairs)} | coding targets: {len(targets)}")
# mapping table
with open(os.path.join(OUT,"antisense_to_coding.tsv"),"w") as fh:
    fh.write("antisense_gene\tcoding_target\tchrom\tantisense_strand\n")
    for p in sorted(set(pairs)): fh.write("\t".join(map(str,p))+"\n")
def emit(d,name,desc,genes,extra,gx=False):
    os.makedirs(d,exist_ok=True); genes=sorted(g for g in genes if g)
    open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(genes)+"\n")
    open(d+"/genesets.gmt","w").write(f"{name}\t{desc}\t"+"\t".join(genes)+"\n")
    cite=f"Derived from GENCODE v{REL} (opposite-strand locus overlap; lncRNA/antisense vs protein_coding), GRCh38, GENCODE/NHGRI, public"+(f"; intersected with GTEx tissue-enrichment (t>={THR}; NIH Common Fund)" if gx else "")+"."
    m={"standard_name":name,"library":"GENCODE_cis_antisense","description":desc,"version":"0.1","file_type":"geneset","n_genes":len(genes),"organism":"human","assembly":"GRCh38","derived_in_this_work":True,"source":cite}; m.update(extra)
    json.dump(m,open(d+"/geneset.meta.json","w"),indent=1)
    json.dump({"focus":name,"operation":"gencode_cis_antisense"+("_x_gtex" if gx else ""),"inputs":[f"GENCODE v{REL} annotation (NHGRI; public)"]+(["GTEx.tstat.hgnc.tsv (NIH Common Fund)"] if gx else []),"source_citation":cite,"public":True,"funding":"NIH/NHGRI (GENCODE) "+("+ NIH Common Fund (GTEx)" if gx else "")},open(d+"/geneset.provenance.json","w"),indent=1)
emit(os.path.join(OUT,"GENCODE_coding_with_cis_antisense"),"GENCODE_coding_with_cis_antisense",
     "Protein-coding genes overlapped by a cis-antisense lncRNA (opposite strand; potential antisense regulation)",targets,{})
# x GTEx
rows=list(csv.reader(open(GTEX),delimiter='\t')); tissues=rows[0][1:]; tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}
def safe(s): return "".join(c if c.isalnum() else "_" for c in s)[:60]
nx=0
for ti,t in enumerate(tissues):
    e=sorted(g for g in targets if g in tmap and tmap[g][ti]>=THR)
    if not e: continue
    emit(os.path.join(OUT,"x_GTEx",f"GENCODE_coding_with_cis_antisense_x_GTEx_enriched_{safe(t)}"),
         f"GENCODE_coding_with_cis_antisense_x_GTEx_enriched_{safe(t)}",f"cis-antisense coding targets GTEx-enriched (t>={THR}) in {t}",e,{"tissue":t},gx=True); nx+=1
print("x_GTEx tissue sets:",nx)
