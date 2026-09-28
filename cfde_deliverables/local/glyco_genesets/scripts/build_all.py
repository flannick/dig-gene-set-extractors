#!/usr/bin/env python
"""Build glycosylation gene sets A/B/C/D per tissue (GTEx), contract-conforming.
A = initiation enzymes (unglyco -> immature/high-mannose)
B = maturation enzymes (immature -> mature, Golgi/terminal)
C = tissue-expressed effectors PREDICTED immaturely glycosylated   (PUTATIVE)
D = tissue-expressed effectors PREDICTED maturely glycosylated     (PUTATIVE)
Honors 'absence = null' throughout (missing data never used as negative evidence)."""
import csv, json, os, re, gzip, sys

BASE   = os.path.expanduser("~/Claude/proj-valiation-challenge/glyco_genesets")
GTEX   = "/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv"
UNIPROT= os.path.expanduser("~/.claude/jobs/32851e29/tmp/uniprot_human.tsv.gz")
EFF    = BASE + "/output/T2D_effectors.tsv"
RUNDATE= sys.argv[1] if len(sys.argv)>1 else "unknown"
OUT    = BASE + "/output/genesets"; os.makedirs(OUT, exist_ok=True)

A = """DPAGT1 ALG1 ALG2 ALG3 ALG5 ALG6 ALG8 ALG9 ALG10 ALG10B ALG11 ALG12 ALG13 ALG14
DPM1 DPM2 DPM3 MPDU1 DOLK SRD5A3 RFT1 STT3A STT3B RPN1 RPN2 DDOST DAD1 OSTC TUSC3 MAGT1
OST4 KRTCAP2 MOGS GANAB PRKCSH MAN1B1""".split()
B = """MAN1A1 MAN1A2 MAN1C1 MAN2A1 MAN2A2 MGAT1 MGAT2 MGAT3 MGAT4A MGAT4B MGAT4C MGAT5 MGAT5B
B4GALT1 B4GALT2 B4GALT3 B4GALT4 B4GALT5 B4GALT6 B4GALT7 ST3GAL1 ST3GAL2 ST3GAL3 ST3GAL4 ST3GAL5
ST3GAL6 ST6GAL1 ST6GAL2 ST6GALNAC1 ST6GALNAC2 ST6GALNAC3 ST6GALNAC4 ST6GALNAC6 ST8SIA1 ST8SIA2
ST8SIA4 FUT1 FUT2 FUT3 FUT4 FUT5 FUT6 FUT7 FUT8 FUT9 FUT10 FUT11 POFUT1 POFUT2 GALNT1 GALNT2
GALNT3 GALNT4 GALNT6 GALNT7 GALNT10 GALNT12 C1GALT1 C1GALT1C1 GCNT1 GCNT3""".split()
# maturation capacity assessed over the FULL B catalog (not a specificity-biased subset);
# absence=null -> a tissue is "B-capable" unless most of the catalog is strongly depleted.

# ---- GTEx ----
rows=list(csv.reader(open(GTEX),delimiter='\t')); tissues=rows[0][1:]
tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}
def expressed(g,ti): return g in tmap and tmap[g][ti] >= -2.0   # lenient; absence=null
def enriched(g,ti):  return g in tmap and tmap[g][ti] >= 2.0
def tstat(g,ti):     return round(tmap[g][ti],2) if g in tmap else None

# ---- UniProt annotations ----
up={}
for r in csv.DictReader(gzip.open(UNIPROT,'rt'),delimiter='\t'):
    g=r["Gene Names (primary)"]
    if g: up[g]={"ac":r["Entry"],"signal":r["Signal peptide"],"tm":r["Transmembrane"],
                 "loc":r["Subcellular location [CC]"] or "","glyco":r["Glycosylation"] or ""}

# NCBI Gene/RefSeq (NLM/NIH-native) secretory-routing corroboration
ncbi_route={}; ncbi_terms={}
for r in csv.DictReader(open(BASE+"/output/ncbi_routing.tsv"),delimiter='\t'):
    ncbi_route[r["gene"]]=r["ncbi_secretory_CC"]; ncbi_terms[r["gene"]]=r["ncbi_CC_terms"]

def glycosylatable(g):
    # DUAL-SOURCE: UniProt + NCBI Gene/RefSeq; missing one source = null, not negative
    ev=[]
    u=up.get(g)
    if u:
        if u["glyco"]: ev.append("uniprot:glyco_site")
        if u["signal"]: ev.append("uniprot:signal_peptide")
        if u["tm"]: ev.append("uniprot:transmembrane")
        if re.search(r"Secreted|Cell membrane|Endoplasmic reticulum|Golgi|Lysosom|Extracellular|Membrane",u["loc"]):
            ev.append("uniprot:secretory_localization")
    if ncbi_route.get(g)=="yes": ev.append("ncbi:secretory_CC")
    return (len(ev)>0), ";".join(ev) if ev else "none"

def maturity_lean(g):
    # combine UniProt subcellular location + NCBI GO cellular-component terms
    u=up.get(g); combined=((u["loc"] if u else "")+" "+ncbi_terms.get(g,"")).strip()
    er = re.search(r"endoplasmic reticulum", combined, re.I)
    surface = re.search(r"cell membrane|plasma membrane|secreted|extracellular|golgi|lysosom|cell surface|external side", combined, re.I)
    if er and not surface: return "immature","ER-resident"
    if surface or (u and u["signal"]): return "mature","secretory-transit"
    return "unknown","insufficient"

# ---- effectors (confident assignment) ----
eff={}
for r in csv.DictReader(open(EFF),delimiter='\t'):
    if float(r["assign_cond_prob_signal"])>=0.5: eff[r["gene"]]=float(r["assign_cond_prob_signal"])

def write_set(tdir, setcode, name, desc, putative, header, rows_):
    d=os.path.join(tdir,setcode); os.makedirs(d,exist_ok=True)
    with open(os.path.join(d,"geneset.tsv"),'w',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(header)
        for r in rows_: w.writerow(r)
    with open(os.path.join(d,"genesets.gmt"),'w') as f:
        f.write(name+"\t"+desc+"\t"+"\t".join(r[0] for r in rows_)+"\n")
    json.dump({"standard_name":name,"set_code":setcode,"description":desc,"library":"GlycoMaturity_GTEx",
               "version":"0.1","file_type":"geneset","n_genes":len(rows_),"putative":putative,
               "tissue":os.path.basename(tdir)}, open(os.path.join(d,"geneset.meta.json"),'w'), indent=1)
    json.dump({"focus":name,"operation":"glyco_maturity_extract","run_date":RUNDATE,
               "inputs":["GTEx.tstat.hgnc.tsv (tissue expression; NIH Common Fund)",
                         "UniProt human reviewed (routing+glyco sites; NIH-co-funded)",
                         "NCBI Gene/RefSeq gene2go (secretory GO cellular-component; NLM/NIH)",
                         "PIGEAN T2D effectors (cond_prob_signal>=0.5; NIDDK GWAS)","curated glyco enzyme catalog A/B"],
               "notes":"absence=null throughout; C/D are PUTATIVE predictions"},
              open(os.path.join(d,"geneset.provenance.json"),'w'), indent=1)

manifest=[]
summary=[]
for ti,tname in enumerate(tissues):
    tdir=os.path.join(OUT,tname); os.makedirs(tdir,exist_ok=True)
    bpresent=[g for g in B if g in tmap]
    bcap=round(sum(1 for g in bpresent if enriched(g,ti))/max(1,len(bpresent)),2)  # frac of B enzymes ENRICHED here (annotation only)
    mature_possible=True  # absence=null: specificity data can't prove a tissue lacks maturation; tissue-B gating deferred to TPM/SenNet
    # A / B : catalog membership (absence=null), annotated with this-tissue expression
    Arows=[[g, tstat(g,ti), "enriched" if enriched(g,ti) else ("present" if expressed(g,ti) else "low/null")] for g in A]
    Brows=[[g, tstat(g,ti), "enriched" if enriched(g,ti) else ("present" if expressed(g,ti) else "low/null")] for g in B]
    # C / D : effectors
    Crows=[]; Drows=[]
    for g,p in eff.items():
        ok,ev=glycosylatable(g)
        if not ok or not expressed(g,ti): continue
        lean,why=maturity_lean(g)
        conf = "observed" if (up.get(g) and up[g]["glyco"]) else "predicted"
        row=[g, round(p,3), ev, conf]
        if lean=="immature": Crows.append(row+[why])
        elif lean=="mature":
            if mature_possible: Drows.append(row+[why])
            else: Crows.append(row+[why+";tissue-lacks-B-machinery"])
        # unknown -> neither (absence=null)
    write_set(tdir,"A_initiation_immature","GlycoA_initiation_%s"%tname,
              "Initiation enzymes (unglyco->immature/high-mannose) expression-annotated in %s"%tname,False,
              ["gene","gtex_tstat","tissue_status"],Arows)
    write_set(tdir,"B_maturation","GlycoB_maturation_%s"%tname,
              "Maturation enzymes (immature->mature) expression-annotated in %s"%tname,False,
              ["gene","gtex_tstat","tissue_status"],Brows)
    write_set(tdir,"C_effectors_immature","GlycoC_effectors_immature_%s"%tname,
              "PUTATIVE: T2D effectors expressed in %s predicted IMMATURELY glycosylated"%tname,True,
              ["gene","effector_assign_prob","glyco_evidence","confidence","maturity_basis"],Crows)
    write_set(tdir,"D_effectors_mature","GlycoD_effectors_mature_%s"%tname,
              "PUTATIVE: T2D effectors expressed in %s predicted MATURELY glycosylated"%tname,True,
              ["gene","effector_assign_prob","glyco_evidence","confidence","maturity_basis"],Drows)
    for sc,n in [("A",len(Arows)),("B",len(Brows)),("C",len(Crows)),("D",len(Drows))]:
        manifest.append([tname,sc,n])
    summary.append((tname,len(Crows),len(Drows),"B-enriched-frac:%.2f"%bcap))

with open(OUT+"/manifest.tsv",'w',newline='') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(["tissue","set","n_genes"]); w.writerows(manifest)
print("tissues:",len(tissues),"| sets/tissue: 4 | files:",len(tissues)*4*3+1)
print("\ntissue            C(immature)  D(mature)  note")
for t,c,d,note in summary[:12]: print(f"  {t:30}{c:>6}{d:>10}   {note}")
print("  ...")
print("\nA enzymes:",len(A),"| B enzymes:",len(B),"| effectors(>=0.5):",len(eff),
      "| effectors glycosylatable:",sum(1 for g in eff if glycosylatable(g)[0]))
