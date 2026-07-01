#!/usr/bin/env python
"""Standalone SenNet tissue-resolved senescence gene set: SenSkin (PUBLISHED, cited) + SenSkin x GTEx.
PROVENANCE: this is a PRE-EXISTING PUBLISHED gene set, transcribed verbatim from the source paper's
Table 1 — NOT derived in this work."""
import csv, json, os
BASE=os.path.expanduser("~/Claude/proj-valiation-challenge/glyco_genesets")

# SenSkin (165 genes) — Wyles SP et al. GeroScience 2025;47(3):2631-2638; doi:10.1007/s11357-025-01568-y, Table 1.
# Transcribed verbatim (symbols as published, incl. legacy aliases). 1 down-regulated: LMNB1.
SENSKIN="""A2M ADAMTS1 ADAMTS4 ADAMTS9 ADAMTSL4 ANGPTL4 ANPEP ANXA2 APOE APOLD1 AQP1 ARID5A B2M BCL3 BCL6
BHLHE40 BTG2 C10orf10 C11orf96 C1QTNF1 C1RL C1S C3 CASP4 CCNL1 CD14 CD59 CD63 CDK2 CEBPB CEBPD CFB CHSY1
CLEC3B CREM CRISPLD2 CSF1 CSRNP1 CTSB CTSL CTSZ DDR1 DEC1 DUSP1 DUSP5 EGFL7 EGR1 EIF4A1 EMP1 ETS2 F3 FOS
FOSL1 FOSL2 GADD45A GADD45B GLIPR1 GPR4 HAPLN3 HILPDA ICAM1 IER2 IER3 IFI16 IGFBP2 IGFBP3 IGFBP4 IGFBP7
IL1R1 IL4R IL6 INHBB ITPKC ITPRIP JUN JUNB KIAA0040 KLF10 KLF6 KRT15 KRT18 LAMA5 LMNB1 LOX LTBP1 MAFF
MAN2B1 MAPK11 MCL1 MIDN MOV10 MT1A MT1M MT1X MYC NAMPT NEDD9 NFIL3 NFKB2 NFKBIZ NR4A1 NUCB1 NXT1 OSMR
PDGFB PDLIM1 PIM1 PLAUR PLSCR1 PLTP PNP PNRC1 PPP1R18 PPRC1 PROS1 RASD1 RGS16 RNASET2 RND3 RNF122 RPS3
SAT1 SBNO2 SEMA3F SEMA4B SERTAD1 SERPINB1 SERPINE1 SERPING1 SFN SHC1 SLC25A25 SLC2A3 SLC39A1 SLC39A14
SLCO4A1 SNAI1 SNRPC SOD2 SOCS3 STAT3 STC1 SUSD6 TGFBI TGIF1 THBD THBS1 TIMP1 TIMP3 TNFRSF10B TNFRSF10D
TNFRSF12A TNFRSF1A TP53 TRIB1 TSKU UBC VCAM1 VEGFA VEGFC YBX3 ZC3H12A ZFP36""".split()
CITE="Wyles SP, Yu GT, Ganier C, Tchkonia T, Lynch MD, Kuchel GA, Kirkland JL. SenSkin: a human skin-specific cellular senescence gene set. GeroScience. 2025;47(3):2631-2638. doi:10.1007/s11357-025-01568-y (Table 1)."

# --- standalone contract gene set ---
d=BASE+"/output/sennet/SenSkin"; os.makedirs(d,exist_ok=True)
open(d+"/geneset.tsv","w").write("gene\n"+"\n".join(SENSKIN)+"\n")
open(d+"/genesets.gmt","w").write("SenSkin_skin_senescence\tSenNet skin-specific senescence gene set (PUBLISHED, Wyles 2025)\t"+"\t".join(SENSKIN)+"\n")
json.dump({"standard_name":"SenSkin_skin_senescence","library":"SenNet","description":
  "Human SKIN-SPECIFIC cellular senescence gene set (SenSkin). PRE-EXISTING PUBLISHED set, transcribed verbatim from source Table 1; NOT derived here.",
  "version":"1.0","file_type":"geneset","n_genes":len(SENSKIN),"organism":"human","putative":False,
  "published_source":CITE,"derived_in_this_work":False}, open(d+"/geneset.meta.json","w"),indent=1)
json.dump({"focus":"SenSkin_skin_senescence","operation":"import_published_geneset",
  "published_source":CITE,"note":"This resource has been published before; transcribed verbatim, not derived here.",
  "funding":"senescence research (NIA-ecosystem; verify exact grant)","public":True}, open(d+"/geneset.provenance.json","w"),indent=1)
print(f"SenSkin standalone written: {len(SENSKIN)} genes")

# --- SenSkin x GTEx (validation: a skin-specific set should enrich in GTEx skin) ---
G="/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv"
rows=list(csv.reader(open(G),delimiter='\t')); tissues=rows[0][1:]
tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}
present=[g for g in SENSKIN if g in tmap]; absent=[g for g in SENSKIN if g not in tmap]
od=BASE+"/output/sennet_gtex"; os.makedirs(od,exist_ok=True)
with open(od+"/SenSkin_x_GTEx.tstat.matrix.tsv","w",newline='') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(["gene"]+tissues)
    for g in present: w.writerow([g]+[f"{v:.3f}" for v in tmap[g]])
counts=sorted([(t,sum(1 for g in present if tmap[g][i]>=2)) for i,t in enumerate(tissues)],key=lambda x:-x[1])
print(f"SenSkin in GTEx: {len(present)}/{len(SENSKIN)} (absent/legacy-symbol=null: {absent})")
print("top GTEx tissues by SenSkin enrichment (validation — expect skin near top):")
for t,n in counts[:8]: print(f"  {t:34} {n}")
