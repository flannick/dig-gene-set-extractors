import csv, os
GTEX="/Users/gage/Codex/PIGEAN_EAGGL/Data/gtex_tstat/GTEx.tstat.hgnc.tsv"

# A: initiation enzymes — unglycosylated -> IMMATURE (high-mannose). Dolichol assembly, OST transfer, ER trim.
A = """DPAGT1 ALG1 ALG2 ALG3 ALG5 ALG6 ALG8 ALG9 ALG10 ALG10B ALG11 ALG12 ALG13 ALG14
DPM1 DPM2 DPM3 MPDU1 DOLK SRD5A3 RFT1
STT3A STT3B RPN1 RPN2 DDOST DAD1 OSTC TUSC3 MAGT1 OST4 KRTCAP2
MOGS GANAB PRKCSH MAN1B1""".split()

# B: maturation enzymes — IMMATURE -> MATURE (Golgi processing + terminal capping).
B = """MAN1A1 MAN1A2 MAN1C1 MAN2A1 MAN2A2
MGAT1 MGAT2 MGAT3 MGAT4A MGAT4B MGAT4C MGAT5 MGAT5B
B4GALT1 B4GALT2 B4GALT3 B4GALT4 B4GALT5 B4GALT6 B4GALT7
ST3GAL1 ST3GAL2 ST3GAL3 ST3GAL4 ST3GAL5 ST3GAL6 ST6GAL1 ST6GAL2
ST6GALNAC1 ST6GALNAC2 ST6GALNAC3 ST6GALNAC4 ST6GALNAC6 ST8SIA1 ST8SIA2 ST8SIA4
FUT1 FUT2 FUT3 FUT4 FUT5 FUT6 FUT7 FUT8 FUT9 FUT10 FUT11 POFUT1 POFUT2
GALNT1 GALNT2 GALNT3 GALNT4 GALNT6 GALNT7 GALNT10 GALNT12 C1GALT1 C1GALT1C1 GCNT1 GCNT3""".split()

rows=list(csv.reader(open(GTEX),delimiter='\t'))
tissues=rows[0][1:]
tmap={r[0]:[float(x) for x in r[1:]] for r in rows[1:]}

def report(name, genes):
    found=[g for g in genes if g in tmap]; missing=[g for g in genes if g not in tmap]
    print(f"\n=== SET {name}: {len(genes)} curated, {len(found)} in GTEx, missing(no GTEx row, =null not absent): {missing}")
    return found

Afound=report("A (initiation→immature)",A)
Bfound=report("B (immature→mature)",B)

# per-tissue annotation; ENRICHED if t>=2 (positive evidence); absence/low = null, NOT exclusion
def enriched(g,ti): return tmap[g][ti]>=2.0
for tname in ["Liver","Brain_Cortex","Pancreas","Muscle_Skeletal"]:
    ti=tissues.index(tname)
    aen=[g for g in Afound if enriched(g,ti)]; ben=[g for g in Bfound if enriched(g,ti)]
    print(f"\n[{tname}]  A enriched(t>=2): {len(aen)}/{len(Afound)} {aen[:8]}")
    print(f"           B enriched(t>=2): {len(ben)}/{len(Bfound)} {sorted(ben, key=lambda g:-tmap[g][ti])[:8]}")
print(f"\n(total GTEx tissues available: {len(tissues)})")
