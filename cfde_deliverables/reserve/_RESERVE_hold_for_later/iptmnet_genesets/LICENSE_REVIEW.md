# iPTMnet-derived gene sets — LICENSE / PROVENANCE REVIEW (HOLD before sharing)

**STATUS: PRODUCED but HELD.** These iPTMnet-derived PTM gene sets are built and ready, but are
**segregated from the shareable handoff** pending team license review. Do **NOT** upload/redistribute
(do not include in `HANDOFF_MANIFEST2.md`) until the group clears them.

## License
iPTMnet (PIR / U. Delaware; NIH/NIGMS-funded), release 6.2, is distributed under **CC BY-NC-SA 4.0** —
**N**on-**C**ommercial, **A**ttribution, **S**hare-**A**like; redistribution requires citing iPTMnet **and**
the per-row source databases. This is more restrictive than the open ENCODE/GTEx/GlyGen sets.

## Provenance — human evidence by SOURCE database (from ptm.txt `source` column)
| source | rows | terms |
|---|---|---|
| **HPRD** | **64,133** | **RESTRICTIVE — academic/non-commercial; redistribution restricted (dominant source)** |
| UniProt | 42,863 | CC BY 4.0 (permissive) |
| GlyGen | 25,012 | free/public |
| IEDB | 9,092 | free |
| RLIMS-P (text-mining) | 7,435 | derived |
| Signor | 5,508 | CC BY-SA / academic |
| neXtProt | 5,107 | CC BY 4.0 |
| PRO | 2,939 | free |
| dbSNO | 2,362 | academic |
| IntAct | 828 | CC BY 4.0 |

➡ The derived sets are **heavily HPRD-influenced**, which is the main redistribution concern.

## Question for the team
The NIH-showcase rule is "NIH-funded + publicly available." iPTMnet qualifies on both, **but**
CC BY-NC-SA + HPRD's restrictive terms may conflict with open sharing/redistribution of derived sets.
**Can we share these, and under what attribution — or must we filter/hold?**

## Mitigation options
1. **HOLD** entirely until cleared (current default).
2. **Permissive-sources-only version** — regenerate excluding HPRD (and any other restrictive sources),
   keeping UniProt/GlyGen/IEDB/neXtProt/IntAct. Clearly shareable. I can generate on request.
3. **Share with full attribution + NC/SA notice** if the team/legal deems it acceptable for an academic showcase.

## What was built (held here, not in handoff)
- `output/by_ptm_type/` — 15 PTM-type substrate sets
- `output/by_enzyme/` — 133 per-enzyme substrate sets
- `output/enzyme_effector_x_gtex/` — per-PTM-type enzyme + substrate sets, standalone + 749 ×GTEx tissue sets
