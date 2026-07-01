# iPTMnet for the gene-set challenge — what it is, sources, and the sharing question

**What it is.** iPTMnet is an integrated resource for protein post-translational modifications (PTMs),
developed by the Protein Information Resource (PIR) at Georgetown University. It links **enzymes →
substrates → sites** across many PTM types (phosphorylation, ubiquitination, acetylation, methylation,
SUMOylation, glycosylation, S-nitrosylation, and more), assembled from curated databases plus text-mining.
We used it to build PTM **enzyme** sets and PTM **substrate (effector)** sets per PTM type, standalone and
intersected with GTEx tissue expression — the PTM analog of our glyco A/B (enzymes) and C/D (effectors).

**Funding.** NIH/NIGMS — grants **U01GM120953**, R01GM080646, P20GM103446. (So it meets "NIH-funded.")

**License.** iPTMnet release 6.2 is distributed under **CC BY-NC-SA 4.0** — Non-Commercial, Attribution,
Share-Alike, and redistribution must **cite iPTMnet and the underlying source databases**. This is more
restrictive than the open ENCODE / GTEx / GlyGen resources we used elsewhere.

**Source-data composition (human evidence in our derived sets).** iPTMnet aggregates many DBs, each with
its own terms. The dominant source is **HPRD**, which is **academic/non-commercial and redistribution-
restricted** — so our derived sets are heavily HPRD-influenced:

| source DB | human rows | terms |
|---|---|---|
| **HPRD** | **64,133** | **restrictive (academic/non-commercial; redistribution limited)** |
| UniProt | 42,863 | CC BY 4.0 (open) |
| GlyGen | 25,012 | open |
| IEDB | 9,092 | open |
| RLIMS-P (text-mining) | 7,435 | derived |
| Signor | 5,508 | CC BY-SA / academic |
| neXtProt | 5,107 | CC BY 4.0 |
| PRO / dbSNO / IntAct | ~6,000 | mostly open |

**The question for the group.** The challenge rule is "NIH-funded + publicly available." iPTMnet is both,
**but** its CC BY-NC-SA license plus HPRD's restrictive terms may conflict with openly redistributing our
*derived* gene sets. Can we share them, under what attribution — or should we filter/hold?

**Options:**
1. **Hold** the iPTMnet sets out of the shared deliverable (current default).
2. **Permissive-sources-only version** — regenerate excluding HPRD (keep UniProt/GlyGen/IEDB/neXtProt/IntAct);
   clearly shareable. (Can be produced quickly.)
3. **Share with full attribution + NC/SA notice** if the team/legal considers it acceptable for an academic showcase.

**Citation if shared.** Huang et al., *iPTMnet: an integrated resource for PTM network discovery*, Nucleic
Acids Research 2018; plus citation of the per-record source databases (see iPTMnet `ptm.txt` `source` column).

*(Status: these iPTMnet sets are BUILT but HELD — segregated from the shareable handoff pending this review.)*
