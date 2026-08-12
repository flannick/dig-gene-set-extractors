# HANDOFF_MANIFEST11 — Batch 11 (cross-resource interoperability gene sets, corrected inputs)
Date: 2026-07-01. All built on corrected/clean NIH-public inputs (Batch7/9/10 + GTEx/ClinVar/ClinGen/MitoCarta).
Contents:
- corrected_regulon_disease  : TF/eCLIP specifically-bound ∩ ClinVar disease genes
- corrected_regulon_mito     : TF/eCLIP specifically-bound ∩ MitoCarta
- corrected_regulon_ClinGenHI: TF/eCLIP specifically-bound ∩ ClinGen haploinsufficient
- corrected_regulon_GlyGenObserved: TF/eCLIP specifically-bound ∩ GlyGen observed glycoproteins
- chromatin_states           : bivalent / active / polycomb / transcribed / heterochromatin (corrected histone combinatorics)
- eqtl_accessibility_convergence, sqtl_accessibility_convergence : GTEx QTL × corrected accessibility (convergent / discordant)
- master_regulator_per_tissue: GTEx tissue-enriched genes bound by top corrected TFs (candidate master regulators)
- triple_convergence         : specifically accessible ∩ expressed ∩ eGene (multi-evidence cis-regulated)
- disease_multievidence      : ClinVar disease genes by # regulatory axes (eQTL + specific binding + specific accessibility)
Uses corrected (background-controlled) binding/accessibility only; excludes superseded raw regulon crosses.
All NIH-funded + public; nothing credential-walled or PIGEAN/EAGGL-derived.
