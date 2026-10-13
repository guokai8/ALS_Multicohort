# Amendment 13: is the ventral-horn rise explained by microglial abundance? (W70)

Date: 2026-10-13
Amends: plan v1.0 (sha256 da63b81a...bfb3bd)
Status: committed **before** W70 was run. Prompted by the internal story evaluation (`09_internal_review/story_evaluation_2026-10-13.md`), which flagged this as the most likely reviewer question.

Inspected before writing (no outcome): the presence of the microglial marker genes below in both spatial platforms. No microglial score has been computed.

## Question

W67 (lumbar) and W68 (cervical) found that the frozen ALS-up module rises more in motor than in sensory grey matter in ALS. Its genes track microglial markers in tissue, so the rise could reflect more microglia (composition) in the ALS ventral horn rather than a change in expression per cell. This analysis tests whether the territory × disease interaction persists after adjusting for microglial content.

## Design

**Inputs and donors.** Identical to the primary analyses:
- Lumbar cord: W67 inputs, ventral vs dorsal horn; independent donors with age; 18 ALS, 5 control.
- Cervical cord: W68 inputs, anterior horn vs other grey matter; 9 ALS, 7 control.

**Microglial score (per pseudobulk unit).**
- Primary: the mean edgeR log2-CPM (TMM, prior count 1, computed over the units of the comparison) of the W66 microglia marker set: P2RY12, CX3CR1, CSF1R, PTPRC, CD68, TREM2, TYROBP, CD163, ITGAX. No marker is a module gene.
- Sensitivity scores:
  - homeostatic only: P2RY12, CX3CR1, CSF1R, TMEM119;
  - activated only: CD68, TREM2, TYROBP, CD163, ITGAX.

**Model.** The W67/W68 model plus the centred unit-level score:
- lumbar: `~ territory * disease + sex + age + microglia_score`;
- cervical: `~ territory * disease + microglia_score`;
- both donor-blocked (duplicateCorrelation, two passes).

The interaction contrast is tested with fry (primary) and camera, for the frozen ALS-up module (primary) and the consensus up arm (secondary).

**Combination.** Both cohort-level results are new, so the combination is prespecified here: one-sided fry P in the Up direction, combined by unweighted Stouffer.

**Attenuation.** Defined as 1 − (adjusted interaction median log2FC / unadjusted interaction median log2FC). The unadjusted values come from W67/W68 primary.

## Reading (frozen ALS-up module, primary score)

| Outcome | Reading |
|---|---|
| Combined one-sided P < 0.05, interaction direction Up in both cohorts, and attenuation < 50% in both | **Not explained by microglial abundance**: a change beyond microglial content |
| Combined P < 0.05 and Up in both, but attenuation ≥ 50% in at least one cohort | **Partly shared** with microglial abundance |
| Combined P ≥ 0.05, or direction not Up in both | **Largely shared** with microglial abundance (consistent with microgliosis in motor territory) |

Every outcome is reported. A "largely shared" result does not contradict the ventral-horn finding; it changes its interpretation from a per-cell expression change to increased microglial content.

**Descriptive.** The microglial score in ALS vs control in each territory (difference of means, Wilcoxon).

## Outputs

- `results/W70_microglia_adjusted_territory/`
- Script: `scripts/run_w70_microglia_adjusted_territory.R`
