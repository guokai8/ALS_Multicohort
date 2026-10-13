# W70: territory interaction adjusted for microglial content (Amendment 13)

Amendment 13 (`edde685`) was committed before this analysis was run.

**Script.** `scripts/run_w70_microglia_adjusted_territory.R`

**Data and donors.** The same as W67 (lumbar, GSE330224: 18 ALS / 5 control) and W68 (cervical, GSE288365: 9 ALS / 7 control).

**Covariate.** A unit-level microglial score (mean log2-CPM of marker genes), centred and added to the W67/W68 interaction model. Three marker sets were used:

| Score | Genes |
|---|---|
| Primary (all microglia) | P2RY12, CX3CR1, CSF1R, PTPRC, CD68, TREM2, TYROBP, CD163, ITGAX |
| Homeostatic | P2RY12, CX3CR1, CSF1R, TMEM119 |
| Activated | CD68, TREM2, TYROBP, CD163, ITGAX |

## Results (`w70_adjusted_interaction_tests.tsv`, `w70_combined_reading.tsv`)

Interaction = motor vs sensory grey matter × ALS. Each cell gives fry P, then the median log2FC of the interaction; the unadjusted values are from W67/W68.

**Frozen ALS-up module**

| Adjustment | Lumbar | Cervical | Attenuation (lumbar / cervical) | Combined one-sided P | Reading (prespecified) |
|---|---|---|---|---|---|
| None (W67 / W68) | 0.012, 0.48 | 0.040, 0.41 | — | 6.1e-4 | — |
| **All microglia (primary)** | 0.061, 0.22 | 0.31, 0.10 | 53% / 75% | **0.020** | **Partly shared with microglial abundance** |
| Activated only | 0.18, 0.08 | 0.45, 0.18 | 84% / 56% | 0.069 | Largely shared |
| Homeostatic only | 0.013, 0.55 | 0.083, 0.20 | -16% / 50% | 0.0015 | Partly shared |

**Consensus up arm**

| Adjustment | Lumbar | Cervical | Attenuation (lumbar / cervical) | Combined one-sided P | Reading |
|---|---|---|---|---|---|
| None | 0.15, 0.086 | 0.0093, 0.19 | — | 0.0022 | — |
| All microglia | 0.041, 0.079 | 0.045, 0.17 | 9% / 14% | 0.0021 | **Not explained by microglial abundance** |
| Activated only | 0.18 | 0.092 | 21% / 24% | 0.016 | Not explained |
| Homeostatic only | 0.048 | 0.0076 | -3% / -25% | 5e-4 | Not explained |

## Microglial score by territory (`w70_microglial_score_by_territory.tsv`; mean ALS minus control, Wilcoxon P)

| Cohort | Score | Motor territory | Sensory territory |
|---|---|---|---|
| Cervical | All microglia | **+1.51 (0.0012)** | +0.71 (0.071) |
| Cervical | Activated | **+1.83 (0.0012)** | +0.81 (0.091) |
| Cervical | Homeostatic | +0.91 (0.042) | +0.40 (0.35) |
| Lumbar | All microglia | -0.41 (0.20) | -0.97 (0.015) |
| Lumbar | Activated | -0.04 (0.86) | -0.66 (0.067) |
| Lumbar | Homeostatic | -1.10 (0.019) | -1.50 (0.0025) |

## Reading

**Frozen module: partly shared with microglial abundance** (prespecified).
- Adjusting for microglial content removes half to three-quarters of the ventral-horn interaction. Most of the removal comes from activated (CD68+, TREM2+) myeloid markers.
- Adjusting for homeostatic markers alone leaves the lumbar interaction intact.
- The direction stays Up in every analysis, and the combined evidence remains nominally significant (P = 0.020).
- In cervical cord, activated myeloid content is itself higher in the ALS anterior horn.
- Interpretation: the response largely reflects an accumulation of activated myeloid cells in motor territory, with a residual component beyond it.

**Consensus signature: not explained by microglial content.** Its astrocyte-compartment change is independent of microglial abundance.

**Caveats.**
- The module contains genes expressed by myeloid cells (GPNMB, C3). Adjusting for myeloid markers is therefore expected to absorb part of the signal by construction, so attenuation is an upper bound on the share explained by composition.
- In lumbar cord, homeostatic marker expression is lower in ALS, consistent with the loss of homeostatic markers in disease-associated microglia. A marker-based score therefore measures microglial state as well as number.
