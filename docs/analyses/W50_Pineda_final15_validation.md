# W50 Pineda BA4 validation of the final W34 MTC astrocyte modules

Date: 2026-10-10
Script: `scripts/run_w50_pineda_final15_validation.py`

## What was tested

The final W34 Synapse/NYGC-resolved MTC astrocyte modules (15 ALS-up genes, 14 ALS-down genes) were projected into Pineda et al. 2024 astrocytes. Scoring followed W08/W09: mean log2(CPM + 1) over module genes per nucleus, averaged within sample and donor with nucleus-count weights. Samples with fewer than 20 astrocytes were dropped. Models were OLS at donor level.

The analysis plan was fixed in the script docstring before results were inspected. The primary test is BA4, ALS vs PN, `module_score ~ ALS + sex`, with BH across the up and down modules.

The fast MatrixMarket reader was checked against the W08 reference reader on sample 191112_ALS_122_snRNA-D2. The maximum absolute difference was 0.0.

## Primary result (BA4; 32 ALS vs 16 PN donors)

| Module | beta | 95% CI | P | BH |
|---|---|---|---|---|
| 15-gene ALS-up | 0.016 | -0.129 to 0.161 | 0.825 | 0.825 |
| 14-gene ALS-down | -0.501 | -0.892 to -0.109 | 0.013 | 0.027 |

**The ALS-up module did not replicate in Pineda BA4 under the prespecified primary test.** The ALS-down module replicated in direction and passed BH.

## Prespecified sensitivities, ALS-up module

- Adding age and PMI (complete case, 30 vs 9 donors): beta 0.117, P 0.33.
- Adding log10 library depth: beta 0.106, P 0.12.
- Donor pseudobulk summarization: beta 0.329, P 0.069.
- sALS vs PN: beta -0.059, P 0.47.
- C9ALS vs PN: beta 0.092, P 0.24.
- BA4/BA9 donor mean: beta 0.034, P 0.48.
- BA9 only: beta 0.035, P 0.50.
- FTLD vs PN, as context: beta 0.181, P 0.028.
- Gene level: 11 of 15 genes had a positive beta. No gene passed BH; GPNMB had nominal P 0.047.
- Leave-one-gene-out: every refit gave beta 0.002 to 0.028, all P > 0.7.

## Post hoc, labelled as such (`w50_post_hoc_depth_adjusted.tsv`)

ALS astrocytes in BA4 have lower library depth than PN astrocytes. The log10 difference is -0.125, about 25% fewer counts, P 0.021. With pseudobulk scoring plus log10 depth adjustment, the ALS-up module gave beta 0.437, 95% CI 0.067 to 0.807, P 0.022. This analysis was chosen after seeing the primary result, so it cannot replace the primary test.

## Deviation from plan

The planned "BA4 plus M1 controls" sensitivity was dropped. The five M1 PN samples (`220616_BL6_MCX_*`, donors MS1-MS5) are C57BL/6 mouse motor cortex and share no human ENSG module genes.

## Claim boundary

Pineda is RNA-only. W50 does not test ATAC, NFI or C3 locus evidence.
