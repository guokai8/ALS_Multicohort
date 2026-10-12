# W52 Ulm/DZNE independent motor-cortex validation of the final W34 astrocyte modules

Date: 2026-10-10

Scripts:
- `scripts/build_w52_ulm_astrocyte_pseudobulk.py`
- `scripts/run_w52_ulm_motor_cortex_validation.R`

Source: Zenodo 10.5281/zenodo.18370646 (CC BY 4.0; Ruf et al. 2026 Nat Commun). The data are human primary motor cortex (M1) 10x multiome from five European brain banks, which makes the cohort independent of NYGC/GSE330130. The local copy is `Data/Other_Datasets/Multiome_Dataset`.

## Design

- Astrocytes are `WNN_L2 == "Astrocytes"` (12,398 nuclei). Raw RNA counts were summed per donor x chip stratum, giving 106 strata. Donors were run on up to two chips.
- Pseudobulk library sizes equal the summed metadata `nCount_RNA` exactly.
- Donors: 30 ALS (all non-C9), 32 HC, 10 non-C9 ALS-FTD and 7 C9 ALS-FTD.
- The well level is mostly single-diagnosis, but chips 3-12 each carry both ALS and HC, so chip is modelled as a block.
- **Age, PMI and RIN are not modelled.** Supplementary Data 2 uses IDs (`ALS_1`, `Ctrl_1`) that do not map to the metadata IDs (`ALS29`, `HC1`): group sizes match, but no ordered alignment reproduces the sexes (best 76%). The NBB controls are old (82-96 years).
- All 29 module symbols are present in the features file. The analysis plan is fixed in the R script header.

## Results

| Analysis | Donors (ALS/HC) | Up: fry P | Up: camera P | Up genes concordant | Down: fry P | Down: camera P |
|---|---|---|---|---|---|---|
| **Primary: donor x chip, `~ disease + sex + chip`, duplicateCorrelation** | 30/32 | **0.619** | n/a | 5/9 | 0.475 | n/a |
| A: one per donor, `+ main_chip` | 30/32 | 0.491 | 0.374 | 6/11 | 0.676 | 0.502 |
| B: primary + QC proxies | 30/32 | 0.323 | n/a | 6/10 | 0.557 | n/a |
| C: one per donor, `~ disease + sex` | 30/32 | 0.457 | 0.300 | 4/7 | 0.365 | 0.201 |
| Context: non-C9 ALS-FTD vs HC | 10/32 | 0.438 (Down) | n/a | 2/7 | 0.311 (Up) | n/a |
| Context: C9 ALS-FTD vs HC | 7/32 | 0.090 | n/a | 7/8 (sign P 0.035) | 0.560 | n/a |

Gene level, primary model:
- GPNMB logFC 0.68, P 0.014 (one per donor: P 0.009). GPNMB is the only up gene with nominal support.
- C3 logFC 0.13, P 0.51.
- KITLG and ANGPT1 had negative logFCs.
- LTF, ANXA1, RAMP3, ANG and two lncRNAs did not pass `filterByExpr`. LTF and ANXA1 have near-zero astrocyte expression in this cohort.

## Power check

The primary genome-wide astrocyte model has 18 genes at FDR < 0.05 and 112 at P < 0.001 (about 12 expected under the null). Top hits include FLT1, PXDN, SLC11A1 and SYTL3. The cohort therefore detects ALS-associated astrocyte signal, so the null module result is not explained by a powerless dataset.

## Post hoc age-bias bound (`scripts/run_w52_post_hoc_age_bias_bound.R`)

Donor ages cannot be linked to metadata IDs: neither Supplementary Data 2 nor the row order of `Multiome_Dataset_Samples_Summary.txt` reproduces the metadata sexes (ordered agreement 16/32 HC, 14/30 ALS, 7/17 FTD, which is chance level). Group ages are known from the summary file: ALS mean 61.2, HC 72.3 (NBB controls mean 87.9). So delta = ALS - HC = -11.1 years.

Per-year age coefficients for the module genes were estimated in Pineda BA4 astrocytes (`~ disease + age + sex`, 32 ALS / 16 PN). The expected Ulm bias for each gene is age coefficient x delta.

| Module | Median discovery logFC | Median expected age bias | Median Ulm observed | Median Ulm age-corrected | Concordant after correction |
|---|---|---|---|---|---|
| Up (9 genes) | +1.95 | -0.30 | +0.11 | +0.41 | 8/9 |
| Down (13 genes) | -1.04 | +0.12 | -0.03 | -0.15 | 11/13 |

Up-module genes rise with age in Pineda (for example GPNMB +0.048/year, P 0.009), so the older Ulm controls bias the ALS-vs-HC contrast toward null. After the correction, GPNMB is +1.21, C3 +0.49, EGF +0.44 and EMP3 +0.48.

Limits of this bound:
- It is post hoc and gives no test statistic.
- It assumes a linear age effect transported from Pineda, which has its own batch structure, and extrapolates to NBB controls aged 82-96.
- Corrected effects remain about 5-fold smaller than discovery estimates.

The Ulm result is therefore reclassified from "null" to "inconclusive, with age-bias-corrected directions concordant for both modules".

## Interpretation boundary

- The 15-gene ALS-up module and the 14-gene ALS-down module do not replicate in this independent 30-vs-32-donor motor-cortex cohort under any prespecified specification.
- GPNMB is the one up gene with consistent direction and nominal support.
- Age is unadjusted, and controls include very old NBB donors. If aged controls carry higher reactive-astrocyte expression, this would bias against ALS-up genes, so the null is not definitive for genes with age-dependent expression.
- RNA only. ATAC (NFI, C3 peak) is not yet tested.
