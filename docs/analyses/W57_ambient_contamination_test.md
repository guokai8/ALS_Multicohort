# W57 Exploratory ambient-RNA / contamination test (Amendment 3)

Date: 2026-10-11
Amendment 3 was committed in `bbd5e98` before any analysis ran. **Exploratory.**

Scripts:
- `scripts/build_w57_ambient_inputs.py`
- `scripts/run_w57_ambient_contamination_test.R`

## Method (as prespecified)

This is a pseudobulk-level, SoupX-style correction without empty droplets:
- The ambient proxy is the pooled counts of all nuclei per library: the GSE330130 library, the Ulm 10x well (about 4 donors) or the Pineda sample.
- rho is estimated from 10 genes that astrocytes do not express (P2RY12, CX3CR1, CSF1R, CLDN5, VWF, MOBP, MOG, SLC17A7, GAD1, NRGN).
- Corrected counts are max(0, c - k x rho x N x ambient proportion).

Variants:
- (a) uncorrected;
- (b) k = 1;
- (c) k = 2;
- (d) rho from microglial markers only (discovery);
- (e) log(rho) as a covariate.

## Integrity

- **Discovery reproduction check passed.** The rebuilt uncorrected astrocyte pseudobulk, under the W33 design, reproduces the W34 seed logFC of all 29 module genes (maximum |diff| 4.9e-15).

## Deviations and caveats found while running

1. One Ulm unit has rho = 0, so log(rho) is undefined. All cohorts' variant (e) uses log(rho + half the minimum positive rho). This deviates from Amendment 3.
2. **The estimator's assumption fails in Ulm.** Median rho is 0.95 and the maximum is 3.1, but a contamination fraction above 1 is impossible. Ulm astrocyte nuclei carry non-astrocyte marker transcripts at about the pooled-well proportion: the multiome data have no ambient removal, and the WNN labels are partly ATAC-driven. Ulm variants (b) and (c) are therefore not interpretable as contamination corrections and are reported only for completeness. Pineda (median rho 0.12) and discovery (0.12-0.14) are in a plausible range.
3. In R, `pmax(0, matrix)` drops dimensions; this was fixed to `pmax(matrix, 0)`. A NULL annotation column silently passed an `anyNA` check; the check now also requires marker counts to match. Both bugs were caught before any reported result was produced.

## Results (`w57_module_tests.tsv`, ALS-up module)

| Cohort and variant | Genes tested | fry P | camera P | Median logFC | Genes Up |
|---|---|---|---|---|---|
| Discovery (a) uncorrected | 15 | 2.5e-7 | 6e-26 | 2.33 | 15/15 |
| Discovery (b) rho | 15 | 3.0e-7 | 2e-22 | 2.81 | 15/15 |
| Discovery (c) 2 x rho | 13 | 3.1e-7 | 2e-20 | 2.83 | 13/13 |
| Discovery (d) microglial rho | 12 | 4.2e-5 | 3e-18 | 3.11 | 12/12 |
| Discovery (e) log(rho) covariate | 15 | 2.0e-7 | 1e-25 | 2.23 | 15/15 |
| Pineda BA9 210526 (a) | 9 | 0.48 | 0.26 | 0.15 | 8/9 |
| Pineda BA9 210526 (b) | 8 | 0.72 | 0.47 | 0.18 | 7/8 |
| Pineda BA9 210526 (c) | 8 | 0.84 | 0.37 | 0.05 | 5/8 |
| Pineda BA9 210526 (e) | 9 | 0.35 | 0.14 | 0.24 | 8/9 |
| Ulm (a) | 9 | 0.62 | — | 0.11 | 5/9 |
| Ulm (e) | 10 | 0.18 | — | 0.11 | 6/10 |

Discovery gene level, uncorrected logFC then corrected (b) and (d):
- C3 3.65 -> 3.89 / 6.20
- GPNMB 2.97 -> 4.45 / not tested
- LTF 4.56 -> 4.67 / 3.76
- ANXA1 3.57 -> 3.81 / 5.49
- RAMP3 2.97 -> 3.28 / not tested
- EMP3 2.13 -> 2.13 / 2.75

rho by disease (`w57_rho_disease_test.tsv`):

| Cohort | ALS median | Control median | Wilcoxon P |
|---|---|---|---|
| Discovery | 0.119 | 0.135 | 0.66 |
| Pineda | 0.195 | 0.152 | 0.37 |
| Ulm | 0.84 | 1.03 | 0.010 |

## Prespecified interpretation (discovery)

**"Robust to contamination correction": MET.**
- (b) up module: fry P 3e-7, direction Up.
- (c) direction Up.
- (b) median logFC retains more than 100% of the uncorrected value.

Contamination does not differ by disease in discovery, so a disease-associated ambient load cannot explain the discovery up signal.

## Interpretation boundary

- Discovery is the selection cohort, so this shows only that the correction does not remove the signal. It is not replication.
- The ambient proxy (pooled nuclei) approximates the true soup and is clearly violated in Ulm.
- This does **not** establish astrocyte origin. W56 showed that most up genes are expressed mainly by microglia and endothelial cells at baseline. W57 argues only against a simple ambient/doublet artefact in the discovery astrocyte compartment.
