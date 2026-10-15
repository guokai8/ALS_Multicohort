# W75: ALS-up module proteins in an independent human spinal-cord proteome (Amendment 17)

Exploratory. Amendment 17 was committed before any expression value was read (ec69e18).

| Item | Details |
|---|---|
| Script | `scripts/run_w75_spinal_cord_proteome.R` |
| Data | PRIDE PXD067060 (Feneberg et al. 2025, Acta Neuropathol Commun 13:175; Oxford Brain Bank), MaxQuant `proteinGroups.txt` for the detergent-soluble (sarkosyl) and detergent-insoluble (urea) fractions |
| Spinal cord, soluble fraction | 14 ALS, 7 control, 6 PD, 9 AD |
| Proteins after filtering | 1,901 (spinal cord, soluble fraction) |

## Implementation details

These were fixed before any result was seen.

- **Fraction labels.** The urea-fraction suffixes come from the sample sheet: 2/4/6/8 = control/PD/AD/ALS. The sarkosyl suffixes 1/3/5/7 were assigned from a one-to-one match of tissue × donor index with the urea fraction. All four groups matched exactly; the sarkosyl fraction has one extra AD donor.
- **Set tests.** fry and camera need complete data, so the set test used only module proteins quantified in every analysed sample. Gene-level tests (limma, moderated t) use all quantified samples.

## A. ALS vs control, spinal cord, detergent-soluble fraction (primary)

Three ALS-up proteins passed the filters: GPNMB, C3 and ANXA1. The other 12 (LTF, EMP3, EGF, ANG, ANGPT1, RAMP3, KITLG, STPG1 and four non-coding genes) were not quantified in ≥ 50% of samples, or were not detected.

| Protein | log₂FC (ALS vs control) | P | BH (leading proteins) | Quantified in |
|---|---|---|---|---|
| ANXA1 | 1.84 | 6.7 × 10⁻⁴ | 0.002 | 21 / 21 |
| C3 | 1.63 | 0.028 | 0.043 | 17 / 21 |
| GPNMB | 1.00 | 0.39 | 0.39 | 14 / 21 |

**Module test: not computable.** Only one of the three proteins was quantified in every sample. Amendment 17 runs the module test when at least 3 proteins pass the filters, which they did, but fry requires complete data. This is recorded as a deviation (Supplementary Table S3).

**Gene-level reading.**
- C3 and ANXA1 protein are higher in ALS spinal cord (BH < 0.05).
- GPNMB is higher in direction but not significant, and is sparsely quantified.
- All three measured proteins are higher in ALS; the median log₂FC is 1.63.

## B. Disease specificity, spinal cord, detergent-soluble fraction (secondary)

| Protein | ALS vs PD | ALS vs AD | PD vs control | AD vs control |
|---|---|---|---|---|
| GPNMB | 2.72 (P = 0.002) | 1.82 (P = 0.061) | −1.72 (P = 0.23) | −0.82 (P = 0.59) |
| C3 | 0.76 (P = 0.17) | 0.95 (P = 0.086) | 0.87 (P = 0.23) | 0.67 (P = 0.35) |
| ANXA1 | 0.60 (P = 0.22) | 0.83 (P = 0.056) | 1.24 (P = 0.030) | 1.01 (P = 0.050) |

**Prespecified reading: undetermined.** The module test could not be computed.

**Descriptive gene-level pattern.**

| Protein | Pattern | Reading |
|---|---|---|
| GPNMB | Higher in ALS than in PD and nominally than in AD; not raised in PD or AD | ALS-selective |
| ANXA1 | Also raised in PD and AD | Shared with other neurodegeneration |
| C3 | Did not separate ALS from PD or AD | Not distinguishable |

## C. Secondary (descriptive)

- **Spinal cord, insoluble fraction.** C3 (log₂FC 1.66, P = 0.019) and ANXA1 (1.63, P < 1 × 10⁻⁴) were higher in ALS than in controls. GPNMB was not quantified. Neither protein separated ALS from PD or AD.
- **Motor cortex.** Only ANXA1 was quantified, and it was not different in the soluble fraction (log₂FC 0.61, P = 0.18).
- **ALS-down module.** Only EPB41L2 and CDH2 were quantified in spinal cord. Neither was lower in ALS vs control.

## Interpretation

In an independent brain bank, two of the three measurable leading proteins of the ALS-up module (C3 and ANXA1) are higher at the protein level in ALS spinal cord, and the third (GPNMB) is higher in direction.
- **GPNMB** separated ALS from PD and, nominally, from AD.
- **ANXA1 and C3** were also raised, or not distinguishable, in PD and AD.

This fits a response that is partly shared with other neurodegenerative diseases and partly ALS-selective. The homogenate cannot test ventral-horn selectivity or cell source. The groups are small, and the models are unadjusted, because covariates are not available in the repository.

## Files

- `w75_sample_counts.tsv`
- `w75_module_tests.tsv`: set tests (not computable) and median log₂FC
- `w75_module_protein_results.tsv`: gene level, all fractions, tissues and contrasts
- `w75_reading.tsv`
- `w75_run.log`
