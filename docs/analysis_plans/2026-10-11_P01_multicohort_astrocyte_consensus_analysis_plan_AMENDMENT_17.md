# Amendment 17: are the ALS-up proteins raised in an independent human spinal-cord proteome, and is this specific to ALS? (W75)

Date: 2026-10-15
Amends: plan v1.0 (sha256 da63b81a...bfb3bd)
Status: committed **before** any expression value of PXD067060 was read. Exploratory.

## Why

All evidence for the ALS-up module so far is RNA, plus a 31-antibody protein panel that does not include the module proteins (W68). Without new tissue, a public human spinal-cord proteome is the most direct protein-level test. PRIDE PXD067060 (Feneberg et al. 2025, Acta Neuropathol Commun, doi:10.1186/s40478-025-02084-y) profiled post-mortem spinal cord and motor cortex from ALS, Parkinson's disease (PD), Alzheimer's disease (AD) and control donors. Because it includes PD and AD, it also allows a first test of disease specificity, which the paper currently lists as untested.

This amendment was prompted by the authors' request for alternatives to immunostaining when no tissue is available.

## Already observed (declared)

Before writing, only the following were read:
- the PRIDE project description and file list;
- the sample key tables and the sample sheet (sample names, tissue, fraction and diagnosis only). Spinal cord ("SPC"): urea fraction ALS 14, control 7, PD 6, AD 8; sarkosyl fraction 36 samples;
- the article's abstract and Methods.

While extracting the Methods, two Results sections were also printed:
- the section on TDP-43, α-synuclein and tau abundance;
- the section on cross-tissue comparison of soluble and insoluble fractions.

Neither names any module gene. No module-protein value from the article or the data has been seen.

## Data

- **Source.** MaxQuant output (`MQ_txt_Sark.zip`, `MQ_txt_Urea.zip`; `proteinGroups.txt`) and the sample key tables from PRIDE PXD067060.
- **Primary fraction.** The detergent-soluble (sarkosyl) fraction is primary, because it is closest to total tissue abundance. The detergent-insoluble (urea) fraction is secondary.
- **Primary tissue.** Spinal cord. Motor cortex is secondary.
- **Samples.** Pooled and quality-control samples are excluded, and so are samples whose key-table entry is missing or flagged.

## Processing (follows the original article where possible)

- **Protein filtering.** Remove reverse, potential-contaminant and only-identified-by-site protein groups. Require at least two unique peptides.
- **Quantity.** LFQ intensity if present, otherwise intensity, log2-transformed with zeros set to missing.
- **Normalisation.** Within each fraction × tissue, subtract each sample's median over the 90% of proteins with the lowest variance across samples.
- **Missing values.** Keep proteins quantified in at least 50% of the samples of that fraction × tissue. Do not impute.
- **Mapping to genes.** Proteins are matched to the frozen module genes through the `Gene names` column. Each gene is matched to the protein group whose gene list contains it with the most unique peptides.
- **Covariates.** The repository does not give age, sex or PMI, so models are unadjusted, as in the original article. This is a stated limitation.

## Tests

### A. ALS vs control, spinal cord, detergent-soluble fraction (primary)

- **Module test.**
  - *Model:* limma (`lmFit`, `eBayes` with trend) on the normalised log2 matrix, design `~ diagnosis`.
  - *Test:* fry of the ALS-up module proteins that pass the filters (requires at least 3), and camera as a secondary test.
- **Reading.**

  | Outcome | Reading |
  |---|---|
  | fry P < 0.05, direction Up, median log2FC > 0 | **Protein-level support** |
  | Median log2FC > 0 but P ≥ 0.05 | **Directionally consistent** |
  | Otherwise | **Not supported** |

  If fewer than 3 module proteins pass, the module test is not run and only gene-level results are reported.
- **Gene level.** Report the moderated t-test for each measured ALS-up protein. BH correction is applied over the measured leading proteins (GPNMB, C3, LTF, ANXA1, EMP3), and separately over all measured module proteins.

### B. Disease specificity, spinal cord, detergent-soluble fraction (secondary)

- **Contrasts.** Same model with four groups (ALS, control, PD, AD), and these contrasts:
  - ALS vs PD;
  - ALS vs AD;
  - PD vs control;
  - AD vs control.
- **Reading.**

  | Outcome | Reading |
  |---|---|
  | ALS higher than both PD and AD (fry P < 0.05, Up) | **ALS-selective among these diseases** |
  | PD or AD also higher than control (P < 0.05, Up) and ALS not higher than that disease | **Shared with other neurodegeneration** |
  | Otherwise | **Undetermined** |

### C. Secondary analyses (descriptive)

- A and B in the detergent-insoluble (urea) fraction.
- A in motor cortex, in both fractions.
- The ALS-down module in A.

## Interpretation fixed in advance

| Outcome | Interpretation |
|---|---|
| A "protein-level support" | The bulk RNA result is reflected in spinal-cord protein, in an independent brain bank |
| A "directionally consistent" | Weak protein-level support, reported as such |
| A "not supported" | The paper states that the protein-level change was not detected in this cohort |

**What this test cannot show.** Whole-tissue homogenate cannot test ventral-horn selectivity or cell source. Those still require immunostaining.

**Reporting.** Every outcome is reported in Supplementary Table S2 and the supplementary information. The main text adds at most two sentences.

## Limits stated in advance

- **Sample size.** Small groups (7 controls).
- **No covariates.** The models are unadjusted.
- **Tissue sampling.** The spinal-cord level was not specified, and fractionation may separate secreted or membrane proteins (C3, GPNMB) between fractions.
- **Missing values.** Missingness in label-free data may not be at random.

## Outputs

- `results/W75_spinal_cord_proteome/`
- `scripts/run_w75_spinal_cord_proteome.R`
