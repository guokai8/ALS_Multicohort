# Amendment 3: exploratory ambient-RNA / contamination test of the astrocyte-compartment signal

Date: 2026-10-11
Amends: plan v1.0 (sha256 da63b81a...bfb3bd), Section 6 (exploratory)
Status: committed **before** any analysis below was run. Known at this point: W55 (bulk up module confirmed) and W56 (most up genes are microglial or endothelial at baseline; no single-nucleus cell type shows a significant up shift). Only metadata (obs fields, cell-type labels, count type) was inspected for this amendment.

All analyses here are **exploratory**.

## Question

Is the ALS-up signal measured inside astrocyte-annotated nuclei explained by contaminating transcripts (ambient RNA or doublets) from microglia, endothelial cells or other lineages?

## Constraint

Only filtered count matrices are available for GSE330130, Ulm and Pineda, with no empty droplets, and DecontX, SoupX and CellBender are not installed. So the correction is applied at the pseudobulk level with a SoupX-style estimator, using no new software.

## Estimator (fixed)

For each 10x library or channel L:

1. **Ambient profile.** a_L,g is the proportion of gene g in the pooled counts of all nuclei in L. Without empty droplets this is the standard proxy for the soup.
2. **Contamination fraction.** Let c_L,g be the astrocyte pseudobulk counts in L and N_L their total. Then
   rho_L = sum over g in M of c_L,g / (N_L x sum over g in M of a_L,g),
   where M is a fixed set of genes not expressed by astrocytes:
   - microglia: P2RY12, CX3CR1, CSF1R
   - endothelial: CLDN5, VWF
   - oligodendrocyte: MOBP, MOG
   - neuron: SLC17A7, GAD1, NRGN

   No marker is a module gene. rho_L captures ambient and doublet contamination together.
3. **Corrected counts.** c*_L,g = max(0, round(c_L,g - rho_L x N_L x a_L,g)) for every gene.

## E4. Discovery cohort (GSE330130 MTC; primary for this amendment)

- **Rebuild.** Rebuild astrocyte pseudobulks (`Cell_Type == "Astrocyte"`) for the 18 W33/W34 model libraries (`used_in_w33`) from `GSE330130_MTC_SNRNASEQ_COUNTS.h5ad`. Ambient comes from all nuclei of each library.
- **Reproduction check.** Uncorrected counts with the W33 design (`~ disease_binary_w33 + age + sex + PMI + tissue_subregion`, filterByExpr, TMM, voom, eBayes trend) must reproduce the W34 seed logFCs for the 29 module genes (maximum absolute difference < 1e-6). If not, stop and report.
- **Tests.** Module fry and camera, and gene-level logFC, on:
  - (a) uncorrected counts;
  - (b) corrected counts with rho;
  - (c) corrected counts with 2 x rho (deliberate over-correction);
  - (d) corrected counts with rho estimated from microglial markers only;
  - (e) uncorrected counts with log(rho) added as a covariate.
- **Descriptive.** rho_L in ALS vs control (Wilcoxon).

Because the modules were selected in this cohort, E4 asks only whether correction removes the signal. It is not a replication.

## E5. Replication in independent cohorts

- **Ulm.** The library is the 10x well, a channel pooling about 4 donors. Astrocyte counts are corrected per donor x well, then summed to donor x chip and analysed with the W52 primary design (`~ disease + sex + chip`, duplicateCorrelation). Variants (a), (b), (c) and (e) are run.
- **Pineda BA9 batch 210526.** The library is the sample (one donor). Analysed with `~ disease + age + sex`. Variants (a), (b), (c) and (e) are run.

## Prespecified interpretation rules

- **"Robust to contamination correction" (discovery):** the E4(b) up module keeps fry P < 0.05 (Up), E4(c) keeps direction Up, and the E4(b) median up-gene logFC keeps at least 50% of the uncorrected value.
- **"Contamination-explained" (discovery):** E4(b) fry P ≥ 0.05, or the median up-gene logFC shrinks by more than 50%.
- Gene level: C3, GPNMB, LTF, ANXA1, RAMP3 and EMP3 are reported individually.
- In E5, the direction and magnitude change from (a) to (b) is reported. Ulm and Pineda were null or weak uncorrected, so they test only whether correction changes the picture.
