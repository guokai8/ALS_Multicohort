# Amendment 2: exploratory cell-source check of the W55 ALS-up signal

Date: 2026-10-11
Amends: plan v1.0 (sha256 da63b81a...bfb3bd), Section 6 (exploratory)
Status: committed **before** any of the analyses below were run. W55 (Section 4) results are known: the ALS-up module was confirmed in NYGC bulk motor cortex and spinal cord, and the ALS-down module was not. Only cell-type label lists and marker/module overlap were inspected for this amendment.

All analyses here are **exploratory** and will be labelled as such. They do not change the W55 confirmatory result.

## Question

Is the W55 bulk ALS-up signal plausibly astrocytic, or is it explained by activated myeloid cells or vascular cells? Module genes such as GPNMB, ANXA1, C3, EMP3, ANGPT1, RAMP3 and EGF are not astrocyte-exclusive.

The W55 microglial covariate uses homeostatic markers (CSF1R, P2RY12, CX3CR1, PTPRC), which fall in disease-associated microglia.

## E1. Cell-type specificity of module genes (single-nucleus)

- **Primary source:** Pineda 2024 (all human BA4 and BA9 samples; mouse M1 samples excluded).
  - Lineages from `CellType`/`CellClass`: Astrocyte, Microglia, Endothelial, Mural/Fibroblast, Oligodendrocyte, OPC, Excitatory neuron, Inhibitory neuron. T cells are dropped as too sparse.
  - Raw counts are summed per sample x lineage and converted to CPM. Per gene, the lineage mean CPM is averaged over samples, with controls and cases pooled.
  - Astrocyte share = astrocyte mean CPM / sum of the eight lineage mean CPMs.
- **Replication:** Ulm WNN_L2 (Astrocytes, Microglia, Oligodendrocytes, OPC, Exc_Neurons, Inh_Neurons). Ulm has no vascular label.
- **Astrocyte-enriched subset rule:** an up-module gene is astrocyte-enriched if its Pineda astrocyte share is at least 0.5. The subset is fixed from E1 before E3 is run.

## E2. Disease shift by cell type (single-nucleus)

The module tests (blocked fry where applicable, plus camera) are repeated in non-astrocyte lineages using the cohort's existing primary design:
- **Ulm**: Microglia, Oligodendrocytes and OPC with the W52 primary design (donor x chip, `~ disease + sex + chip`, duplicateCorrelation). The age caveat from W52 applies.
- **Pineda BA9 batch 210526** (the only batch-balanced Pineda contrast): Astrocyte, Microglia and Endothelial with `~ disease + age + sex` and one profile per donor.

E2 asks whether the ALS-up shift is present in astrocytes and how it compares with microglia and endothelial cells.

## E3. NYGC bulk re-test with added composition covariates

These use the W55 inputs, the primary complete-case design and blocked fry.

- **E3a:** add an activated-myeloid score (CD68, TREM2, TYROBP, SPP1, CD163, ITGAX) and a vascular score (CLDN5, PECAM1, VWF, PDGFRB). Both are computed as the mean z-scored TMM log-CPM, as in W55. None of these markers is a module gene.
- **E3b:** test the astrocyte-enriched subset of the up module from E1, both with the W55 primary design and with the E3a design.
- Both regions are tested. The down module is reported for completeness.

## Prespecified interpretation rules

- **"Not explained by activated-myeloid or vascular composition":** E3a up module keeps blocked fry P < 0.05, direction Up, in both regions.
- **"Astrocyte-attributable":** the above holds, the E3b astrocyte-enriched subset keeps P < 0.05 (Up) in both regions under the E3a design, and E2 shows an Up shift in astrocytes (Pineda BA9 batch 210526, fry or camera P < 0.05).
- Otherwise the manuscript will describe the confirmed signal as a **"glial/tissue ALS-up program, with astrocyte contribution not isolated"**.
