# Amendment 9: spatial context of the reproducible programme in GSE330224 Visium (W66)

Date: 2026-10-13
Amends: plan v1.0 (sha256 da63b81a...bfb3bd); adds an exploratory, review-prompted analysis
Status: committed **before** any analysis below was run.

## What was inspected before writing

Only the following, none of which involves an outcome:
- the obs and var layout of the two GSE330224 files (raw integer counts; no deconvolution or cell-type field);
- the anatomical annotation labels and their spot counts;
- the donor, section and diagnosis counts, and age and PMI missingness;
- whether marker genes are present;
- donor overlap with GSE330130 (exact Subject ID match).

No expression summary by diagnosis has been computed. Scripts were written before this amendment, but none has been run.

| Tissue | Spots | Sections | Donors (ALS / control) | Donors also in GSE330130 (ALS / control) | Independent donors (ALS / control) |
|---|---|---|---|---|---|
| LSC (lumbar spinal cord) | 590,029 | 137 | 27 / 8 | 9 / 3 | 18 / 5 |
| MTC (motor cortex) | 746,829 | 183 | 38 / 12 | 14 / 6 | 24 / 6 |

Notes:
- "ALS" is the file's Level 1 field. It includes ALS-FTD and "ALS Spectrum MND, Other Neurological Disorders"; control is "Non-neurological control".
- 18 donors contribute to both tissues.
- In MTC, 2 donors lack age and 12 lack PMI.

## Why

External review asked for orthogonal, spatial evidence on where the reproducible programme sits, and on whether the candidate communication axes from W65 co-occur in tissue. GSE330224 comes from the same NYGC resource as the discovery cohort. Donors shared with GSE330130 are therefore excluded from every disease test (the primary analysis) and kept only in one sensitivity analysis.

## Inputs

Built by `scripts/build_w66_visium_pseudobulk.py`.

**Spot-level lineage marker scores.** Each score is log1p(10^4 × summed marker counts / spot total). No module gene is a marker, and SPP1 is excluded.
- astrocyte: AQP4, GJA1, SLC1A3, SLC1A2, ALDH1L1, SLC4A4
- microglia: P2RY12, CX3CR1, CSF1R, PTPRC, CD68, TREM2, TYROBP, CD163, ITGAX
- endothelial: CLDN5, PECAM1, VWF, FLT1
- mural/fibroblast: PDGFRB, DCN, COL1A2, LUM
- vascular: endothelial plus mural/fibroblast
- oligodendrocyte: MBP, PLP1, MOBP, MOG, ST18
- neuron: SNAP25, RBFOX3, SYT1, NRGN, SLC17A7, GAD1

**Module scores.** The same transform, applied to the frozen ALS-up and ALS-down modules and to the W59 consensus arms.

**Pseudobulk units.** Raw counts summed per donor over all of that donor's sections, giving one profile per donor in each unit:
- **whole**: all spots;
- **region**: anatomical region groups.
  - LSC: grey (Vent_Horn, Dors_Horn, Med_Grey); white (the four *_White labels); edge (Vent_Edge, Dors_Edge, Lat_Edge). Cent_Can enters the whole unit only.
  - MTC: L1 (Layer_1); L2_3 (Layer_2, Layer_3); L5_6 (Layer_5, Layer_6); WM (White_matter).
- **niche**: lineage niches, defined within each section.
  - astrocyte_rich: astrocyte score >= the section's 75th percentile.
  - vascular_rich, microglia_rich: defined the same way from the vascular and microglia scores.
  - astrocyte_only: astrocyte score >= 75th percentile, and both vascular and microglia scores below the section's median.

## W66 A: location of the module genes (descriptive; independent donors, both diagnoses)

**A1. Region share.** For each ALS-up gene: mean CPM in each region group divided by the sum across groups, computed per donor and summarised as the median over donors.

**A2. Spot co-localisation.**
- Within each section, take spots with total counts >= 500.
- Compute the Spearman correlation between the ALS-up module score (and each ALS-up gene, log1p CPM × 10^4) and each lineage score.
- Report the median correlation over sections and the fraction of sections with a positive correlation.

## W66 B: disease tests in donor pseudobulk (exploratory)

**Model.** ALS vs control; edgeR filterByExpr and TMM; limma-voom `~ disease + sex + age`, complete cases. There is one profile per donor and unit, so no donor blocking is needed. fry is the primary set test, and camera is reported.

**Gene sets.** Frozen ALS-up and ALS-down modules; W59 consensus up and down arms.

**Primary family.** Frozen up and down modules, whole unit, independent donors, both tissues: 4 tests, BH within the family.

**Sensitivity analyses** (whole unit):
- ALS (Diagnosis = "ALS") vs control only;
- `+ PMI`, complete cases;
- all donors, including GSE330130 donors (not independent).

**Region units.** Same model, independent donors; nominal P values.

**Reading for the ALS-up module, whole unit, primary:**

| Outcome | Reading |
|---|---|
| fry P (BH) < 0.05, direction Up, in a tissue | Spatial pseudobulk supports the programme in that tissue |
| Direction Up, but not significant | Directionally consistent; underpowered (5-6 controls) |
| Direction Down | Not supported in that tissue |

## W66 C: niche-resolved disease shift (exploratory)

Run the W66 B model in each niche unit (independent donors) and report fry P and the median log2FC of the ALS-up module per niche.

**Reading (per tissue):**

| Outcome | Reading |
|---|---|
| ALS-up fry P < 0.05 (Up) in astrocyte_only, and its median log2FC at least that in both vascular_rich and microglia_rich | Shift concentrated in astrocyte-dominated tissue |
| ALS-up fry P < 0.05 (Up) in vascular_rich or microglia_rich, and astrocyte_only median log2FC below the larger of the two | Shift carried by vascular or microglial niches (multicellular) |
| Any other pattern | Niche not resolved |

## W66 D: ligand-receptor co-detection (descriptive)

**Pairs** (from W65): C3-C3AR1, C3-CR1, C3-LRP1, ANGPT1-TEK, KITLG-KIT, EFNA5-EPHA4, EGF-ERBB4, ANXA1-EGFR.

**Statistic.** Within each section (spots with total counts >= 500): the observed fraction of spots in which both genes are detected, over the expected fraction (product of the two marginal detection fractions). This is summarised as log2((observed + 0.001) / (expected + 0.001)), averaged per donor.

**Reported.**
- The median over donors (co-occurrence above chance when > 0).
- ALS vs control on independent donors (two-sided Wilcoxon; nominal, exploratory).

**Reading.** Co-detection in the same 55-µm spot is consistent with spatial proximity of the two genes. It does not show signalling.

## Outputs

`results/W66_spatial_context/`, built by these scripts:
- `scripts/build_w66_visium_pseudobulk.py`
- `scripts/run_w66_spatial_disease_tests.R`
- `scripts/run_w66_spatial_colocalisation.py`

Everything is reported whatever the outcome.
