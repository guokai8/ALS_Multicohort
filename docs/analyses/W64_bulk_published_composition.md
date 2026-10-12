# W64: bulk hold-out against the consortium's published cell-type composition

This analysis is **exploratory**. It was prespecified in Amendment 8 (committed before it was run) and answers internal review item R-3.

- Script: `scripts/run_w64_bulk_published_composition.R` (sources the W55 helpers).
- Composition covariates: NYGC Flagship ST7 deconvolution proportions, with one reference type per tissue dropped (Excitatory in cortex, Neurons in spinal cord).
  - Motor cortex: 7 covariates.
  - Spinal cord: 10 covariates.
- Samples: W55 hold-out samples with an ST7 estimate (complete case).
  - Motor cortex: 174 samples (89 ALS, 12 control donors).
  - Spinal cord: 243 samples (106 ALS, 31 control donors).

**Designs**
- **(a)** ST7 proportions replace the four marker-based scores.
- **(b)** ST7 proportions are added to the marker-based scores.

**Implementation fixes made while running (no outcome involved)**
1. Some ST7 estimates and ST4 values are stored as text in the workbook; they are coerced to numeric.
2. Subregion levels left empty by the ST7 restriction are dropped, so the design is of full rank.

## Set tests (`w64_set_tests.tsv`; blocked fry, camera)

| Design | Region | Gene set | fry P | camera P | Concordant genes | Median log2FC |
|---|---|---|---|---|---|---|
| (a) ST7 replaces markers | Motor cortex | ALS-up | **1.1e-4** | 6.3e-5 | 13/14 | 0.31 |
| (a) | Spinal cord | ALS-up | **0.047** | 0.14 | 9/13 | 0.06 |
| (a) | Motor cortex | ALS-down | 0.74 | 0.39 | 8/14 | -0.02 |
| (a) | Spinal cord | ALS-down | 0.99 | 0.18 | 6/14 | 0.01 |
| (a) | Motor cortex | Consensus up | 0.039 | 0.024 | 51/81 | 0.05 |
| (a) | Motor cortex | Consensus down | 0.096 | 1.2e-4 | 53/73 | -0.05 |
| (a) | Spinal cord | Consensus up | 0.66 | 0.13 | 43/81 | 0.02 |
| (a) | Spinal cord | Consensus down | 0.086 | 0.021 | 44/72 | -0.04 |
| (b) markers + ST7 | Motor cortex | ALS-up | 6.7e-4 | 1.9e-4 | 13/14 | 0.19 |
| (b) | Spinal cord | ALS-up | 0.19 | 0.39 | 7/13 | 0.01 |
| (b) | Spinal cord | Consensus down | 0.025 | 0.010 | 51/72 | -0.04 |

## Reading (prespecified, ALS-up module, design (a))

**"Not explained by published composition":** fry P < 0.05 with the Up direction in both regions.

The two regions differ:
- **Motor cortex:** the signal is robust; the median fold change of 0.31 is similar to the primary W55 analysis (0.25).
- **Spinal cord:** the result is borderline (P = 0.047), and the median fold change falls from 0.32 to 0.06. It is no longer significant when both composition sets are included (design (b), P = 0.19).

So in spinal cord, most of the ALS-up signal is shared with composition captured by the consortium's deconvolution.

## Overlap with the consortium's own DE (`w64_module_genes_in_nygc_st4.tsv`; context)

Frozen module genes with padj < 0.05 and the hypothesised direction in NYGC Flagship ST4:

| Module | Motor cortex (MCX) | Cervical spinal cord (CSC) |
|---|---|---|
| ALS-up (15 genes) | 8 | 10 |
| ALS-down (14 genes) | 2 | 2 |

The up genes include LTF, KITLG, ANXA1, GPNMB, EGF, EMP3, ANGPT1 and LINC01094 in both regions, plus C3 and ANG in CSC. The down genes are ADAMTS17 and ZNF462 in MCX, and ZNF462 and ENOX1 in CSC.

ST4 used all NYGC donors, including the discovery donors, so this is not an independent test. It shows that most frozen ALS-up genes are already ALS-associated in the consortium's own bulk analysis. The W55 confirmation therefore establishes that a set derived from single-nucleus astrocyte-compartment data recovers this known bulk response. It does not establish a new bulk finding.
