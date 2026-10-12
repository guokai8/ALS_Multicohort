# W66: spatial context of the programme in GSE330224 Visium

This analysis is **exploratory**. It was prespecified in Amendment 9 (`0576362`), which was committed before any analysis was run.

## Data and scripts

**Data.** GSE330224 Visium raw counts, from the same NYGC resource as the discovery cohort.
- LSC (lumbar spinal cord): 35 donors, 137 sections.
- MTC (motor cortex): 50 donors, 183 sections.

**Primary disease tests** use donors not in GSE330130:
- LSC: 18 ALS, 5 control;
- MTC: 24 ALS, 6 control (23 ALS and 5 controls with age recorded).

**Scripts.**
- `scripts/build_w66_visium_pseudobulk.py`: spot table and donor pseudobulk.
- `scripts/run_w66_spatial_disease_tests.R`: B and C.
- `scripts/run_w66_spatial_colocalisation.py`: A and D.

## Deviations, fixed before any disease test (Supplementary Table S3)

1. **Lineage niche definitions.** Marker scores are zero in many low-count spots: the median vascular and microglia scores were 0 in 95 of 137 LSC sections. As written in the amendment, the definition would have left the astrocyte_only niche nearly empty, and the "rich" niches could have taken in whole sections. Niches were therefore redefined, before any disease test:
   - only spots with total counts >= 500 are used, the same threshold as A2;
   - "rich" also requires a score > 0;
   - "below the median" became "at or below the median".
2. **Complete cases.** In MTC, two donors without age are dropped from models with age, as the amendment specifies; the script now does this explicitly.

## B: disease tests in donor pseudobulk (`w66_set_tests.tsv`; fry; camera in brackets)

| Tissue | Unit | ALS-up module | ALS-down module | Consensus up | Consensus down |
|---|---|---|---|---|---|
| LSC | whole (primary) | P 0.61 (0.47), median log2FC -0.02, 7/14 up | 0.60 (Up) | **0.019 (2e-4)**, 62/82 up | 0.054 (**0.005**), 53/72 down |
| LSC | whole, ALS only | 0.47 | 0.72 | **0.005** | **0.026** |
| LSC | whole, all donors (incl. GSE330130) | **0.029 (0.007)**, median 0.30, 10/14 | 0.19 | **8e-4** | **0.015** |
| MTC | whole (primary) | 0.92 (Down), median 0.05, 10/15 | 0.46 (Up) | 0.31 | 0.57 (Up) |
| MTC | whole, all donors | 0.14 (camera **0.041**), median 0.27, 13/15 | 0.21 | 0.42 | 0.21 |

**Primary family** (frozen up and down modules, whole unit, both tissues): no test passes (BH 0.82-0.92).

**Region units** (LSC grey and white; MTC L1): consensus up and down are often nominally significant, for example LSC grey up P 0.006 and MTC L1 up 0.033, down 0.027. The frozen modules are not significant in any region.

**Reading (prespecified).**
- The ALS-up module is **not supported** in independent-donor spatial pseudobulk of either tissue: the direction is Up but non-significant in LSC, and Down in MTC.
- The frozen module becomes nominally significant only when GSE330130 donors are included. Part of its spatial signal therefore depends on the donors in whom it was selected, even in a different tissue.
- The discovery-agnostic consensus signature is supported in independent LSC donors, but not in MTC.
- With 5 or 6 controls, a null result has low power.

## C: niche-resolved shift (independent donors)

**LSC:**

| Niche | ALS-up fry P (median log2FC) | Consensus up P | Consensus down P |
|---|---|---|---|
| astrocyte_only | 0.28 (0.27) | **0.003** | **0.004** |
| astrocyte_rich | 0.21 (0.25) | **0.019** | **0.005** |
| vascular_rich | 0.97 (-0.02) | 0.11 (camera **0.028**) | 0.25 |
| microglia_rich | 0.87 (-0.08) | 0.077 (camera **0.007**) | 0.17 |

**MTC:** nothing is significant in any niche.

**Reading (prespecified, ALS-up module): niche not resolved** in either tissue.

**Descriptive pattern (not a prespecified reading).** In LSC, the consensus signature shifts most strongly in astrocyte-dominated spots, and the ALS-up module's median fold change is also largest there (0.27 vs -0.02 to -0.08). Both point to an astrocyte-associated component that the vascular- and microglia-rich niches do not show.

## A: location (`w66_region_share.tsv`, `w66_colocalisation_summary.tsv`)

**A1. Region share** (median over independent donors).
- LSC: most ALS-up genes have their largest share in the edge (pial/meningeal) group, 0.47-0.60 for ANG, ANGPT1, ANXA1, C3, EGF, EMP3, KITLG, LINC01094 and RAMP3. GPNMB and LINC02668 are highest in white matter.
- MTC: the largest shares are in white matter for LTF (0.84), AL133482.1, C3, RAMP3, ANG and ANGPT1, and in layer 1 (pial surface) for LINC01094 (0.58).

**A2. Spot co-localisation** (median Spearman rho over sections, spots with >= 500 counts; correlations are weak because Visium counts are sparse).
- ALS-up module score: microglia 0.26 (LSC) and 0.13 (MTC); astrocyte 0.11 and 0.03; mural/fibroblast 0.10 and 0.04; neuron -0.09 and -0.01.
- C3 with microglia: 0.22 and 0.17. GPNMB with microglia: 0.17 (LSC).
- RAMP3 is highest with endothelial (0.06 and 0.04); ANXA1 with mural/fibroblast (0.10, LSC).

## D: ligand-receptor co-detection (`w66_lr_codetection_summary.tsv`; median donor log2 observed/expected)

| Pair | LSC | MTC | ALS vs control (Wilcoxon, nominal) |
|---|---|---|---|
| C3-C3AR1 | 0.58 | 0.63 | LSC 0.094 (0.63 vs 0.26) |
| C3-LRP1 | 0.41 | 0.16 | |
| C3-CR1 | 0.21 | 0.06 | |
| KITLG-KIT | 0.12 | 0.62 | |
| EFNA5-EPHA4 | 0.17 | 0.65 | |
| ANXA1-EGFR | 0.26 | 0.17 | LSC **0.030** (0.32 vs -0.01) |
| EGF-ERBB4 | 0.06 | 0.07 | |
| ANGPT1-TEK | 0.03 | 0.09 | |

- C3 with C3AR1 and with LRP1 are detected together more often than chance in both tissues; KITLG-KIT and EFNA5-EPHA4 are as well, in cortex.
- ANGPT1-TEK and EGF-ERBB4 are not.
- Same-spot co-detection means spatial proximity at 55-µm resolution, not signalling.

## Overall reading

- The spatial data place the module genes at the **pial/meningeal surface and in white matter**, and their spot-level expression tracks **microglial** markers more than astrocytic ones. This is consistent with a multicellular glial–vascular location.
- The disease tests are mixed:
  - the frozen ALS-up module is **not reproduced** in independent spatial donors (5-6 controls);
  - the discovery-agnostic consensus signature **is** reproduced in independent spinal-cord donors, with its strongest shift in astrocyte-dominated spots.
- C3-C3AR1 is the candidate axis with the most consistent spatial co-detection.
