# W54 Li et al. 2023 independent motor-cortex validation

Date: 2026-10-10

Scripts:
- `scripts/w54_marker_annotation.py`: fixed astrocyte rule
- `scripts/run_w54_annotation_rule_check_ulm.py`: rule check
- `scripts/build_w54_li2023_astrocyte_pseudobulk.py`
- `scripts/run_w54_li2023_motor_cortex_validation.R`
- `scripts/run_w54_post_hoc_depth.R`

## Data

The source is GSE219280 (Li et al. 2023 *Nat Commun*, PMID 37714849): BA4 motor cortex snRNA from 6 C9-ALS, 5 C9-FTD and 6 controls. The 17 per-sample CellBender-corrected filtered h5 files (2.66 GB) were size-verified against the GEO filelist and stored in `Data/GSE219280/motor_cortex_h5/`. PMI is missing for every control, so the W33 rule drops it and the primary model is `~ disease + age + sex`.

GEO has no cell annotation. Astrocytes were called with a marker-dominance rule fixed before testing, and no module gene is used as a marker. Checked against Ulm WNN labels, the rule has precision 0.936 and recall 0.624 (`w54_annotation_rule_check_ulm.tsv`). Each sample yields 198-1,088 astrocytes.

## Prespecified results (`w54_module_set_tests.tsv`)

| Analysis | Donors | Up: fry P | Up: camera P | Up genes concordant | Down: fry P (BH) | Down: camera P | Down genes concordant |
|---|---|---|---|---|---|---|---|
| **Primary C9-ALS vs control** | 6/6 | 0.063 | 0.0001 | 11/11 (sign P 0.0005; Spearman 0.62) | 0.013 (0.026) | <0.0001 | 11/12 |
| `~ disease + sex` | 6/6 | 0.21 | 0.003 | 9/11 | 0.018 | <0.0001 | 11/12 |
| Context C9-FTD vs control | 5/6 | 0.072 | 0.002 | 10/10 | 0.052 | 0.021 | 9/12 |

No gene reached genome-wide FDR < 0.05 in the 12-donor primary model, which has low power.

Up-gene logFCs: C3 +1.51 (P 0.027), EMP3 +1.01 (P 0.038), STPG1 +0.73 (P 0.013), LTF +2.14, ANXA1 +2.31 and ANGPT1 +0.97. RAMP3, ANG, LINC02668 and AL133482.1 did not pass `filterByExpr`.

## Post hoc depth sensitivity (`w54_post_hoc_depth_sensitivity.tsv`)

Median astrocyte UMI is 5,670 in C9-ALS vs 10,271 in controls. The imbalance was noted before the primary model was run, but this check was not prespecified. Adding log10 depth:
- Up module: fry P 0.39, camera P 0.039.
- Down module: fry P 0.065, camera P 1.1e-5.

With 12 donors, depth and diagnosis are nearly collinear, so this adjustment cannot separate a depth artefact from disease biology.

## Interpretation boundary

- In an independent motor-cortex cohort, both modules are directionally concordant under the prespecified model.
- The down module is robust to depth adjustment. The up module is depth-sensitive.
- The cohort is C9-ALS only and very small (6 vs 6), and the annotation is marker-based. Treat this as supportive, not as definitive replication.
