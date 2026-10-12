# W56 Exploratory cell-source check of the W55 ALS-up signal (Amendment 2)

Date: 2026-10-11
Amendment 2 (`docs/analysis_plans/..._AMENDMENT_2.md`) was committed in `dfa6954` before any of these analyses ran. **Everything here is exploratory.**

Scripts:
- `scripts/build_w56_lineage_pseudobulk.py`
- `scripts/run_w56_e1_celltype_specificity.py`
- `scripts/run_w56_e2_celltype_disease_shift.R`
- `scripts/run_w56_e3_nygc_composition_retest.R` (sources the refactored W55 helpers)

Integrity checks:
- The W56 astrocyte strata reproduce the W51 (Pineda) and W52 (Ulm) library sizes exactly.
- E2 Ulm astrocytes reproduce the W52 primary fry P values exactly (0.6191 up, 0.4753 down).
- The refactored W55 script reproduces all four committed W55 result tables byte-for-byte.

## E1: cell-type specificity (`w56_e1_module_gene_specificity.tsv`)

Astrocyte share is the astrocyte mean CPM divided by the sum of mean CPMs across 8 Pineda lineages (human BA4 and BA9). Ulm (6 lineages, no vascular label) is the replication.

| Up gene | Pineda astrocyte share | Pineda top lineage |
|---|---|---|
| AL133482.1 | 0.82 | Astrocyte |
| EGF | 0.53 | Astrocyte |
| ANGPT1 | 0.34 | Vascular/Meningeal (0.35) |
| LTF | 0.32 | Endothelial (0.48); Ulm astrocyte share 0.68 |
| LINC01094 | 0.23 | Microglia (0.75) |
| EMP3 | 0.11 | Vascular/Meningeal (0.45) |
| KITLG | 0.08 | Vascular/Meningeal (0.42) |
| ANG | 0.05 | Endothelial (0.44) |
| ANXA1 | 0.03 | Endothelial (0.64) |
| C3 | 0.025 | **Microglia (0.95)** |
| GPNMB | 0.009 | OPC/Microglia |
| RAMP3 | 0.006 | **Endothelial (0.93)** |

Only 2 of 15 up genes are astrocyte-enriched by the prespecified rule (share ≥ 0.5): EGF and AL133482.1. No down gene is astrocyte-enriched.

## E2: disease shift by cell type (`w56_e2_celltype_disease_shift.tsv`)

| Analysis | Donors (ALS/ctrl) | Up: fry P (direction) | Up genes with logFC > 0 / tested |
|---|---|---|---|
| Ulm Astrocytes (= W52) | 30/32 | 0.62 (Up) | 5/9 |
| Ulm Microglia | 30/32 | 0.37 (Up) | 2/3 |
| Ulm Oligodendrocytes | 30/32 | 0.73 (Down) | 4/5 |
| Ulm OPC | 30/32 | 0.096 (Down) | 0/3 |
| Pineda BA9 210526 Astrocyte | 23/15 | 0.48 (Up); camera 0.26 | 8/9 |
| Pineda BA9 210526 Microglia | 22/16 | 0.36 (Down) | 0/2 |
| Pineda BA9 210526 Endothelial | 23/16 | 0.36 (Down) | 0/2 |

No cell type shows a significant single-nucleus ALS-up shift. In Pineda, microglia and endothelial cells are sparse (median 166 and 26 nuclei per sample), so very few module genes pass `filterByExpr`.

## E3: NYGC bulk re-test (`w56_e3_nygc_composition_retest.tsv`)

| Analysis | Motor cortex fry P (camera) | Spinal cord fry P (camera) | Median logFC MC / SC |
|---|---|---|---|
| W55 primary (reference) | 0.0040 (1e-4) | <1e-4 (1e-4) | 0.25 / 0.32 |
| **E3a: + activated-myeloid + vascular scores** | **0.022 (0.011)** | **0.0042 (0.023)** | 0.11 / 0.15 |
| E3b: astrocyte-enriched subset, primary design | 0.019 | 0.11 | — |
| E3b: astrocyte-enriched subset, E3a design | 0.055 | 0.81 | — |

E3b tested 1 gene (EGF); AL133482.1 did not pass filtering in NYGC.

## Prespecified interpretation (Amendment 2)

- **"Not explained by activated-myeloid or vascular composition": MET.** The E3a up module stays P < 0.05, Up, in both regions, with the effect about halved.
- **"Astrocyte-attributable": NOT MET.** The E3b subset is not significant in both regions under the E3a design, and the E2 Pineda astrocyte test is not significant.
- Manuscript wording: **a tissue-level ALS-up program (C3/GPNMB/ANXA1/LTF/EMP3, ...) not explained by activated-myeloid or vascular composition scores, whose astrocyte contribution is not isolated.**

## Context for discussion (not part of the rules)

- Most up genes are expressed mainly by microglia or endothelial cells at baseline. The discovery astrocyte-compartment signal may therefore partly reflect ambient RNA or multi-lineage transcripts.
- Against a purely ambient origin: Li 2023 (W54) used CellBender ambient-corrected counts and still showed C3 +1.51 (P 0.027) in marker-dominant astrocytes, though that cohort was depth-confounded.
- The two lineage count matrices (`w56_*_lineage_counts.tsv.gz`) are not committed. `build_w56_lineage_pseudobulk.py` regenerates them.
