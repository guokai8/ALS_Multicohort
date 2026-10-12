# W62: discovery confounder audit and tissue-programme sensitivity

Governed by Amendment 7 (commit 284d22c), committed before this analysis was run. This is a **post hoc, review-prompted** robustness check, not a replication. The ALS libraries are the same as in discovery, and the module genes were selected on the full cohort.

- Script: `scripts/run_w62_discovery_confounder_audit.R`
- Inputs: W57 rebuilt discovery astrocyte pseudobulk counts, W33 library audit, `metadata/gse330130_geo_sample_subject_bridge_2026-08-12.tsv`, NYGC ST1 (`site_sample_collected`).
- Self-check: the full W34 design reproduces the W34 seed logFCs (max |diff| 6.7e-15).
  - Camera P values equal W57's.
  - The full-design fry P is 2.2e-7 here vs 2.5e-7 in W57 for the up module (same genes, same logFC). The cause was not traced. No reading depends on it.

## Balance by diagnosis (`w62_balance_by_diagnosis.tsv`, `w62_collection_site_by_diagnosis.tsv`)

| Variable | ALS (n = 11) | Control (n = 7) | P |
|---|---|---|---|
| **Tissue programme (geo_library_id)** | CGND 11, U54 0 | **CGND 2, U54 5** | **0.0025 (Fisher)** |
| Subregion BA4 / lateral / medial | 3 / 5 / 3 | 5 / 1 / 1 | 0.23 |
| Age, median | 62 | 71 | 0.27 |
| PMI (h), median | 15.8 | 25 | 0.24 |
| Nuclei per library, median | 872 | 1,114 | 0.42 |
| Library size, median | 4.51 M | 4.67 M | 0.42 |
| Sex F / M | 5 / 6 | 5 / 2 | 0.37 |

- The five U54 libraries are all controls, all Cortex_Motor_BA4, with donors NY7xxx. None of them is listed in NYGC ST1, so their collection site is unknown.
- Listed collection sites:
  - ALS: Columbia (5), Edinburgh (2), Johns Hopkins (1);
  - CGND controls: Georgetown (1); the other is not listed.
- W38 (commit 841cd1a): the 29 module genes reach FDR < 0.05 only in the full model with `tissue_subregion`. That term largely absorbs the programme difference. Without it, no gene reaches FDR < 0.05, but 29/29 keep their direction and 27-28/29 stay nominally significant.

## Programme sensitivity (`w62_module_set_tests.tsv`, `w62_module_gene_results.tsv`, `w62_summary.tsv`)

| Model | ALS / ctrl | Genes keeping W34 direction | Spearman vs W34 | Up median logFC, fry P | Down median logFC, fry P |
|---|---|---|---|---|---|
| W34 primary (reproduction) | 11 / 7 | 29/29 | 1 | 2.33, 2.2e-7 | -1.05, 1.8e-7 |
| **CGND only, `~ disease + age + sex`** | 11 / 2 | **29/29** | **0.95** | **2.28, 0.013** | **-1.04, 0.006** |
| Programme-adjusted, all libraries | 11 / 7 | 29/29 | 0.97 | 2.75, 0.0034 | -1.04, 0.0009 |

In the CGND-only model, 23 of 29 genes have nominal P < 0.05 (up 10/15, down 13/14).

## Control-source comparison (`w62_module_scores_by_control_source.tsv`)

| Group | n | Up score, median [range] | Down score, median [range] |
|---|---|---|---|
| ALS | 11 | 2.81 [1.92, 4.80] | 6.76 [5.67, 6.99] |
| CGND controls | 2 | 0.61 [0.50, 0.72] | 7.88 [7.79, 7.97] |
| U54 controls | 5 | 1.82 [0.96, 2.53] | 7.36 [7.28, 7.45] |

The two same-programme (CGND) controls lie further from ALS than the U54 controls on both modules. So the U54 controls dilute the discovery contrast rather than create it.

## Reading (prespecified in Amendment 7)

- **Not explained by tissue programme**: 29/29 genes keep their direction, and both module median logFCs keep their sign, in the CGND-only model.
- Limits:
  - Only two same-programme controls exist.
  - The ALS libraries and the gene selection are shared with discovery, so this is a robustness check, not independent evidence.
  - Collection site is unknown for 9 of 18 libraries.
