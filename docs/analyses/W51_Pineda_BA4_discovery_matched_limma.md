# W51 Discovery-matched limma-voom replication in Pineda BA4

Date: 2026-10-10

Scripts:
- `scripts/build_w51_pineda_astrocyte_pseudobulk.py`
- `scripts/run_w51_pineda_ba4_discovery_matched_limma.R`

W51 was added after the W50 cell-level module score was null in BA4. Its own analysis plan was fixed in the R script header before results were inspected.

## Method

Astrocyte (`CellType == "Astro"`) raw counts were summed per Pineda sample over all 36,503 genes, matching the discovery pseudobulk unit. No minimum-cell threshold was applied, matching discovery. Mouse M1 samples (`220616_BL6_MCX_*`) were excluded.

The W51 module-gene count sums equal the W50 values exactly (maximum difference 0.0 over 4,205 sample-gene rows).

The model reuses the W33 pipeline. Candidate covariates are age, sex and PMI; tissue subregion does not apply. Missing covariates are dropped only if a group falls below 5 complete-case donors. Each model then runs `filterByExpr(design)`, TMM, `voom` and `eBayes(trend = TRUE)`.

Module tests use `limma::fry` (self-contained, directional) as the primary test and `camera` (competitive) as a sensitivity, with BH across the up and down sets. Four of the 15 up genes did not pass `filterByExpr` in Pineda: RAMP3, ANG, LINC02668 and AL133482.1. The up-module tests therefore use 11 genes.

## Results

| Analysis | Donors (ALS/PN) | Up: fry P | Up: camera P | Up genes concordant | Down: fry P | Down: camera P |
|---|---|---|---|---|---|---|
| **Primary, discovery-matched `~ ALS + age + sex + PMI`** | 30/9 | **0.864** | 0.992 | 5/11 | 0.078 | 0.0004 |
| No PMI `~ ALS + age + sex` | 32/16 | 0.023 | 0.0017 | 11/11 (sign test P 0.0005) | 0.014 | 0.0001 |
| At least 20 astrocytes, discovery-matched | 30/9 | 0.864 | 0.992 | 5/11 | 0.078 | 0.0004 |

In the no-PMI model, GPNMB was the strongest single gene (logFC 1.66, P < 0.001).

## Post hoc diagnostics (`w51_post_hoc_pmi_vs_composition.tsv`)

| Analysis | Donors (ALS/PN) | Up: fry P | Up: camera P | Down: camera P |
|---|---|---|---|---|
| Same 30/9 donors, no PMI | 30/9 | 0.840 | 0.832 | 0.0007 |
| All donors, sex only | 32/16 | 0.310 | 0.092 | 0.029 |

The positive no-PMI result therefore depends on including the 7 PN donors that lack PMI, and on age adjustment. Adjusting for PMI itself does not change the result.

Disease group is nearly confounded with processing batch, inferred from the sample-folder date prefix:

| Batch | ALS | PN |
|---|---|---|
| 191112 | 16 | 0 |
| 191114 | 0 | 14 |
| 200721 | 0 | 1 |
| 201019 | 6 | 1 |
| 210910 | 10 | 0 |

Batch cannot be adjusted for, so any Pineda BA4 ALS-vs-PN estimate, positive or negative, is batch-confounded.

### Batch verification (`scripts/run_w51_post_hoc_batch_verification.py`)

Neither Pineda S01 (`donor_information`, `sample_mapping`) nor the Cell 2024 STAR Methods records a batch field. The Methods report 10x v3.1, NovaSeq sequencing, ACTIONet/batchelor batch correction and limma pseudobulk differential expression with SVA. The sample-ID date prefix was therefore tested empirically. Within a fixed diagnosis group, prefix cannot reflect disease biology, so a strong prefix effect indicates technical batch. The statistic is the variance of the top five astrocyte pseudobulk PCs (2,000 most variable genes) explained by prefix, with 5,000 label permutations.

| Within-group test | n | Prefix R2 | Perm P |
|---|---|---|---|
| BA4 ALS, all prefixes | 32 | 0.349 | 0.0002 |
| BA4 C9ALS, 191112 vs 210910 | 15 | 0.505 | 0.0004 |
| BA4 SALS, 191112 vs 201019 | 16 | 0.169 | 0.026 |
| BA4 FTLD, all prefixes | 25 | 0.153 | 0.076 |
| BA9 ALS, 210526 vs 210910 | 33 | 0.042 | 0.233 |

The prefix acts as a technical batch in BA4, and the effect holds within a fixed ALS genotype.

The sample-ID suffix (`snRNA-B9`, ...) is a 96-well index-plate position. Index-well reuse shows the campaign structure:
- The 2019 BA4 dates use disjoint wells of one plate: 191112 ALS A1-D2, 191113 FTLD D3-H5, 191114 PN E7-F9. This fits one library campaign, possibly sequenced in one pool, with each diagnosis processed on its own day.
- 200721 and 201019 share wells B2 and B4, so they are separate pools.
- 210910 reuses wells C2-C4 from 191112, so it is separate as well.

BA4 therefore spans three campaigns (2019, 2020, 2021). The within-ALS prefix effect above compares campaigns. The ALS-vs-PN contrast is mainly 191112 vs 191114, two days within the 2019 campaign. That day-level confound cannot be measured from these data. PN donors come almost entirely from 2019 (14/16), whereas 16 of the 32 ALS donors come from the 2020 and 2021 campaigns, which have no controls processed in the same campaign. BA9 is different: ALS (23), FTLD (24) and PN (15) were co-processed in batch 210526 (`w51_pineda_prefix_by_condition.tsv`), so BA9 is the only batch-balanced Pineda contrast, but it is prefrontal rather than motor cortex.

## Post hoc 2019-campaign test (`scripts/run_w51_post_hoc_campaign2019.R`)

The test was restricted to 191112 (ALS) vs 191114 (PN), which share one index plate, so cross-year campaign heterogeneity is removed. The processing-day confound remains. The plan was fixed in the script header before results were inspected; the primary model is `~ disease + age + sex`. Outputs are `w51_campaign2019_*`.

| Model | Donors (ALS/PN) | Up genes tested | Up: fry P | Up: camera P | Up genes concordant | Down: fry P | Down: camera P | Down genes concordant |
|---|---|---|---|---|---|---|---|---|
| **Primary `~ ALS + age + sex`** | 16/14 | 9 | **0.454** | 0.358 | 7/9 | 0.057 | 0.0003 | 11/13 |
| With PMI | 16/9 | 11 | 0.944 | 0.981 | 7/11 | 0.031 | 0.0001 | 12/13 |

GPNMB and LTF did not pass `filterByExpr` in the 30-donor primary model. GPNMB was the strongest up gene in the full-donor no-PMI W51 model, and that model's ALS group includes the 2020 and 2021 campaigns, which have no same-campaign controls. In the down module, CDH2 (P 0.001, BH within 29 = 0.025), SLC35F1, NOL4 and EFNA5 were the most concordant genes.

So in the cleanest available Pineda BA4 contrast, the 15-gene up module is not supported. The single positive specification (full donors, no PMI) mixed later-campaign ALS libraries with 2019 controls. The 14-gene down module is concordant by camera in every specification. The down-module result is still within a processing-day confound and should not be presented as clean replication.

## Pineda author DE lookup (`w51_pineda_S06_author_DE_lookup.tsv`)

The full Synapse download was inventoried:
- `Processed Data`: 152 complete sample folders matching S01.
- `Supplemental Materials`: S01-S11.

No file contains a batch, library-prep or sequencing-run field.

S06 holds the authors' limma + SVA astrocyte DE results, split into GFAP-neg and GFAP-pos subtypes for BA4/BA9 x sALS/C9ALS. Each table reports only about 6,100 expression-filtered genes, so S06 contains just 1 (GFAP-neg) or 4-5 (GFAP-pos) of the 15 up genes. LTF, GPNMB, RAMP3 and ANG are absent. S06 therefore cannot test the up module. Where reported, the BA4 GFAP-pos logFCs for C3, EGF and ANGPT1 were small and positive, none at padj < 0.05.

## Interpretation boundary

- The prespecified discovery-matched test does not replicate the 15-gene up module in Pineda BA4.
- The full-donor no-PMI model does, with all 11 tested genes concordant.
- The 14-gene down module is concordant in every model by camera.
- Because disease and batch are nearly confounded in this cohort, Pineda BA4 should be reported as an inconclusive external motor-cortex check, not as replication and not as refutation.
