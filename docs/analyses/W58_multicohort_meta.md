# W58 Multi-cohort meta-analysis (plan v1.0 Section 5; Amendment 4)

Date: 2026-10-11
Amendment 4 was committed in `4624cc1` before any pooled estimate was computed. Section 5 is a **consolidation** analysis: per-cohort module results were already known (plan Section 2).

Scripts:
- `scripts/build_w58_zelic_astrocyte_pseudobulk.py`: W54 marker rule on all 114,161 Zelic nuclei, 616-1,922 astrocytes per donor.
- `scripts/run_w58_multicohort_meta.R`

The six cohorts are Li 2023 BA4 (C9), Ulm M1, Pineda BA4 2019 campaign, Pineda BA9 batch 210526, Zelic cervical SC, and NYGC MTC discovery (rebuilt in W57; used in H3 only).

## A4.2 Module level (`w58_module_effects_by_cohort.tsv`, `w58_module_meta.tsv`)

Standardised effects d per cohort, up module:

| Cohort | d (SE) |
|---|---|
| Li | 1.17 (0.58) |
| Ulm, uncorrected | 0.11 (0.27) |
| Ulm, age-corrected | 0.39 (0.30) |
| Pineda BA4 2019 | 0.33 (0.40) |
| Pineda BA9 | 0.32 (0.39) |
| Zelic | 2.22 (1.16) |
| Discovery (selection cohort) | 5.39 |

Pooled age effect (5 cohorts): up module +0.025 SD per year (SE 0.0105); down module -0.009 SD per year.

| Pooled analysis | Module | k | d [95% CI] | P | I^2 | LOO max P | A4.2 criterion |
|---|---|---|---|---|---|---|---|
| **Primary: motor-cortex independent** (Li, Ulm age-corrected, Pineda BA4) | up | 3 | **0.48 [0.05, 0.92]** | **0.028** | 0% | 0.12 | **MET** |
| Primary | down | 3 | -0.75 [-1.63, 0.13] | 0.094 | 72% | 0.26 | not met |
| **Secondary: all independent** (+ Pineda BA9, Zelic) | up | 5 | **0.49 [0.12, 0.86]** | **0.009** | 0% | 0.040 | **MET** |
| Secondary | down | 5 | -1.24 [-2.50, 0.02] | 0.053 | 90% | 0.17 | not met |
| Sensitivity: motor, Ulm uncorrected | up | 3 | 0.33 [-0.12, 0.78] | 0.15 | 12% | — | not met |
| Sensitivity: all, Ulm uncorrected | up | 5 | 0.36 [-0.0001, 0.71] | 0.050 | 0% | — | not met (borderline) |

Reading:
- The up module is consistent across independent single-nucleus astrocyte cohorts (I^2 = 0), with a modest pooled effect.
- Significance depends on the Ulm age correction, which is model-based.
- The primary motor-cortex result is not robust to leaving one cohort out (max P 0.12); the 5-cohort result is (max P 0.040).
- The down module is highly heterogeneous.

## A4.4 Genome-wide consensus (`w58_genome_wide_meta.tsv.gz`, `w58_consensus_genes.tsv`)

| Run | Genes meta-analysed | Meta BH < 0.05 | Consensus (BH < 0.05, at least 2/3 same sign, LOO BH < 0.10) |
|---|---|---|---|
| (i) all six | 14,742 | 825 | 316 (186 up / 130 down) |
| **(ii) without discovery** | 12,784 | 556 | **159 (86 up / 73 down)** |
| (iii) motor cortex only | 10,279 | 452 | 111 (72 up / 39 down) |

Top discovery-agnostic consensus genes in run (ii):
- **Up:** PXDN, TNS1, SPP1, MAN1C1, CORO1C, MYRIP, SNED1, MYOF, MYLK, BCL6, ALPK1, PDXK, GALNT15, CAST
- **Down:** FKBP8, FAM171B, SF3B4, SLITRK2, SLC1A4, SLC20A1, ID4, PI4K2A, CEND1, RHBDD2

Run (i) adds, among others, ALDOC, NDRG4, DBN1, HAMP, GBP2 and ASAH1.

Frozen module genes at gene level:
- CDH2 is a consensus gene in runs (i) and (ii). ZNF462 is a consensus gene in run (ii). SLC35F1 is close.
- No up-module gene is individually consensus. For example, C3 has meta BH 0.18 in run (i) and 0.34 in run (ii).

## A4.5 Reproducibility matrix (`w58_reproducibility_matrix.tsv`)

The statistic is the fraction of a source cohort's top 200 genes that replicate (same sign, P < 0.05) in a target cohort, against an expression-matched null.
- Li, Ulm, Zelic and discovery replicate each other above the null in most pairs (empirical P ≤ 0.005 in 10 of 12 directed pairs). Replicated fractions are 12-43%.
- **Pineda cohorts do not replicate the other cohorts** (most empirical P > 0.2), but Pineda BA4 and BA9 replicate each other strongly (48-59%, P 0.001). This mainly reflects shared donors (see Correction below) as well as the Pineda-specific processing structure documented in W51.

## Correction (2026-10-12): Pineda BA4 and BA9 share donors

29 of the 30 Pineda BA4 2019-campaign donors are also in the Pineda BA9 batch-210526 contrast (16/16 ALS, 13/14 PN), because each Pineda donor contributed both regions. Consequences:
- The **secondary** all-independent module pool and the run (ii) consensus definition treat Pineda BA4 and BA9 as independent, so these donors count twice. They should be described as "five contrasts from four studies", and the secondary pooled estimate is partially non-independent.
- The **primary** motor-cortex pool (Li, Ulm, Pineda BA4) and run (iii) are unaffected.
- The high BA4-BA9 replication in the reproducibility matrix mainly reflects shared donors (paired regions), not only shared processing. This replaces the earlier reading above.
- W59 is unaffected as a test, because the signature is a fixed definition. Its derivation over-weights Pineda donors.

## Notes

- `metafor` emitted more than 50 warnings in each of the two gene-wise passes. These are typical for per-gene REML on small k: extreme variance ratios or non-convergence, with DL used as a fallback in `rma_safe`. The number of fallbacks was not tracked.
- The H3 consensus genes are **not** cell-source-validated. For example, SPP1 is a microglial or disease-associated marker, and W56 showed that compartment signals can be multi-lineage.
