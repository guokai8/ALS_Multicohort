# P01 multi-cohort astrocyte consensus: analysis plan

Status: v1.0 FROZEN 2026-10-11. The co-author supports the change of study question. Open decisions 2-5 (Section 10) are adopted as proposed; decision 4 is checked at the start of Phase 1. Nothing in Sections 4-7 has been run at freeze, except where Section 2 lists it as already observed. A frozen copy is committed in the analysis repository under `docs/analysis_plans/`.
Date: 2026-10-11
Target journal: Molecular Neurodegeneration (Research Article)
Repository: `Data/P01_ALS_CNS_regulatory_map_release_candidate_2026-08-15` (W50-W54 not yet committed)

## 1. Why the study changes

The submitted-package claim ("donor-reproducible ALS motor-cortex astrocyte program with spinal RNA recurrence") rested on one 18-donor discovery cohort and one 12-donor spinal cohort. External motor-cortex checks (W50-W54) showed:

- No motor-cortex cohort cleanly refutes the modules, but every cohort carries a cohort-specific confounder: processing batch (Pineda), age (Ulm), depth (Li).
- Effect sizes are several-fold smaller than the discovery estimates (winner's curse).
- The 14-gene ALS-down module is directionally the most reproducible component. The 15-gene ALS-up module (C3, LTF, ANXA1, GPNMB, ...) is concordant only where confounders are handled.

New study question: **Which human ALS astrocyte transcriptional changes reproduce across independent cohorts once cohort-specific confounding is handled, and which features of cohort design make single-atlas astrocyte signatures unreliable?**

Working title: *Cross-cohort donor-level analysis resolves a reproducible loss of homeostatic astrocyte programs in human ALS.* The title is provisional and depends on Section 8.

## 2. What has already been observed (transparency statement)

These results are known to the analysts and will be reported in full in the Supplement. They cannot serve as confirmatory evidence for hypotheses they suggested.

| Work package | Cohort | Observed |
|---|---|---|
| W33/W34 | NYGC MTC discovery (GSE330130; 11 ALS / 7 ctrl) | 29 FDR genes in the focused model; modules defined here |
| W09/W44/W49 | Zelic 2025, GSE287257 cervical SC (8/4) | Up module replicated |
| W08/W50/W51 | Pineda 2024 BA4 | Diagnosis confounded with processing day/campaign (index-plate evidence). Up module null in every interpretable specification; down module concordant by camera |
| W50 | Pineda 2024 BA9 | Cell-level module score only: up null (P 0.50). Discovery-matched limma not run |
| W52 | Ulm/DZNE, Ruf 2026 M1 (30/32) | Both modules null raw. Controls 11.1 y older. Post hoc age-bias bound: up 8/9 and down 11/13 concordant |
| W53 | Masrori 2025 MC + SC (5/5/5) | Summary cell-level DE: up not concordant; down 13/14 MC, 12/14 SC |
| W54 | Li 2023 BA4 C9-ALS (6/6) | Up 11/11 concordant (camera 1e-4; fry 0.063); down fry 0.013. Up is depth-sensitive (post hoc) |

Consequence: the down-module hypothesis (H1) arose from these data. Confirmatory status is given only to data not yet examined (Section 4).

## 3. Hypotheses

- **H1 (primary):** The 14-gene ALS-down module is lower in ALS than in controls.
- **H2 (primary):** The 15-gene ALS-up module is higher in ALS than in controls.
- **H3 (secondary):** A discovery-agnostic genome-wide meta-analysis across cohorts identifies consensus astrocyte genes that survive leave-one-cohort-out.
- **H4 (secondary):** Within ALS, module scores are associated with disease duration. This test is two-sided, since no direction is prespecified.
- **H5 (secondary, design):** Diagnosis-batch confounding, age imbalance and depth imbalance measurably distort per-cohort astrocyte DE.

Module definitions are frozen as the W34 seed lists:
- `results/W34_synapse_resolved_regulatory_reanalysis/inputs/rna_seed_lists/GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_up.tsv`
- `results/W34_synapse_resolved_regulatory_reanalysis/inputs/rna_seed_lists/GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_down.tsv`

No gene is added or removed.

## 4. Confirmatory analysis: NYGC bulk hold-out (never examined)

### Data

NYGC ALS Consortium bulk RNA-seq is stored locally in `Other_Datasets/NYGC_ALS_Consortium`. It holds 2,574 samples (ST3), including about 568 motor-cortex (BA4, lateral, medial, unspecified) and about 903 spinal-cord (cervical, lumbar, thoracic) samples. ST1 gives donor-level age at death, sex, C9orf72 status, site of onset and disease duration; ST3 gives RIN and PMI.

### Exclusions (fixed before any outcome is computed)

1. All samples from every donor in GSE330130, MTC or LSC, identified via the W30/W33 `CGND-HRA` and NYGC donor bridge. Unmatched GSE330130 libraries are resolved by donor ID; if overlap cannot be ruled out, that donor is excluded.
2. Non-ALS/non-control groups (FTD-only, other neurological). ALS-FTD is kept as a sensitivity analysis only.
3. Samples failing the consortium QC flags in ST3.
4. Data-use check: confirm that the `data_access_restriction` field in ST1 permits this analysis before running.

### Composition handling

Bulk tissue confounds cell composition with disease. An astrocyte-abundance covariate is computed per sample as the mean z-scored log-CPM of the W54 astrocyte markers (AQP4, GJA1, SLC1A3, SLC1A2, ALDH1L1, SLC4A4). These markers contain no module genes. Neuronal, oligodendrocyte and microglial marker scores from the same list are included as covariates. Sensitivity: NNLS deconvolution against an Ulm WNN-labelled pseudobulk reference.

### Model

For each region (motor cortex, spinal cord):
- limma-voom with duplicateCorrelation, with donor as block (multiple subregions per donor).
- Design: `~ disease + age_at_death + sex + PMI + RIN + subregion + astro_score + neuron_score + oligo_score + micro_score + sequencing_platform/prep`.
- Module tests: fry with block and correlation (primary statistic) and camera, for H1 and H2.
- Multiplicity: BH over the 4 tests (2 modules x 2 regions).
- Success criterion: BH < 0.05 in the hypothesised direction.

### Clinical association (H4)

Within ALS donors who pass the exclusions:
- Response: donor-mean module score (composition-residualised).
- Model: `log(disease_duration_in_months)` ~ module score + age at onset + sex + site of onset + C9orf72 status, per region.
- Multiplicity: BH over 4 tests.
- Sensitivity: Cox model treating duration as survival time, with all donors deceased.

## 5. Consolidation analysis: all-cohort meta-analysis (partly observed; reported as such)

### 5.1 Cohort-specific primary specifications (fixed now)

| Cohort | Region | Astrocyte definition | Unit | Design | Notes |
|---|---|---|---|---|---|
| Zelic GSE287257 | Cervical SC | Author annotated subset (W09) plus W54 marker rule | Donor | `~ disease + age + sex` (PMI/RIN by W33 rule) | Age, sex, PMI, RIN from `GSE287257.xlsx` |
| Li 2023 GSE219280 | BA4 | W54 marker rule | Donor | `~ disease + age + sex` | C9 only |
| Ulm, Zenodo 18370646 | M1 | Author WNN_L2 | Donor x chip, blocked | `~ disease + sex + chip` | Age not linkable; see 5.3 |
| Pineda BA4 | BA4 | Author CellType | Donor | 2019 campaign only (191112 vs 191114), `~ disease + age + sex` | Full-cohort result as sensitivity |
| Pineda BA9 | BA9 (cross-region) | Author CellType | Donor | Batch 210526 only, `~ disease + age + sex (+PMI rule)` | The only batch-balanced Pineda contrast |
| Masrori 2025 | MC + SC | Author | Summary | Direction only, not pooled | Cell-level supplementary DE |
| NYGC MTC (discovery) | MTC | Author | Donor | W34 | Excluded from pooled module estimates; included in H3 with and without |

Common pipeline: raw counts summed per unit, `filterByExpr(design)`, TMM, voom, `eBayes(trend = TRUE)`.

### 5.2 Module-level meta-analysis (H1, H2)

- Per-cohort effect: the donor-level module score is the mean voom log-CPM over module genes present, residualised on the cohort's non-disease covariates. The standardised ALS effect (Hedges' g) and its SE come from the cohort design.
- Pooling: random-effects (REML, `metafor`) over independent motor-cortex cohorts (Li, Ulm, Pineda BA4) and separately over all independent cohorts (adding Zelic SC and Pineda BA9). Report tau^2, I^2 and a leave-one-cohort-out analysis.
- Ulm enters twice, uncorrected and age-bias-corrected (5.3). The primary pooled estimate uses the corrected value with inflated SE.

### 5.3 Age-bias correction for Ulm (K2)

Per-gene age coefficients are estimated by random-effects meta across cohorts with donor ages: NYGC MTC, Pineda BA4 and BA9, Li and Zelic.

- Corrected logFC = observed - beta_age x delta_age, with delta_age = -11.1 years from `Multiome_Dataset_Samples_Summary.txt`.
- Variance = SE_obs^2 + delta_age^2 x SE_beta^2.
- Sensitivity: beta_age estimated in donors aged 70 or older only, which probes non-linearity at the NBB control ages of 82-96.

### 5.4 Discovery-agnostic genome-wide consensus (H3)

- Per-cohort limma logFC and SE (SE = logFC / t) for all genes.
- Random-effects meta for genes tested in at least 4 cohorts.
- Consensus gene: meta BH < 0.05, the same direction in at least two-thirds of cohorts, and BH < 0.10 in every leave-one-cohort-out refit.
- Run with and without the discovery cohort, and with motor cortex only.

### 5.5 Cross-cohort reproducibility matrix (H5)

For each cohort, take the top 200 genes by P and compute the fraction replicated in every other cohort (same sign, P < 0.05). Compare with 1,000 expression- and detection-matched random gene sets.

### 5.6 Confounder audit (H5)

For each cohort, tabulate diagnosis x batch, the age difference, the depth ratio, the PMI and C9 composition, and the annotation method. Demonstrate how the effect changes:
- Pineda full cohort vs within-campaign;
- Ulm uncorrected vs age-corrected;
- Li with vs without depth.

## 6. Exploratory analyses (labelled exploratory in the manuscript)

- E1 ATAC: for consensus genes and modules, test NFI motif chromVAR and the C3-locus peak in Ulm ATAC (Zenodo raw ATAC counts) and Li snATAC (GSE219279; fragments partly local).
- E2 Spatial: astrocyte-rich spot colocalisation of the modules in GSE330224 (Barretto, same resource; an orthogonal modality, not independent donors).
- E3 Biology of the down module: pathway context (adhesion, synaptic support), relation to GFAP-pos vs GFAP-neg astrocytes (Pineda subtypes) and to reactive-state scores.
- E4 ALS-FTD and C9 moderators: genotype as a meta-regression moderator.

## 7. Outputs

Figures:
1. Cohort map and confounder audit.
2. Reproducibility matrix and winner's-curse shrinkage.
3. Module meta-analysis forest plots.
4. NYGC bulk hold-out and clinical association.
5. Consensus genes and down-module biology.

Supplement:
- every W08-W54 result, including the superseded W08 seed-list tests;
- per-cohort design tables;
- per-gene meta table.

## 8. Decision rules for framing and venue (agreed before Section 4 is run)

| Outcome | Framing | Venue |
|---|---|---|
| NYGC hold-out confirms H1 (BH < 0.05) in at least one region, and the 5.2 pooled H1 estimate is significant | Reproducible loss of homeostatic astrocyte program; reactive/complement program smaller and confounder-sensitive | Molecular Neurodegeneration |
| Hold-out confirms H2 but not H1 | Reactive/complement program reproducible after confounder control | Molecular Neurodegeneration |
| Hold-out null for both, but 5.4 yields robust consensus genes | Consensus genes as the main result, modules as a cautionary example | Molecular Neurodegeneration or Genome Medicine |
| No hold-out support and no consensus | Reproducibility and design-pitfall paper | Genome Biology, eLife or Cell Genomics |

H4 (clinical) is reported whatever the outcome. A positive H4 strengthens any row.

## 9. Work plan

| Phase | Task | Estimate |
|---|---|---|
| 0 | Author sign-off on this plan; commit W50-W54; freeze this file (hash in commit) | 1 day |
| 1 | NYGC bulk: data-use check, donor-overlap exclusion, composition scores. **Run Section 4 first, before any further cohort-level work** | 3-4 days |
| 2 | Harmonise cohorts per 5.1 (Zelic re-run with W54 rule; Pineda BA9 limma) | 3 days |
| 3 | Sections 5.2-5.6 meta-analyses | 3-4 days |
| 4 | Exploratory E1-E4 | 1-2 weeks |
| 5 | Manuscript rewrite (MN format), new figures, supplement | 2 weeks |

## 10. Open decisions for the authors

1. Agreement from the co-author (Zhihan Wang) on the change of study question and on this plan.
2. Whether C9-only cohorts (Li) enter the primary pool or only a sensitivity analysis. Proposal: primary pool, with genotype as moderator.
3. Whether ALS-FTD donors are excluded from the primary analyses. Proposal: yes, as in discovery.
4. NYGC data-use permission for the bulk expression analysis.
5. Whether the success criterion in Section 4 should require both regions instead of at least one.

## 11. Main risks

- NYGC bulk composition adjustment may be imperfect; the marker-score and deconvolution approaches must agree.
- Residual donor overlap between NYGC bulk and GSE330130 would invalidate the hold-out, so it is resolved by donor ID before outcomes are computed.
- Age-bias correction for Ulm is model-based and cannot be verified at the donor level.
- Pooled cohorts differ in nuclei isolation, chemistry and annotation. Heterogeneity is reported, not hidden.
