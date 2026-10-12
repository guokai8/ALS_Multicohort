# Amendment 1: Section 4 implementation details (pre-outcome)

Date: 2026-10-11
Amends: `2026-10-11_P01_multicohort_astrocyte_consensus_analysis_plan.md` v1.0 (sha256 da63b81a...bfb3bd)
Status: written and committed **before** any NYGC bulk module or clinical outcome was computed. Only metadata fields (IDs, covariate types, missingness) were inspected.

The frozen plan does not specify the items below. They are fixed here so that Section 4 involves no outcome-dependent choices.

## A1. Discovery-donor exclusion

The GSE330130 subjects are the 27 `subject_id` values in `results/W32_takeuchi_gse330130_library_reconciliation/w32_gse330130_unique_public_library_inventory.tsv`.

1. **Exact:** exclude every NYGC sample whose ST3 `donor_id` equals a GSE330130 subject (13 subjects).
2. **Fingerprint:** for the 14 subjects with no exact ID match, exclude every NYGC donor (ST1) that matches on all three of:
   - diagnosis class (ALS or ALS-FTD vs control; "ALS Spectrum MND, Other Neurological Disorders" counts as ALS class);
   - sex;
   - age at death within ±1 year ("90 or older" matches any subject aged 89 or above).

   This is deliberately conservative: it removes some non-overlapping donors in exchange for removing all plausible overlaps.
3. The exclusion table is written to `results/W55_NYGC_bulk_holdout/w55_discovery_overlap_exclusions.tsv` before any model is fitted.

## A2. Samples, groups and regions

- QC: keep ST3 `used_in_analysis == TRUE`.
- Groups: ST3 `disease_group` ALS vs Control is primary. ALS-FTD vs Control is a sensitivity analysis.
- Motor cortex: `tissue_clean` in {Unspecified_Motor_Cortex, Lateral_Motor_Cortex, Medial_Motor_Cortex}, plus `tissue_full == Cortex_Motor_BA4`.
- Spinal cord: `tissue_clean` in {Cervical, Lumbar, Thoracic, Unspecified}_Spinal_Cord.
- Subregion covariate: `tissue_clean`.

## A3. Covariates and missing data

- `age_at_death` is numeric, with "90 or older" set to 90. Sex uses ST1 `sex_call`.
- PMI and RIN from ST3.
- Platform is ST3 `rna_platform`. `rna_prep` is collinear with platform and is not added.
- Primary analysis: complete case on all design covariates, consistent with the W33 rule because both groups keep at least 5 donors.
- Prespecified sensitivity: drop PMI (about 37% missing in ST3) and use complete case on the remaining covariates.

## A4. Composition scores

For each region:
1. Compute TMM log2-CPM (`edgeR::cpm(log = TRUE, prior.count = 2)`) on `filterByExpr` genes.
2. z-score each marker gene across samples.
3. Average the z-scores per lineage, using the W54 marker lists (astrocyte, neuron, oligodendrocyte, microglia). These lists contain no module gene.

## A5. Model and tests

- ENSG version suffixes are stripped. Duplicate ENSG rows are summed.
- Per region: `voom`, then `duplicateCorrelation(block = donor_id)`, then `voom` again with block and correlation, then `duplicateCorrelation` again, then `lmFit` with block and `eBayes(trend = TRUE)`.
- **Primary statistic:** `fry` with block and correlation, giving 4 tests (2 modules x 2 regions) with BH correction.
- `camera` is reported as secondary. It does not accept a block, so it is anticonservative with repeated donor samples.

## A6. Clinical association (H4)

- Sample module score: mean z-scored voom log-CPM of the module genes present, residualised by `lm` on subregion, RIN, platform and the four composition scores. Averaged per donor within region.
- **Primary:** among ALS donors with numeric `disease_duration_in_months`, fit `log(duration) ~ score + age_at_symptom_onset + sex + onset_site + C9`.
  - `onset_site` is coded Limb, Bulbar, or Other/Unknown.
  - C9 positive means `c9orf72_rp_pcr == "Yes"` or numeric `c9orf72_max_call >= 30`.
  - Entries of the form ">x (from onset to date last seen)" are excluded from the primary model.
- **Sensitivity:** a Cox model on duration, with ">x" entries right-censored at x.
- Multiplicity: BH over 4 tests (2 modules x 2 regions).
