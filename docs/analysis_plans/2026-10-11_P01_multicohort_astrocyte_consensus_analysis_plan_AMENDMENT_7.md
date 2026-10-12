# Amendment 7: review-prompted checks of the meta-analysis inference and of discovery-cohort confounding

Date: 2026-10-12
Amends: plan v1.0 (sha256 da63b81a...bfb3bd); adds post hoc, review-prompted analyses
Status: committed **before** W62 was run. W61 had already been computed during an internal simulated review and is reported here as **already observed**. For W62, only library identifiers, subregion, diagnosis, age, sex, PMI and nucleus counts were inspected. No module score, fold change or set test was computed in any discovery subset.

Neither analysis is confirmatory. Neither changes W55 or W59.

## Why

An internal simulated peer review raised two issues.

1. **Small-k inference.** The W58 pools use REML with a normal (z) test. With k = 3 this is anti-conservative.
2. **Asymmetric confounder audit.** Discordant external cohorts were audited for batch, age and depth, but the discovery cohort was not audited the same way.

While preparing (2), the following was found in the library metadata (`metadata/gse330130_geo_sample_subject_bridge_2026-08-12.tsv`):
- the 18 discovery RNA libraries come from two tissue programmes, labelled by `geo_library_id`:
  - **U54-HRA**: 5 libraries, all controls, all Cortex_Motor_BA4, donors NY7xxx;
  - **CGND-HRA**: 13 libraries, 11 ALS and 2 controls.
- W38 (committed in 841cd1a) had already shown that the 29 module genes reach FDR < 0.05 only in the model that includes `tissue_subregion`. That term is largely aliased with tissue programme.

## W61: small-k inference (already observed)

- For every W58 module pool and the W60 primary age >= 70 pool, report:
  - REML with the z test (reproducing W58 and W60);
  - REML with the Knapp-Hartung (HKSJ) adjustment (`test = "knha"`);
  - the 95% prediction interval.
- The A4.2 criterion is re-evaluated with HKSJ intervals. This is reported alongside the prespecified result, not in place of it.

## W62: discovery confounder audit and tissue-programme sensitivity

**1. Balance table** (descriptive). Age, sex, PMI, nuclei per library, library size, subregion, tissue programme, and the collection site where NYGC ST1 lists one, all by diagnosis.
- Continuous variables: Wilcoxon test.
- Categorical variables: Fisher test.

**2. Programme-restricted model** (primary sensitivity).
- Data: CGND-HRA libraries only (11 ALS, 2 controls), using the W57 rebuilt discovery pseudobulk counts.
- Model: `~ disease + age + sex`. PMI and subregion are dropped because 2 controls cannot support them.
- Statistics:
  - for each of the 29 module genes, logFC and its direction;
  - for each module, fry and camera P;
  - Spearman correlation of the 29 logFCs with the W34 primary logFCs.

**3. Programme-adjusted model** (secondary sensitivity).
- Data: all 18 libraries.
- Model: `~ disease + age + sex + PMI + programme`. Disease is then identified from CGND libraries only, while the covariates use all libraries.
- Statistics: the same as in step 2.

**4. Control-source comparison** (descriptive).
- Compute per-library module scores (mean voom log-CPM over module genes) from the full-cohort voom.
- Compare CGND controls with U54 controls.

## Reading (fixed now)

Based on the programme-restricted model (step 2):

| Outcome | Reading |
|---|---|
| >= 24 of 29 genes keep the W34 direction, and both module median logFCs keep their sign | Discovery modules are not explained by tissue programme |
| < 15 of 29 genes keep their direction | Tissue programme is a plausible explanation for the discovery modules |
| Anything else | Inconclusive |

With 2 controls, set-test P values are reported but not used for this reading. The manuscript reports the outcome and the programme imbalance whatever the result.

## Outputs

- `results/W61_small_k_meta_inference/`
- `results/W62_discovery_confounder_audit/`
