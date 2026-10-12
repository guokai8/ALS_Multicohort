# Amendment 6: restore the plan 5.3 age-correction sensitivity (donors aged 70 or older)

Date: 2026-10-12
Amends: plan v1.0 (sha256 da63b81a...bfb3bd), Section 5.3; supplements Amendment 4 (A4.3)
Status: committed **before** the analysis below was run. W58 results are known, including the primary Ulm age correction and the pooled module effects. For this amendment, only donor ages and diagnoses were inspected (counts below); no module score or age coefficient was computed in the restricted donors.

## Why

Plan 5.3 specified a sensitivity analysis in which beta_age is estimated in donors aged 70 or older only, to probe non-linearity at the NBB control ages of 82-96. Amendment 4 did not carry it forward and W58 did not run it (Supplementary Table S3). This amendment restores it as a prespecified **sensitivity analysis** of the W58 consolidation result. It does not change W58 or any confirmatory test.

## Eligible cohorts (from metadata only)

| Cohort | Donors aged >= 70 (ALS / control) | Eligible |
|---|---|---|
| Pineda BA4 2019 campaign | 7 / 9 | Yes |
| Pineda BA9 batch 210526 | 7 / 10 | Yes |
| NYGC MTC discovery | 2 / 5 | No |
| Li 2023 | 2 / 2 | No |
| Zelic | 0 / 0 | No |

Eligibility rule: at least 8 donors aged >= 70 and at least 3 per diagnosis group.

## Analysis (W60)

1. **Module scores.** Use the W58 donor-level module scores unchanged: the mean voom log-CPM of module genes, with voom fitted on the full cohort.
2. **Restricted age model.** In each eligible cohort, restrict to donors aged >= 70 and fit `score ~ disease + age + sex`, which is the W58 design for both Pineda cohorts. The standardised age effect is age_std = beta_age / sigma_residual, with SE = SE(beta_age) / sigma_residual, as in A4.3.
3. **Pooling.** Pool age_std over the eligible cohorts with REML random effects (`metafor::rma`), per module.
   - **Primary sensitivity:** Pineda BA4 2019 and Pineda BA9 210526.
   - **Secondary:** Pineda BA4 2019 only, because BA4 and BA9 share 29 of 30 donors.
   - **Reference:** the same two cohorts with all ages, to separate the effect of the age restriction from the change of cohorts.
4. **Ulm correction.** For each beta_age estimate: d_Ulm,corr = d_Ulm - beta_age x (-11.1), var_corr = var + 11.1^2 x var(beta_age), as in A4.3.
5. **Re-pooling.** Repeat the W58 primary (Li, Ulm corrected, Pineda BA4 2019) and secondary (adding Pineda BA9 and Zelic) module pools with each corrected Ulm value. Report pooled d, 95% CI, P, I^2, leave-one-out and the A4.2 criterion.

## Reading (fixed now)

- If the primary pooled ALS-up d keeps its hypothesised sign with a 95% CI excluding 0 under the primary sensitivity, the W58 primary result is reported as robust to estimating the age effect in donors aged 70 or older.
- Otherwise, the manuscript states that the primary pooled result depends on a linear age effect extrapolated to the oldest controls.
- The result is reported whatever the outcome.

## Not done

- Gene-level (H3) age correction: A4.3 pools gene-level age coefficients over at least three cohorts. Only two cohorts are eligible, so the genome-wide consensus is not recomputed. This is reported as a limitation.

## Outputs

`results/W60_age70_sensitivity/` with the restricted age coefficients, the pooled beta_age estimates, corrected Ulm effects and re-pooled module meta-analyses.
