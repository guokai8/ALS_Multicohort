# W60: age >= 70 sensitivity of the W58 Ulm age correction

Governed by Amendment 6 (`docs/analysis_plans/2026-10-11_P01_multicohort_astrocyte_consensus_analysis_plan_AMENDMENT_6.md`, commit 896ce41), which restores the plan 5.3 sensitivity. This is a **sensitivity analysis** of the W58 consolidation result, run after W58 results were known.

- Script: `scripts/run_w60_age70_sensitivity.R` (sources `scripts/w58_cohort_fits.R`)
- Self-checks: the full-age refit reproduces the committed W58 Pineda age coefficients (max |diff| 4.9e-17), and the W58 five-cohort age estimate reproduces the committed W58 pools.

## Age coefficients (`w60_age_coefficients.tsv`)

Model `score ~ disease + age + sex` on the W58 donor-level module scores. age_std is per year, in residual-SD units.

| Cohort | Ages | Donors (ALS/ctrl) | Up age_std (SE) | Down age_std (SE) |
|---|---|---|---|---|
| Pineda BA4 2019 | all | 16/14 | 0.038 (0.019) | -0.024 (0.019) |
| Pineda BA4 2019 | >= 70 | 7/9 | 0.070 (0.060) | -0.007 (0.060) |
| Pineda BA9 210526 | all | 23/15 | 0.034 (0.017) | -0.015 (0.017) |
| Pineda BA9 210526 | >= 70 | 7/10 | 0.081 (0.057) | 0.074 (0.057) |

## Pooled age effect and corrected Ulm effect (`w60_age_effect_pool.tsv`, `w60_ulm_corrected_effects.tsv`)

Raw Ulm d: up 0.111, down -0.068.

| beta_age estimate | Up beta_age (SE) | Up Ulm d corrected (SE) | Down beta_age (SE) | Down Ulm d corrected (SE) |
|---|---|---|---|---|
| W58 primary (five cohorts, all ages) | 0.025 (0.010) | 0.39 (0.30) | -0.009 (0.010) | -0.17 (0.30) |
| Reference: Pineda BA4 + BA9, all ages | 0.035 (0.013) | 0.51 (0.30) | -0.019 (0.013) | -0.28 (0.30) |
| **Primary sensitivity: Pineda BA4 + BA9, >= 70** | 0.076 (0.041) | **0.96 (0.53)** | 0.036 (0.041) | 0.33 (0.53) |
| Secondary: Pineda BA4 only, >= 70 | 0.070 (0.060) | 0.89 (0.72) | -0.007 (0.060) | -0.15 (0.72) |

## Re-pooled module meta-analyses (`w60_module_meta.tsv`)

| beta_age estimate | Pool | Module | Pooled d [95% CI] | P | I² | Max LOO P | Criterion |
|---|---|---|---|---|---|---|---|
| W58 primary | Primary motor | up | 0.48 [0.05, 0.92] | 0.028 | 0% | 0.12 | met |
| **>= 70, BA4 + BA9** | **Primary motor** | **up** | **0.70 [0.15, 1.26]** | **0.013** | 1% | 0.12 | **met** |
| >= 70, BA4 + BA9 | Secondary all | up | 0.63 [0.19, 1.07] | 0.005 | 0% | 0.027 | met |
| >= 70, BA4 only | Primary motor | up | 0.65 [0.06, 1.25] | 0.031 | 1% | 0.19 | met |
| Reference, all ages | Primary motor | up | 0.55 [0.11, 0.99] | 0.014 | 0% | 0.12 | met |
| >= 70, BA4 + BA9 | Primary motor | down | -0.67 [-1.83, 0.48] | 0.25 | 76% | 0.67 | not met |
| >= 70, BA4 + BA9 | Secondary all | down | -1.20 [-2.61, 0.21] | 0.096 | 89% | 0.25 | not met |

## Reading (prespecified in Amendment 6)

- **ALS-up module:** the primary pooled effect keeps its sign with a 95% CI excluding 0 when beta_age is estimated in donors aged 70 or older. By the Amendment 6 rule, the W58 primary result is robust to estimating the age effect in the oldest donors.
- The age slope in donors aged 70 or older is steeper than the all-age estimate, so the correction becomes larger, not smaller. The data give no sign that a linear extrapolation over-corrects the Ulm controls. The restricted slopes are imprecise (16-17 donors per cohort), and the corrected Ulm effect has a wide SE (0.53).
- What does not change: the primary pool is still not robust to leaving out single cohorts (max LOO P 0.12). The conclusion still depends on correcting the Ulm age imbalance, now under either age model.
- **ALS-down module:** still heterogeneous and the criterion is not met under the primary sensitivity. The down-module secondary pool meets the criterion only under two non-primary estimates (BA4 only, and the all-age reference); this is not interpreted.

## Limitations

- Only Pineda BA4 and BA9 were eligible, and they share 29 of 30 donors. The secondary BA4-only estimate addresses this.
- The gene-level (H3) age correction was not recomputed: A4.3 needs at least three cohorts per gene and only two are eligible (Amendment 6).
