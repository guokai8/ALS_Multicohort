# W61: small-k inference for the module meta-analyses

Governed by Amendment 7 (commit 284d22c). This is a **post hoc, review-prompted** analysis. The HKSJ results were already observed during an internal simulated review, before Amendment 7 was written.

- Script: `scripts/run_w61_small_k_meta_inference.R`
- Self-check: the REML z-test rows reproduce every committed W58 pool and the W60 primary pool (max |diff| < 1e-8).

## Results (`w61_small_k_meta_inference.tsv`)

| Pool | Module | k | d | z-test 95% CI, P | HKSJ 95% CI, P | 95% PI |
|---|---|---|---|---|---|---|
| Primary motor (W58) | up | 3 | 0.48 | [0.05, 0.92], 0.028 | **[-0.38, 1.35], 0.137** | [0.05, 0.92] |
| Primary motor (W58) | down | 3 | -0.75 | [-1.63, 0.13], 0.094 | [-2.76, 1.26], 0.248 | [-2.31, 0.80] |
| Secondary all (W58) | up | 5 | 0.49 | [0.12, 0.86], 0.009 | [-0.04, 1.02], 0.062 | [0.12, 0.86] |
| Secondary all (W58) | down | 5 | -1.24 | [-2.50, 0.02], 0.053 | [-3.24, 0.76], 0.160 | [-4.12, 1.63] |
| Primary motor, Ulm uncorrected | up | 3 | 0.33 | [-0.12, 0.78], 0.152 | [-0.77, 1.44], 0.327 | [-0.20, 0.86] |
| Primary motor, Ulm >= 70 slope (W60) | up | 3 | 0.70 | [0.15, 1.26], 0.013 | [-0.43, 1.83], 0.117 | [0.14, 1.26] |

## Reading

- **Under HKSJ, no pool meets the A4.2 criterion** (the 95% CI excludes 0 in the hypothesised direction). This includes the W58 primary ALS-up pool and its W60 age >= 70 sensitivity.
- The ALS-up effect is **directionally consistent** (every independent contrast has d > 0, and I² ≈ 0). With three or five contrasts it is **not statistically established**.
- With τ² ≈ 0, the prediction interval equals the z-test CI and should not be read as extra support.
- The prespecified W58 result (z test) stays on record. The manuscript reports both and bases its wording on HKSJ.
