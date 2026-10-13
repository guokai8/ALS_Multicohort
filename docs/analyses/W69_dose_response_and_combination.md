# W69: dose-response with motor-neuron loss, and cross-cohort combination (Amendment 12)

Amendment 12 (`1e4c470`) was committed before W69 A was run. W69 B combines results that were already known and is a post hoc summary.

- Script: `scripts/run_w69_dose_response_and_combination.R` (seed 20261013).

## A. Dose-response with motor-neuron loss (`w69a_dose_response_by_cohort.tsv`, `w69a_dose_response_combined.tsv`)

**Measures.**
- Motor-neuron (MN) index: mean motor-territory log2-CPM of the MN genes passing filterByExpr (6 of 10 in lumbar cord, 5 of 10 in cervical cord).
- Module response: mean log2-CPM of the frozen ALS-up genes, motor minus sensory territory.

**The index reflects loss.** It is lower in ALS than in controls in both cohorts:
- lumbar: median 7.88 vs 9.17, Wilcoxon P 0.037;
- cervical: median 8.62 vs 10.33, P 0.0033.

**Results, within ALS donors:**

| Cohort | Gene set | n ALS | Spearman rho (predicted < 0) | One-sided P | Random-set median rho | Percentile vs random |
|---|---|---|---|---|---|---|
| Lumbar | Frozen ALS-up | 18 | **+0.25** | 0.85 | 0.34 | 0.40 |
| Cervical | Frozen ALS-up | 9 | **+0.20** | 0.71 | 0.17 | 0.52 |
| Lumbar | Consensus up | 18 | +0.30 | 0.89 | 0.63 | 0.11 |
| Cervical | Consensus up | 9 | +0.17 | 0.68 | 0.47 | 0.16 |

**Reading (prespecified): not supported** for either gene set (combined Stouffer one-sided P 0.87 and 0.88).
- Within ALS donors, the motor-territory response does not grow with MN loss as indexed by these transcripts.
- Expression-matched random gene sets show similar or larger positive correlations. A positive coupling between the MN index and the response is therefore expected from composition alone.

**Descriptive observation (not a prespecified reading; `w69a_donor_values.tsv`).**
- In both cohorts, ALS donors with a reduced MN index show a positive module response, and controls show responses near zero.
- In cervical cord, the one ALS case with a control-like MN index (NEUPG876HJT, whose primary neuropathological diagnosis in the source table is Alzheimer's disease) has no response (-0.33).
- The pattern looks more like presence or absence of motor-neuron pathology than a graded relationship. With 9 to 18 ALS donors, these data cannot distinguish the two.

## B. Cross-cohort combination of the territory interaction (`w69b_interaction_combination.tsv`; post hoc summary)

| Gene set | Lumbar one-sided P (W67) | Cervical one-sided P (W68) | Stouffer Z | Combined one-sided P |
|---|---|---|---|---|
| **Frozen ALS-up** (primary) | 0.0058 | 0.020 | 3.23 | **6.1 × 10⁻⁴** |
| Consensus up (secondary) | 0.077 | 0.0047 | 2.85 | 0.0022 |

The two independent, prespecified tests point the same way, and together they give strong evidence that the frozen module rises more in spinal motor than in sensory grey matter. This summarises known results; it is not a new test.
