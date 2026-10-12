# W66 post hoc sensitivity analyses (Amendment 10 A)

These analyses were **review-prompted and already observed**: they were computed during the internal simulated review of 2026-10-13, before Amendment 10 was written. Amendment 10 (`ed60da9`) fixes their reproducible implementation.

- Script: `scripts/run_w66_post_hoc_sensitivity.py` (seed 20261013).
- Donors: independent only (not in GSE330130).
- Spots: total counts >= 500.

## A1. Ligand-receptor co-detection against detection-matched random gene pairs (`w66_lr_codetection_empirical_null.tsv`)

The observed values reproduce W66 D exactly. Each pair is compared with 200 random gene pairs, matched to the ligand and the receptor on detection-fraction decile.

| Pair | LSC observed (random median) | LSC empirical P | MTC observed (random median) | MTC empirical P |
|---|---|---|---|---|
| C3-C3AR1 | 0.58 (0.32) | **0.015** | 0.63 (0.33) | **0.005** |
| KITLG-KIT | 0.12 (0.08) | 0.25 | 0.62 (0.37) | **0.005** |
| EFNA5-EPHA4 | 0.17 (0.08) | 0.18 | 0.65 (0.38) | **0.005** |
| C3-LRP1 | 0.41 (0.30) | 0.18 | 0.16 (0.27) | 0.80 |
| C3-CR1 | 0.21 (0.21) | 0.52 | 0.06 (0.27) | 0.96 |
| ANXA1-EGFR | 0.26 (0.37) | 0.81 | 0.17 (0.37) | 0.94 |
| ANGPT1-TEK | 0.03 (0.11) | 0.91 | 0.09 (0.07) | 0.37 |
| EGF-ERBB4 | 0.06 (0.12) | 0.81 | 0.07 (0.25) | 0.93 |

**Reading.**
- The independence expectation used in W66 D is lenient: most random pairs also exceed it.
- Against matched random pairs, only C3-C3AR1 is co-detected above chance in both tissues. KITLG-KIT and EFNA5-EPHA4 are above chance in motor cortex only.
- C3-LRP1 (LSC) no longer qualifies.

## A2. Depth-partial co-localisation (`w66_colocalisation_depth_partial.tsv`; median over sections)

| Lineage | LSC raw / depth-partial | MTC raw / depth-partial |
|---|---|---|
| Microglia | 0.26 / 0.25 | 0.13 / 0.11 |
| Astrocyte | 0.11 / 0.09 | 0.03 / 0.03 |
| Mural/fibroblast | 0.10 / 0.09 | 0.04 / 0.02 |
| Neuron | -0.09 / -0.13 | -0.01 / -0.08 |

**Reading.** Spot depth does not explain the co-localisation of the ALS-up module with microglial markers.
