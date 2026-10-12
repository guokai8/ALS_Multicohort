# W59 Consensus-signature confirmatory test in the NYGC bulk hold-out (Amendment 5)

Date: 2026-10-11
Amendment 5 and the frozen signature (sha256 `9a4bd65b...1480ce`) were committed in `5605120` before this test. The script verifies the checksum before running.
Script: `scripts/run_w59_nygc_consensus_signature_holdout.R` (sources the W55 helpers).

The signature is the W58 run (ii) discovery-agnostic consensus from five single-nucleus cohorts, with SPP1 removed as a composition marker: 85 up and 73 down genes. NYGC bulk took no part in defining it.

## Results (`w59_signature_tests.tsv`)

| Test | Region | Set | Donors (ALS/ctrl) | Genes | fry P | fry BH | camera P | Median logFC | Confirmed |
|---|---|---|---|---|---|---|---|---|---|
| **1 Primary (W55 design)** | Motor cortex | up | 90/12 | 81 | 0.096 | 0.096 | 0.0099 | +0.034 | no |
| | Motor cortex | **down** | 90/12 | 73 | 0.014 | **0.025** | 9e-9 | -0.070 | **yes** |
| | Spinal cord | **up** | 107/32 | 81 | 1.2e-4 | **4.6e-4** | 2e-6 | +0.081 | **yes** |
| | Spinal cord | **down** | 107/32 | 72 | 0.019 | **0.025** | 0.0044 | -0.051 | **yes** |
| 2 Secondary (+ myeloid/vascular) | Motor cortex | up | | | 0.25 | 0.25 | 0.21 | +0.009 | no |
| | Motor cortex | down | | | 0.064 | 0.13 | 6e-6 | -0.049 | no |
| | Spinal cord | up | | | 0.21 | 0.25 | 0.19 | +0.034 | no |
| | Spinal cord | down | | | 0.014 | 0.057 | 0.0013 | -0.062 | no |
| 3 Sensitivity (158 overlap-candidate donors removed) | Motor cortex | up | 32/8 | | 0.89 | 0.89 | 0.72 | -0.017 | no |
| | Motor cortex | down | 32/8 | | 0.54 | 0.72 | 0.062 | -0.023 | no |
| | Spinal cord | up | 36/19 | | 0.0086 | **0.017** | 1.8e-4 | +0.089 | **yes** |
| | Spinal cord | down | 36/19 | | 0.0063 | **0.017** | 0.0029 | -0.099 | **yes** |

Gene-level sign concordance in the primary test (`w59_signature_sign_concordance.tsv`):

| Set | Motor cortex | Spinal cord |
|---|---|---|
| Down | 58/73 = 79.5% (binomial P 2e-7) | 53/72 = 73.6% (P 4e-5) |
| Up | 47/81 = 58% (P 0.091) | 62/81 = 76.5% (P 9e-7) |

## Prespecified reading (Amendment 5)

- **Primary criterion: MET for both sets.** The up set is confirmed in spinal cord. The down set is confirmed in both regions.
- The statement "a single-nucleus-derived cross-cohort consensus signature was validated in an independent bulk hold-out" is supported.
- **Qualifications:**
  - (2) After adding activated-myeloid and vascular scores, no set reaches BH < 0.05, although the down set keeps strong camera support. The bulk signal is therefore composition-sensitive.
  - (3) Removing every NYGC donor that could overlap with signature-defining cohorts leaves small motor-cortex groups (8 controls), where the result is null. Spinal-cord confirmation survives for both sets.
  - Bulk effect sizes are small (median |logFC| 0.03-0.10).
