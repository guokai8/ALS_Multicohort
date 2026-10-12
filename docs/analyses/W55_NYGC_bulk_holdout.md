# W55 NYGC bulk hold-out (confirmatory, plan v1.0 Section 4 + Amendment 1)

Date: 2026-10-11

Scripts:
- `scripts/build_w55_nygc_holdout_inputs.py`
- `scripts/run_w55_nygc_bulk_holdout.R`

Plan and amendment: `docs/analysis_plans/`. Amendment 1 was committed in `02f1136` before any outcome was computed.

## Inputs

- NYGC ALS Consortium bulk RNA-seq (local `Other_Datasets/NYGC_ALS_Consortium`). The modules had never been tested in these data.
- Data use: ST1 consent codes (GRU, HMB, neurologic DS) all permit ALS research.
- Discovery overlap (`w55_discovery_overlap_exclusions.tsv`): 13 GSE330130 subjects matched exactly; the 14 unmatched subjects were handled with a conservative fingerprint (diagnosis class, sex, age at death ±1). In total 344 NYGC donors were excluded.
- Retained after QC (`used_in_analysis`) and exclusion: motor cortex 201 ALS / 23 control donors (306 samples); spinal cord 218 ALS / 46 control donors (514 samples).
- No ALS-FTD sample survived QC and exclusion, so the planned ALS-FTD sensitivity could not be run.

## Primary confirmatory results (`w55_primary_module_tests.tsv`)

The model is complete case: `~ disease + age + sex + rin + platform + subregion + pmi + astro/neuron/oligo/micro scores`, with donor-blocked voom and blocked fry. BH is applied over 4 tests.

| Region | Module | Samples | Donors (ALS/ctrl) | fry P | fry BH | Direction | camera P | Genes concordant | Median logFC | Confirmed |
|---|---|---|---|---|---|---|---|---|---|---|
| Motor cortex | up | 175 | 90/12 | 0.0040 | 0.0079 | Up | 1e-4 | 12/14 | +0.25 | **Yes** |
| Motor cortex | down | 175 | 90/12 | 0.95 | 0.95 | Up (wrong) | 0.83 | 7/14 | -0.003 | No |
| Spinal cord | up | 279 | 107/32 | <1e-4 | <1e-4 | Up | 1e-4 | 10/13 | +0.32 | **Yes** |
| Spinal cord | down | 279 | 107/32 | 0.30 | 0.40 | Down | 0.22 | 9/14 | -0.03 | No |

**H2 (ALS-up module) is confirmed in both regions. H1 (ALS-down module) is not confirmed.**

The prespecified sensitivity without PMI (`w55_sensitivity_no_pmi_module_tests.tsv`; motor cortex 174/19 donors, spinal cord 200/46) agrees:
- up: fry P 1e-4 (motor cortex) and <1e-4 (spinal cord);
- down: fry P 0.77 and 0.20.

Up-module drivers, logFC in motor cortex / spinal cord:
- GPNMB +0.69 / +1.98
- LTF +0.53 / +1.48
- EMP3 +0.45 / +0.52
- ANXA1 +0.39 / +0.46
- C3 +0.23 / +0.45
- ANGPT1 +0.71 / +0.06
- EGF +0.39 / +0.17

## Clinical association, H4 (`w55_clinical_association.tsv`)

Within ALS donors. The response is the donor-mean, composition-residualised module score.

| Region | Module | n | log-duration beta per SD | P | BH | Cox HR per SD | Cox P |
|---|---|---|---|---|---|---|---|
| Motor cortex | up | 87 | +0.158 | 0.046 | 0.093 | 0.80 | 0.081 |
| Motor cortex | down | 87 | -0.015 | 0.84 | 0.84 | 0.97 | 0.79 |
| Spinal cord | up | 91 | -0.032 | 0.66 | 0.84 | 1.08 | 0.50 |
| Spinal cord | down | 91 | -0.203 | 0.0067 | **0.027** | 1.25 | 0.036 |

In spinal cord, a higher down-module score goes with shorter disease duration (BH 0.027). In motor cortex, a higher up-module score goes with longer duration, but BH is 0.093. Tissue is sampled at death, so duration-dependent accumulation (reverse causation) cannot be separated from prognostic effects.

## Interpretation boundary

- Bulk tissue cannot assign the signal to astrocytes. GPNMB, ANXA1, C3 and EMP3 are also expressed by activated microglia and macrophages.
- The microglial composition score uses homeostatic markers (CSF1R, P2RY12, CX3CR1, PTPRC), which fall in disease-associated microglia, so residual myeloid confounding is possible. An exploratory check with an activated-myeloid score is needed.
- Motor-cortex control donors are few after complete case (12).
- The plan's decision table (Section 8, row 2) now applies: H2 confirmed, H1 not.
