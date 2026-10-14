# W74: is the ALS-up signal astrocyte-intrinsic or shared with the microglial compartment? (Amendment 16)

Exploratory. Amendment 16 was committed before this analysis (ebd506f). Its "already observed" section lists the earlier results it declared.

Script: `scripts/run_w74_within_lineage_source.R`.

**Check against earlier results.** The within-lineage tests in Ulm and Pineda reproduce W56 E2 exactly:
- Ulm: astrocytes fry P = 0.62, microglia 0.37.
- Pineda BA9: astrocytes 0.48, microglia 0.36.

## A. Discovery microglial compartment (primary)

**Data and model.** GSE330130 motor-cortex microglia, using the same libraries as the astrocyte compartment: 11 ALS and 7 control donors. W33 diagnoses and covariate rules.

| Compartment | Module | Genes tested | fry P (direction) | camera P | Median log₂FC |
|---|---|---|---|---|---|
| Astrocyte (reference; module selected here) | up | 15 | 2.5 × 10⁻⁷ (Up) | 5.8 × 10⁻²⁶ | 2.33 |
| Microglia | up | 7 | **0.021 (Up)** | 0.18 | **0.17** |
| Microglia | down | 13 | 0.11 (Down) | 0.047 | −0.31 |

**Prespecified reading: "ALS-up signal shared with the microglial compartment".** The criteria were fry P < 0.05, direction Up and median log₂FC > 0.

The microglial effect is small, and the shared part is mainly one gene.
- Gene-level values are descriptive, not prespecified (`w74_discovery_lineage_module_genes_descriptive.tsv`).
- GPNMB is up in microglia (log₂FC 1.92; astrocytes 2.97).
- C3 (0.42), ANXA1 (0.19), EMP3 (0.11) and KITLG (0.17) are close to zero in microglia, while in astrocytes they are 1.7-3.7.
- LTF, ANG and RAMP3 do not pass expression filtering in microglia.

## B. Discovery astrocyte vs microglia, same libraries (descriptive)

**Model.** `~ lineage * disease + age + sex + PMI + tissue_subregion`, library-blocked (consensus correlation 0.20). There are 36 units.

| Module | Genes | Interaction fry P | Median interaction log₂FC | Microglia / astrocyte median log₂FC |
|---|---|---|---|---|
| up | 12 | 1.6 × 10⁻⁵ (larger in astrocytes) | 1.37 | 0.07 |
| down | 14 | 0.0067 (larger in astrocytes) | −0.58 | — |

The prespecified caveat applies: the modules were selected in the astrocyte compartment, so B is biased toward astrocytes. The microglial effect is far smaller than the astrocytic one. That does not argue against an astrocyte-intrinsic signal. On its own it does not establish one either.

## C. Independent cohorts, same donors (secondary)

| Cohort | Astrocytes, fry P (median log₂FC) | Microglia, fry P (median log₂FC) | Interaction, fry P (genes) | Reading |
|---|---|---|---|---|
| Ulm (cortex) | 0.62 (0.11) | 0.37 (0.09; 3 genes) | 0.89, Up (4) | Not distinguishable |
| Pineda BA9, batch 210526 | 0.48 (0.15) | 0.36 Down (−0.13; 2 genes) | 0.29, Up (9) | Not distinguishable |
| Zelic (spinal cord, 8 ALS / 4 control) | 0.012 (0.49) | 0.020 (0.22) | 0.048, Up (15) | Larger in astrocytes |

**Combined.** Stouffer of signed z gives z = 1.84, two-sided P = 0.066. The sign is positive in all three cohorts, i.e. toward astrocytes.

**Zelic.** This is the only independent cohort in which the module shifts within lineages. It shifts in both astrocytes and microglia, and more in astrocytes. The ALS-down module also falls in both (P = 0.00016 and 0.019).

## D. Discovery, all lineages (descriptive)

| Lineage | ALS-up fry P (direction) | Median log₂FC | ALS-down fry P (direction) |
|---|---|---|---|
| Astrocyte | 2.5 × 10⁻⁷ (Up) | 2.33 | 3.8 × 10⁻⁷ (Down) |
| Microglia | 0.021 (Up) | 0.17 | 0.11 (Down) |
| Oligodendrocyte | 0.055 (Up) | 0.44 | 0.0087 (Down) |
| OPC | 0.076 (Up) | 0.68 | 0.0047 (Down) |
| Excitatory neuron | 0.044 (Up) | 0.24 | 4.7 × 10⁻⁴ (Down) |
| Inhibitory neuron | 0.0040 (Up) | 0.37 | 8.1 × 10⁻⁴ (Down) |

In the discovery cortex the module moves in the same direction in every lineage. The shifts outside astrocytes are smaller (median 0.17-0.68, against 2.33). This was not expected. Two explanations cannot be separated here:
- ambient RNA from astrocytes or other cells (the stated limit; W57 corrected the astrocytes only);
- a donor-level effect shared across compartments.

C3, GPNMB and EMP3 carry most of the non-astrocytic rise (descriptive gene table); in inhibitory neurons ANXA1 also rises (1.62).

## Overall reading against the prespecified table

**Mixed.**
- **Row 1 (shared) is not cleanly met.** It reads "A positive, and B or C not larger in astrocytes".
  - A is positive.
  - B is larger in astrocytes, but B is biased by selection.
  - C combined is not significant (P = 0.066), yet the only distinguishable cohort, Zelic, is larger in astrocytes.
  - The outcome falls on both sides of the row, so it is reported as mixed. Neither interpretation is chosen.
- **Row 2 (astrocyte-associated) is not met,** because A is detected.

**What the data say.**
- In single-nucleus data the module is mainly an astrocyte-compartment signal.
- The microglial compartment carries a smaller, same-direction shift, driven mainly by GPNMB.
- Where the two can be compared in an independent cohort (Zelic spinal cord), the ALS effect is larger in astrocytes.
- In the discovery cortex every lineage shifts slightly. Ambient RNA or a donor-level effect therefore cannot be excluded as the source of part of the non-astrocytic shift.

**Consequence for the manuscript.** The tissue-level myeloid association (W55-W73) cannot be explained by the snRNA signal being microglial in origin. The module's astrocyte-compartment genes are mostly astrocytic in single-nucleus data, except GPNMB, which is up in both compartments. In tissue they track a myeloid-led state.

## Limits

- **Few genes in microglia.** Only 2-9 module genes pass filtering in microglia, and in Ulm only 4 pass in the interaction model.
- **Small Zelic cohort.** Zelic has 4 controls.
- **Ambient RNA.** Ambient RNA is not corrected outside the discovery astrocytes.

## Files

- `w74_discovery_lineage_tests.tsv` (A, D)
- `w74_discovery_astrocyte_vs_microglia.tsv` (B)
- `w74_independent_cohort_tests.tsv` (C, within-lineage and interaction)
- `w74_independent_cohort_interaction_reading.tsv` (C, signed z and per-cohort reading)
- `w74_reading.tsv`
- `w74_discovery_lineage_module_genes_descriptive.tsv` (gene-level, descriptive, not prespecified)
- `w74_run.log`
