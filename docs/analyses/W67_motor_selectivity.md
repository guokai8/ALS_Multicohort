# W67: motor-system selectivity in GSE330224 Visium (Amendment 10 B)

**Status.** Exploratory, but the hypothesis and reading rule were fixed in Amendment 10 (`ed60da9`) before the analysis was run.

**Scripts.**
- `scripts/build_w67_annotation_pseudobulk.py`
- `scripts/run_w67_motor_selectivity.R`

**Design.**
- Donor-level pseudobulk per territory, from independent donors (LSC 18 ALS / 5 control; MTC 22 / 5 with age).
- Model: limma-voom `~ territory * disease + sex + age`, blocked on donor.
- The interaction asks whether the ALS change is larger in motor than in sensory territory.

## Set tests (`w67_motor_selectivity_tests.tsv`; fry P, median log2FC)

| Comparison | Gene set | Interaction (motor − sensory) | ALS in motor territory | ALS in sensory territory | Reading |
|---|---|---|---|---|---|
| LSC ventral vs dorsal horn | **Frozen ALS-up** | **P 0.012 (BH 0.047), Up, 0.48** | P 0.0495, Up, 0.33 | P 0.84, 0.02 | **Motor-system selective** |
| | Consensus up | P 0.15 (BH 0.31), 0.09 | P 4e-4, 0.21 | P 0.011, 0.14 | Not motor-selective (diffuse) |
| | Consensus down | P 0.25 | P 0.029, -0.11 | P 0.20, -0.02 | (secondary) |
| LSC lateral funiculus vs dorsal columns | Frozen ALS-up | P 0.42 (BH 0.42) | P 0.74, -0.03 | P 0.40, -0.12 | Not resolved (no ALS increase in either) |
| | Consensus up | P 0.24 (BH 0.32) | P 0.0015, 0.18 | P 0.032, 0.11 | Not motor-selective (diffuse) |
| | Consensus down | P 0.50 | P 0.041, -0.12 | P 0.013, -0.12 | (secondary) |
| MTC layer 5 vs layers 2-3 (secondary) | Frozen ALS-up | P 0.77 | P 0.36 | P 0.29 | Not resolved |
| | Consensus up | P 0.70 | P 0.26 | P 0.20 | Not resolved |

## Gene drivers, ventral vs dorsal horn (`w67_module_gene_results.tsv`; log2FC ALS vs control)

| Gene | Ventral horn | Dorsal horn | Interaction P |
|---|---|---|---|
| LTF | **2.54** (P 0.005) | 0.74 | 0.008 |
| ANXA1 | **1.10** (P 0.021) | 0.21 | 0.024 |
| GPNMB | 0.68 (P 0.077) | 0.06 | 0.050 |
| EMP3 | 0.61 (P 0.039) | 0.39 | 0.35 |
| C3 | -0.10 | -0.94 (P 0.023) | 0.018 |
| ANG | 0.18 | -0.82 (P 0.030) | 0.001 |
| KITLG | 0.50 | 0.54 | 0.84 |

The interaction for C3 and ANG comes from lower expression in the ALS dorsal horn, not from higher expression in the ventral horn.

## Reading

**Frozen ALS-up module.** It rises in the ventral horn, where lower motor neurons degenerate, and not in the dorsal horn of the same sections. The increase is carried by LTF, ANXA1, GPNMB and EMP3, the same genes that drove the bulk spinal-cord signal (W55). This explains the null whole-section test (W66): the change is focal and diluted by the larger dorsal-horn and white-matter areas. The lateral funiculus shows no increase.

**Consensus signature.** It shifts in both motor and sensory grey and white matter, a diffuse response of astrocyte-compartment genes.

**Motor cortex.** No laminar selectivity is detected (5 controls).

**Caveat.** With 5 controls, these are exploratory results.
