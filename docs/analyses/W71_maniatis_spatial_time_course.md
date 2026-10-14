# W71: SOD1-G93A spatial time course and a third human spatial cohort (Amendment 14)

Exploratory analysis. Amendment 14 was committed before this analysis (dc0e51d).

Data: Maniatis et al. 2019 (Science 364:89), ALS-ST portal. Archive md5 values were checked.

Scripts:
- `scripts/build_w71_maniatis_pseudobulk.py`
- `scripts/run_w71_mouse_time_course.R`
- `scripts/run_w71_human_colocalisation.py`

## A. Mouse: SOD1-G93A vs SOD1-WT lumbar cord, ventral vs dorsal horn

The model was fitted at each age: `~ territory * genotype + sex`, with animal blocking and the fry test. The primary family was the frozen-module grey-matter interaction at the four ages (BH over four).

**Module genes tested.** Most mouse orthologs of the frozen ALS-up genes were too sparse in 100-µm spots to pass filterByExpr. The gene set tested at each age was:
- P30: C3, Kitl, Angpt1;
- P70 and P100: C3, Kitl;
- P120: C3, Kitl, Emp3, Gpnmb.

The module test in mouse is therefore driven mainly by C3. This follows from the prespecified pipeline and is not a deviation.

| Age | G93A / WT | Ventral-horn effect, fry P (median log₂FC) | Dorsal-horn effect | Interaction P (BH) |
|---|---|---|---|---|
| P30 | 10 / 6 | 0.94 (0.10) | 0.43 (−0.13) | 0.47 (0.63) |
| P70 | 8 / 4 | 0.42 (0.32) | 0.73 (0.17) | 0.75 (0.75) |
| P100 | 5 / 8 | 0.0012 (1.83) | 0.25 (0.76) | **0.043 (0.17)** |
| P120 | 9 / 6 | 5.2 × 10⁻⁵ (3.52) | 0.0049 (2.62) | 0.32 (0.63) |

**Prespecified readings.**
- Cross-species replication: **nominal**. At P100 the interaction is Up with P < 0.05 but BH ≥ 0.05.
- Timing: **symptomatic or late**. The first Up ventral-horn effect is at P100; there is none at P30 or P70.

**Gene level (ventral-horn G93A vs WT log₂FC).**
- C3: 0.10 (P30), 0.74 (P70), 3.82 (P100) and 3.55 (P120). Its interaction at P100 is 1.69 (P = 0.034).
- Gpnmb and Emp3 pass the filter only at P120: Gpnmb 5.26 in the ventral horn and 5.26 in the dorsal horn (interaction 0.00); Emp3 3.50 in the ventral horn and 2.21 in the dorsal horn (interaction P = 0.024).
- Ltf and Anxa1 were not testable.

**Microglia adjustment at P100.** The interaction disappeared (median log₂FC 1.08 to 0.003; fry P = 0.81; attenuation 100%).

**Consensus up arm.** Ventral-horn effect at P100 P = 1.6 × 10⁻⁴, but no interaction at any age (BH 0.98).

**White matter** (ventral-lateral vs dorsal-medial). Both territories rose from P100 onwards, with no interaction.

**Reading.** In SOD1-G93A mice the response appears at the symptomatic stage, not before onset. It first favours the ventral horn (P100) and then extends to the dorsal horn and white matter (P120). In mouse it is carried by C3 and is entirely shared with microglial content. SOD1-G93A has no TDP-43 pathology, so this is evidence of conservation, not a model of sporadic ALS.

## B. Human: seven sporadic ALS patients, no controls (descriptive)

58,087 spots with at least 500 counts, in 77 sections. The table gives median Spearman ρ over sections between the frozen ALS-up module score and lineage scores.

| Lineage | ρ | Depth-partial ρ | Sections with ρ > 0 |
|---|---|---|---|
| Microglia | 0.062 | **0.055** | 90% |
| Astrocyte | 0.045 | 0.034 | 78% |
| Mural/fibroblast | 0.034 | 0.025 | 73% |
| Endothelial | −0.001 | −0.004 | 47% |
| Oligodendrocyte | −0.035 | −0.010 | 36% |
| Neuron | −0.078 | −0.074 | 0% |

- Prespecified reading: **replicates W66** (microglia highest), although every correlation is weak at 100-µm resolution.
- Module CPM was higher in the ventral than the dorsal horn in all seven patients (log₂ ratio 0.10-0.56). Without controls, this cannot be separated from baseline anatomy.

## Files

- `w71_mouse_tests.tsv`, `w71_mouse_gene_effects.tsv`, `w71_mouse_microglia_adjusted.tsv`, `w71_mouse_reading.tsv`
- `w71_human_colocalisation_summary.tsv`, `w71_human_colocalisation_sections.tsv`, `w71_human_colocalisation_by_patient.tsv`, `w71_human_ventral_dorsal_ratio.tsv`
- `inputs/` (unit counts, spot table, section log), plus run logs
