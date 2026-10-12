# Amendment 10: review-prompted W66 sensitivity analyses and an anatomical test of motor-system selectivity (W67)

Date: 2026-10-13
Amends: plan v1.0 (sha256 da63b81a...bfb3bd)
Status: committed before W67 B was run. Part A was observed during an internal simulated review and is declared as such.

## A. Review-prompted W66 sensitivity analyses (post hoc; already observed)

During the internal simulated review of 2026-10-13 (`09_internal_review/simulated_peer_review_spatial_2026-10-13.md`), two checks were computed on the committed W66 inputs. Their outcomes were seen before this amendment was written. This amendment fixes their reproducible implementation and reporting, and reports them as post hoc.

**A1. Empirical null for ligand-receptor co-detection.**
- For each W66 D pair, draw 200 random gene pairs from genes detected in >= 1% of eligible spots (total counts >= 500). The first gene of each random pair is matched to the ligand and the second to the receptor, each within the same overall detection-fraction decile.
- Compute the W66 D statistic for each random pair: per section, log2((both + 0.001) / (p1 × p2 + 0.001)); donor mean; median over independent donors.
- Report the empirical percentile and an empirical one-sided P: (1 + number of random pairs >= observed) / 201. Seed 20261013.
- *Wording rule:* co-detection is called "above matched random pairs" only if the empirical P < 0.05.

**A2. Depth-partial co-localisation.** Repeat W66 A2 with Spearman ranks of both variables residualised on the rank of log total counts within each section. Report the median over sections.

## B. W67: motor-system selectivity in lumbar spinal cord (new; fixed before running)

**Hypothesis.** If the bulk-confirmed programme reflects the motor-system pathology of ALS, its ALS-related change should be larger in motor than in sensory territory of the same sections:
- grey matter: ventral horn (lower motor neurons) vs dorsal horn (sensory);
- white matter: lateral funiculus, which contains the lateral corticospinal tract (`Med_Lat_White` + `Vent_Lat_White`), vs dorsal columns (`Dors_Med_White`), which are largely spared in ALS.

Spot coordinates were inspected only to confirm the dorsoventral and mediolateral orientation of the labels.

**Units.** Donor-level pseudobulk of raw counts per anatomical label group, summed over sections. Independent donors only: 18 ALS and 5 control (W66 rule).

**Model.** For each comparison (grey: ventral vs dorsal horn; white: lateral funiculus vs dorsal columns):
- Units of the two territories; edgeR filterByExpr and TMM normalisation.
- limma-voom with `~ territory * disease + sex + age`, blocked on donor (duplicateCorrelation, two passes).
- Interaction contrast `territory(motor):diseaseALS`, tested with fry (primary) and camera.
- Gene sets: frozen ALS-up and ALS-down modules; consensus up and down arms.
- Also reported: the ALS vs control contrast within each territory (same fit).

**Primary family.** ALS-up module and consensus up arm, both comparisons: 4 interaction tests, BH within the family.

**Reading (per gene set):**

| Outcome | Reading |
|---|---|
| Interaction fry BH < 0.05, direction Up (Down for down sets) | ALS change larger in motor territory: motor-system selective |
| Interaction not significant, same-direction ALS change in both territories | Not motor-selective (diffuse tissue response) |
| Otherwise | Not resolved |

There are only 5 controls, so a non-significant interaction is weak evidence against selectivity.

**Secondary (motor cortex).** The same design in MTC, comparing Layer_5 (Betz cells, upper motor neurons) with Layer_2 + Layer_3. Independent donors with age. Nominal P values.

## Outputs

- `results/W66_spatial_context/post_hoc/`: A1, A2.
- `results/W67_motor_selectivity/`: B.

Scripts:
- `scripts/run_w66_post_hoc_sensitivity.py`
- `scripts/build_w67_annotation_pseudobulk.py`
- `scripts/run_w67_motor_selectivity.R`
