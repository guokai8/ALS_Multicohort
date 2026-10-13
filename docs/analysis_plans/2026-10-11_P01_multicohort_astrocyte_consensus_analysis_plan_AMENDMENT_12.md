# Amendment 12: dose-response with motor-neuron loss, and cross-cohort combination of the territory tests (W69)

Date: 2026-10-13
Amends: plan v1.0 (sha256 da63b81a...bfb3bd)
Status: committed **before** W69 A was run. W69 B combines results already known (W67, W68) and is declared as a post hoc summary.

## Why

The W68 replication of the motor-territory interaction was nominal (P 0.040, BH 0.081). Two biologically motivated questions follow:
1. Does the response scale with motor-neuron loss, as a reaction to degeneration would?
2. What is the combined evidence of the two independent, prespecified territory tests?

**Adding donor-level covariates is not used here.** Only age and sex are available, they are incomplete in GSE288365, and the interaction is a within-donor contrast. The primary W67 and W68 results stay as reported.

Inspected before writing (no outcomes): the presence of the motor-neuron genes below in both platforms. No motor-neuron score has been computed.

## W69 A: dose-response with motor-neuron loss (new; fixed before running)

**Cohorts.**
- Lumbar cord: W67 inputs (GSE330224); independent donors with both territories (18 ALS, 5 control); ventral vs dorsal horn.
- Cervical cord: W68 inputs (GSE288365); primary donors (9 ALS, 7 control); anterior horn vs other grey matter.

**Normalisation.** Within each cohort, the motor- and sensory-territory units of all included donors are TMM-normalised together and converted to edgeR log2-CPM (prior count 1).

**Measures, per donor.**
- **Motor-neuron index:** mean log2-CPM in the motor territory of 10 motor-neuron genes: CHAT, SLC5A7, SLC18A3, MNX1, ISL1, PRPH, NEFH, NEFM, NEFL, STMN2. Lower values mean more motor-neuron loss.
  - Panel choice: the 17-gene panel of the source study is not listed in text, so this panel uses established cholinergic and large-motor-neuron markers. STMN2 is also a TDP-43 target, so its loss reflects both neuron loss and TDP-43 dysfunction.
- **Module response:** mean log2-CPM of the frozen ALS-up genes in the motor territory minus that in the sensory territory. This is the donor-level counterpart of the W67/W68 interaction.

**Test.** Within ALS donors only, the Spearman correlation between the motor-neuron index and the module response.
- *Prediction:* rho < 0 (more motor-neuron loss, larger motor-territory response).
- One-sided P per cohort, from the exact permutation distribution (Spearman, `cor.test`, alternative "less").
- Combined over the two cohorts by unweighted Stouffer.

**Specificity control.** Draw 1,000 random gene sets with the same number of genes as the measured frozen ALS-up set, matched on motor-territory mean log2-CPM decile within the cohort. Compute the same rho for each. Report the empirical percentile of the module rho (seed 20261013). This controls for compositional coupling: when neurons are lost, other genes gain relative share.

**Reading (frozen ALS-up module):**

| Outcome | Reading |
|---|---|
| Combined one-sided P < 0.05, both rhos < 0, and module rho below the 5th percentile of random sets in both cohorts | **Dose-response, specific** |
| Combined P < 0.05 and both rhos < 0, but not below the 5th percentile in both cohorts | Dose-response, not distinguishable from compositional coupling |
| Otherwise | Not supported |

**Secondary.**
- The same analysis for the consensus up arm.
- Descriptive: motor-neuron index in ALS vs control (Wilcoxon), to confirm that the index reflects loss.

## W69 B: cross-cohort combination of the territory interaction (post hoc summary of known results)

- **Inputs:** the primary grey-matter interaction fry tests of W67 (lumbar) and W68 (cervical) for the frozen ALS-up module (primary) and the consensus up arm (secondary).
- **Transformation:** two-sided fry P is converted to one-sided in the hypothesised (Up) direction: P/2 if the direction is Up, 1 − P/2 otherwise.
- **Combination:** unweighted Stouffer, Z = (z1 + z2)/√2.
- **Reporting:** the result is reported as a summary of two independent, prespecified tests. It is not a new test, and it does not replace either primary result.

## Outputs

- `results/W69_dose_response_and_combination/`
- Script: `scripts/run_w69_dose_response_and_combination.R`
