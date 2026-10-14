# Amendment 16: is the ALS-up signal astrocyte-intrinsic or shared with the microglial compartment? (W74)

Date: 2026-10-14
Amends: plan v1.0 (sha256 da63b81a...bfb3bd)
Status: committed **before** W74 was run. Exploratory.

## Why

The tissue-level analyses (W55-W73) show which cell types are associated with the module in tissue. They do not show where the astrocyte-compartment signal itself comes from. Two readings are possible:
- **(i) Astrocyte-intrinsic.** ALS changes these genes inside astrocytes, as W72 suggests for GPNMB and C3 in spinal cord.
- **(ii) Shared with myeloid cells.** The same ALS change is present, or larger, in the microglial compartment of the same donors.

This amendment tests the two directly, within donors. It was prompted by the pre-submission story review.

## Already observed (declared)

- **W56 E2.** Module-level ALS effects within single lineages of Ulm and Pineda BA9 (batch 210526). No lineage was significant:
  - Ulm astrocytes fry P = 0.62, microglia 0.37;
  - Pineda BA9 astrocytes P = 0.48, microglia 0.36.
  - Only 2-3 module genes passed filtering in microglia.
- **W72.** Gene-level within-lineage ALS effects in GSE287257 (Zelic). Among them, GPNMB and C3 were higher in ALS astrocytes. No module-level within-lineage test was run.
- **W01 (August 2026, before the plan).** Per-lineage differential expression in GSE330130, including the motor-cortex microglia compartment, used pre-reconciliation diagnosis labels. Module-gene values in the W01 microglia table were not used in any governed analysis. They may have been seen when W01 was run, so the discovery microglia result below is declared as "possibly seen in an earlier, differently labelled run".
- **Before writing,** only the sample tables of the GSE330130 astrocyte and microglia pseudobulk were read: 22 libraries each, the same donors and libraries. No expression value was read.

## Design

**A. Discovery microglial compartment (primary).**
- *Data:* GSE330130 motor-cortex microglia pseudobulk, from the same libraries as the discovery astrocyte compartment, with W33 synapse-resolved diagnoses.
- *Model:* edgeR filterByExpr / TMM / limma-voom, `~ disease + age + sex + PMI + tissue_subregion`, with the W33 covariate-drop rules.
- *Test:* the frozen ALS-up module (and, secondarily, the ALS-down module) by fry (primary) and camera.
- *Reading:* "**ALS-up signal shared with the microglial compartment**" if fry P < 0.05, direction Up and median log2FC > 0. Otherwise "not detected in the microglial compartment".

**B. Discovery astrocyte vs microglia, same libraries (descriptive).**
- *Model:* `~ lineage * disease + age + sex + PMI + tissue_subregion`, with library blocking (duplicateCorrelation, two passes).
- *Report:*
  - the lineage × disease interaction for the module (fry);
  - the ratio of microglial to astrocytic median log2FC.
- *Caveat:* the modules were selected in the astrocyte compartment, so the astrocyte effect is inflated by selection. A microglial effect close to or above the astrocytic one therefore argues against reading (i). A smaller one does not establish it.

**C. Independent cohorts, astrocyte vs microglia in the same donors (secondary).**
- *Cohorts:*
  - Ulm: W56 lineage pseudobulk, donor × chip, donor-blocked, `~ lineage * disease + sex + chip`;
  - Pineda BA9 batch 210526: W56 lineage pseudobulk, `~ lineage * disease + age + sex`;
  - Zelic: W72 lineage pseudobulk, `~ lineage * disease + age + sex`.
  - Li et al. is not included, because no lineage pseudobulk exists for it.
- *Reported per cohort:*
  - module ALS effect within astrocytes and within microglia (fry; the Ulm and Pineda values are already known from W56 E2);
  - the lineage × disease interaction.
- *Reading per cohort:*
  - "larger in astrocytes" if the astrocyte-minus-microglia interaction is Up with fry P < 0.05;
  - "larger in microglia" if it is Down with P < 0.05;
  - otherwise "not distinguishable".
- *Combination:* each cohort's interaction fry P is converted to a signed z (positive when the ALS effect is larger in astrocytes than in microglia). The three z values are combined by unweighted Stouffer and reported with a two-sided P.

**D. Discovery, all lineages (descriptive).** The ALS-up module effect in each GSE330130 motor-cortex lineage (astrocyte, microglia, oligodendrocyte, OPC, excitatory and inhibitory neurons). Same model as A.

## Interpretation fixed in advance

| Outcome | Interpretation |
|---|---|
| A positive, and B or C not "larger in astrocytes" | The discovery signal is shared with the microglial compartment. This supports a multicellular, myeloid-associated programme that astrocyte-compartment profiling also captures |
| A not detected, and C "larger in astrocytes" or "not distinguishable" | The signal is mainly astrocyte-compartment-associated in single-nucleus data. The tissue-level myeloid association then reflects co-regulation or composition in tissue, not a microglial origin of the snRNA signal |
| Mixed | Reported as such |

Every outcome is reported, and the main text is revised to match whichever holds.

## Limits stated in advance

- **Low power.** Microglia are sparse in some cohorts (Pineda median about 166 nuclei per sample), and few module genes may pass filtering.
- **Selection bias.** Selection in the discovery astrocyte compartment biases B toward astrocytes.
- **Ambient RNA.** Ambient RNA can move transcripts between compartments in either direction. W57 addressed this for the discovery astrocytes only.

## Outputs

- `results/W74_within_lineage_source/`
- `scripts/run_w74_within_lineage_source.R`
