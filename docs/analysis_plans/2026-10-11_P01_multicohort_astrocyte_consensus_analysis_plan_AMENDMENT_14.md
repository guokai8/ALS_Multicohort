# Amendment 14: spatial time course in SOD1-G93A mice and a third human spatial cohort (W71), and disease-state cell source in human spinal cord (W72)

Date: 2026-10-14
Amends: plan v1.0 (sha256 da63b81a...bfb3bd)
Status: committed **before** W71 and W72 were run. All analyses are exploratory.

## Why

Two questions remain open after W67-W70:
- When does the ventral-horn response appear relative to motor neuron loss? W69 found no dose-response within ALS donors.
- Which cells express the module genes in diseased spinal cord? Tanaka et al. (Sci Rep 2012;2:573) reported GPNMB immunoreactivity in motor neurons and GFAP-positive astrocytes, but not in IBA1-positive microglia, in ALS patients and SOD1-G93A mice. This conflicts with a myeloid reading of GPNMB. Our own baseline lineage shares (W56) come from unaffected cortex.

## Inspected before writing (no outcome)

- **Maniatis et al. 2019 (Science 364:89; ALS-ST portal).** The md5 checksums of the count-matrix and metadata archives matched. The following were read:
  - the metadata tables (design counts below);
  - the annotation column names;
  - the count-file structure: genes × spots; spot identifiers are `x_y` array coordinates matching the annotation `xPos`/`yPos`; arrays carrying two sections are split into `_1`/`_2` files, one per animal; a gene absent from a file has zero counts.
- **Number of mouse sections in which each module ortholog is detected**, with all genotypes and ages pooled. Of 327 mouse count files: Kitl 327, Angpt1 325, C3 325, Emp3 314, Anxa1 301, Ang 280, Ramp3 276, Gpnmb 259, Egf 212, Stpg1 171, Ltf 112.
- **Donor overlap.** Of the seven human patients, one (CGND-HRA-00041) is in the W55 bulk hold-out. None is in GSE330130, GSE330224 or GSE288365 under the identifiers available.
- **GSE287257 (Zelic).** The obs columns of the processed object (donor, disease) were read. The W58 builder already applies the W54 marker rule to this object.

No expression summary stratified by genotype, age, territory or diagnosis has been computed.

## W71 A: SOD1-G93A spatial time course (primary)

**Data.**
- Lumbar spinal cord (L3-L5), B6SJL SOD1-G93A vs SOD1-WT, at P30, P70, P100 and P120.
- Animals (G93A / WT): P30 10 / 6; P70 8 / 4; P100 5 / 8; P120 9 / 6.
- The ATG7(fl/fl);ChAT-Cre;SOD line is excluded, because it has no matched control.

**Territories.**
- Spots with exactly one annotation region are used.
- Motor territory: `Vent_Horn`. Sensory territory: `Dors_Horn`.
- Secondary white-matter comparison: `Vent_Lat_White` vs `Dors_Med_White`.

**Units.**
- Raw counts are summed per animal and territory over all of that animal's sections.
- A unit needs at least 5 spots and a library size above 0. Animals need both territories of a comparison.

**Gene sets.** Mouse orthologs from babelgene 22.9, with support from at least 3 sources:
- Frozen ALS-up module (11 genes): Ang, Angpt1, Anxa1, C3, Egf, Emp3, Gpnmb, Kitl, Ltf, Ramp3, Stpg1. The four human lncRNAs have no ortholog.
- Consensus up arm: the 74 of 85 genes with an ortholog.

**Model.** Fitted separately at each age:
- edgeR filterByExpr and TMM normalisation;
- limma-voom `~ territory * genotype + sex`, with the sex term dropped if an age has one sex only;
- animal blocking with duplicateCorrelation (two passes).

Contrasts:
- the territory × genotype interaction;
- the genotype effect within the ventral horn;
- the genotype effect within the dorsal horn.

Tests: fry (primary) and camera.

**Primary family.** The frozen ALS-up interaction at the four ages (BH over four tests).

**Readings.**

| Question | Rule | Reading |
|---|---|---|
| Cross-species replication | Interaction Up with BH < 0.05 at P100 or P120 | Replicated in SOD1-G93A |
| | Interaction Up with fry P < 0.05 at P100 or P120, BH ≥ 0.05 | Nominal |
| | Otherwise | Not replicated |
| Timing (descriptive) | Earliest age with ventral-horn genotype effect Up at fry P < 0.05 | P30: presymptomatic; P70: at onset; P100/P120: symptomatic or late; none: no ventral-horn response |

**Secondary analyses.**
- **Consensus up arm.** The same tests, with BH over its four interaction tests.
- **Microglia adjustment.** At ages where the unadjusted frozen-module interaction has fry P < 0.05, the W70 adjustment is applied with the mouse orthologs of the all-microglia score (P2ry12, Cx3cr1, Csf1r, Ptprc, Cd68, Trem2, Tyrobp, Cd163, Itgax). Attenuation is defined as in W70.
- **Gene level.** The ventral-horn genotype log₂ fold change of Gpnmb, Anxa1, Ltf, Emp3 and C3 at each age. The prediction, carried over from the human cohorts, is that each is higher in G93A at P100 and P120.
- **White matter.** The interaction in the white-matter comparison, reported without a reading.

**Caveat.** SOD1-G93A is an overexpression model of SOD1 familial ALS, without TDP-43 pathology. Agreement supports conservation of the response. Disagreement does not contradict the human data.

## W71 B: human ALS spinal cord, descriptive (secondary)

**Data.** Seven sporadic ALS patients (lumbar or bulbar onset), with 80 cervical and lumbar sections. There are no control donors, so **no disease test is made**.

**Spot scores.** The W66 transforms and marker sets are used: lineage scores and module score as log1p(10⁴ × summed counts / spot total). Spots need at least 500 counts, and sections need at least 30 qualifying spots.

**Co-localisation.** Spearman correlation between the frozen ALS-up module score and each lineage score, computed within each section, raw and depth-partial (W66 method). Medians are reported over sections and per patient.

**Reading.** The result "replicates W66" if the depth-partial median correlation with microglia exceeds those with astrocyte, endothelial and mural/fibroblast scores.

**Descriptive.** The ventral-horn to dorsal-horn ratio of module CPM per patient, without a test.

## W72: disease-state cell source in human spinal cord (secondary)

**Data.** GSE287257 (Zelic et al.): cervical spinal cord snRNA-seq, raw counts, 8 ALS and 4 control donors.

**Lineage calls.** W54 marker sets and margin rule, applied to every lineage:
- A nucleus is assigned to the lineage with the highest score if that score exceeds the second highest by at least 0.5 and the nucleus has at least 500 UMIs. Other nuclei are unassigned.
- The lineages are astrocyte, neuron (including motor neurons), oligodendrocyte, OPC, microglia and vascular.
- No module gene is a marker.

**Disease-state lineage share.**
- Counts are summed per donor and lineage, and CPM is computed within each unit.
- For each frozen ALS-up gene, the share is the mean CPM over ALS donors in a lineage divided by the sum over lineages. The same is computed in controls.

**Reading (descriptive, no P value).**

| GPNMB share in ALS donors | Reading |
|---|---|
| Microglia ≥ 0.5 | Myeloid in ALS spinal cord |
| Astrocyte or neuron ≥ 0.5 | Astrocytic or neuronal, consistent with Tanaka et al. |
| Otherwise | Mixed |

The same rule is reported for C3, ANXA1, LTF and EMP3.

**Secondary.**
- The share shift, ALS minus control, for each gene and lineage.
- Within-lineage ALS vs control log₂ fold change: limma-voom `~ disease`, for lineages with at least 3 donors per group having 20 or more nuclei. With 4 controls this is descriptive.

**Caveat.** Nuclear RNA may under-represent transcripts of secreted proteins and their cell of origin. Protein localisation needs immunostaining.

## Outputs

- `results/W71_maniatis_spatial_time_course/`
- `results/W72_zelic_disease_cell_source/`
- Scripts:
  - `scripts/build_w71_maniatis_pseudobulk.py`
  - `scripts/run_w71_mouse_time_course.R`
  - `scripts/run_w71_human_colocalisation.py`
  - `scripts/run_w72_zelic_disease_cell_source.py`

Every outcome is reported.
