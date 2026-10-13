# Amendment 11: replication of motor-territory selectivity in an independent spatial cohort (GSE288365; W68)

Date: 2026-10-13
Amends: plan v1.0 (sha256 da63b81a...bfb3bd)
Status: committed **before** any expression summary of GSE288365 was computed.

## What was inspected before writing

None of the following involves an outcome:
- the GEO file list;
- the Supplementary Table 1 of Zhang et al., *Nat Neurosci* 2026 (sample-to-donor map, diagnosis, age, sex);
- the Seurat metadata of all 18 samples: author-made spot annotations (`manual`), condition, motor-neuron distance group and IBA1 staining flag;
- the gene and protein feature lists;
- donor identifiers matched against our earlier cohorts.

No expression value has been compared between groups.

**Cohort.** GSE288365 is Visium CytAssist (FFPE) cervical spinal cord with a 35-antibody protein panel: 18 donors, 9 non-neurological controls, 4 C9orf72 ALS and 5 sporadic ALS. It comes from a different laboratory and different brain banks than GSE330224, but its donors carry NYGC/Target ALS identifiers. Donor overlap:

| Overlap | Donors |
|---|---|
| Discovery donor (GSE330130) | NEUVM032BF2 (control) |
| Also in W67 (GSE330224 lumbar cord) | NEUWD570BTK (control) |
| Also in the W55 bulk hold-out | NEUHD481VCL, NEUEP010ZFE (controls); NEUBB078VEG (C9-ALS) |

**Gene coverage.** The probe panel contains 11 of 15 frozen ALS-up genes (the four lncRNAs are absent), 14 of 14 ALS-down genes, and 71/85 and 69/73 consensus up and down genes.

## W68 A: primary replication of W67 (motor vs sensory territory)

**Donors.** All except NEUVM032BF2 and NEUWD570BTK: 9 ALS, 7 control.

**Territories.** From the authors' `manual` annotation:
- Grey matter: `Anterior_horns` (motor) vs `Grey_matter` (the remaining grey matter: dorsal horn and intermediate zone).
- White matter: `White_matter` (lateral and ventral funiculi, including the corticospinal tracts) vs `Posterior_white_matter` (dorsal columns).
- Meninges, Border and unannotated spots are excluded.

**Units.** Raw Spatial-assay counts summed per donor and territory. The two capture areas of a section are summed.

**Model.** Same as W67, but without covariates in the primary analysis:
- `~ territory * disease`, limma-voom, donor-blocked (duplicateCorrelation, two passes).
- The interaction contrast is tested with fry (primary) and camera.
- *Why no covariates:* age or sex is missing for 4 of 18 donors, and the interaction is a within-donor contrast, so donor-level covariates largely cancel.

**Gene sets.** Frozen ALS-up (primary) and ALS-down modules; consensus up and down arms.

**Primary family.** Frozen ALS-up and consensus up interactions, grey and white comparisons: 4 tests, BH.

**Reading for the frozen ALS-up module, grey-matter interaction:**

| Outcome | Reading |
|---|---|
| fry BH < 0.05, Up | **Replicated** |
| fry P < 0.05, Up (BH ≥ 0.05) | Replicated at nominal level |
| Up, P ≥ 0.05 | Directionally consistent, not replicated |
| Down | Not replicated |

**Consensus up arm.** W67 found a diffuse shift. A non-significant interaction with ALS increases in both territories would replicate that pattern.

**Sensitivity analyses** (grey-matter comparison, all four gene sets):
1. Also exclude the W55 hold-out donors (8 ALS, 5 control).
2. All 18 donors.
3. Sporadic ALS only vs control (5 vs 7).
4. Add sex and age, complete cases.

**Gene level.** Per-gene log2FC in each territory and the interaction, for the frozen ALS-up genes. W67 predicts that LTF, ANXA1, GPNMB and EMP3 are higher in ALS in the anterior horn.

## W68 B: motor-neuron proximity (secondary, exploratory)

Within `Anterior_horns`, compare spots flagged `MN_adjacent` with the remaining anterior-horn spots, using the same interaction model (`~ proximity * disease`, donor-blocked) on the frozen ALS-up module. Nominal P values.

## W68 C: protein-level co-localisation (secondary, descriptive)

Within each section, compute the Spearman correlation across spots between the frozen ALS-up module score (log1p(10^4 × module counts / spot Spatial counts)) and each protein. Protein values are raw counts with a centred log-ratio transform within each spot over the 31 non-isotype antibodies.

| Protein | Marks |
|---|---|
| CD68, CD163, ITGAM, HLA-DRA, PTPRC | Myeloid and immune cells |
| PECAM1 | Endothelium |
| ACTA2 | Mural cells |
| VIM | Reactive glia and mesenchyme |

Correlations are computed over all annotated spots and over anterior-horn spots, and also repeated with ranks residualised on log Spatial counts. The median over sections is reported. W66 predicts that the myeloid proteins correlate positively with the module.

## Outputs

`results/W68_GSE288365_replication/`, from:
- `scripts/extract_w68_gse288365_metadata.R`
- `scripts/build_w68_gse288365_pseudobulk.R`
- `scripts/run_w68_replication.R`
- `scripts/run_w68_protein_colocalisation.R`
