# W68: replication of motor-territory selectivity in an independent spatial cohort (GSE288365)

This analysis is exploratory. Amendment 11 (`c5af36b`) fixed the territories, model, reading rule, sensitivity analyses and gene-level prediction before any expression summary was computed.

**Data.** GSE288365 (Zhang et al., *Nat Neurosci* 2026): Visium CytAssist FFPE cervical spinal cord, with a 35-antibody protein panel; 18 donors. The authors' spot annotation (`manual`) defines the territories. The sample-to-donor map comes from their Supplementary Table 1 and is hard-coded in `build_w68_gse288365_pseudobulk.R`.

**Donors.** The primary analysis excludes NEUVM032BF2 (a discovery donor) and NEUWD570BTK (also in W67), leaving 9 ALS (4 C9orf72, 5 sporadic) and 7 controls.

**Scripts.**
- `extract_w68_gse288365_metadata.R`
- `build_w68_gse288365_pseudobulk.R`
- `run_w68_replication.R`
- `run_w68_protein_colocalisation.R`

## A. Motor vs sensory territory (`w68_replication_tests.tsv`; fry P, median log2FC)

| Comparison | Gene set | Interaction (motor − sensory) | ALS in motor | ALS in sensory |
|---|---|---|---|---|
| Grey: anterior horn vs other grey | **Frozen ALS-up** (8 genes tested) | **P 0.040 (BH 0.081), Up, 0.41** | **P 1e-4, 1.20** | P 0.047, 0.39 |
| | Consensus up | **P 0.0093 (BH 0.037), Up, 0.19** | P 1e-4, 0.21 | P 0.059, 0.04 |
| | Consensus down | P 4e-4, Down | P < 1e-4, -0.18 | P 0.27, -0.09 |
| White: lateral/ventral vs dorsal columns | Frozen ALS-up | P 0.14 (BH 0.19), Up, 0.24 | P 5e-4, 0.87 | P 0.072, 0.40 |
| | Consensus up | P 0.20 (BH 0.20) | P 2e-4, 0.28 | P 0.025, 0.14 |

**Reading (prespecified, frozen ALS-up module, grey-matter interaction): replicated at nominal level** (P < 0.05, Up; BH 0.081).

**Sensitivity analyses** (grey interaction, frozen ALS-up):

| Analysis | Interaction P | ALS in anterior horn |
|---|---|---|
| Excluding W55 hold-out donors (8/5) | 0.059 | P 0.0046, median 1.01 |
| All 18 donors | 0.034 | P 1e-4, 1.15 |
| Sporadic ALS only (5/7) | 0.081 | P 0.0025, 0.95 |
| + sex and age (6/6) | 0.11 | P 0.0012, 0.95 |

All analyses agree in direction. The anterior-horn increase is significant in every one; the interaction is significant in two.

**Consensus signature.** Its grey-matter interaction was significant here (BH 0.037), unlike in lumbar cord (W67, BH 0.31). In both cohorts the anterior-horn increase is similar (median 0.21). The cohorts differ in the sensory-territory change: 0.14 in lumbar cord, 0.04 in cervical cord.

## Gene level (`w68_module_gene_results.tsv`; anterior horn vs other grey, ALS vs control log2FC)

| Gene | Anterior horn | Other grey | Interaction P |
|---|---|---|---|
| GPNMB | **2.83** (P 5e-4) | 1.52 | 0.073 |
| LTF | 2.24 (P 0.057) | 1.81 | 0.69 |
| C3 | **1.97** (P < 1e-4) | 0.89 | 0.0027 |
| ANXA1 | **1.44** (P 1e-4) | 0.10 | 0.0002 |
| EMP3 | **0.95** (P 6e-4) | 0.57 | 0.12 |
| ANGPT1 | 0.08 | 0.21 | 0.60 |
| RAMP3 | -0.33 | -0.19 | 0.52 |

**Prediction from W67:** LTF, ANXA1, GPNMB and EMP3 are higher in the ALS anterior horn. All four are in the predicted direction, and three are nominally significant. In contrast to lumbar cord, C3 also rises in the cervical anterior horn.

## B. Motor-neuron proximity (secondary)

Within the anterior horn, the ALS increase of the frozen module is similar in spots next to motor neurons and in other spots (median 1.69 and 1.59; interaction P 0.92). The response spans the anterior horn and is not confined to the immediate neighbourhood of surviving motor neurons.

## C. Protein co-localisation (`w68_protein_colocalisation_summary.tsv`; median Spearman over 18 sections)

| Protein | All annotated spots: raw / depth-partial (fraction of sections > 0) | Anterior horn: raw (fraction > 0) |
|---|---|---|
| VIM | 0.21 / 0.14 (1.00) | 0.06 (0.78) |
| CD68 | 0.12 / 0.11 (1.00) | 0.07 (0.89) |
| HLA-DRA | 0.07 / 0.08 (0.78) | 0.07 (0.89) |
| PECAM1 | -0.06 / -0.03 (0.28) | 0.00 |
| ACTA2 | 0.01 / 0.03 (0.50) | -0.01 |
| ITGAM | -0.18 / -0.12 (0.06) | -0.04 |

**Reading.** At the protein level, the module tracks the phagocytic myeloid marker CD68, MHC class II (HLA-DR) and vimentin. It does not track endothelial (PECAM1) or mural (ACTA2) proteins. The negative ITGAM correlation is unexplained; antibody performance in FFPE tissue is one possibility.

## Overall

- In an independent cohort, the frozen module's ALS increase is concentrated in the anterior horn (nominal replication), carried by GPNMB, C3, ANXA1, EMP3 and LTF.
- At the protein level, it co-localises with CD68 and HLA-DR, consistent with a phagocytic myeloid response in motor-neuron territory.
- The consensus signature also rises most in the anterior horn here, so the focal-versus-diffuse contrast seen in lumbar cord is weaker in cervical cord.
