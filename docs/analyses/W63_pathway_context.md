# W63: pathway context of the reproducible programme

This is an **exploratory** analysis, prespecified in Amendment 8 (committed before it was run).

- Script: `scripts/run_w63_pathway_context.R`
- Gene sets: msigdbr 26.1.0.
- Ranking statistic: the W58 run (ii) meta statistic (meta_logFC / meta_se). Run (ii) is the discovery-agnostic meta-analysis over five independent single-nucleus contrasts; the universe is its 12,784 tested genes.

## Themed preranked GSEA (`w63_themed_gsea.tsv`; primary family, 14 sets, BH across them)

| Theme | Gene set | NES | FDR | Leading-edge genes (first) |
|---|---|---|---|---|
| Extracellular matrix | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 1.88 | 0.0002 | SPP1, MYLK, SPOCK1, CALD1, COL4A2, TPM4, ECM2, CD44, EMP3 |
| Inflammation | HALLMARK_INTERFERON_GAMMA_RESPONSE | 1.68 | 0.002 | PLSCR1, IFNAR2, NOD1, CIITA, CD40, C1R, STAT3, B2M, CCL2 |
| Extracellular matrix | GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION | 1.55 | 0.004 | PXDN, COL4A2, THSD4, ECM2, COLEC12, COL4A1, TGFB2 |
| Inflammation | HALLMARK_TNFA_SIGNALING_VIA_NFKB | 1.63 | 0.004 | BCL6, SPSB1, SLC2A3, CD44, REL, SMAD3, PLAUR, CCL2, TLR2 |
| Inflammation | HALLMARK_INFLAMMATORY_RESPONSE | 1.63 | 0.008 | IFNGR2, AHR, EMP3, CD40, OSMR, P2RX7, PLAUR, CCL2, TLR2 |
| Complement | GOBP_COMPLEMENT_ACTIVATION | 1.73 | 0.022 | MASP1, C1RL, C1R, A2M, CD55, C1S, C3, C1QB, C2, CFI |
| Vascular | HALLMARK_ANGIOGENESIS | 1.59 | 0.048 | SPP1, PDGFA, VCAN, VAV2 |
| Vascular | GOBP_BLOOD_VESSEL_MORPHOGENESIS | 1.21 | 0.074 | |
| Neuronal support | GOBP_SYNAPSE_ORGANIZATION | 1.12 | 0.17 | |
| Complement | HALLMARK_COMPLEMENT | 1.09 | 0.38 | |
| Blood-brain barrier | GOBP_MAINTENANCE_OF_BLOOD_BRAIN_BARRIER | -1.09 | 0.42 | |
| Neuronal support | GOBP_AMINO_ACID_TRANSPORT | 1.02 | 0.47 | |
| Glial activation | GOBP_ASTROCYTE_ACTIVATION | -0.96 | 0.56 | |
| Glial activation | GOBP_GLIAL_CELL_ACTIVATION | 0.82 | 0.74 | |

## Hallmark GSEA (`w63_hallmark_gsea.tsv`; 50 sets)

- **Up in ALS:** epithelial-mesenchymal transition (ECM), allograft rejection, IL6-JAK-STAT3, interferon-gamma response, hypoxia and TNFA-NF-kB (FDR <= 0.004).
- **Down in ALS:** oxidative phosphorylation, DNA repair, MYC targets V1 and fatty-acid metabolism (FDR <= 0.004).

## Consensus-arm over-representation (`w63_consensus_ora.tsv`)

- Consensus-up (86 genes): only epithelial-mesenchymal transition is enriched (5 genes: CALD1, COL4A2, MYLK, SPOCK1, SPP1; 5.5-fold; FDR 0.031).
- Consensus-down (73 genes): no themed set reached FDR < 0.05. Amino-acid transport was nominal (GABBR1, SLC1A4, SLC25A22; P 0.023, FDR 0.32).

## Frozen module genes (`w63_module_gene_themes.tsv`)

- 8 of the 15 ALS-up genes fall in at least one themed set:
  - vascular: EGF, GPNMB, ANGPT1, ANG, ANXA1, C3;
  - complement: LTF, ANG, C3;
  - ECM and inflammation: EMP3;
  - neuronal support: C3, as a complement-mediated synapse gene.
- 6 of the 14 ALS-down genes fall in a themed set: ECM (ADAMTS17, CDH2), vascular (SPRED1, CDH2), neuronal support (EFNA5, CDH2, LRRC8C) and inflammation (FOSB).

## NYGC Flagship ST5 context (`w63_nygc_st5_themed_go_context.tsv`)

In the consortium's own bulk GSEA, complement activation is up in ALS motor cortex (NES 1.88) and in cervical and lumbar spinal cord. Glial cell activation and blood vessel morphogenesis are also up in spinal cord. These are published results, shown only as context.

## Reading

- The discovery-agnostic single-nucleus statistic shows coordinated **ECM remodelling, interferon and NF-kB inflammatory signalling, complement activation and angiogenesis** up in ALS, and **oxidative phosphorylation and DNA repair** down.
- Astrocyte- and glial-activation GO terms are **not** enriched, which is consistent with a multicellular glial-vascular programme rather than an astrocyte-specific reactive state.
- All of this is exploratory; it describes biological context and does not establish mechanism.
