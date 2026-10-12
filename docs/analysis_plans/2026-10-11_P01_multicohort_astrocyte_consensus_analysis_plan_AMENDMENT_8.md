# Amendment 8: exploratory biological context of the reproducible programme (pathways, bulk specificity, ligand-receptor)

Date: 2026-10-12
Amends: plan v1.0 (sha256 da63b81a...bfb3bd); adds exploratory, review-prompted analyses
Status: committed **before** any analysis below was run.

For this amendment, only the following were inspected, and none of them involves an outcome:
- gene-set names in msigdbr 26.1.0;
- the column layout and cell-type labels of the NYGC Flagship ST7 deconvolution table (763 of the 820 W55 hold-out samples have estimates);
- the column layout of the liana Consensus ligand-receptor resource (4,701 pairs);
- the W56 lineage pseudobulk files.

No enrichment, fold change, set test or interaction score has been computed.

All three analyses are **exploratory**. They do not change W55 or W59, and they are reported whatever the outcome.

## Why

External review asked two questions. Why does the reproducible signal matter biologically? And does the bulk hold-out confirm something beyond the glial activation already reported in NYGC data (Humphrey et al. 2023; Tam et al. 2019)?

## W63: pathway context

**Prespecified themed gene sets** (msigdbr 26.1.0, *Homo sapiens*), forming one family:

| Theme | Gene sets |
|---|---|
| Complement | HALLMARK_COMPLEMENT; GOBP_COMPLEMENT_ACTIVATION |
| Extracellular matrix | GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION; HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION |
| Vascular and blood-brain barrier | HALLMARK_ANGIOGENESIS; GOBP_BLOOD_VESSEL_MORPHOGENESIS; GOBP_MAINTENANCE_OF_BLOOD_BRAIN_BARRIER |
| Inflammation and glial activation | HALLMARK_INFLAMMATORY_RESPONSE; HALLMARK_TNFA_SIGNALING_VIA_NFKB; HALLMARK_INTERFERON_GAMMA_RESPONSE; GOBP_GLIAL_CELL_ACTIVATION; GOBP_ASTROCYTE_ACTIVATION |
| Homeostatic and neuronal support | GOBP_SYNAPSE_ORGANIZATION; GOBP_AMINO_ACID_TRANSPORT |

**Analyses:**

1. **Preranked GSEA** (fgsea, 10,000 permutations, seed 20261012).
   - Ranking: the W58 run (ii) genome-wide meta-analysis statistic (meta_logFC / meta_se), over its tested universe, using gene symbols.
   - Primary family: the 14 themed sets, with BH across them.
   - Also reported: all 50 Hallmark sets (BH within Hallmark).
2. **Over-representation** (one-sided Fisher test) of the consensus-up (86) and consensus-down (73) genes.
   - Background: the run (ii) tested universe.
   - Gene sets: the 14 themed sets, with BH within each arm.
3. **Descriptive**: theme membership of each frozen module gene.
4. **Descriptive context**: the NYGC Flagship ST5 GSEA results for the themed GO terms in MCX and CSC.

## W64: bulk specificity relative to published composition (R-3)

**Samples.** W55 hold-out samples with an ST7 deconvolution estimate (complete case).

**Designs.**
- **(a) Primary for W64.** The W55 primary design, with the four marker-based composition scores replaced by the ST7 cell-type proportions for that tissue. One proportion per tissue is dropped as the reference to avoid collinearity: Excitatory in cortex, Neurons in spinal cord.
- **(b)** The W55 primary design plus the ST7 proportions, keeping the marker scores.

**Model and gene sets.** Donor-blocked limma-voom; blocked fry primary, camera reported. Gene sets: the frozen ALS-up and ALS-down modules, and the W59 consensus up and down arms.

**Reading for the ALS-up module, design (a):**

| Outcome | Reading |
|---|---|
| fry P < 0.05 with the hypothesised direction in both regions | Signal not explained by the composition captured in the consortium's own deconvolution |
| fry P < 0.05 in one region only | Partly explained |
| Neither region | Not separable from published composition |

**Descriptive overlap.** Frozen module genes with padj < 0.05 and the hypothesised direction in NYGC Flagship ST4 for MCX and for CSC. ST4 includes all NYGC donors, including discovery donors, so this is context only.

## W65: ligand-receptor context

**Resource.** liana Consensus. Keep pairs in which a frozen module gene is the ligand or the receptor, or a subunit of either; complexes are split on "_".

**Cohorts.** Pineda BA9 batch 210526 and Ulm/DZNE, both from the W56 lineage pseudobulk. Pineda BA4 is excluded because of its processing confound.

**Expression.**
- Per-lineage mean CPM over all donors.
- A pair is expressed for a sender-receiver combination if the ligand has mean CPM >= 1 in the sender and the receptor has mean CPM >= 1 in the receiver. For complexes, the minimum over subunits is used.

**Disease change.**
- Per-lineage limma-voom `~ disease + sex` (Ulm: donor-level sum over chips) gives log2 fold changes.
- The pair score change is the sum of the ligand logFC in the sender and the receptor logFC in the receiver.

**Reported.**
- All expressed pairs that involve module genes, with their sender and receiver lineages.
- The pairs whose score change has the same sign in both cohorts.

**Reading.** Hypothesis-generating only. No inference on mechanism, and no P values are used for selection.

## Outputs

- `results/W63_pathway_context/`
- `results/W64_bulk_published_composition/`
- `results/W65_ligand_receptor_context/`
