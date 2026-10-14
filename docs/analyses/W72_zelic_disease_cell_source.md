# W72: disease-state cell source of the frozen ALS-up genes in human spinal cord (Amendment 14)

Exploratory analysis. Amendment 14 was committed before this analysis (dc0e51d).

Data: GSE287257 (Zelic et al. 2025), cervical spinal cord snRNA-seq, raw counts, 8 ALS and 4 control donors, 114,161 nuclei.

Scripts:
- `scripts/run_w72_zelic_disease_cell_source.py`: lineage calls and shares;
- `scripts/run_w72_lineage_de.R`: within-lineage ALS vs control model, split out because it uses limma.

**Lineage calls.** The W54 marker-dominance rule, applied to all six lineages:

| Lineage | Nuclei |
|---|---|
| Oligodendrocyte | 59,550 |
| Astrocyte | 14,858 |
| Microglia | 8,396 |
| OPC | 6,437 |
| Neuron | 1,673 |
| Vascular | 1,437 |
| Unassigned | 21,810 |

## Lineage share in ALS donors (mean CPM per lineage / sum over lineages)

| Gene | Top lineage (share) | Prespecified reading |
|---|---|---|
| GPNMB | Microglia (0.74) | **Myeloid in ALS spinal cord** |
| C3 | Microglia (0.81) | **Myeloid** |
| LTF | Astrocyte (0.74) | Astrocytic |
| ANXA1 | Vascular (0.54); astrocyte 0.27 | Mixed |
| EMP3 | Vascular (0.55); microglia 0.18 | Mixed |
| ANG | Microglia (0.42) | Mixed |
| EGF, ANGPT1, AL133482.1, LINC02668, LINC01094 | Astrocyte (0.71-0.95) | Astrocytic |
| RAMP3, KITLG | Vascular | Mixed |

## ALS induces the genes in non-myeloid lineages as well

Within-lineage ALS vs control log₂FC (limma-voom; 8 vs 4; descriptive). P values are in brackets.

| Gene | Astrocyte | Oligodendrocyte | OPC | Microglia | Vascular |
|---|---|---|---|---|---|
| GPNMB | 2.62 (2.8 × 10⁻⁴) | 2.66 (6.6 × 10⁻⁶) | 1.71 (0.0028) | 1.17 (0.023) | 1.68 (0.24) |
| C3 | 2.60 (1.2 × 10⁻⁴) | 1.16 (0.013) | 1.28 (0.0034) | 0.08 (0.71) | – |
| ANXA1 | 0.99 (0.05) | 1.52 (2.7 × 10⁻⁴) | 1.14 (0.014) | 0.32 (0.33) | 0.27 (0.52) |
| EMP3 | 0.64 (0.043) | 0.82 (0.0039) | 1.81 (0.0012) | 0.28 (0.17) | 0.88 (0.03) |
| LTF | 1.84 (0.47) | 1.97 (0.073) | – | – | – |

As a result, the microglial share of GPNMB and C3 is lower in ALS than in control (shift −0.12 and −0.13), and the astrocyte share of LTF is higher (+0.35).

## Reading

- **Main source.** In ALS spinal cord, GPNMB and C3 transcripts come mainly from microglia. This answers the Tanaka et al. (2012) discrepancy at the RNA level for the main source.
- **ALS-induced expression elsewhere.** ALS also raises GPNMB and C3 in astrocytes and in oligodendrocyte-lineage nuclei, which is consistent with Tanaka's astrocytic GPNMB immunoreactivity.
- **LTF.** LTF is astrocytic in ALS cord.
- **ANXA1 and EMP3.** Both are mainly vascular.
- **Overall.** The module is multicellular at the gene level: its strongest ventral-horn genes come from microglia (GPNMB, C3), astrocytes (LTF) and vascular cells (ANXA1, EMP3).

## Caveats

- **Ambient RNA.** Ambient RNA from expanded microglia could inflate microglial transcripts in other nuclei. The non-myeloid gains in GPNMB and C3 should therefore be confirmed at the protein level.
- **Nuclear RNA.** Nuclear RNA may under-represent transcripts of secreted proteins.
- **Small control group.** There are only 4 controls.

## Files

- `w72_lineage_share.tsv`: ALS and control shares and shifts
- `w72_within_lineage_de.tsv`
- `w72_lineage_calls_by_donor.tsv`
- `inputs/`: lineage pseudobulk counts and units
- `w72_run.log`
