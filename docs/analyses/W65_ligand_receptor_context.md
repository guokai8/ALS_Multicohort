# W65: ligand-receptor context of the frozen module genes

This analysis is **exploratory and hypothesis-generating**. It was prespecified in Amendment 8 (committed before it was run). No P values are used, and no mechanism is inferred.

- Script: `scripts/run_w65_ligand_receptor_context.R`
- Resource: liana Consensus. 95 ligand-receptor pairs involve a frozen module gene as ligand, receptor or subunit (`w65_module_gene_lr_pairs.tsv`).
- Cohorts: W56 lineage pseudobulk from Pineda BA9 batch 210526 (8 lineages) and Ulm/DZNE (6 lineages; donor × chip profiles summed to donor level).
- Expression threshold: mean CPM >= 1 for the ligand in the sender lineage and for the receptor in the receiver lineage. For complexes, the minimum over subunits is used.
- Disease change: the ligand log2FC in the sender plus the receptor log2FC in the receiver, from per-lineage limma-voom `~ disease + sex`.
- Implementation detail not fixed in Amendment 8: for complexes, the log2FC is the mean over subunits.

## Results

**Coverage.**
- Expressed pairs: Pineda 76 of 95 (2,646 sender-receiver rows); Ulm 61 of 95 (1,346 rows).
- Module genes with at least one expressed pair: ANG, ANGPT1, ANXA1, C3, CDH2, EFNA5, EGF, GPNMB, KITLG, LTF and RAMP3.
- Sender-receiver rows evaluable in both cohorts: 1,214 (Ulm has no vascular lineages).

**Direction of the pair score change** (`w65_pairs_in_both_cohorts.tsv`):
- 269 of the 1,214 rows have the same sign in both cohorts: 224 lower in ALS and 45 higher.
- **Only ever higher in ALS** (whenever the sign agrees between cohorts): C3-CR1 (complement receptor 1), the neuropeptide-receptor pairs onto RAMP3 (CRH, VIP, CALCB, PTHLH), ADAM10-GPNMB, ANXA1-EGFR, ANXA1-DYSF and KITLG-EPOR.
- **Only ever lower in ALS:**
  - C3 onto C3AR1, CD81, LRP1, NRP1, ITGAM, GRM7 and ADRA2A;
  - KITLG-KIT and KITLG-LRRN3;
  - ANGPT1 onto TEK, ITGA5 and ITGB1;
  - NCAN-CDH2 and CDH1-CDH2.

  Many of these receptors are expressed by neurons or oligodendrocytes.
- **Both directions** across sender-receiver combinations: EFNA5 with most EPH receptors, EGF with ERBB3/ERBB4, LDLR, LRP2 and AXL, and C3 with CD46 and ITGAX.

## Reading

- The module genes take part in complement (C3 with its receptors), ephrin (EFNA5-EPH), stem-cell factor (KITLG-KIT), angiopoietin (ANGPT1-TEK and integrins) and EGF-ERBB communication across glial, vascular and neuronal lineages.
- The pair score mixes ligand and receptor changes in different lineages, and in disease tissue it can reflect changes in lineage composition or state. It does not show signalling.
- These pairs are candidates for spatial or functional follow-up only.
