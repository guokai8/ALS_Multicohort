# Amendment 4: Section 5 implementation details (pre-meta-analysis)

Date: 2026-10-11
Amends: plan v1.0 (sha256 da63b81a...bfb3bd), Section 5
Status: committed **before** any pooled (meta-analytic) estimate, genome-wide meta-analysis or reproducibility matrix was computed. Per-cohort module results from W51-W57 are already known (plan Section 2 and W55-W57 READMEs), so Section 5 remains a consolidation analysis, as the plan states.

## A4.1 Cohorts, units and designs

| Cohort | Counts source | Unit | Design | Role |
|---|---|---|---|---|
| Li 2023 BA4 | W54 marker-rule astrocyte pseudobulk | Donor | `~ disease + age + sex` | Independent motor cortex |
| Ulm M1 | W52 donor x chip pseudobulk | Donor x chip (blocked) for genes; donor for module | `~ disease + sex + chip` (genes); `~ disease + sex + main_chip` (module) | Independent motor cortex |
| Pineda BA4 2019 campaign (191112 vs 191114) | W51 pseudobulk | Donor | `~ disease + age + sex` | Independent motor cortex |
| Pineda BA9 batch 210526 | W56 Astrocyte lineage counts | Donor | `~ disease + age + sex` | Independent, cross-region |
| Zelic 2025 cervical SC (GSE287257) | **New:** W54 marker rule applied to all 114,161 nuclei of `GSE287257_als_spinal.h5ad`, one pseudobulk per donor | Donor | `~ disease + age + sex`; PMI/RIN added by the W33 rule from `GSE287257.xlsx` | Independent, spinal |
| NYGC MTC discovery | W57 rebuilt astrocyte pseudobulk | Library = donor | W33 design | H3 only (with and without) |
| Masrori 2025 | Summary DE | — | — | Direction only, not pooled |

All cohorts use filterByExpr(design), TMM, voom and eBayes(trend). Disease means ALS vs control in each cohort; Li uses C9-ALS.

## A4.2 Module-level effect per cohort (H1, H2)

- **Module score.** On the cohort's voom log-CPM (unblocked voom; Ulm units averaged to donor), the module score is the mean over module genes present.
- **Model.** A donor-level linear model, `score ~ disease + covariates`, with the cohort covariates above (Ulm: sex + main chip).
- **Effect.** Standardised effect d = beta_disease / sigma_residual, with SE_d = SE(beta_disease) / sigma_residual.
- **Pooling.** REML random effects with `metafor::rma`:
  - **Primary:** independent motor cortex (Li, Ulm age-corrected, Pineda BA4 2019).
  - **Secondary:** all independent cohorts (adding Pineda BA9 and Zelic).
  - Report tau^2, I^2 and leave-one-cohort-out.
  - Ulm uncorrected is a sensitivity analysis.
- **Primary criterion for "cross-cohort single-nucleus consistency" of a module:** the pooled primary d has the hypothesised sign and a 95% CI excluding 0.

## A4.3 Ulm age correction

- **Module level.** beta_age (per year, standardised) comes from the same donor-level model in every age-available cohort (Li, Pineda BA4 2019, Pineda BA9, Zelic, discovery), pooled with `rma`.
  - d_Ulm,corr = d_Ulm - beta_age x (-11.1)
  - var_corr = var + 11.1^2 x var(beta_age)
- **Gene level (H3).** The same procedure is applied per gene, with logFC and SE from the limma age coefficient. The age-corrected Ulm estimates are primary in H3; uncorrected estimates are a sensitivity analysis.

## A4.4 Genome-wide consensus (H3)

- **ID harmonisation.** Ulm and Zelic symbols are mapped to ENSG with the Pineda GENCODE row metadata. Unmapped or ambiguous symbols are dropped.
- **Per-gene meta-analysis.** `metafor::rma` (REML) over per-cohort limma logFC with SE = logFC / t. Only genes tested in at least 4 cohorts are included. BH is applied across genes.
- **Consensus gene:** meta BH < 0.05, the same sign in at least two-thirds of contributing cohorts, and meta BH < 0.10 in every leave-one-cohort-out refit.
- **Runs:**
  - (i) all six cohorts, primary;
  - (ii) without discovery;
  - (iii) motor cortex only (discovery, Li, Ulm, Pineda BA4 2019).

## A4.5 Cross-cohort reproducibility matrix (H5)

- For each source cohort, take the top 200 genes by P among genes tested in both source and target.
- Replication in a target cohort means the same sign and P < 0.05.
- **Null:** 1,000 random gene sets of the same size, matched on source-cohort AveExpr deciles. The empirical P is (exceedances + 1) / 1001.
- Seed: 20261011.

## A4.6 Reporting

All pooled estimates, forest-plot inputs, per-gene meta tables and the matrix go to `results/W58_multicohort_meta/`. The criterion in A4.2 is reported whatever the outcome.
