#!/usr/bin/env Rscript

# W58 multi-cohort meta-analysis (plan v1.0 Section 5; Amendment 4).
# Per-cohort astrocyte pseudobulk limma-voom fits, module-level standardised
# effects pooled with metafor REML (Ulm age-corrected with propagated variance),
# genome-wide per-gene random-effects consensus with leave-one-cohort-out, and
# a cross-cohort reproducibility matrix with expression-matched nulls.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "w58_cohort_fits.R"))

cohorts <- list(load_li(), load_ulm(), load_pineda_ba4(), load_pineda_ba9(), load_zelic(), load_discovery())
fits <- lapply(cohorts, fit_cohort)
names(fits) <- vapply(fits, `[[`, "", "name")
mod <- do.call(rbind, lapply(fits, `[[`, "module"))

# ---- A4.3 module-level Ulm age correction ----
age_cohorts <- c("Li2023_BA4", "Pineda_BA4_2019", "Pineda_BA9_210526", "Zelic_cervicalSC", "NYGC_MTC_discovery")
age_pool <- do.call(rbind, lapply(names(seeds), function(m) {
  a <- mod[mod$module == m & mod$cohort %in% age_cohorts, ]
  r <- rma(yi = a$age_std, sei = a$se_age_std, method = "REML")
  data.frame(module = m, beta_age_std = as.numeric(r$b), se_beta_age_std = r$se, k = r$k)
}))
ulm_corr <- mod[mod$cohort == "Ulm_M1", ]
ulm_corr <- merge(ulm_corr, age_pool, by = "module")
ulm_corr$d <- ulm_corr$d - ulm_corr$beta_age_std * DELTA_AGE_ULM
ulm_corr$se_d <- sqrt(ulm_corr$se_d^2 + DELTA_AGE_ULM^2 * ulm_corr$se_beta_age_std^2)
ulm_corr$cohort <- "Ulm_M1_agecorrected"
mod_all <- rbind(mod, ulm_corr[, names(mod)])

# ---- A4.2 module-level pooling ----
pool <- function(m, set, label) {
  a <- mod_all[mod_all$module == m & mod_all$cohort %in% set, ]
  r <- rma(yi = a$d, sei = a$se_d, method = "REML")
  loo <- if (nrow(a) > 2) leave1out(r) else NULL
  data.frame(analysis = label, module = m, k = r$k, cohorts = paste(a$cohort, collapse = ";"), pooled_d = as.numeric(r$b),
             ci_low = r$ci.lb, ci_high = r$ci.ub, p = r$pval, tau2 = r$tau2, I2 = r$I2,
             loo_min_d = if (is.null(loo)) NA else min(loo$estimate), loo_max_p = if (is.null(loo)) NA else max(loo$pval),
             criterion_met = (if (m == "up") r$b > 0 else r$b < 0) & (r$ci.lb > 0 | r$ci.ub < 0))
}
sets <- list(
  primary_motor_independent = c("Li2023_BA4", "Ulm_M1_agecorrected", "Pineda_BA4_2019"),
  secondary_all_independent = c("Li2023_BA4", "Ulm_M1_agecorrected", "Pineda_BA4_2019", "Pineda_BA9_210526", "Zelic_cervicalSC"),
  sens_motor_ulm_uncorrected = c("Li2023_BA4", "Ulm_M1", "Pineda_BA4_2019"),
  sens_all_ulm_uncorrected = c("Li2023_BA4", "Ulm_M1", "Pineda_BA4_2019", "Pineda_BA9_210526", "Zelic_cervicalSC")
)
module_meta <- do.call(rbind, lapply(names(sets), function(l) do.call(rbind, lapply(names(seeds), pool, set = sets[[l]], label = l))))

# ---- A4.3/A4.4 gene-level: Ulm age correction and genome-wide meta ----
rma_safe <- function(yi, sei) {
  r <- tryCatch(rma(yi = yi, sei = sei, method = "REML"), error = function(e) NULL)
  if (is.null(r)) r <- tryCatch(rma(yi = yi, sei = sei, method = "DL"), error = function(e) NULL)
  r
}
age_long <- do.call(rbind, lapply(age_cohorts, function(n) cbind(cohort = n, fits[[n]]$age)))
age_gene <- do.call(rbind, lapply(split(age_long, age_long$gene), function(a) {
  if (nrow(a) < 3) return(NULL)
  r <- rma_safe(a$logFC, a$SE)
  if (is.null(r)) return(NULL)
  data.frame(gene = a$gene[1], beta_age = as.numeric(r$b), se_beta_age = r$se)
}))
ulm_gene <- merge(fits$Ulm_M1$disease, age_gene, by = "gene")
ulm_gene$logFC <- ulm_gene$logFC - ulm_gene$beta_age * DELTA_AGE_ULM
ulm_gene$SE <- sqrt(ulm_gene$SE^2 + DELTA_AGE_ULM^2 * ulm_gene$se_beta_age^2)
disease_long <- function(names_vec, ulm_corrected = TRUE) {
  do.call(rbind, lapply(names_vec, function(n) {
    tab <- if (n == "Ulm_M1" && ulm_corrected) ulm_gene[, c("gene", "logFC", "SE", "P", "AveExpr")] else fits[[n]]$disease
    cbind(cohort = n, tab)
  }))
}
gene_meta <- function(long, min_k) {
  rows <- lapply(split(long, long$gene), function(a) {
    if (nrow(a) < min_k) return(NULL)
    r <- rma_safe(a$logFC, a$SE)
    if (is.null(r)) return(NULL)
    data.frame(gene = a$gene[1], k = nrow(a), meta_logFC = as.numeric(r$b), meta_se = r$se, meta_p = r$pval,
               I2 = r$I2, frac_same_sign = mean(sign(a$logFC) == sign(as.numeric(r$b))))
  })
  res <- do.call(rbind, rows)
  res$meta_BH <- p.adjust(res$meta_p, method = "BH")
  res
}
consensus_run <- function(label, cohort_names, min_k) {
  long <- disease_long(cohort_names)
  main <- gene_meta(long, min_k)
  loo_bh <- sapply(cohort_names, function(drop) {
    g <- gene_meta(long[long$cohort != drop, ], max(3, min_k - 1))
    setNames(g$meta_BH, g$gene)[main$gene]
  })
  main$loo_max_BH <- apply(loo_bh, 1, function(x) if (all(is.na(x))) NA else max(x, na.rm = FALSE))
  main$consensus <- main$meta_BH < 0.05 & main$frac_same_sign >= 2 / 3 & !is.na(main$loo_max_BH) & main$loo_max_BH < 0.10
  main$analysis <- label
  main
}
all6 <- names(fits)
runs <- list(
  i_all_six = list(cohorts = all6, min_k = 4),
  ii_without_discovery = list(cohorts = setdiff(all6, "NYGC_MTC_discovery"), min_k = 4),
  iii_motor_cortex = list(cohorts = c("NYGC_MTC_discovery", "Li2023_BA4", "Ulm_M1", "Pineda_BA4_2019"), min_k = 4)
)
gm <- do.call(rbind, lapply(names(runs), function(l) consensus_run(l, runs[[l]]$cohorts, runs[[l]]$min_k)))
gm$symbol <- gene_map$Gene[match(gm$gene, gene_map$ENSEMBL)]
gm$module <- ifelse(gm$gene %in% seeds$up$gene_id, "up", ifelse(gm$gene %in% seeds$down$gene_id, "down", ""))

# ---- A4.5 reproducibility matrix ----
repro <- do.call(rbind, lapply(all6, function(src) do.call(rbind, lapply(setdiff(all6, src), function(tgt) {
  a <- fits[[src]]$disease; b <- fits[[tgt]]$disease
  common <- intersect(a$gene, b$gene)
  a <- a[match(common, a$gene), ]; b <- b[match(common, b$gene), ]
  top <- order(a$P)[1:min(200, length(common))]
  rep_flag <- sign(a$logFC) == sign(b$logFC) & b$P < 0.05
  obs <- mean(rep_flag[top])
  dec <- cut(a$AveExpr, quantile(a$AveExpr, 0:10 / 10), include.lowest = TRUE, labels = FALSE)
  null <- replicate(1000, {
    pick <- unlist(lapply(split(top, dec[top]), function(ix) {
      pool_ix <- which(dec == dec[ix[1]])
      sample(pool_ix, length(ix))
    }))
    mean(rep_flag[pick])
  })
  data.frame(source = src, target = tgt, n_common = length(common), n_top = length(top), replicated_fraction = obs,
             null_mean = mean(null), empirical_p = (sum(null >= obs) + 1) / 1001)
}))))

write.table(mod_all, file.path(out, "w58_module_effects_by_cohort.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(age_pool, file.path(out, "w58_module_age_effect_pool.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(module_meta, file.path(out, "w58_module_meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gm, gzfile(file.path(out, "w58_genome_wide_meta.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gm[gm$consensus, ], file.path(out, "w58_consensus_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(repro, file.path(out, "w58_reproducibility_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(mod_all[, c("cohort", "module", "n_als", "n_control", "n_genes", "d", "se_d")], digits = 3, row.names = FALSE)
print(module_meta[, c("analysis", "module", "k", "pooled_d", "ci_low", "ci_high", "p", "I2", "loo_max_p", "criterion_met")], digits = 3, row.names = FALSE)
print(table(gm$analysis, gm$consensus))
print(aggregate(cbind(replicated_fraction, null_mean, empirical_p) ~ source, data = repro, FUN = median), digits = 3)
