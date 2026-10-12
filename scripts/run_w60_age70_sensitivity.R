#!/usr/bin/env Rscript

# W60 age >= 70 sensitivity of the W58 Ulm age correction (plan 5.3; Amendment 6).
# beta_age is re-estimated in donors aged 70 or older in the eligible cohorts
# (Pineda BA4 2019, Pineda BA9 210526), Ulm is re-corrected, and the W58
# primary and secondary module pools are re-run.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "w58_cohort_fits.R"))

w58 <- file.path(res_dir, "W58_multicohort_meta")
out <- file.path(res_dir, "W60_age70_sensitivity")
dir.create(out, showWarnings = FALSE)
AGE_MIN <- 70
ELIGIBLE <- c("Pineda_BA4_2019", "Pineda_BA9_210526")

fits <- lapply(list(load_pineda_ba4(), load_pineda_ba9()), fit_cohort)
names(fits) <- vapply(fits, `[[`, "", "name")
committed <- read.delim(file.path(w58, "w58_module_effects_by_cohort.tsv"))

# ---- restricted (and full-age reference) age models on the W58 donor scores ----
age_fit <- function(cohort, m, min_age) {
  dd <- fits[[cohort]]$scores
  dd <- dd[dd$module == m & dd$age >= min_age, ]
  dd$disease <- factor(dd$disease, levels = c("Control", "ALS"))
  lmf <- lm(score ~ disease + age + sex, data = dd)
  sig <- summary(lmf)$sigma
  cf <- summary(lmf)$coefficients
  data.frame(cohort = cohort, module = m, min_age = min_age, n_als = sum(dd$disease == "ALS"),
             n_control = sum(dd$disease == "Control"), residual_df = lmf$df.residual,
             age_std = cf["age", 1] / sig, se_age_std = cf["age", 2] / sig)
}
grid <- expand.grid(cohort = ELIGIBLE, module = c("up", "down"), min_age = c(-Inf, AGE_MIN), stringsAsFactors = FALSE)
age_coef <- do.call(rbind, Map(age_fit, grid$cohort, grid$module, grid$min_age))
age_coef$age_set <- ifelse(is.finite(age_coef$min_age), "age_70_plus", "all_ages")

# The full-age refit must reproduce the committed W58 age coefficients.
ref <- merge(age_coef[age_coef$age_set == "all_ages", ], committed, by = c("cohort", "module"), suffixes = c("", "_w58"))
max_diff <- max(abs(ref$age_std - ref$age_std_w58), abs(ref$se_age_std - ref$se_age_std_w58))
if (max_diff > 1e-10) stop("Full-age refit does not reproduce W58 age coefficients: ", max_diff)

# ---- pooled beta_age per estimate ----
estimates <- list(
  primary_age70_BA4_BA9 = list(cohorts = ELIGIBLE, age_set = "age_70_plus"),
  secondary_age70_BA4_only = list(cohorts = "Pineda_BA4_2019", age_set = "age_70_plus"),
  reference_all_ages_BA4_BA9 = list(cohorts = ELIGIBLE, age_set = "all_ages")
)
w58_pool <- read.delim(file.path(w58, "w58_module_age_effect_pool.tsv"))
age_pool <- do.call(rbind, c(lapply(names(estimates), function(e) do.call(rbind, lapply(c("up", "down"), function(m) {
  a <- age_coef[age_coef$module == m & age_coef$cohort %in% estimates[[e]]$cohorts & age_coef$age_set == estimates[[e]]$age_set, ]
  if (nrow(a) == 1) return(data.frame(estimate = e, module = m, k = 1, beta_age_std = a$age_std, se_beta_age_std = a$se_age_std))
  r <- rma(yi = a$age_std, sei = a$se_age_std, method = "REML")
  data.frame(estimate = e, module = m, k = r$k, beta_age_std = as.numeric(r$b), se_beta_age_std = r$se)
}))), list(data.frame(estimate = "w58_primary_five_cohorts", module = w58_pool$module, k = w58_pool$k,
                      beta_age_std = w58_pool$beta_age_std, se_beta_age_std = w58_pool$se_beta_age_std))))

# ---- Ulm re-correction and re-pooling (A4.2 pools, Amendment 6 step 5) ----
ulm <- committed[committed$cohort == "Ulm_M1", ]
others <- committed[committed$cohort %in% c("Li2023_BA4", "Pineda_BA4_2019", "Pineda_BA9_210526", "Zelic_cervicalSC"), ]
pool_sets <- list(primary_motor_independent = c("Li2023_BA4", "Ulm_corr", "Pineda_BA4_2019"),
                  secondary_all_independent = c("Li2023_BA4", "Ulm_corr", "Pineda_BA4_2019", "Pineda_BA9_210526", "Zelic_cervicalSC"))
ulm_corr <- list(); meta <- list()
for (e in unique(age_pool$estimate)) {
  for (m in c("up", "down")) {
    b <- age_pool[age_pool$estimate == e & age_pool$module == m, ]
    u <- ulm[ulm$module == m, ]
    d_corr <- u$d - b$beta_age_std * DELTA_AGE_ULM
    se_corr <- sqrt(u$se_d^2 + DELTA_AGE_ULM^2 * b$se_beta_age_std^2)
    ulm_corr[[paste(e, m)]] <- data.frame(estimate = e, module = m, d_raw = u$d, d_corrected = d_corr, se_corrected = se_corr)
    for (ps in names(pool_sets)) {
      a <- others[others$module == m, c("cohort", "d", "se_d")]
      a <- rbind(a, data.frame(cohort = "Ulm_corr", d = d_corr, se_d = se_corr))
      a <- a[a$cohort %in% pool_sets[[ps]], ]
      r <- rma(yi = a$d, sei = a$se_d, method = "REML")
      loo <- leave1out(r)
      meta[[paste(e, m, ps)]] <- data.frame(estimate = e, pool = ps, module = m, k = r$k, pooled_d = as.numeric(r$b),
        ci_low = r$ci.lb, ci_high = r$ci.ub, p = r$pval, I2 = r$I2, loo_min_d = min(loo$estimate), loo_max_p = max(loo$pval),
        criterion_met = (if (m == "up") r$b > 0 else r$b < 0) & (r$ci.lb > 0 | r$ci.ub < 0))
    }
  }
}
ulm_corr <- do.call(rbind, ulm_corr)
meta <- do.call(rbind, meta)

# The w58_primary_five_cohorts row must reproduce the committed W58 pools.
w58_meta <- read.delim(file.path(w58, "w58_module_meta.tsv"))
chk <- merge(meta[meta$estimate == "w58_primary_five_cohorts", ], w58_meta, by.x = c("pool", "module"), by.y = c("analysis", "module"))
if (nrow(chk) != 4 || max(abs(chk$pooled_d.x - chk$pooled_d.y)) > 1e-8) stop("W58 pools not reproduced")

write.table(age_coef, file.path(out, "w60_age_coefficients.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(age_pool, file.path(out, "w60_age_effect_pool.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ulm_corr, file.path(out, "w60_ulm_corrected_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta, file.path(out, "w60_module_meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("full-age refit vs W58 max diff:", max_diff, "\n")
print(age_coef, digits = 3, row.names = FALSE)
print(age_pool, digits = 3, row.names = FALSE)
print(ulm_corr, digits = 3, row.names = FALSE)
print(meta[, c("estimate", "pool", "module", "k", "pooled_d", "ci_low", "ci_high", "p", "I2", "loo_max_p", "criterion_met")], digits = 3, row.names = FALSE)
