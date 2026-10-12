#!/usr/bin/env Rscript

# Figure source data for the W55 clinical association (Fig. 5c). Re-runs the
# W55 primary fits and the W55 clinical models unchanged, then saves the
# donor-level scores and coefficient confidence intervals that W55 did not write.
# The beta and P values must equal the committed w55_clinical_association.tsv.

suppressPackageStartupMessages(library(survival))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "run_w55_nygc_bulk_holdout.R"))

rows <- list(); coefs <- list()
for (region in c("motor_cortex", "spinal_cord")) {
  res <- fit_region(region)
  v <- res$v; s <- res$s
  for (m in names(modules)) {
    g <- intersect(modules[[m]], rownames(v))
    raw <- colMeans(t(scale(t(v$E[g, , drop = FALSE]))))
    resid <- residuals(lm(raw ~ subregion + rin + platform + astro_score + neuron_score + oligo_score + micro_score, data = s))
    d <- aggregate(data.frame(score = resid), by = list(donor_id = s$donor_id), FUN = mean)
    d <- merge(d, s[!duplicated(s$donor_id), ], by = "donor_id")
    d <- d[d$disease == "ALS", ]
    dur_txt <- as.character(d$disease_duration_in_months)
    d$censored <- grepl("^>", dur_txt)
    d$duration <- suppressWarnings(as.numeric(ifelse(d$censored, sub("^>([0-9.]+).*", "\\1", dur_txt), dur_txt)))
    d$onset_age <- suppressWarnings(as.numeric(sub("90 or [Oo]lder", "90", d$age_at_symptom_onset)))
    d$onset_site <- factor(ifelse(d$site_of_motor_onset %in% c("Limb", "Bulbar"), d$site_of_motor_onset, "Other_Unknown"))
    maxcall <- suppressWarnings(as.numeric(d$c9orf72_max_call))
    d$c9 <- factor(ifelse(d$c9orf72_rp_pcr %in% "Yes" | (!is.na(maxcall) & maxcall >= 30), "C9", "nonC9"))
    prim <- d[!d$censored & !is.na(d$duration) & d$duration > 0 & !is.na(d$onset_age), ]
    fit <- lm(log(duration) ~ score + onset_age + sex + onset_site + c9, data = prim)
    ci <- confint(fit)["score", ] * sd(prim$score)
    coefs[[paste(region, m)]] <- data.frame(region = region, module = m, n = nrow(prim),
      beta_per_sd = coef(fit)["score"] * sd(prim$score), ci_low = ci[1], ci_high = ci[2],
      p = summary(fit)$coefficients["score", 4])
    rows[[paste(region, m)]] <- data.frame(region = region, module = m, donor_id = prim$donor_id,
      score_sd = prim$score / sd(prim$score), duration_months = prim$duration)
  }
}
coef_tab <- do.call(rbind, coefs)
committed <- read.delim(file.path(out, "w55_clinical_association.tsv"))
chk <- merge(coef_tab, committed, by = c("region", "module"))
max_diff <- max(abs(chk$beta_per_sd - chk$lm_beta_log_duration_per_sd), abs(chk$p - chk$lm_p))
if (max_diff > 1e-10) stop("Clinical source data do not reproduce committed W55 values: ", max_diff)
write.table(coef_tab, file.path(out, "w55_clinical_coefficients_with_ci.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, rows), file.path(out, "w55_clinical_donor_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("clinical source data written; max diff vs committed =", max_diff, "\n")
