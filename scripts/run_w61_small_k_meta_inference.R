#!/usr/bin/env Rscript

# W61 small-k inference for the W58/W60 module meta-analyses (Amendment 7; post hoc,
# review-prompted, already observed). REML with the z test (reproducing W58/W60),
# REML with the Knapp-Hartung adjustment, and 95% prediction intervals.

suppressPackageStartupMessages(library(metafor))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
res_dir <- file.path(p01_project_root(), "results")
out <- file.path(res_dir, "W61_small_k_meta_inference")
dir.create(out, showWarnings = FALSE)
DELTA_AGE_ULM <- -11.1

eff <- read.delim(file.path(res_dir, "W58_multicohort_meta", "w58_module_effects_by_cohort.tsv"))
w58_meta <- read.delim(file.path(res_dir, "W58_multicohort_meta", "w58_module_meta.tsv"))
u70 <- read.delim(file.path(res_dir, "W60_age70_sensitivity", "w60_ulm_corrected_effects.tsv"))
w60_meta <- read.delim(file.path(res_dir, "W60_age70_sensitivity", "w60_module_meta.tsv"))

# W60 primary Ulm estimate (beta_age from donors aged >= 70, Pineda BA4 + BA9)
u70 <- u70[u70$estimate == "primary_age70_BA4_BA9", ]
eff <- rbind(eff[, c("cohort", "module", "d", "se_d")],
             data.frame(cohort = "Ulm_M1_age70", module = u70$module, d = u70$d_corrected, se_d = u70$se_corrected))

pools <- list(
  primary_motor_independent = c("Li2023_BA4", "Ulm_M1_agecorrected", "Pineda_BA4_2019"),
  secondary_all_independent = c("Li2023_BA4", "Ulm_M1_agecorrected", "Pineda_BA4_2019", "Pineda_BA9_210526", "Zelic_cervicalSC"),
  sens_motor_ulm_uncorrected = c("Li2023_BA4", "Ulm_M1", "Pineda_BA4_2019"),
  sens_all_ulm_uncorrected = c("Li2023_BA4", "Ulm_M1", "Pineda_BA4_2019", "Pineda_BA9_210526", "Zelic_cervicalSC"),
  w60_primary_motor_age70 = c("Li2023_BA4", "Ulm_M1_age70", "Pineda_BA4_2019")
)

rows <- list()
for (p in names(pools)) for (m in c("up", "down")) {
  a <- eff[eff$module == m & eff$cohort %in% pools[[p]], ]
  if (nrow(a) != length(pools[[p]])) stop("missing cohort in pool ", p)
  z <- rma(yi = a$d, sei = a$se_d, method = "REML")
  hk <- rma(yi = a$d, sei = a$se_d, method = "REML", test = "knha")
  pr <- predict(z)
  hit <- function(r) (if (m == "up") r$b > 0 else r$b < 0) & (r$ci.lb > 0 | r$ci.ub < 0)
  rows[[paste(p, m)]] <- data.frame(pool = p, module = m, k = z$k, pooled_d = as.numeric(z$b),
    z_ci_low = z$ci.lb, z_ci_high = z$ci.ub, z_p = z$pval, hksj_ci_low = hk$ci.lb, hksj_ci_high = hk$ci.ub,
    hksj_p = hk$pval, pi_low = pr$pi.lb, pi_high = pr$pi.ub, tau2 = z$tau2, I2 = z$I2,
    criterion_met_z = hit(z), criterion_met_hksj = hit(hk))
}
tab <- do.call(rbind, rows)

# The z-test rows must reproduce the committed W58 and W60 pools.
ref <- rbind(data.frame(pool = w58_meta$analysis, module = w58_meta$module, d_ref = w58_meta$pooled_d, p_ref = w58_meta$p),
             data.frame(pool = "w60_primary_motor_age70", module = w60_meta$module, d_ref = w60_meta$pooled_d, p_ref = w60_meta$p)[
               w60_meta$estimate == "primary_age70_BA4_BA9" & w60_meta$pool == "primary_motor_independent", ])
chk <- merge(tab, ref, by = c("pool", "module"))
if (nrow(chk) != nrow(tab)) stop("reference pools not all found")
max_diff <- max(abs(chk$pooled_d - chk$d_ref), abs(chk$z_p - chk$p_ref))
if (max_diff > 1e-8) stop("z-test pools do not reproduce W58/W60: ", max_diff)

write.table(tab, file.path(out, "w61_small_k_meta_inference.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("reproduction max diff vs W58/W60:", max_diff, "\n")
print(tab[, c("pool", "module", "k", "pooled_d", "z_p", "hksj_ci_low", "hksj_ci_high", "hksj_p", "pi_low", "pi_high", "criterion_met_z", "criterion_met_hksj")],
      digits = 3, row.names = FALSE)
