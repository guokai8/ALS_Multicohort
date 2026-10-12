#!/usr/bin/env Rscript

# W51 post hoc: discovery-matched module test restricted to the 2019 BA4
# campaign, 191112 (ALS) vs 191114 (PN). Index wells show these dates share one
# plate, so cross-year campaign heterogeneity is removed; the processing-day
# confound between diagnosis groups remains and cannot be estimated.
#
# Plan fixed before results were inspected:
# - Primary: ~ disease + age + sex (5 of 14 PN lack PMI).
# - Sensitivity: W33 covariate rule with PMI (complete case).
# Tests as in W51: fry (BH across up/down), camera, gene-level sign test.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "run_w51_pineda_ba4_discovery_matched_limma.R"))

ba4 <- ba4_samples()
campaign <- ba4[substr(ba4$sample_folder, 1, 6) %in% c("191112", "191114"), ]
write_results(list(
  campaign2019_primary_no_PMI = run_analysis("campaign2019_primary_no_PMI", campaign, c("age", "sex")),
  campaign2019_with_PMI = run_analysis("campaign2019_with_PMI", campaign, c("age", "sex", "PMI"))
), "w51_campaign2019")
