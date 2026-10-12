#!/usr/bin/env Rscript

# W52 post hoc: bound the bias that the unadjustable Ulm age imbalance could
# introduce. Ulm donor ages cannot be linked to metadata IDs, but the released
# Samples_Summary gives group ages: mean ALS 61.2 vs HC 72.3 years, so
# delta = ALS - HC = -11.1 years. Per-year age coefficients for the module
# genes are estimated in Pineda BA4 astrocytes (~ disease + age + sex, the
# W51 no-PMI donor set), and the expected Ulm bias is age_coef * delta.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "run_w51_pineda_ba4_discovery_matched_limma.R"))

delta_age <- -11.1
ulm_out <- file.path(root, "results", "W52_Ulm_motor_cortex_validation")

s <- ba4_samples()
s$disease <- factor(ifelse(s$Condition == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
s$sex <- factor(s$Sex)
s <- s[!is.na(s$age), ]
design <- model.matrix(~ disease + age + sex, data = s)
y <- DGEList(counts = as.matrix(counts_all[, s$sample_folder]))
y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
fit <- eBayes(lmFit(voom(y, design), design), trend = TRUE)
age_tab <- topTable(fit, coef = "age", number = Inf, sort.by = "none")
age_tab$gene_id <- rownames(age_tab)

ulm <- read.delim(file.path(ulm_out, "w52_module_gene_results.tsv"), stringsAsFactors = FALSE)
ulm <- ulm[ulm$analysis_id == "primary_donor_chip_blocked", c("gene", "logFC")]
names(ulm)[2] <- "ulm_observed_logFC"

bound <- merge(modules, age_tab[, c("gene_id", "logFC", "P.Value")], by = "gene_id", all.x = TRUE)
names(bound)[names(bound) == "logFC"] <- "pineda_age_logFC_per_year"
names(bound)[names(bound) == "P.Value"] <- "pineda_age_P"
bound$expected_ulm_age_bias <- bound$pineda_age_logFC_per_year * delta_age
bound <- merge(bound, ulm, by = "gene", all.x = TRUE)
bound$ulm_age_corrected_logFC <- bound$ulm_observed_logFC - bound$expected_ulm_age_bias
bound <- bound[order(bound$module, -bound$discovery_logFC), ]
write.table(bound, file.path(ulm_out, "w52_post_hoc_age_bias_bound.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summ <- do.call(rbind, lapply(c("up", "down"), function(m) {
  b <- bound[bound$module == m & !is.na(bound$expected_ulm_age_bias) & !is.na(bound$ulm_observed_logFC), ]
  expected_sign <- if (m == "up") 1 else -1
  data.frame(module = m, n_genes = nrow(b),
             median_discovery_logFC = median(b$discovery_logFC),
             median_expected_age_bias = median(b$expected_ulm_age_bias),
             median_ulm_observed = median(b$ulm_observed_logFC),
             median_ulm_age_corrected = median(b$ulm_age_corrected_logFC),
             n_concordant_after_correction = sum(sign(b$ulm_age_corrected_logFC) == expected_sign))
}))
write.table(summ, file.path(ulm_out, "w52_post_hoc_age_bias_bound_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, digits = 3)
print(bound[, c("gene", "module", "discovery_logFC", "pineda_age_logFC_per_year", "expected_ulm_age_bias", "ulm_observed_logFC", "ulm_age_corrected_logFC")], digits = 3, row.names = FALSE)
