#!/usr/bin/env Rscript

# W59 (Amendment 5): confirmatory test of the frozen W58 discovery-agnostic
# consensus signature in the NYGC bulk hold-out, using the W55 inputs/helpers.
# 1 primary (W55 design), 2 secondary (+ myeloid/vascular scores),
# 3 sensitivity excluding NYGC donors fingerprint-matching signature-defining
# donors (Li, Pineda, Zelic), 4 gene-level sign concordance.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "run_w55_nygc_bulk_holdout.R"))

w59 <- file.path(root, "results", "W59_NYGC_consensus_signature_holdout")
sig_file <- file.path(w59, "w59_frozen_consensus_signature.tsv")
expected_sha <- "9a4bd65b0d17995081f5ff7f4071522b447d4aea3f103258eced23ef0c1480ce"
observed_sha <- sub(" .*$", "", system2("shasum", c("-a", "256", shQuote(sig_file)), stdout = TRUE))
if (!identical(observed_sha, expected_sha)) stop("Frozen signature checksum mismatch: ", observed_sha)
sig <- read.delim(sig_file, stringsAsFactors = FALSE)
sig_sets <- list(up_consensus = sig$gene[sig$direction == "up"], down_consensus = sig$gene[sig$direction == "down"])
extra <- list(
  myeloid_act = symbols_to_ensg(c("CD68", "TREM2", "TYROBP", "SPP1", "CD163", "ITGAX")),
  vascular = symbols_to_ensg(c("CLDN5", "PECAM1", "VWF", "PDGFRB"))
)
stopifnot(!any(unlist(c(marker_ensg, extra)) %in% unlist(sig_sets)))

run_both <- function(label, extra_markers = list()) {
  rows <- lapply(c("motor_cortex", "spinal_cord"), function(region) {
    r <- fit_region(region, extra_markers = extra_markers, gene_sets = sig_sets)
    list(sets = cbind(analysis = label, r$sets), genes = cbind(analysis = label, r$genes))
  })
  sets <- do.call(rbind, lapply(rows, `[[`, "sets"))
  sets$fry_BH <- p.adjust(sets$fry_p, method = "BH")
  sets$direction_met <- (startsWith(sets$module, "up") & sets$fry_direction == "Up") | (startsWith(sets$module, "down") & sets$fry_direction == "Down")
  sets$confirmed_region <- sets$direction_met & sets$fry_BH < 0.05
  list(sets = sets, genes = do.call(rbind, lapply(rows, `[[`, "genes")))
}

# Test 3 fingerprints: signature-defining donors (diagnosis class, sex, age).
norm_sex <- function(x) toupper(substr(as.character(x), 1, 1))
fp <- function(dx, sex, age) data.frame(dx = dx, sex = norm_sex(sex), age = as.numeric(age))
li <- read.delim(file.path(root, "results/W54_Li2023_motor_cortex_validation/w54_li2023_astrocyte_pseudobulk_samples.tsv"))
li <- li[li$disease %in% c("C9-ALS", "Control"), ]
pin <- read.delim(file.path(root, "results/W51_Pineda_BA4_discovery_matched_limma/w51_pineda_astrocyte_pseudobulk_samples.tsv"))
ba4 <- pin[pin$Region == "BA4" & substr(pin$sample_folder, 1, 6) %in% c("191112", "191114") & pin$Condition %in% c("ALS", "PN"), ]
ba9 <- pin[pin$Region == "BA9" & startsWith(pin$sample_folder, "210526") & pin$Condition %in% c("ALS", "PN"), ]
zel <- read.delim(file.path(root, "results/W58_multicohort_meta/w58_zelic_astrocyte_pseudobulk_samples.tsv"))
defining <- rbind(
  fp(ifelse(li$disease == "Control", "Control", "ALS"), li$sex, li$age),
  fp(ifelse(ba4$Condition == "PN", "Control", "ALS"), ba4$Sex, ba4$age),
  fp(ifelse(ba9$Condition == "PN", "Control", "ALS"), ba9$Sex, ba9$age),
  fp(ifelse(zel$disease == "Control", "Control", "ALS"), zel$sex, zel$age)
)
defining <- defining[!is.na(defining$age), ]
nygc_donors <- samples_all[!duplicated(samples_all$donor_id), ]
nygc_fp <- fp(ifelse(nygc_donors$disease_group == "Control", "Control", "ALS"), nygc_donors$sex_call, nygc_donors$age)
match_any <- vapply(seq_len(nrow(nygc_fp)), function(i) {
  x <- nygc_fp[i, ]
  if (is.na(x$age)) return(FALSE)
  age_ok <- abs(defining$age - x$age) <= 1 | (x$age == 90 & defining$age >= 89)
  any(defining$dx == x$dx & defining$sex == x$sex & age_ok)
}, logical(1))
overlap_excluded <- nygc_donors$donor_id[match_any]
write.table(data.frame(nygc_donor_id = overlap_excluded), file.path(w59, "w59_signature_cohort_overlap_exclusions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

primary <- run_both("1_primary_W55_design")
secondary <- run_both("2_secondary_myeloid_vascular", extra)
samples_full <- samples_all
samples_all <- samples_full[!samples_full$donor_id %in% overlap_excluded, ]
sensitivity <- run_both("3_sensitivity_signature_cohort_overlap_excluded")
samples_all <- samples_full

sets <- rbind(primary$sets, secondary$sets, sensitivity$sets)
genes <- rbind(primary$genes, secondary$genes, sensitivity$genes)
concord <- do.call(rbind, lapply(split(genes, list(genes$analysis, genes$region, genes$module), drop = TRUE), function(g) {
  expected <- if (startsWith(g$module[1], "up")) 1 else -1
  k <- sum(sign(g$logFC) == expected)
  data.frame(analysis = g$analysis[1], region = g$region[1], module = g$module[1], n_genes = nrow(g), n_expected_sign = k,
             frac_expected_sign = k / nrow(g), binom_p_one_sided = binom.test(k, nrow(g), 0.5, alternative = "greater")$p.value)
}))
summary_rows <- do.call(rbind, lapply(split(sets, sets$analysis), function(s) {
  data.frame(analysis = s$analysis[1],
             up_confirmed_any_region = any(s$confirmed_region[startsWith(s$module, "up")]),
             up_confirmed_both_regions = all(s$confirmed_region[startsWith(s$module, "up")]),
             down_confirmed_any_region = any(s$confirmed_region[startsWith(s$module, "down")]),
             down_confirmed_both_regions = all(s$confirmed_region[startsWith(s$module, "down")]))
}))

write.table(sets, file.path(w59, "w59_signature_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes, file.path(w59, "w59_signature_gene_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(concord, file.path(w59, "w59_signature_sign_concordance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_rows, file.path(w59, "w59_confirmation_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("signature sha256 verified; overlap-excluded NYGC donors:", length(overlap_excluded), "\n")
print(sets[, c("analysis", "region", "module", "n_samples", "n_als_donors", "n_control_donors", "n_module_genes_tested", "fry_direction", "fry_p", "fry_BH", "camera_p", "n_genes_concordant", "median_logFC", "confirmed_region")], digits = 3, row.names = FALSE)
print(concord, digits = 3, row.names = FALSE)
print(summary_rows, row.names = FALSE)
