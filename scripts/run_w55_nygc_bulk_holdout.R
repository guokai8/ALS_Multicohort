#!/usr/bin/env Rscript

# W55 NYGC bulk hold-out (plan v1.0 Section 4; Amendment 1 A3-A6).
# Confirmatory tests of the frozen W34 modules in NYGC bulk motor cortex and
# spinal cord after excluding GSE330130 donors (build_w55_nygc_holdout_inputs.py).
# Primary: blocked fry, 2 modules x 2 regions, BH over 4. Secondary: camera.
# Sensitivity: drop PMI. Clinical (H4): within ALS donors, log duration and Cox.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(survival)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
out <- file.path(root, "results", "W55_NYGC_bulk_holdout")
seed_dir <- file.path(root, "results", "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")

counts_all <- as.matrix(read.delim(gzfile(file.path(out, "w55_holdout_counts.tsv.gz")), row.names = 1, check.names = FALSE))
samples_all <- read.delim(file.path(out, "w55_holdout_samples.tsv"), stringsAsFactors = FALSE)
modules <- lapply(c(up = "up", down = "down"), function(d) {
  read.delim(file.path(seed_dir, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", d, ".tsv")))$gene_id
})
# W54 marker lists (symbols), mapped to ENSG with a local GENCODE-based table.
markers_symbol <- list(
  astro = c("AQP4", "GJA1", "SLC1A3", "SLC1A2", "ALDH1L1", "SLC4A4"),
  neuron = c("SNAP25", "RBFOX3", "SYT1", "GAD1", "SLC17A7", "NRGN"),
  oligo = c("MBP", "PLP1", "MOBP", "MOG", "ST18"),
  micro = c("CSF1R", "P2RY12", "CX3CR1", "PTPRC")
)
gene_map <- read.delim(file.path(p01_source_data_root(), "Other_Datasets/Pineda_snRNAseq/synapse_download/Processed Data/191112_ALS_110_snRNA-B9/row_metadata.tsv"), stringsAsFactors = FALSE)
symbols_to_ensg <- function(sym) {
  ids <- gene_map$ENSEMBL[match(sym, gene_map$Gene)]
  if (anyNA(ids)) stop("Unmapped marker symbols: ", paste(sym[is.na(ids)], collapse = ", "))
  ids
}
marker_ensg <- lapply(markers_symbol, symbols_to_ensg)
stopifnot(!any(unlist(marker_ensg) %in% unlist(modules)))

prepare <- function(region, drop_pmi = FALSE, markers = marker_ensg) {
  s <- samples_all[samples_all$region == region & samples_all$disease_group %in% c("ALS", "Control"), ]
  s$disease <- factor(ifelse(s$disease_group == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
  s$sex <- factor(s$sex_call)
  s$subregion <- factor(s$tissue_clean)
  s$platform <- factor(s$rna_platform)
  needed <- c("age", "sex", "rin", "platform", "subregion", if (!drop_pmi) "pmi")
  s <- s[complete.cases(s[, needed]), ]
  y <- DGEList(counts = counts_all[, s$rna_id])
  y <- calcNormFactors(y[filterByExpr(y, group = s$disease), , keep.lib.sizes = FALSE], method = "TMM")
  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  for (lin in names(markers)) {
    g <- intersect(markers[[lin]], rownames(logcpm))
    s[[paste0(lin, "_score")]] <- colMeans(t(scale(t(logcpm[g, , drop = FALSE]))))
  }
  list(s = s, y = y, needed = needed)
}

# extra_markers adds composition scores (named ENSG lists); gene_sets defaults to
# the frozen modules. Sets whose name starts with "up" are expected Up.
fit_region <- function(region, drop_pmi = FALSE, extra_markers = list(), gene_sets = modules) {
  markers <- c(marker_ensg, extra_markers)
  p <- prepare(region, drop_pmi, markers)
  s <- p$s
  rhs <- c(p$needed, paste0(names(markers), "_score"))
  design <- model.matrix(as.formula(paste("~ disease +", paste(rhs, collapse = " + "))), data = s)
  v <- voom(p$y, design)
  cor1 <- duplicateCorrelation(v, design, block = s$donor_id)$consensus.correlation
  v <- voom(p$y, design, block = s$donor_id, correlation = cor1)
  cons <- duplicateCorrelation(v, design, block = s$donor_id)$consensus.correlation
  fit <- eBayes(lmFit(v, design, block = s$donor_id, correlation = cons), trend = TRUE)
  index <- lapply(gene_sets, intersect, rownames(v))
  fr <- fry(v, index = index, design = design, contrast = "diseaseALS", block = s$donor_id, correlation = cons, sort = "none")
  ca <- camera(v, index = index, design = design, contrast = "diseaseALS", sort = FALSE)
  tab <- topTable(fit, coef = "diseaseALS", number = Inf, sort.by = "none")
  sets <- do.call(rbind, lapply(names(index), function(m) {
    lfc <- tab[index[[m]], "logFC"]
    expected <- if (startsWith(m, "up")) 1 else -1
    data.frame(region = region, module = m, pmi_in_model = !drop_pmi,
               n_samples = ncol(v), n_als_donors = length(unique(s$donor_id[s$disease == "ALS"])),
               n_control_donors = length(unique(s$donor_id[s$disease == "Control"])),
               consensus_correlation = cons, n_module_genes_tested = length(index[[m]]),
               fry_direction = fr[m, "Direction"], fry_p = fr[m, "PValue"],
               camera_direction = ca[m, "Direction"], camera_p = ca[m, "PValue"],
               n_genes_concordant = sum(sign(lfc) == expected), median_logFC = median(lfc))
  }))
  genes <- data.frame(region = region, pmi_in_model = !drop_pmi, gene_id = unlist(index),
                      module = rep(names(index), lengths(index)), tab[unlist(index), c("logFC", "AveExpr", "t", "P.Value")])
  list(sets = sets, genes = genes, v = v, s = s)
}

clinical <- function(res, region) {
  v <- res$v
  s <- res$s
  rows <- list()
  for (m in names(modules)) {
    g <- intersect(modules[[m]], rownames(v))
    raw <- colMeans(t(scale(t(v$E[g, , drop = FALSE]))))
    resid <- residuals(lm(raw ~ subregion + rin + platform + astro_score + neuron_score + oligo_score + micro_score, data = s))
    d <- aggregate(data.frame(score = resid), by = list(donor_id = s$donor_id), FUN = mean)
    meta <- s[!duplicated(s$donor_id), ]
    d <- merge(d, meta, by = "donor_id")
    d <- d[d$disease == "ALS", ]
    dur_txt <- as.character(d$disease_duration_in_months)
    d$censored <- grepl("^>", dur_txt)
    d$duration <- suppressWarnings(as.numeric(ifelse(d$censored, sub("^>([0-9.]+).*", "\\1", dur_txt), dur_txt)))
    d$onset_age <- suppressWarnings(as.numeric(sub("90 or [Oo]lder", "90", d$age_at_symptom_onset)))
    d$onset_site <- factor(ifelse(d$site_of_motor_onset %in% c("Limb", "Bulbar"), d$site_of_motor_onset, "Other_Unknown"))
    maxcall <- suppressWarnings(as.numeric(d$c9orf72_max_call))
    d$c9 <- factor(ifelse(d$c9orf72_rp_pcr %in% "Yes" | (!is.na(maxcall) & maxcall >= 30), "C9", "nonC9"))
    prim <- d[!d$censored & !is.na(d$duration) & d$duration > 0 & !is.na(d$onset_age), ]
    lmfit <- lm(log(duration) ~ score + onset_age + sex + onset_site + c9, data = prim)
    cx <- d[!is.na(d$duration) & d$duration > 0 & !is.na(d$onset_age), ]
    coxfit <- coxph(Surv(duration, !censored) ~ score + onset_age + sex + onset_site + c9, data = cx)
    rows[[m]] <- data.frame(region = region, module = m, n_lm = nrow(prim),
                            lm_beta_log_duration_per_sd = coef(lmfit)["score"] * sd(prim$score),
                            lm_p = summary(lmfit)$coefficients["score", 4],
                            n_cox = nrow(cx), n_events = sum(!cx$censored),
                            cox_hr_per_sd = exp(coef(coxfit)["score"] * sd(cx$score)),
                            cox_p = summary(coxfit)$coefficients["score", "Pr(>|z|)"])
  }
  do.call(rbind, rows)
}

main <- function() {
  primary <- lapply(c("motor_cortex", "spinal_cord"), fit_region)
  sens <- lapply(c("motor_cortex", "spinal_cord"), fit_region, drop_pmi = TRUE)
  sets <- do.call(rbind, lapply(primary, `[[`, "sets"))
  sets$fry_BH <- p.adjust(sets$fry_p, method = "BH")
  sets$hypothesis_direction_met <- (sets$module == "up" & sets$fry_direction == "Up") | (sets$module == "down" & sets$fry_direction == "Down")
  sets$confirmed_BH05 <- sets$hypothesis_direction_met & sets$fry_BH < 0.05
  sens_sets <- do.call(rbind, lapply(sens, `[[`, "sets"))
  sens_sets$fry_BH <- p.adjust(sens_sets$fry_p, method = "BH")
  clin <- do.call(rbind, mapply(clinical, primary, c("motor_cortex", "spinal_cord"), SIMPLIFY = FALSE))
  clin$lm_BH <- p.adjust(clin$lm_p, method = "BH")
  clin$cox_BH <- p.adjust(clin$cox_p, method = "BH")

  write.table(sets, file.path(out, "w55_primary_module_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(sens_sets, file.path(out, "w55_sensitivity_no_pmi_module_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(do.call(rbind, c(lapply(primary, `[[`, "genes"), lapply(sens, `[[`, "genes"))), file.path(out, "w55_module_gene_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(clin, file.path(out, "w55_clinical_association.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(sets, digits = 3)
  print(sens_sets[, c("region", "module", "n_samples", "fry_direction", "fry_p", "fry_BH", "camera_p")], digits = 3)
  print(clin, digits = 3)
}

# Run only when executed directly, so W56 can source the helpers.
if (sys.nframe() == 0) main()
