#!/usr/bin/env Rscript

# W70 (Amendment 13): territory x disease interaction of the frozen ALS-up module (and consensus up
# arm) adjusted for a unit-level microglial score, in lumbar (W67 inputs) and cervical (W68 inputs)
# spinal cord. Primary score: W66 microglia marker set; sensitivity: homeostatic-only and
# activated-only sets. Stouffer combination of one-sided fry P; attenuation versus W67/W68 primary.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(Matrix)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
res <- file.path(root, "results")
out <- file.path(res, "W70_microglia_adjusted_territory")
dir.create(out, showWarnings = FALSE)

SCORES <- list(
  microglia = c("P2RY12", "CX3CR1", "CSF1R", "PTPRC", "CD68", "TREM2", "TYROBP", "CD163", "ITGAX"),
  homeostatic = c("P2RY12", "CX3CR1", "CSF1R", "TMEM119"),
  activated = c("CD68", "TREM2", "TYROBP", "CD163", "ITGAX"))
seed <- file.path(res, "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
up <- read.delim(file.path(seed, "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_up.tsv"))$common_name
sig <- read.delim(file.path(res, "W59_NYGC_consensus_signature_holdout", "w59_frozen_consensus_signature.tsv"))
sets <- list(up = up, up_consensus = sig$symbol[sig$direction == "up"])

# ---- cohort inputs (same donors and territories as W67 / W68 primary), symbol row names ----
lumbar <- function() {
  inp <- file.path(res, "W67_motor_selectivity", "inputs")
  u <- read.delim(file.path(inp, "LSC_territory_units.tsv"), stringsAsFactors = FALSE)
  u$gse330130_donor <- u$gse330130_donor %in% c(TRUE, "True", "TRUE")
  u <- u[u$group %in% c("ALS", "Control") & !u$gse330130_donor & u$lib_size > 0 & !is.na(u$age) &
           u$territory %in% c("ventral_horn", "dorsal_horn"), ]
  u <- u[u$donor %in% names(which(table(u$donor) == 2)), ]
  x <- as.matrix(read.delim(gzfile(file.path(inp, "LSC_territory_counts.tsv.gz")), row.names = 1, check.names = FALSE))[, u$unit]
  sym <- read.delim(pipe(paste("gzip -dc", shQuote(file.path(res, "W66_spatial_context", "inputs", "LSC_pseudobulk_counts.tsv.gz")), "| cut -f1,2")))
  s <- setNames(sym[[2]], sym[[1]])[rownames(x)]
  keep <- !is.na(s) & !duplicated(s)
  x <- x[keep, ]; rownames(x) <- s[keep]
  meta <- data.frame(unit = u$unit, donor = u$donor, territory = ifelse(u$territory == "ventral_horn", "motor", "sensory"),
                     disease = ifelse(u$group == "ALS", "ALS", "Control"), sex = u$sex, age = u$age)
  list(name = "lumbar_GSE330224", counts = x, meta = meta, rhs = "+ sex + age")
}
cervical <- function() {
  inp <- file.path(res, "W68_GSE288365_replication", "inputs")
  x <- as.matrix(readMM(gzfile(file.path(inp, "w68_pseudobulk_counts.mtx.gz"))))
  rownames(x) <- readLines(gzfile(file.path(inp, "w68_genes.txt.gz"))); colnames(x) <- readLines(file.path(inp, "w68_unit_order.txt"))
  u <- read.delim(file.path(inp, "w68_units.tsv"), stringsAsFactors = FALSE)
  u <- u[u$unit_type == "territory" & u$stratum %in% c("anterior_horn", "other_grey") & !u$donor %in% c("NEUVM032BF2", "NEUWD570BTK") & u$lib_size > 0, ]
  meta <- data.frame(unit = u$unit, donor = u$donor, territory = ifelse(u$stratum == "anterior_horn", "motor", "sensory"),
                     disease = ifelse(u$diagnosis == "Control", "Control", "ALS"), sex = u$sex, age = u$age)
  list(name = "cervical_GSE288365", counts = x[, u$unit], meta = meta, rhs = "")
}

fit_adjusted <- function(co, score_name) {
  meta <- co$meta
  meta$disease <- factor(meta$disease, levels = c("Control", "ALS"))
  meta$territory <- factor(meta$territory, levels = c("sensory", "motor"))
  meta$sex <- factor(meta$sex)
  y0 <- calcNormFactors(DGEList(counts = co$counts), method = "TMM")
  lc <- cpm(y0, log = TRUE, prior.count = 1)
  mk <- intersect(SCORES[[score_name]], rownames(lc))
  meta$mg <- colMeans(lc[mk, , drop = FALSE]); meta$mg <- meta$mg - mean(meta$mg)
  design <- model.matrix(as.formula(paste("~ territory * disease", co$rhs, "+ mg")), data = meta)
  colnames(design) <- make.names(colnames(design))
  y <- DGEList(counts = co$counts)
  y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(y, design)
  cor1 <- duplicateCorrelation(v, design, block = meta$donor)$consensus.correlation
  v <- voom(y, design, block = meta$donor, correlation = cor1)
  cons <- duplicateCorrelation(v, design, block = meta$donor)$consensus.correlation
  fit <- eBayes(lmFit(v, design, block = meta$donor, correlation = cons), trend = TRUE)
  index <- lapply(sets, intersect, rownames(v))
  fr <- fry(v, index = index, design = design, contrast = "territorymotor.diseaseALS", block = meta$donor, correlation = cons, sort = "none")
  ca <- camera(v, index = index, design = design, contrast = "territorymotor.diseaseALS", sort = FALSE)
  tab <- topTable(fit, coef = "territorymotor.diseaseALS", number = Inf, sort.by = "none")
  sets_res <- do.call(rbind, lapply(names(index), function(m) data.frame(
    cohort = co$name, score = score_name, n_markers = length(mk), gene_set = m, n_genes_tested = length(index[[m]]),
    n_als = length(unique(meta$donor[meta$disease == "ALS"])), n_control = length(unique(meta$donor[meta$disease == "Control"])),
    fry_direction = fr[m, "Direction"], fry_p = fr[m, "PValue"], camera_p = ca[m, "PValue"],
    median_logFC_interaction = median(tab[index[[m]], "logFC"]))))
  desc <- do.call(rbind, lapply(c("motor", "sensory"), function(t) {
    a <- meta$mg[meta$territory == t & meta$disease == "ALS"]; c0 <- meta$mg[meta$territory == t & meta$disease == "Control"]
    data.frame(cohort = co$name, score = score_name, territory = t, mean_ALS_minus_control = mean(a) - mean(c0),
               wilcoxon_p = wilcox.test(a, c0)$p.value)
  }))
  list(sets = sets_res, desc = desc)
}

cohorts <- list(lumbar(), cervical())
runs <- unlist(lapply(cohorts, function(co) lapply(names(SCORES), function(sc) fit_adjusted(co, sc))), recursive = FALSE)
tab <- do.call(rbind, lapply(runs, `[[`, "sets"))
desc <- do.call(rbind, lapply(runs, `[[`, "desc"))

# unadjusted reference (W67 / W68 primary interaction)
w67 <- read.delim(file.path(res, "W67_motor_selectivity", "w67_motor_selectivity_tests.tsv"))
w67 <- w67[w67$comparison == "grey_ventral_vs_dorsal_horn" & w67$contrast == "interaction", ]
w68 <- read.delim(file.path(res, "W68_GSE288365_replication", "w68_replication_tests.tsv"))
w68 <- w68[w68$analysis == "primary" & w68$comparison == "grey" & w68$contrast == "interaction", ]
ref <- rbind(data.frame(cohort = "lumbar_GSE330224", gene_set = w67$gene_set, unadjusted_fry_p = w67$fry_p, unadjusted_median_logFC = w67$median_logFC),
             data.frame(cohort = "cervical_GSE288365", gene_set = w68$gene_set, unadjusted_fry_p = w68$fry_p, unadjusted_median_logFC = w68$median_logFC))
tab <- merge(tab, ref, by = c("cohort", "gene_set"), all.x = TRUE)
tab$attenuation <- 1 - tab$median_logFC_interaction / tab$unadjusted_median_logFC
tab$one_sided_p <- ifelse(tab$fry_direction == "Up", tab$fry_p / 2, 1 - tab$fry_p / 2)

stouffer <- function(p) pnorm(sum(qnorm(1 - p)) / sqrt(length(p)), lower.tail = FALSE)
comb <- do.call(rbind, lapply(split(tab, list(tab$score, tab$gene_set)), function(d) data.frame(
  score = d$score[1], gene_set = d$gene_set[1], both_up = all(d$fry_direction == "Up"),
  stouffer_p_one_sided = stouffer(d$one_sided_p), max_attenuation = max(d$attenuation))))
comb$reading <- with(comb, ifelse(stouffer_p_one_sided < 0.05 & both_up,
  ifelse(max_attenuation < 0.5, "not explained by microglial abundance", "partly shared with microglial abundance"),
  "largely shared with microglial abundance"))

write.table(tab, file.path(out, "w70_adjusted_interaction_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(comb, file.path(out, "w70_combined_reading.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(desc, file.path(out, "w70_microglial_score_by_territory.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(tab[, c("cohort", "score", "gene_set", "n_als", "n_control", "fry_direction", "fry_p", "camera_p", "median_logFC_interaction",
              "unadjusted_median_logFC", "attenuation")], digits = 3, row.names = FALSE)
print(comb, digits = 3, row.names = FALSE); print(desc, digits = 3, row.names = FALSE)
