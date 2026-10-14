#!/usr/bin/env Rscript

# W73 (Amendment 15): does the ventral-horn response reduce to generic inflammation, microglial
# (TYROBP-network) activation or ageing? A: territory x disease interaction of six generic msigdbr
# sets (module and consensus-up genes removed) in lumbar (W67 inputs) and cervical (W68 inputs) cord.
# B: the frozen ALS-up module (and consensus up arm) refitted with each generic unit score as a
# covariate, with attenuation against the W67/W68 primary interaction and a Stouffer combination.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(Matrix)
  library(msigdbr)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
res <- file.path(root, "results")
out <- file.path(res, "W73_generic_signature_comparison")
dir.create(out, showWarnings = FALSE)

seed <- file.path(res, "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
up <- read.delim(file.path(seed, "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_up.tsv"))$common_name
sig <- read.delim(file.path(res, "W59_NYGC_consensus_signature_holdout", "w59_frozen_consensus_signature.tsv"))
modules <- list(up = up, up_consensus = sig$symbol[sig$direction == "up"])
GENERIC <- c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
             "HALLMARK_COMPLEMENT", "WP_TYROBP_CAUSAL_NETWORK_IN_MICROGLIA", "LU_AGING_BRAIN_UP")
ms <- rbind(msigdbr(species = "Homo sapiens", collection = "H")[, c("gs_name", "gene_symbol")],
            msigdbr(species = "Homo sapiens", collection = "C2")[, c("gs_name", "gene_symbol")])
generic <- lapply(setNames(GENERIC, GENERIC), function(s) setdiff(unique(ms$gene_symbol[ms$gs_name == s]), unlist(modules)))

# ---- cohort inputs: identical to W67 / W68 primary (as in W70) ----
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

fit <- function(co, index, covariate = NULL) {
  meta <- co$meta
  meta$disease <- factor(meta$disease, levels = c("Control", "ALS"))
  meta$territory <- factor(meta$territory, levels = c("sensory", "motor"))
  meta$sex <- factor(meta$sex)
  rhs <- co$rhs
  n_cov <- NA
  if (!is.null(covariate)) {
    lc <- cpm(calcNormFactors(DGEList(counts = co$counts), method = "TMM"), log = TRUE, prior.count = 1)
    g <- intersect(generic[[covariate]], rownames(lc)); n_cov <- length(g)
    meta$gs <- colMeans(lc[g, , drop = FALSE]); meta$gs <- meta$gs - mean(meta$gs)
    rhs <- paste(rhs, "+ gs")
  }
  design <- model.matrix(as.formula(paste("~ territory * disease", rhs)), data = meta)
  colnames(design) <- make.names(colnames(design))
  y <- DGEList(counts = co$counts)
  y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(y, design)
  cor1 <- duplicateCorrelation(v, design, block = meta$donor)$consensus.correlation
  v <- voom(y, design, block = meta$donor, correlation = cor1)
  cons <- duplicateCorrelation(v, design, block = meta$donor)$consensus.correlation
  ft <- eBayes(lmFit(v, design, block = meta$donor, correlation = cons), trend = TRUE)
  idx <- lapply(index, intersect, rownames(v))
  fr <- fry(v, index = idx, design = design, contrast = "territorymotor.diseaseALS", block = meta$donor, correlation = cons, sort = "none")
  ca <- camera(v, index = idx, design = design, contrast = "territorymotor.diseaseALS", sort = FALSE)
  tab <- topTable(ft, coef = "territorymotor.diseaseALS", number = Inf, sort.by = "none")
  do.call(rbind, lapply(names(idx), function(m) data.frame(
    cohort = co$name, gene_set = m, adjusted_for = ifelse(is.null(covariate), "none", covariate), n_cov_genes = n_cov,
    n_genes_tested = length(idx[[m]]), n_als = length(unique(meta$donor[meta$disease == "ALS"])),
    n_control = length(unique(meta$donor[meta$disease == "Control"])),
    fry_direction = fr[m, "Direction"], fry_p = fr[m, "PValue"], camera_p = ca[m, "PValue"],
    median_logFC_interaction = median(tab[idx[[m]], "logFC"]))))
}
one_sided <- function(d) ifelse(d$fry_direction == "Up", d$fry_p / 2, 1 - d$fry_p / 2)
stouffer <- function(p) pnorm(sum(qnorm(1 - p)) / sqrt(length(p)), lower.tail = FALSE)

cohorts <- list(lumbar(), cervical())
# ---- A: generic sets in the territory test ----
A <- do.call(rbind, lapply(cohorts, function(co) fit(co, generic)))
A$fry_BH_within_cohort <- ave(A$fry_p, A$cohort, FUN = function(p) p.adjust(p, "BH"))
A$one_sided_p <- one_sided(A)
A_comb <- do.call(rbind, lapply(split(A, A$gene_set), function(d) data.frame(
  gene_set = d$gene_set[1], both_up = all(d$fry_direction == "Up"), stouffer_p_one_sided = stouffer(d$one_sided_p))))
A_comb$reading <- ifelse(A_comb$stouffer_p_one_sided < 0.05 & A_comb$both_up, "shares the ventral-horn pattern", "does not share it")

# ---- B: module adjusted for each generic score ----
B <- do.call(rbind, lapply(cohorts, function(co) do.call(rbind, lapply(GENERIC, function(g) fit(co, modules, covariate = g)))))
w67 <- read.delim(file.path(res, "W67_motor_selectivity", "w67_motor_selectivity_tests.tsv"))
w67 <- w67[w67$comparison == "grey_ventral_vs_dorsal_horn" & w67$contrast == "interaction", ]
w68 <- read.delim(file.path(res, "W68_GSE288365_replication", "w68_replication_tests.tsv"))
w68 <- w68[w68$analysis == "primary" & w68$comparison == "grey" & w68$contrast == "interaction", ]
ref <- rbind(data.frame(cohort = "lumbar_GSE330224", gene_set = w67$gene_set, unadjusted_median_logFC = w67$median_logFC),
             data.frame(cohort = "cervical_GSE288365", gene_set = w68$gene_set, unadjusted_median_logFC = w68$median_logFC))
B <- merge(B, ref, by = c("cohort", "gene_set"), all.x = TRUE)
B$attenuation <- 1 - B$median_logFC_interaction / B$unadjusted_median_logFC
B$one_sided_p <- one_sided(B)
B_comb <- do.call(rbind, lapply(split(B, list(B$gene_set, B$adjusted_for)), function(d) data.frame(
  gene_set = d$gene_set[1], adjusted_for = d$adjusted_for[1], both_up = all(d$fry_direction == "Up"),
  stouffer_p_one_sided = stouffer(d$one_sided_p), max_attenuation = max(d$attenuation))))
B_comb$reading <- with(B_comb, ifelse(stouffer_p_one_sided < 0.05 & both_up,
  ifelse(max_attenuation < 0.5, "beyond this programme", "partly shared"), "largely shared"))
core <- B_comb[B_comb$gene_set == "up" & B_comb$adjusted_for != "WP_TYROBP_CAUSAL_NETWORK_IN_MICROGLIA", ]
overall <- data.frame(question = "frozen module not reducible to generic inflammation or ageing",
                      reading = ifelse(all(core$reading == "beyond this programme"), "yes", "no"),
                      basis = paste(core$adjusted_for, core$reading, sep = ": ", collapse = "; "))

write.table(A, file.path(out, "w73_generic_set_territory_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(A_comb, file.path(out, "w73_generic_set_combined.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(B, file.path(out, "w73_module_adjusted_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(B_comb, file.path(out, "w73_module_adjusted_combined.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(overall, file.path(out, "w73_reading.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene_set = names(generic), n_genes_after_removing_module = lengths(generic)), file.path(out, "w73_gene_sets.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(A[, c("cohort", "gene_set", "n_genes_tested", "fry_direction", "fry_p", "fry_BH_within_cohort", "median_logFC_interaction")], digits = 3, row.names = FALSE)
print(A_comb, digits = 3, row.names = FALSE)
print(B[, c("cohort", "gene_set", "adjusted_for", "fry_direction", "fry_p", "median_logFC_interaction", "attenuation")], digits = 3, row.names = FALSE)
print(B_comb, digits = 3, row.names = FALSE)
print(overall, row.names = FALSE)
