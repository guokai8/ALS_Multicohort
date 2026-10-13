#!/usr/bin/env Rscript

# W69 (Amendment 12).
# A: within ALS donors, Spearman correlation between a motor-neuron index (motor-territory mean
#    log2-CPM of 10 motor-neuron genes) and the frozen ALS-up module response (motor minus sensory
#    territory mean log2-CPM), in lumbar (W67 inputs) and cervical (W68 inputs) cord; one-sided
#    (rho < 0), Stouffer-combined; specificity against 1,000 expression-matched random gene sets.
# B: Stouffer combination of the W67 and W68 primary grey-matter interaction tests (post hoc summary).

suppressPackageStartupMessages({
  library(edgeR)
  library(Matrix)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
res <- file.path(root, "results")
out <- file.path(res, "W69_dose_response_and_combination")
dir.create(out, showWarnings = FALSE)
set.seed(20261013)
N_RANDOM <- 1000
MN_GENES <- c("CHAT", "SLC5A7", "SLC18A3", "MNX1", "ISL1", "PRPH", "NEFH", "NEFM", "NEFL", "STMN2")

seed <- file.path(res, "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
up <- read.delim(file.path(seed, "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_up.tsv"))
sig <- read.delim(file.path(res, "W59_NYGC_consensus_signature_holdout", "w59_frozen_consensus_signature.tsv"))

# ---- cohort inputs, returned with symbol row names ----
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
  list(name = "lumbar_GSE330224", counts = x, donor = u$donor, motor = u$territory == "ventral_horn", als = u$group == "ALS")
}
cervical <- function() {
  inp <- file.path(res, "W68_GSE288365_replication", "inputs")
  x <- as.matrix(readMM(gzfile(file.path(inp, "w68_pseudobulk_counts.mtx.gz"))))
  rownames(x) <- readLines(gzfile(file.path(inp, "w68_genes.txt.gz"))); colnames(x) <- readLines(file.path(inp, "w68_unit_order.txt"))
  u <- read.delim(file.path(inp, "w68_units.tsv"), stringsAsFactors = FALSE)
  u <- u[u$unit_type == "territory" & u$stratum %in% c("anterior_horn", "other_grey") & !u$donor %in% c("NEUVM032BF2", "NEUWD570BTK") & u$lib_size > 0, ]
  list(name = "cervical_GSE288365", counts = x[, u$unit], donor = u$donor, motor = u$stratum == "anterior_horn", als = u$diagnosis != "Control")
}

dose_response <- function(co, set_genes, set_name) {
  y <- DGEList(counts = co$counts)
  y <- calcNormFactors(y[filterByExpr(y, group = co$motor), , keep.lib.sizes = FALSE], method = "TMM")
  lc <- cpm(y, log = TRUE, prior.count = 1)
  mn <- intersect(MN_GENES, rownames(lc)); gs <- intersect(set_genes, rownames(lc))
  donors <- unique(co$donor)
  mi <- function(d) which(co$donor == d & co$motor); si <- function(d) which(co$donor == d & !co$motor)
  mn_index <- sapply(donors, function(d) mean(lc[mn, mi(d)]))
  resp <- function(g) sapply(donors, function(d) mean(lc[g, mi(d)]) - mean(lc[g, si(d)]))
  is_als <- sapply(donors, function(d) co$als[co$donor == d][1])
  r <- resp(gs)
  ct <- suppressWarnings(cor.test(mn_index[is_als], r[is_als], method = "spearman", alternative = "less", exact = TRUE))
  # matched random sets: one gene from the motor-territory expression decile of each set gene
  motor_mean <- rowMeans(lc[, co$motor, drop = FALSE])
  pool <- setdiff(rownames(lc), c(gs, mn, up$common_name))
  dec <- cut(motor_mean, quantile(motor_mean, 0:10 / 10), include.lowest = TRUE, labels = FALSE); names(dec) <- rownames(lc)
  rnd <- replicate(N_RANDOM, {
    g <- vapply(gs, function(x) sample(pool[dec[pool] == dec[x]], 1), "")
    suppressWarnings(cor(mn_index[is_als], resp(g)[is_als], method = "spearman"))
  })
  list(summary = data.frame(cohort = co$name, gene_set = set_name, n_als = sum(is_als), n_control = sum(!is_als),
                            n_genes = length(gs), n_mn_genes = length(mn), rho = unname(ct$estimate), p_one_sided = ct$p.value,
                            random_median_rho = median(rnd), percentile_vs_random = mean(rnd < ct$estimate),
                            mn_index_wilcox_p = wilcox.test(mn_index[is_als], mn_index[!is_als])$p.value,
                            mn_index_median_als = median(mn_index[is_als]), mn_index_median_control = median(mn_index[!is_als])),
       donors = data.frame(cohort = co$name, gene_set = set_name, donor = donors, als = is_als, mn_index = mn_index, module_response = r))
}

cohorts <- list(lumbar = lumbar(), cervical = cervical())
sets <- list(up = up$common_name, up_consensus = sig$symbol[sig$direction == "up"])
a <- unlist(lapply(cohorts, function(co) lapply(names(sets), function(m) dose_response(co, sets[[m]], m))), recursive = FALSE)
summ <- do.call(rbind, lapply(a, `[[`, "summary"))
donors <- do.call(rbind, lapply(a, `[[`, "donors"))
stouffer <- function(p) { z <- qnorm(1 - p); pnorm(sum(z) / sqrt(length(z)), lower.tail = FALSE) }
comb_a <- do.call(rbind, lapply(split(summ, summ$gene_set), function(d) data.frame(
  gene_set = d$gene_set[1], both_rho_negative = all(d$rho < 0), stouffer_p_one_sided = stouffer(d$p_one_sided),
  below_5th_percentile_both = all(d$percentile_vs_random < 0.05))))
comb_a$reading <- ifelse(comb_a$stouffer_p_one_sided < 0.05 & comb_a$both_rho_negative,
                         ifelse(comb_a$below_5th_percentile_both, "dose-response, specific", "dose-response, not distinguishable from compositional coupling"),
                         "not supported")

# ---- B: cross-cohort combination of the primary grey-matter interaction tests ----
w67 <- read.delim(file.path(res, "W67_motor_selectivity", "w67_motor_selectivity_tests.tsv"))
w67 <- w67[w67$comparison == "grey_ventral_vs_dorsal_horn" & w67$contrast == "interaction", ]
w68 <- read.delim(file.path(res, "W68_GSE288365_replication", "w68_replication_tests.tsv"))
w68 <- w68[w68$analysis == "primary" & w68$comparison == "grey" & w68$contrast == "interaction", ]
one_sided <- function(p, dir) ifelse(dir == "Up", p / 2, 1 - p / 2)
comb_b <- do.call(rbind, lapply(c("up", "up_consensus"), function(m) {
  p <- c(lumbar = one_sided(w67$fry_p[w67$gene_set == m], w67$fry_direction[w67$gene_set == m]),
         cervical = one_sided(w68$fry_p[w68$gene_set == m], w68$fry_direction[w68$gene_set == m]))
  data.frame(gene_set = m, lumbar_one_sided_p = p[["lumbar"]], cervical_one_sided_p = p[["cervical"]],
             stouffer_z = sum(qnorm(1 - p)) / sqrt(2), stouffer_p_one_sided = stouffer(p),
             role = if (m == "up") "primary" else "secondary")
}))

write.table(summ, file.path(out, "w69a_dose_response_by_cohort.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(comb_a, file.path(out, "w69a_dose_response_combined.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(donors, file.path(out, "w69a_donor_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(comb_b, file.path(out, "w69b_interaction_combination.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, digits = 3, row.names = FALSE); print(comb_a, digits = 3, row.names = FALSE); print(comb_b, digits = 3, row.names = FALSE)
