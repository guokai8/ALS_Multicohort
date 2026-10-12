#!/usr/bin/env Rscript

# W51 discovery-matched replication of the final W34 MTC astrocyte modules in
# Pineda BA4 astrocytes. Reuses the W33 pipeline: covariate candidates
# age/sex/PMI (tissue_subregion does not apply), complete-case covariates are
# dropped only if a group falls below 5 donors, then filterByExpr(design), TMM,
# voom and eBayes(trend = TRUE).
#
# Plan fixed before results were inspected:
# - Primary: discovery-matched design, ALS vs PN; module tests with limma::fry
#   (self-contained, directional) on the 15-gene up and 14-gene down sets, BH
#   across the two sets.
# - Sensitivities: drop PMI (all donors); restrict to samples with >=20
#   astrocytes; camera competitive test; gene-level directions with a sign test.
# W51 was added after the W50 cell-level module score was null in BA4.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "p01_paths.R"))

root <- p01_project_root()
out <- file.path(root, "results", "W51_Pineda_BA4_discovery_matched_limma")
seed_dir <- file.path(root, "results", "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")

counts_all <- read.delim(gzfile(file.path(out, "w51_pineda_astrocyte_pseudobulk_counts.tsv.gz")), row.names = 1, check.names = FALSE)
samples_all <- read.delim(file.path(out, "w51_pineda_astrocyte_pseudobulk_samples.tsv"), stringsAsFactors = FALSE)

read_seed <- function(direction) {
  seed <- read.delim(file.path(seed_dir, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", direction, ".tsv")), stringsAsFactors = FALSE)
  data.frame(gene_id = seed$gene_id, gene = seed$common_name, discovery_logFC = seed$logFC, module = direction, stringsAsFactors = FALSE)
}
modules <- rbind(read_seed("up"), read_seed("down"))

# W33 rule: keep candidates with >=2 levels; drop missing covariates in order
# PMI, sex, age only if a group falls below 5 complete-case donors.
build_design <- function(samples, candidates) {
  samples$disease <- factor(ifelse(samples$Condition == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
  samples$sex <- factor(samples$Sex)
  used <- candidates[vapply(candidates, function(cv) length(unique(na.omit(samples[[cv]]))) >= 2, logical(1))]
  repeat {
    model_samples <- samples[complete.cases(samples[, c("disease", used), drop = FALSE]), , drop = FALSE]
    n_group <- table(model_samples$disease)
    if (all(n_group >= 5)) break
    removable <- intersect(c("PMI", "sex", "age"), used[vapply(used, function(cv) anyNA(samples[[cv]]), logical(1))])
    if (length(removable) == 0) break
    used <- setdiff(used, removable[1])
  }
  formula_text <- paste("~ disease", if (length(used)) paste("+", paste(used, collapse = " + ")) else "")
  list(samples = model_samples, design = model.matrix(as.formula(formula_text), data = model_samples), formula_text = formula_text)
}

run_analysis <- function(analysis_id, samples, candidates) {
  info <- build_design(samples, candidates)
  y <- DGEList(counts = as.matrix(counts_all[, info$samples$sample_folder, drop = FALSE]))
  keep <- filterByExpr(y, design = info$design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(y, info$design, plot = FALSE)
  fit <- eBayes(lmFit(v, info$design), trend = TRUE)
  table_all <- topTable(fit, coef = "diseaseALS", number = Inf, sort.by = "none")
  table_all$gene_id <- rownames(table_all)

  index <- list(up = intersect(modules$gene_id[modules$module == "up"], rownames(v)),
                down = intersect(modules$gene_id[modules$module == "down"], rownames(v)))
  fry_res <- fry(v, index = index, design = info$design, contrast = "diseaseALS", sort = "none")
  camera_res <- camera(v, index = index, design = info$design, contrast = "diseaseALS", sort = FALSE)
  set_rows <- lapply(names(index), function(set) {
    expected <- if (set == "up") "Up" else "Down"
    data.frame(
      analysis_id = analysis_id, module = set, formula = info$formula_text,
      n_als = sum(info$samples$disease == "ALS"), n_control = sum(info$samples$disease == "Control"),
      n_genes_tested_total = nrow(v), n_module_genes_seed = sum(modules$module == set), n_module_genes_tested = length(index[[set]]),
      fry_direction = fry_res[set, "Direction"], fry_p_directional = fry_res[set, "PValue"],
      fry_direction_matches_discovery = fry_res[set, "Direction"] == expected,
      camera_direction = camera_res[set, "Direction"], camera_p = camera_res[set, "PValue"],
      stringsAsFactors = FALSE
    )
  })
  sets <- do.call(rbind, set_rows)
  sets$fry_BH <- p.adjust(sets$fry_p_directional, method = "BH")
  sets$camera_BH <- p.adjust(sets$camera_p, method = "BH")

  genes <- merge(modules, table_all[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], by = "gene_id", all.x = TRUE)
  genes$analysis_id <- analysis_id
  genes$tested_after_filterByExpr <- !is.na(genes$logFC)
  genes$direction_concordant <- sign(genes$logFC) == sign(genes$discovery_logFC)
  genes$BH_within_29 <- NA_real_
  genes$BH_within_29[genes$tested_after_filterByExpr] <- p.adjust(genes$P.Value[genes$tested_after_filterByExpr], method = "BH")
  names(genes)[names(genes) == "adj.P.Val"] <- "BH_genome_wide"

  for (set in c("up", "down")) {
    g <- genes[genes$module == set & genes$tested_after_filterByExpr, ]
    k <- sum(g$direction_concordant)
    sets[sets$module == set, "n_genes_direction_concordant"] <- k
    sets[sets$module == set, "sign_test_p_one_sided"] <- binom.test(k, nrow(g), 0.5, alternative = "greater")$p.value
    sets[sets$module == set, "spearman_logFC_vs_discovery"] <- suppressWarnings(cor(g$logFC, g$discovery_logFC, method = "spearman"))
  }
  write.table(table_all, gzfile(file.path(out, paste0("w51_", analysis_id, "_genome_wide.tsv.gz"))), sep = "\t", quote = FALSE, row.names = FALSE)
  list(sets = sets, genes = genes)
}

ba4_samples <- function() {
  samples_all[samples_all$Region == "BA4" & samples_all$Condition %in% c("ALS", "PN") & samples_all$lib_size > 0, ]
}

write_results <- function(analyses, prefix) {
  sets <- do.call(rbind, lapply(analyses, `[[`, "sets"))
  genes <- do.call(rbind, lapply(analyses, `[[`, "genes"))
  write.table(sets, file.path(out, paste0(prefix, "_module_set_tests.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(genes[order(genes$analysis_id, genes$module, genes$P.Value), ], file.path(out, paste0(prefix, "_module_gene_results.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  print(sets, digits = 3)
}

main <- function() {
  ba4 <- ba4_samples()
  write_results(list(
    primary_discovery_matched = run_analysis("primary_discovery_matched", ba4, c("age", "sex", "PMI")),
    sensitivity_no_PMI = run_analysis("sensitivity_no_PMI", ba4, c("age", "sex")),
    sensitivity_min20_cells = run_analysis("sensitivity_min20_cells", ba4[ba4$n_cells >= 20, ], c("age", "sex", "PMI"))
  ), "w51")
}

# Run only when executed directly, so post hoc scripts can source the helpers.
if (sys.nframe() == 0) main()
