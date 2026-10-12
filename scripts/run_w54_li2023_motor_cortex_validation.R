#!/usr/bin/env Rscript

# W54 independent motor-cortex validation of the final W34 astrocyte modules in
# Li et al. 2023 (GSE219280, BA4; C9-ALS, C9-FTD, controls). One pseudobulk per
# donor from marker-dominant astrocytes (see build_w54_li2023_astrocyte_pseudobulk.py).
#
# Plan fixed before results were inspected:
# - Primary: C9-ALS vs control; W33 covariate rule over age/sex/PMI (PMI is
#   missing for every control, so the rule drops it); filterByExpr(design),
#   TMM, voom, eBayes(trend); fry on the 15-gene up and 14-gene down modules
#   with BH across the two sets; camera as the competitive test.
# - Sensitivity: ~ disease + sex.
# - Context: C9-FTD vs control with the primary rule.
# - Gene level: directions, sign test, Spearman vs discovery logFC; genome-wide
#   FDR counts as a power check.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))

root <- p01_project_root()
out <- file.path(root, "results", "W54_Li2023_motor_cortex_validation")
seed_dir <- file.path(root, "results", "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")

counts <- as.matrix(read.delim(gzfile(file.path(out, "w54_li2023_astrocyte_pseudobulk_counts.tsv.gz")), row.names = 1, check.names = FALSE))
samples <- read.delim(file.path(out, "w54_li2023_astrocyte_pseudobulk_samples.tsv"), stringsAsFactors = FALSE)

read_seed <- function(direction) {
  seed <- read.delim(file.path(seed_dir, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", direction, ".tsv")), stringsAsFactors = FALSE)
  data.frame(gene_id = seed$gene_id, gene = seed$common_name, discovery_logFC = seed$logFC, module = direction, stringsAsFactors = FALSE)
}
modules <- rbind(read_seed("up"), read_seed("down"))

# W33 rule: drop a missing covariate (order PMI, sex, age) only if a group
# falls below 5 complete-case donors.
build_design <- function(s, candidates) {
  used <- candidates[vapply(candidates, function(cv) length(unique(na.omit(s[[cv]]))) >= 2, logical(1))]
  repeat {
    m <- s[complete.cases(s[, c("disease", used), drop = FALSE]), , drop = FALSE]
    if (all(table(m$disease) >= 5)) break
    removable <- intersect(c("PMI", "sex", "age"), used[vapply(used, function(cv) anyNA(s[[cv]]), logical(1))])
    if (length(removable) == 0) break
    used <- setdiff(used, removable[1])
  }
  f <- paste("~ disease", if (length(used)) paste("+", paste(used, collapse = " + ")) else "")
  list(samples = m, design = model.matrix(as.formula(f), data = m), formula = f)
}

run_analysis <- function(analysis_id, case, candidates) {
  s <- samples[samples$disease %in% c(case, "Control"), ]
  s$disease <- factor(ifelse(s$disease == case, "Case", "Control"), levels = c("Control", "Case"))
  s$sex <- factor(s$sex)
  info <- build_design(s, candidates)
  y <- DGEList(counts = counts[, info$samples$GSM, drop = FALSE])
  keep <- filterByExpr(y, design = info$design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(y, info$design)
  fit <- eBayes(lmFit(v, info$design), trend = TRUE)
  tab <- topTable(fit, coef = "diseaseCase", number = Inf, sort.by = "none")
  tab$gene_id <- rownames(tab)
  index <- list(up = intersect(modules$gene_id[modules$module == "up"], rownames(v)),
                down = intersect(modules$gene_id[modules$module == "down"], rownames(v)))
  fr <- fry(v, index = index, design = info$design, contrast = "diseaseCase", sort = "none")
  ca <- camera(v, index = index, design = info$design, contrast = "diseaseCase", sort = FALSE)
  sets <- do.call(rbind, lapply(names(index), function(set) {
    g <- merge(modules[modules$module == set, ], tab[, c("gene_id", "logFC")], by = "gene_id")
    k <- sum(sign(g$logFC) == sign(g$discovery_logFC))
    data.frame(analysis_id = analysis_id, contrast = paste(case, "vs Control"), module = set, formula = info$formula,
               n_case = sum(info$samples$disease == "Case"), n_control = sum(info$samples$disease == "Control"),
               n_genes_tested_total = nrow(v), n_genome_wide_FDR05 = sum(tab$adj.P.Val < 0.05),
               n_module_genes_tested = length(index[[set]]),
               fry_direction = fr[set, "Direction"], fry_p = fr[set, "PValue"],
               camera_direction = ca[set, "Direction"], camera_p = ca[set, "PValue"],
               n_genes_concordant = k, sign_test_p_one_sided = binom.test(k, nrow(g), 0.5, alternative = "greater")$p.value,
               spearman_logFC_vs_discovery = suppressWarnings(cor(g$logFC, g$discovery_logFC, method = "spearman")),
               stringsAsFactors = FALSE)
  }))
  sets$fry_BH <- p.adjust(sets$fry_p, method = "BH")
  sets$camera_BH <- p.adjust(sets$camera_p, method = "BH")
  genes <- merge(modules, tab[, c("gene_id", "logFC", "AveExpr", "P.Value", "adj.P.Val")], by = "gene_id", all.x = TRUE)
  genes$analysis_id <- analysis_id
  genes$tested_after_filterByExpr <- !is.na(genes$logFC)
  write.table(tab, gzfile(file.path(out, paste0("w54_", analysis_id, "_genome_wide.tsv.gz"))), sep = "\t", quote = FALSE, row.names = FALSE)
  list(sets = sets, genes = genes)
}

analyses <- list(
  run_analysis("primary_C9ALS_rule", "C9-ALS", c("age", "sex", "PMI")),
  run_analysis("sens_C9ALS_sex_only", "C9-ALS", c("sex")),
  run_analysis("context_C9FTD_rule", "C9-FTD", c("age", "sex", "PMI"))
)
sets <- do.call(rbind, lapply(analyses, `[[`, "sets"))
genes <- do.call(rbind, lapply(analyses, `[[`, "genes"))
write.table(sets, file.path(out, "w54_module_set_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes[order(genes$analysis_id, genes$module, genes$P.Value), ], file.path(out, "w54_module_gene_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(sets[, c("analysis_id", "module", "formula", "n_case", "n_control", "n_genome_wide_FDR05", "n_module_genes_tested", "fry_direction", "fry_p", "fry_BH", "camera_p", "n_genes_concordant", "sign_test_p_one_sided")], digits = 3)
