#!/usr/bin/env Rscript

# W52 independent motor-cortex validation of the final W34 MTC astrocyte
# modules in the Ulm/DZNE multiome cohort (Zenodo 18370646; Ruf et al. 2026).
# Pipeline as in discovery: filterByExpr(design), TMM, voom, eBayes(trend).
# Age, PMI and RIN cannot be linked to the released metadata IDs, so they are
# not modelled; chip is modelled as the processing block.
#
# Plan fixed before results were inspected:
# - Primary: ALS vs HC on donor x chip pseudobulk, ~ disease + sex + chip,
#   duplicateCorrelation(block = donor); fry with block/correlation on the
#   15-gene up and 14-gene down modules, BH across the two sets.
# - Sensitivity A: one profile per donor, ~ disease + sex + main_chip; fry, camera.
# - Sensitivity B: primary + astrocyte QC proxies (mean percent_mito,
#   log10 median nCount_RNA).
# - Sensitivity C: one profile per donor, ~ disease + sex; fry, camera.
# - Context: non-C9 ALS-FTD vs HC and C9 ALS-FTD vs HC with the primary design.
# - Gene level: primary-model directions, sign test, Spearman vs discovery logFC.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))

root <- p01_project_root()
out <- file.path(root, "results", "W52_Ulm_motor_cortex_validation")
seed_dir <- file.path(root, "results", "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")

counts_strata <- as.matrix(read.delim(gzfile(file.path(out, "w52_ulm_astrocyte_pseudobulk_counts.tsv.gz")), row.names = 1, check.names = FALSE))
strata <- read.delim(file.path(out, "w52_ulm_astrocyte_pseudobulk_samples.tsv"), stringsAsFactors = FALSE)
counts_strata <- counts_strata[, strata$stratum]

read_seed <- function(direction) {
  seed <- read.delim(file.path(seed_dir, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", direction, ".tsv")), stringsAsFactors = FALSE)
  data.frame(gene = seed$common_name, discovery_logFC = seed$logFC, module = direction, stringsAsFactors = FALSE)
}
modules <- rbind(read_seed("up"), read_seed("down"))

# One profile per donor: sum strata; main chip = chip with most astrocytes.
donor_ids <- unique(strata$ID)
counts_donor <- sapply(donor_ids, function(id) rowSums(counts_strata[, strata$ID == id, drop = FALSE]))
donors <- do.call(rbind, lapply(donor_ids, function(id) {
  s <- strata[strata$ID == id, ]
  data.frame(ID = id, Sex = s$Sex[1], Case = s$Case[1], Case_Type = s$Case_Type[1], C9ORF72 = s$C9ORF72[1],
             main_chip = s$Chip[which.max(s$n_cells)], n_cells = sum(s$n_cells), stringsAsFactors = FALSE)
}))

select_contrast <- function(samples, case_type) {
  keep <- samples$Case_Type %in% c(case_type, "HC")
  s <- samples[keep, , drop = FALSE]
  s$disease <- factor(ifelse(s$Case_Type == case_type, "ALS", "Control"), levels = c("Control", "ALS"))
  s$sex <- factor(tolower(s$Sex))
  s
}

run_analysis <- function(analysis_id, samples, counts, rhs, block = FALSE, case_type = "ALS") {
  s <- select_contrast(samples, case_type)
  if ("Chip" %in% names(s)) s$chip <- factor(s$Chip)
  if ("main_chip" %in% names(s)) s$main_chip <- factor(s$main_chip)
  if ("median_nCount_RNA" %in% names(s)) s$log10_nCount <- log10(s$median_nCount_RNA)
  design <- model.matrix(as.formula(paste("~ disease +", rhs)), data = s)
  if (qr(design)$rank < ncol(design)) stop("Rank-deficient design in ", analysis_id)
  sample_ids <- if (block) s$stratum else s$ID
  y <- DGEList(counts = counts[, sample_ids, drop = FALSE])
  keep <- filterByExpr(y, design = design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")

  consensus <- NA_real_
  if (block) {
    v <- voom(y, design)
    consensus <- duplicateCorrelation(v, design, block = s$ID)$consensus.correlation
    v <- voom(y, design, block = s$ID, correlation = consensus)
    consensus <- duplicateCorrelation(v, design, block = s$ID)$consensus.correlation
    fit <- lmFit(v, design, block = s$ID, correlation = consensus)
  } else {
    v <- voom(y, design)
    fit <- lmFit(v, design)
  }
  fit <- eBayes(fit, trend = TRUE)
  tab <- topTable(fit, coef = "diseaseALS", number = Inf, sort.by = "none")
  tab$gene <- rownames(tab)

  index <- list(up = intersect(modules$gene[modules$module == "up"], rownames(v)),
                down = intersect(modules$gene[modules$module == "down"], rownames(v)))
  fry_res <- if (block) {
    fry(v, index = index, design = design, contrast = "diseaseALS", block = s$ID, correlation = consensus, sort = "none")
  } else {
    fry(v, index = index, design = design, contrast = "diseaseALS", sort = "none")
  }
  camera_res <- if (block) NULL else camera(v, index = index, design = design, contrast = "diseaseALS", sort = FALSE)

  sets <- do.call(rbind, lapply(names(index), function(set) {
    g <- merge(modules[modules$module == set, ], tab[, c("gene", "logFC")], by = "gene")
    concordant <- sum(sign(g$logFC) == sign(g$discovery_logFC))
    data.frame(
      analysis_id = analysis_id, contrast = paste(case_type, "vs HC"), module = set, formula = paste("~ disease +", rhs),
      unit = if (block) "donor x chip (duplicateCorrelation)" else "donor",
      n_case_donors = length(unique(s$ID[s$disease == "ALS"])), n_control_donors = length(unique(s$ID[s$disease == "Control"])),
      n_profiles = ncol(y), consensus_correlation = consensus, n_genes_tested_total = nrow(v),
      n_module_genes_tested = length(index[[set]]),
      fry_direction = fry_res[set, "Direction"], fry_p = fry_res[set, "PValue"],
      camera_direction = if (is.null(camera_res)) NA else camera_res[set, "Direction"],
      camera_p = if (is.null(camera_res)) NA else camera_res[set, "PValue"],
      n_genes_concordant = concordant,
      sign_test_p_one_sided = binom.test(concordant, nrow(g), 0.5, alternative = "greater")$p.value,
      spearman_logFC_vs_discovery = suppressWarnings(cor(g$logFC, g$discovery_logFC, method = "spearman")),
      stringsAsFactors = FALSE
    )
  }))
  sets$fry_BH <- p.adjust(sets$fry_p, method = "BH")
  sets$camera_BH <- if (all(is.na(sets$camera_p))) NA else p.adjust(sets$camera_p, method = "BH")

  genes <- merge(modules, tab[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], by = "gene", all.x = TRUE)
  genes$analysis_id <- analysis_id
  genes$tested_after_filterByExpr <- !is.na(genes$logFC)
  genes$direction_concordant <- sign(genes$logFC) == sign(genes$discovery_logFC)
  genes$BH_within_29 <- NA_real_
  genes$BH_within_29[genes$tested_after_filterByExpr] <- p.adjust(genes$P.Value[genes$tested_after_filterByExpr], method = "BH")
  names(genes)[names(genes) == "adj.P.Val"] <- "BH_genome_wide"
  write.table(tab, gzfile(file.path(out, paste0("w52_", analysis_id, "_genome_wide.tsv.gz"))), sep = "\t", quote = FALSE, row.names = FALSE)
  list(sets = sets, genes = genes)
}

analyses <- list(
  run_analysis("primary_donor_chip_blocked", strata, counts_strata, "sex + chip", block = TRUE),
  run_analysis("sensA_donor_main_chip", donors, counts_donor, "sex + main_chip"),
  run_analysis("sensB_primary_plus_QC", strata, counts_strata, "sex + chip + mean_percent_mito + log10_nCount", block = TRUE),
  run_analysis("sensC_donor_no_chip", donors, counts_donor, "sex"),
  run_analysis("context_ALS_FTD_nonC9", strata, counts_strata, "sex + chip", block = TRUE, case_type = "ALS_FTD"),
  run_analysis("context_C9_ALS_FTD", strata, counts_strata, "sex + chip", block = TRUE, case_type = "C9_ALS_FTD")
)
sets <- do.call(rbind, lapply(analyses, `[[`, "sets"))
genes <- do.call(rbind, lapply(analyses, `[[`, "genes"))
write.table(sets, file.path(out, "w52_module_set_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes[order(genes$analysis_id, genes$module, genes$P.Value), ], file.path(out, "w52_module_gene_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(sets[, c("analysis_id", "module", "n_case_donors", "n_control_donors", "n_module_genes_tested", "fry_direction", "fry_p", "fry_BH", "camera_p", "n_genes_concordant", "sign_test_p_one_sided")], digits = 3)
