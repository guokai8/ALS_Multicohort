#!/usr/bin/env Rscript

# W33 retains the frozen W01/W03 methods but uses Synapse/NYGC diagnosis
# resolution for exact matched libraries, leaving unmatched GEO labels unchanged.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "p01_paths.R"))

root <- p01_project_root()
out <- file.path(root, "results", "W33_synapse_resolved_diagnosis_sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)

bridge_path <- file.path(
  root,
  "results",
  "Supplemental_QC_covariate_audit",
  "gse330130_nygc_latest_zenodo_bridge_2026-08-12.tsv"
)
if (!file.exists(bridge_path)) stop("Missing NYGC bridge: ", bridge_path)
bridge <- read.delim(bridge_path, stringsAsFactors = FALSE, check.names = FALSE)
bridge <- bridge[!duplicated(bridge$mo_id), , drop = FALSE]

empty_to_na <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "<NA>")] <- NA_character_
  x
}

read_counts <- function(path, id_column) {
  tab <- read.delim(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE)
  counts <- as.matrix(tab[, setdiff(colnames(tab), id_column), drop = FALSE])
  rownames(counts) <- tab[[id_column]]
  storage.mode(counts) <- "integer"
  counts
}

resolve_diagnosis <- function(samples) {
  merged <- merge(
    samples,
    bridge[, c("mo_id", "nygc_latest_direct_rna_id_match", "resolved_disease_group", "resolved_conflict_note")],
    by.x = "library",
    by.y = "mo_id",
    all.x = TRUE,
    sort = FALSE
  )
  merged <- merged[match(samples$sample_col, merged$sample_col), , drop = FALSE]
  direct <- tolower(empty_to_na(merged$nygc_latest_direct_rna_id_match)) %in% c("yes", "true", "1")
  resolved_primary <- ifelse(merged$resolved_disease_group %in% c("ALS", "Control"), merged$resolved_disease_group, NA_character_)
  merged$disease_binary_w33 <- merged$disease_binary
  merged$disease_binary_w33[direct] <- resolved_primary[direct]
  merged$label_source_w33 <- ifelse(
    direct,
    "synapse_nygc_resolved_exact_match",
    "geo_fallback_no_exact_synapse_match"
  )
  merged$label_decision_w33 <- ifelse(
    direct & is.na(resolved_primary),
    "exclude_nonprimary_resolved_diagnosis",
    ifelse(
      direct & merged$disease_binary != merged$disease_binary_w33,
      "reclassified_by_synapse_nygc",
      "retained"
    )
  )
  merged
}

build_design <- function(samples) {
  samples$disease_binary_w33 <- factor(samples$disease_binary_w33, levels = c("Control", "ALS"))
  samples$sex <- factor(empty_to_na(samples$sex))
  samples$tissue_subregion <- factor(empty_to_na(samples$tissue_subregion))
  samples$age <- suppressWarnings(as.numeric(samples$age))
  samples$PMI <- suppressWarnings(as.numeric(samples$PMI))

  candidates <- c("age", "sex", "PMI", "tissue_subregion")
  used <- c()
  dropped <- c()
  reasons <- c()
  for (covar in candidates) {
    value <- samples[[covar]]
    non_missing <- sum(!is.na(value))
    n_levels <- length(unique(value[!is.na(value)]))
    if (non_missing == 0 || n_levels < 2) {
      dropped <- c(dropped, covar)
      reason <- if (non_missing == 0) "all_missing" else "constant_after_filtering"
      if (n_levels < 2) reason <- paste(reason, "levels_lt_2", sep = ";")
      reasons <- c(reasons, paste0(covar, ":", reason))
    } else {
      used <- c(used, covar)
    }
  }

  keep <- complete.cases(samples[, c("disease_binary_w33", used), drop = FALSE])
  samples_model <- samples[keep, , drop = FALSE]
  low_power_drop_order <- c("tissue_subregion", "PMI", "sex", "age")
  repeat {
    n_als <- length(unique(samples_model$donor[samples_model$disease_binary_w33 == "ALS"]))
    n_control <- length(unique(samples_model$donor[samples_model$disease_binary_w33 == "Control"]))
    if (n_als >= 5 && n_control >= 5) break
    missing_used <- used[vapply(used, function(covar) any(is.na(samples[[covar]])), logical(1))]
    removable <- intersect(low_power_drop_order, missing_used)
    if (length(removable) == 0) break
    drop_one <- removable[1]
    used <- setdiff(used, drop_one)
    dropped <- c(dropped, drop_one)
    reasons <- c(reasons, paste0(drop_one, ":complete_case_donor_count_below_threshold"))
    keep <- complete.cases(samples[, c("disease_binary_w33", used), drop = FALSE])
    samples_model <- samples[keep, , drop = FALSE]
  }

  rank_drop_order <- c("tissue_subregion", "PMI", "sex", "age")
  repeat {
    formula_text <- paste("~ disease_binary_w33", if (length(used) > 0) paste("+", paste(used, collapse = " + ")) else "")
    design <- model.matrix(as.formula(formula_text), data = samples_model)
    if (qr(design)$rank == ncol(design)) break
    removable <- intersect(rank_drop_order, used)
    if (length(removable) == 0) break
    drop_one <- removable[1]
    used <- setdiff(used, drop_one)
    dropped <- c(dropped, drop_one)
    reasons <- c(reasons, paste0(drop_one, ":rank_deficient"))
  }
  formula_text <- paste("~ disease_binary_w33", if (length(used) > 0) paste("+", paste(used, collapse = " + ")) else "")
  design <- model.matrix(as.formula(formula_text), data = samples_model)
  list(
    samples = samples_model,
    design = design,
    formula_text = formula_text,
    used = used,
    dropped = unique(dropped),
    reasons = unique(reasons),
    rank = qr(design)$rank
  )
}

run_context <- function(context) {
  samples <- read.delim(context$samples_path, stringsAsFactors = FALSE, check.names = FALSE)
  counts <- read_counts(context$counts_path, context$id_column)
  samples <- samples[match(colnames(counts), samples$sample_col), , drop = FALSE]
  samples <- resolve_diagnosis(samples)
  samples$context_id <- context$context_id
  samples$used_in_w33 <- FALSE
  base_primary <- samples$disease_binary %in% c("ALS", "Control")
  model_candidates <- samples[base_primary & samples$disease_binary_w33 %in% c("ALS", "Control"), , drop = FALSE]
  design_info <- build_design(model_candidates)
  samples$used_in_w33[match(design_info$samples$sample_col, samples$sample_col)] <- TRUE
  counts_model <- counts[, design_info$samples$sample_col, drop = FALSE]
  n_als <- length(unique(design_info$samples$donor[design_info$samples$disease_binary_w33 == "ALS"]))
  n_control <- length(unique(design_info$samples$donor[design_info$samples$disease_binary_w33 == "Control"]))

  audit <- data.frame(
    context_id = context$context_id,
    modality = context$modality,
    n_libraries_raw_label_primary = sum(base_primary),
    n_libraries_resolved_label_primary = nrow(model_candidates),
    n_libraries_model = ncol(counts_model),
    n_als_donors_model = n_als,
    n_control_donors_model = n_control,
    n_label_changed_or_excluded = sum(samples$label_decision_w33 != "retained"),
    model_formula = design_info$formula_text,
    covariates_used = paste(design_info$used, collapse = ";"),
    covariates_dropped = paste(design_info$dropped, collapse = ";"),
    covariate_drop_reasons = paste(design_info$reasons, collapse = "|"),
    design_ncol = ncol(design_info$design),
    design_rank = design_info$rank,
    residual_df = ncol(counts_model) - design_info$rank,
    status = "not_run",
    message = "",
    stringsAsFactors = FALSE
  )
  write.table(samples, file.path(out, paste0(context$context_id, ".library_label_audit.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  if (n_als < 5 || n_control < 5 || design_info$rank < ncol(design_info$design) || !("disease_binary_w33ALS" %in% colnames(design_info$design))) {
    audit$status <- "skipped"
    audit$message <- "Insufficient resolved-label donor coverage or rank-deficient design"
    return(list(audit = audit, comparison = NULL))
  }

  y <- DGEList(counts = counts_model)
  keep_feature <- filterByExpr(y, design = design_info$design)
  y <- y[keep_feature, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  v <- voom(y, design_info$design, plot = FALSE)
  fit <- eBayes(lmFit(v, design_info$design), trend = TRUE)
  current <- topTable(fit, coef = "disease_binary_w33ALS", number = Inf, sort.by = "P")
  current$feature_id <- rownames(current)
  current$context_id <- context$context_id
  current$model_formula <- design_info$formula_text
  current$n_libraries_model <- ncol(counts_model)
  current$n_als_donors_model <- n_als
  current$n_control_donors_model <- n_control
  write.table(current, gzfile(file.path(out, "tables", paste0(context$context_id, ".synapse_resolved.tsv.gz"))), sep = "\t", quote = FALSE, row.names = FALSE)

  frozen <- read.delim(gzfile(context$frozen_path), stringsAsFactors = FALSE, check.names = FALSE)
  frozen$feature_id <- frozen[[context$id_column]]
  compare <- merge(
    frozen[, c("feature_id", "logFC", "P.Value", "adj.P.Val")],
    current[, c("feature_id", "logFC", "P.Value", "adj.P.Val")],
    by = "feature_id",
    suffixes = c("_frozen", "_w33")
  )
  compare$frozen_fdr05 <- compare$adj.P.Val_frozen < 0.05
  compare$w33_fdr05 <- compare$adj.P.Val_w33 < 0.05
  compare$direction_retained <- sign(compare$logFC_frozen) == sign(compare$logFC_w33)
  compare$context_id <- context$context_id
  write.table(compare, gzfile(file.path(out, "tables", paste0(context$context_id, ".frozen_vs_w33_feature_comparison.tsv.gz"))), sep = "\t", quote = FALSE, row.names = FALSE)

  frozen_set <- compare[compare$frozen_fdr05, , drop = FALSE]
  summary <- data.frame(
    context_id = context$context_id,
    modality = context$modality,
    frozen_features_tested = nrow(frozen),
    w33_features_tested = nrow(current),
    frozen_fdr05 = sum(compare$frozen_fdr05),
    w33_fdr05 = sum(compare$w33_fdr05),
    frozen_fdr05_direction_retained = sum(frozen_set$direction_retained),
    frozen_fdr05_direction_retained_fraction = mean(frozen_set$direction_retained),
    frozen_fdr05_retained_at_w33_fdr05 = sum(frozen_set$w33_fdr05),
    frozen_fdr05_retained_at_w33_nominal_p05 = sum(frozen_set$P.Value_w33 < 0.05),
    all_feature_logFC_pearson_r = cor(compare$logFC_frozen, compare$logFC_w33),
    all_feature_logFC_spearman_rho = cor(compare$logFC_frozen, compare$logFC_w33, method = "spearman"),
    stringsAsFactors = FALSE
  )
  audit$status <- "completed"
  audit$message <- paste0("Completed; W33 FDR05 features=", sum(compare$w33_fdr05))
  list(audit = audit, comparison = summary)
}

input_root <- file.path(root, "intermediate")
contexts <- list(
  list(
    context_id = "MTC_Astrocyte_RNA",
    modality = "RNA",
    samples_path = file.path(input_root, "W01_RNA_DE", "GSE330130_MTC_RNA", "Astrocyte", "samples.tsv"),
    counts_path = file.path(input_root, "W01_RNA_DE", "GSE330130_MTC_RNA", "Astrocyte", "counts.tsv.gz"),
    frozen_path = file.path(root, "results", "W01_RNA_DE", "tables", "GSE330130_MTC_RNA__Astrocyte.limma_voom.tsv.gz"),
    id_column = "gene_id"
  ),
  list(
    context_id = "MTC_Astrocyte_ATAC",
    modality = "ATAC",
    samples_path = file.path(input_root, "W03_ATAC_DA", "GSE330130_MTC_ATAC", "Astrocyte", "samples.tsv"),
    counts_path = file.path(input_root, "W03_ATAC_DA", "GSE330130_MTC_ATAC", "Astrocyte", "counts.tsv.gz"),
    frozen_path = file.path(root, "results", "W03_ATAC_DA", "tables", "GSE330130_MTC_ATAC__Astrocyte.DA_limma_voom.tsv.gz"),
    id_column = "peak_id"
  ),
  list(
    context_id = "LSC_Microglia_RNA",
    modality = "RNA",
    samples_path = file.path(input_root, "W01_RNA_DE", "GSE330130_LSC_RNA", "Microglia", "samples.tsv"),
    counts_path = file.path(input_root, "W01_RNA_DE", "GSE330130_LSC_RNA", "Microglia", "counts.tsv.gz"),
    frozen_path = file.path(root, "results", "W01_RNA_DE", "tables", "GSE330130_LSC_RNA__Microglia.limma_voom.tsv.gz"),
    id_column = "gene_id"
  ),
  list(
    context_id = "LSC_Microglia_ATAC",
    modality = "ATAC",
    samples_path = file.path(input_root, "W03_ATAC_DA", "GSE330130_LSC_ATAC", "Microglia", "samples.tsv"),
    counts_path = file.path(input_root, "W03_ATAC_DA", "GSE330130_LSC_ATAC", "Microglia", "counts.tsv.gz"),
    frozen_path = file.path(root, "results", "W03_ATAC_DA", "tables", "GSE330130_LSC_ATAC__Microglia.DA_limma_voom.tsv.gz"),
    id_column = "peak_id"
  )
)

for (context in contexts) {
  for (path in unlist(context[c("samples_path", "counts_path", "frozen_path")])) {
    if (!file.exists(path)) stop("Missing W33 input: ", path)
  }
}

results <- lapply(contexts, run_context)
audit <- do.call(rbind, lapply(results, `[[`, "audit"))
comparison <- do.call(rbind, lapply(results, `[[`, "comparison"))
write.table(audit, file.path(out, "w33_model_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(comparison, file.path(out, "w33_frozen_feature_concordance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  "# W33 Synapse/NYGC-resolved diagnosis sensitivity",
  "",
  "W33 does not replace frozen W01/W03. It reuses their lead-context counts, feature filtering, TMM, limma-voom and covariate-selection rules.",
  "For exact NYGC ST3 RNA-ID matches, resolved diagnosis is authoritative: ALS and Control remain primary; ALS-FTD and Other are excluded from an ALS-versus-Control sensitivity model.",
  "For libraries without an exact match, original GEO labels are retained. All changed/excluded libraries, model membership and frozen-versus-W33 feature comparisons are written here.",
  "This is required because the frozen analyses used source-object GEO labels while later clinical reconciliation identified exact-match diagnosis conflicts."
), file.path(out, "README.md"))
cat("[W33] wrote ", out, "\n", sep = "")
print(audit)
print(comparison)
