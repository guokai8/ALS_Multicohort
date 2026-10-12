#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, mustWork = FALSE)) else getwd()
source(file.path(script_dir, "p01_paths.R"))

root <- p01_project_root()
input_root <- file.path(root, "intermediate/W01_RNA_DE")
output_root <- file.path(root, "results/W01_RNA_DE")

resolve_project_path <- function(path) {
  if (grepl("^/|^[A-Za-z]:", path)) {
    return(path)
  }
  file.path(root, path)
}

project_rel <- function(path) {
  if (is.na(path) || !nzchar(path)) {
    return(path)
  }
  root_norm <- normalizePath(root, mustWork = FALSE)
  path_norm <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(root_norm, .Platform$file.sep)
  if (startsWith(path_norm, prefix)) {
    return(substring(path_norm, nchar(prefix) + 1))
  }
  path
}

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "module_seed_lists"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "logs"), recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(input_root, "w01_export_manifest.tsv")
if (!file.exists(manifest_path)) {
  stop("Missing export manifest: ", manifest_path)
}

manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "NA", x)
}

empty_to_na <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "NA" | x == "NaN" | x == "<NA>"] <- NA
  x
}

build_design <- function(samples) {
  samples$disease_binary <- factor(samples$disease_binary, levels = c("Control", "ALS"))
  samples$sex <- factor(empty_to_na(samples$sex))
  samples$tissue_subregion <- factor(empty_to_na(samples$tissue_subregion))
  samples$age <- suppressWarnings(as.numeric(samples$age))
  samples$PMI <- suppressWarnings(as.numeric(samples$PMI))

  candidates <- c("age", "sex", "PMI", "tissue_subregion")
  used <- c()
  dropped <- c()
  drop_reasons <- c()

  for (covar in candidates) {
    value <- samples[[covar]]
    non_missing <- sum(!is.na(value))
    n_levels <- if (is.factor(value)) length(unique(value[!is.na(value)])) else length(unique(value[!is.na(value)]))
    if (non_missing == 0 || n_levels < 2) {
      dropped <- c(dropped, covar)
      reason <- if (non_missing == 0) {
        "all_missing"
      } else {
        "constant_after_filtering"
      }
      if (n_levels < 2) {
        reason <- paste(reason, "levels_lt_2", sep = ";")
      }
      drop_reasons <- c(drop_reasons, paste0(covar, ":", reason))
    } else {
      used <- c(used, covar)
    }
  }

  keep <- complete.cases(samples[, c("disease_binary", used), drop = FALSE])
  samples_model <- samples[keep, , drop = FALSE]

  # If complete-case filtering removes too many biological replicates, drop
  # missingness-driving covariates before looking at any disease effect.
  low_power_drop_order <- c("tissue_subregion", "PMI", "sex", "age")
  repeat {
    n_als <- length(unique(samples_model$donor[samples_model$disease_binary == "ALS"]))
    n_control <- length(unique(samples_model$donor[samples_model$disease_binary == "Control"]))
    if (n_als >= 5 && n_control >= 5) {
      break
    }
    missing_used <- used[vapply(used, function(covar) any(is.na(samples[[covar]])), logical(1))]
    removable <- intersect(low_power_drop_order, missing_used)
    if (length(removable) == 0) {
      break
    }
    drop_one <- removable[1]
    used <- setdiff(used, drop_one)
    dropped <- c(dropped, drop_one)
    drop_reasons <- c(drop_reasons, paste0(drop_one, ":complete_case_donor_count_below_threshold"))
    keep <- complete.cases(samples[, c("disease_binary", used), drop = FALSE])
    samples_model <- samples[keep, , drop = FALSE]
  }

  # Build rank-stable design. Drop least-central technical covariates first.
  rank_drop_order <- c("tissue_subregion", "PMI", "sex", "age")
  repeat {
    formula_text <- paste("~ disease_binary", if (length(used) > 0) paste("+", paste(used, collapse = " + ")) else "")
    design <- model.matrix(as.formula(formula_text), data = samples_model)
    qr_rank <- qr(design)$rank
    if (qr_rank == ncol(design)) {
      break
    }
    removable <- intersect(rank_drop_order, used)
    if (length(removable) == 0) {
      break
    }
    drop_one <- removable[1]
    used <- setdiff(used, drop_one)
    dropped <- c(dropped, drop_one)
    drop_reasons <- c(drop_reasons, paste0(drop_one, ":rank_deficient"))
  }

  formula_text <- paste("~ disease_binary", if (length(used) > 0) paste("+", paste(used, collapse = " + ")) else "")
  design <- model.matrix(as.formula(formula_text), data = samples_model)

  list(
    samples = samples_model,
    design = design,
    formula_text = formula_text,
    used = used,
    dropped = unique(dropped),
    drop_reasons = unique(drop_reasons),
    rank = qr(design)$rank
  )
}

run_one <- function(row) {
  dataset_id <- row[["dataset_id"]]
  cell_type <- row[["cell_type"]]
  safe_id <- paste(sanitize(dataset_id), sanitize(cell_type), sep = "__")

  counts_df <- read.delim(gzfile(resolve_project_path(row[["counts_path"]])), stringsAsFactors = FALSE, check.names = FALSE)
  genes <- read.delim(resolve_project_path(row[["genes_path"]]), stringsAsFactors = FALSE, check.names = FALSE)
  samples <- read.delim(resolve_project_path(row[["samples_path"]]), stringsAsFactors = FALSE, check.names = FALSE)

  stopifnot("gene_id" %in% colnames(counts_df))
  stopifnot(all(samples$sample_col %in% colnames(counts_df)))

  counts <- as.matrix(counts_df[, samples$sample_col, drop = FALSE])
  rownames(counts) <- counts_df$gene_id
  storage.mode(counts) <- "integer"

  samples <- samples[match(colnames(counts), samples$sample_col), , drop = FALSE]
  samples <- samples[samples$disease_binary %in% c("ALS", "Control"), , drop = FALSE]
  counts <- counts[, samples$sample_col, drop = FALSE]

  n_als_before <- length(unique(samples$donor[samples$disease_binary == "ALS"]))
  n_control_before <- length(unique(samples$donor[samples$disease_binary == "Control"]))

  design_info <- build_design(samples)
  samples_model <- design_info$samples
  counts_model <- counts[, samples_model$sample_col, drop = FALSE]
  sample_log_path <- file.path(output_root, "logs", paste0(safe_id, ".model_samples.tsv"))
  sample_log <- samples
  sample_log$used_in_model <- sample_log$sample_col %in% samples_model$sample_col
  write.table(sample_log, sample_log_path, sep = "\t", quote = FALSE, row.names = FALSE)

  n_als <- length(unique(samples_model$donor[samples_model$disease_binary == "ALS"]))
  n_control <- length(unique(samples_model$donor[samples_model$disease_binary == "Control"]))

  audit <- data.frame(
    dataset_id = dataset_id,
    cell_type = cell_type,
    method = "edgeR_TMM_limma_voom_eBayes_trend",
    model_formula = design_info$formula_text,
    covariates_used = paste(design_info$used, collapse = ";"),
    covariates_dropped = paste(design_info$dropped, collapse = ";"),
    covariate_drop_reasons = paste(design_info$drop_reasons, collapse = "|"),
    n_samples_before_covariate_filter = ncol(counts),
    n_samples_model = ncol(counts_model),
    n_als_donors_before_covariate_filter = n_als_before,
    n_control_donors_before_covariate_filter = n_control_before,
    n_als_donors_model = n_als,
    n_control_donors_model = n_control,
    design_ncol = ncol(design_info$design),
    design_rank = design_info$rank,
    residual_df = ncol(counts_model) - design_info$rank,
    coefficient = "disease_binaryALS",
    n_genes_input = nrow(counts_model),
    n_genes_tested = NA_integer_,
    sample_log_path = project_rel(sample_log_path),
    status = "not_run",
    message = "",
    stringsAsFactors = FALSE
  )

  if (n_als < 5 || n_control < 5) {
    audit$status <- "skipped"
    audit$message <- "Insufficient donor counts after covariate/model filtering"
    return(list(audit = audit, result_path = NA_character_, seed_up = NA_character_, seed_down = NA_character_))
  }
  if (design_info$rank < ncol(design_info$design) || !("disease_binaryALS" %in% colnames(design_info$design))) {
    audit$status <- "skipped"
    audit$message <- "Rank-deficient design or missing ALS coefficient"
    return(list(audit = audit, result_path = NA_character_, seed_up = NA_character_, seed_down = NA_character_))
  }

  y <- DGEList(counts = counts_model)
  keep_gene <- filterByExpr(y, design = design_info$design)
  y <- y[keep_gene, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  v <- voom(y, design_info$design, plot = FALSE)
  fit <- lmFit(v, design_info$design)
  fit <- eBayes(fit, trend = TRUE)
  tt <- topTable(fit, coef = "disease_binaryALS", number = Inf, sort.by = "P")
  tt$gene_id <- rownames(tt)

  gene_cols <- c("gene_id")
  for (col in c("common_name", "gene_name", "symbol", "chrom", "id")) {
    if (col %in% colnames(genes)) {
      gene_cols <- unique(c(gene_cols, col))
    }
  }
  tt <- merge(tt, genes[, gene_cols, drop = FALSE], by = "gene_id", all.x = TRUE, sort = FALSE)
  tt <- tt[order(tt$P.Value), , drop = FALSE]
  tt$dataset_id <- dataset_id
  tt$cell_type <- cell_type
  tt$method <- audit$method
  tt$model_formula <- audit$model_formula
  tt$covariates_used <- audit$covariates_used
  tt$n_samples_model <- audit$n_samples_model
  tt$n_als_donors_model <- audit$n_als_donors_model
  tt$n_control_donors_model <- audit$n_control_donors_model
  tt$n_genes_tested <- nrow(tt)

  front <- c(
    "dataset_id",
    "cell_type",
    "gene_id",
    intersect(c("common_name", "gene_name", "symbol", "chrom", "id"), colnames(tt)),
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B",
    "method",
    "model_formula",
    "covariates_used",
    "n_samples_model",
    "n_als_donors_model",
    "n_control_donors_model",
    "n_genes_tested"
  )
  tt <- tt[, unique(c(front, setdiff(colnames(tt), front))), drop = FALSE]

  out_path <- file.path(output_root, "tables", paste0(safe_id, ".limma_voom.tsv.gz"))
  write.table(tt, gzfile(out_path), sep = "\t", quote = FALSE, row.names = FALSE)

  # Primary seed modules use donor-level FDR only.  The historical
  # FDR05_logFC025 lists are retained below for reproducibility audits.
  seed_up <- tt[!is.na(tt$adj.P.Val) & tt$adj.P.Val < 0.05 & tt$logFC > 0, , drop = FALSE]
  seed_down <- tt[!is.na(tt$adj.P.Val) & tt$adj.P.Val < 0.05 & tt$logFC < 0, , drop = FALSE]
  seed_up_path <- file.path(output_root, "module_seed_lists", paste0(safe_id, ".FDR05_no_abslogFC_up.tsv"))
  seed_down_path <- file.path(output_root, "module_seed_lists", paste0(safe_id, ".FDR05_no_abslogFC_down.tsv"))
  historical_up <- tt[!is.na(tt$adj.P.Val) & tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.25 & tt$logFC > 0, , drop = FALSE]
  historical_down <- tt[!is.na(tt$adj.P.Val) & tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.25 & tt$logFC < 0, , drop = FALSE]
  historical_up_path <- file.path(output_root, "module_seed_lists", paste0(safe_id, ".FDR05_logFC025_up.tsv"))
  historical_down_path <- file.path(output_root, "module_seed_lists", paste0(safe_id, ".FDR05_logFC025_down.tsv"))
  write.table(seed_up, seed_up_path, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(seed_down, seed_down_path, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(historical_up, historical_up_path, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(historical_down, historical_down_path, sep = "\t", quote = FALSE, row.names = FALSE)

  audit$n_genes_tested <- nrow(tt)
  audit$status <- "completed"
  audit$message <- paste0(
    "Completed; up_seeds=", nrow(seed_up),
    "; down_seeds=", nrow(seed_down),
    "; historical_up_seeds=", nrow(historical_up),
    "; historical_down_seeds=", nrow(historical_down)
  )

  list(
    audit = audit,
    result_path = project_rel(out_path),
    seed_up = project_rel(seed_up_path),
    seed_down = project_rel(seed_down_path),
    historical_seed_up = project_rel(historical_up_path),
    historical_seed_down = project_rel(historical_down_path)
  )
}

audits <- list()
result_manifest <- list()

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]
  message("Running W01: ", row$dataset_id, " / ", row$cell_type)
  result <- tryCatch(
    run_one(row),
    error = function(e) {
      data <- data.frame(
        dataset_id = row$dataset_id,
        cell_type = row$cell_type,
        method = "edgeR_TMM_limma_voom_eBayes_trend",
        model_formula = NA_character_,
        covariates_used = NA_character_,
        covariates_dropped = NA_character_,
        covariate_drop_reasons = NA_character_,
        n_samples_before_covariate_filter = as.integer(row$n_samples),
        n_samples_model = NA_integer_,
        n_als_donors_before_covariate_filter = as.integer(row$n_als_donors),
        n_control_donors_before_covariate_filter = as.integer(row$n_control_donors),
        n_als_donors_model = NA_integer_,
        n_control_donors_model = NA_integer_,
        design_ncol = NA_integer_,
        design_rank = NA_integer_,
        residual_df = NA_integer_,
        coefficient = "disease_binaryALS",
        n_genes_input = as.integer(row$n_genes),
        n_genes_tested = NA_integer_,
        sample_log_path = NA_character_,
        status = "error",
        message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      list(audit = data, result_path = NA_character_, seed_up = NA_character_, seed_down = NA_character_, historical_seed_up = NA_character_, historical_seed_down = NA_character_)
    }
  )
  audits[[length(audits) + 1]] <- result$audit
  result_manifest[[length(result_manifest) + 1]] <- data.frame(
    dataset_id = row$dataset_id,
    cell_type = row$cell_type,
    result_path = result$result_path,
    seed_up_path = result$seed_up,
    seed_down_path = result$seed_down,
    historical_seed_up_path = result$historical_seed_up,
    historical_seed_down_path = result$historical_seed_down,
    stringsAsFactors = FALSE
  )
}

audit_df <- do.call(rbind, audits)
manifest_df <- do.call(rbind, result_manifest)

write.table(
  audit_df,
  file.path(output_root, "w01_model_design_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  manifest_df,
  file.path(output_root, "w01_de_result_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_df <- audit_df[, c(
  "dataset_id",
  "cell_type",
  "status",
  "model_formula",
  "covariates_used",
  "covariates_dropped",
  "n_samples_model",
  "n_als_donors_model",
  "n_control_donors_model",
  "residual_df",
  "n_genes_tested",
  "message"
)]
write.table(
  summary_df,
  file.path(output_root, "w01_run_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
