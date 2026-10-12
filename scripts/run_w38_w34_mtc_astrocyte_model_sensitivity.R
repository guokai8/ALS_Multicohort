#!/usr/bin/env Rscript

# W38 rebuilds the MTC astrocyte model-sensitivity analyses on the exact
# Synapse/NYGC-resolved W33/W34 donor sets. The full model must reproduce W33
# before nested-model or BA4-only outputs are accepted.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "p01_paths.R"))

root <- p01_project_root()
w33_root <- file.path(root, "results", "W33_synapse_resolved_diagnosis_sensitivity")
out <- file.path(root, "results", "W38_W34_MTC_astrocyte_model_sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)

model_specs <- list(
  list(model_id = "D0_disease_only", formula = "~ disease_binary_w33"),
  list(model_id = "D1_disease_sex", formula = "~ disease_binary_w33 + sex"),
  list(model_id = "D2_disease_age_sex", formula = "~ disease_binary_w33 + age + sex"),
  list(model_id = "D3_disease_age_sex_PMI", formula = "~ disease_binary_w33 + age + sex + PMI"),
  list(
    model_id = "D4_W34_primary_full",
    formula = "~ disease_binary_w33 + age + sex + PMI + tissue_subregion"
  )
)

specs <- list(
  list(
    modality = "RNA",
    context_id = "MTC_Astrocyte_RNA",
    dataset_id = "GSE330130_MTC_RNA",
    input_dir = file.path(
      root, "intermediate", "W01_RNA_DE", "GSE330130_MTC_RNA", "Astrocyte"
    ),
    w33_table = file.path(
      w33_root, "tables", "MTC_Astrocyte_RNA.synapse_resolved.tsv.gz"
    ),
    label_audit = file.path(w33_root, "MTC_Astrocyte_RNA.library_label_audit.tsv"),
    feature_id = "gene_id",
    expected_n_primary = 29L
  ),
  list(
    modality = "ATAC",
    context_id = "MTC_Astrocyte_ATAC",
    dataset_id = "GSE330130_MTC_ATAC",
    input_dir = file.path(
      root, "intermediate", "W03_ATAC_DA", "GSE330130_MTC_ATAC", "Astrocyte"
    ),
    w33_table = file.path(
      w33_root, "tables", "MTC_Astrocyte_ATAC.synapse_resolved.tsv.gz"
    ),
    label_audit = file.path(w33_root, "MTC_Astrocyte_ATAC.library_label_audit.tsv"),
    feature_id = "peak_id",
    expected_n_primary = 1705L
  )
)

read_counts <- function(path, feature_id) {
  tab <- read.delim(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE)
  if (!feature_id %in% colnames(tab)) stop("Missing ", feature_id, " in ", path)
  counts <- as.matrix(tab[, setdiff(colnames(tab), feature_id), drop = FALSE])
  rownames(counts) <- tab[[feature_id]]
  storage.mode(counts) <- "integer"
  counts
}

prepare_samples <- function(path) {
  samples <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  samples <- samples[samples$used_in_w33 %in% c(TRUE, "TRUE", "true", 1, "1"), , drop = FALSE]
  samples$disease_binary_w33 <- factor(samples$disease_binary_w33, levels = c("Control", "ALS"))
  samples$sex <- factor(samples$sex)
  samples$tissue_subregion <- factor(samples$tissue_subregion)
  samples$age <- suppressWarnings(as.numeric(samples$age))
  samples$PMI <- suppressWarnings(as.numeric(samples$PMI))
  required <- c("disease_binary_w33", "age", "sex", "PMI", "tissue_subregion")
  if (!all(complete.cases(samples[, required, drop = FALSE]))) {
    stop("W38 label audit contains incomplete full-model samples: ", path)
  }
  if (length(unique(samples$donor)) != nrow(samples)) {
    stop("W38 requires one library per donor: ", path)
  }
  samples
}

fit_model <- function(counts, samples, keep_feature, model_spec) {
  design <- model.matrix(as.formula(model_spec$formula), data = samples)
  if (qr(design)$rank != ncol(design)) {
    stop("Rank-deficient W38 design: ", model_spec$model_id)
  }
  if (!"disease_binary_w33ALS" %in% colnames(design)) {
    stop("Missing W38 disease coefficient: ", model_spec$model_id)
  }
  y <- DGEList(counts = counts[keep_feature, , drop = FALSE])
  y <- calcNormFactors(y, method = "TMM")
  v <- voom(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design), trend = TRUE)
  tt <- topTable(
    fit,
    coef = "disease_binary_w33ALS",
    number = Inf,
    sort.by = "P"
  )
  tt
}

nested_audit_rows <- list()
nested_feature_rows <- list()
ba4_summary_rows <- list()
ba4_feature_rows <- list()
reproduction_rows <- list()
manifest_rows <- list()

for (spec in specs) {
  samples <- prepare_samples(spec$label_audit)
  counts <- read_counts(file.path(spec$input_dir, "counts.tsv.gz"), spec$feature_id)
  if (!all(samples$sample_col %in% colnames(counts))) {
    stop("W38 sample columns are missing from counts for ", spec$modality)
  }
  counts <- counts[, samples$sample_col, drop = FALSE]

  full_design <- model.matrix(
    ~ disease_binary_w33 + age + sex + PMI + tissue_subregion,
    data = samples
  )
  if (qr(full_design)$rank != ncol(full_design)) {
    stop("Rank-deficient W38 full design for ", spec$modality)
  }
  keep_feature <- filterByExpr(DGEList(counts = counts), design = full_design)

  fitted <- list()
  for (model_spec in model_specs) {
    tt <- fit_model(counts, samples, keep_feature, model_spec)
    tt[[spec$feature_id]] <- rownames(tt)
    tt$feature_id <- rownames(tt)
    tt$modality <- spec$modality
    tt$model_id <- model_spec$model_id
    tt$model_formula <- model_spec$formula
    fitted[[model_spec$model_id]] <- tt
  }

  w33 <- read.delim(gzfile(spec$w33_table), stringsAsFactors = FALSE, check.names = FALSE)
  w33$feature_id <- w33$feature_id
  primary <- fitted[["D4_W34_primary_full"]]
  reproduced <- merge(
    w33[, c("feature_id", "logFC", "P.Value", "adj.P.Val")],
    primary[, c("feature_id", "logFC", "P.Value", "adj.P.Val")],
    by = "feature_id",
    suffixes = c("_w33", "_w38")
  )
  max_logfc_diff <- max(abs(reproduced$logFC_w33 - reproduced$logFC_w38))
  max_p_diff <- max(abs(reproduced$P.Value_w33 - reproduced$P.Value_w38))
  max_bh_diff <- max(abs(reproduced$adj.P.Val_w33 - reproduced$adj.P.Val_w38))
  n_primary <- sum(primary$adj.P.Val < 0.05)
  reproduction_pass <- (
    nrow(reproduced) == nrow(w33) &&
      n_primary == spec$expected_n_primary &&
      max_logfc_diff < 1e-10 &&
      max_p_diff < 1e-10 &&
      max_bh_diff < 1e-10
  )
  reproduction_rows[[length(reproduction_rows) + 1]] <- data.frame(
    modality = spec$modality,
    n_w33_features = nrow(w33),
    n_w38_features = nrow(primary),
    n_common_features = nrow(reproduced),
    expected_primary_fdr05 = spec$expected_n_primary,
    observed_primary_fdr05 = n_primary,
    max_abs_logFC_difference = max_logfc_diff,
    max_abs_p_value_difference = max_p_diff,
    max_abs_BH_difference = max_bh_diff,
    reproduction_pass = reproduction_pass,
    stringsAsFactors = FALSE
  )
  if (!reproduction_pass) {
    stop("W38 failed to reproduce W33/W34 primary model for ", spec$modality)
  }

  primary_small <- primary[, c("feature_id", "logFC", "P.Value", "adj.P.Val")]
  colnames(primary_small)[2:4] <- c(
    "primary_logFC", "primary_P_Value", "primary_adj_P_Val"
  )
  primary_ids <- primary_small$feature_id[primary_small$primary_adj_P_Val < 0.05]

  for (model_spec in model_specs) {
    current <- fitted[[model_spec$model_id]]
    current_small <- current[, c("feature_id", "logFC", "P.Value", "adj.P.Val")]
    colnames(current_small)[2:4] <- c(
      "model_logFC", "model_P_Value", "model_adj_P_Val"
    )
    comparison <- merge(primary_small, current_small, by = "feature_id", all = FALSE)
    primary_features <- comparison[comparison$feature_id %in% primary_ids, , drop = FALSE]
    primary_features$direction_retained <- (
      sign(primary_features$primary_logFC) == sign(primary_features$model_logFC)
    )
    primary_features$nominal_directional <- (
      primary_features$direction_retained & primary_features$model_P_Value < 0.05
    )
    primary_features$fdr05_directional <- (
      primary_features$direction_retained & primary_features$model_adj_P_Val < 0.05
    )
    primary_features$modality <- spec$modality
    primary_features$model_id <- model_spec$model_id
    primary_features$model_formula <- model_spec$formula
    nested_feature_rows[[length(nested_feature_rows) + 1]] <- primary_features

    design <- model.matrix(as.formula(model_spec$formula), data = samples)
    nested_audit_rows[[length(nested_audit_rows) + 1]] <- data.frame(
      modality = spec$modality,
      model_id = model_spec$model_id,
      model_formula = model_spec$formula,
      n_libraries = nrow(samples),
      n_als_donors = sum(samples$disease_binary_w33 == "ALS"),
      n_control_donors = sum(samples$disease_binary_w33 == "Control"),
      design_ncol = ncol(design),
      residual_df = nrow(samples) - ncol(design),
      n_features_tested = nrow(current),
      n_model_fdr05 = sum(current$adj.P.Val < 0.05),
      n_primary_fdr05 = length(primary_ids),
      n_primary_direction_retained = sum(primary_features$direction_retained),
      n_primary_nominal_directional = sum(primary_features$nominal_directional),
      n_primary_fdr05_directional = sum(primary_features$fdr05_directional),
      all_feature_logFC_pearson_r_vs_primary = cor(
        comparison$primary_logFC, comparison$model_logFC
      ),
      all_feature_logFC_spearman_rho_vs_primary = cor(
        comparison$primary_logFC, comparison$model_logFC, method = "spearman"
      ),
      stringsAsFactors = FALSE
    )
  }

  ba4_samples <- samples[samples$tissue_subregion == "Cortex_Motor_BA4", , drop = FALSE]
  if (
    sum(ba4_samples$disease_binary_w33 == "ALS") < 3 ||
      sum(ba4_samples$disease_binary_w33 == "Control") < 3
  ) {
    stop("W38 BA4 restriction lacks at least three donors per group for ", spec$modality)
  }
  ba4_counts <- counts[, ba4_samples$sample_col, drop = FALSE]
  ba4_design <- model.matrix(~ disease_binary_w33, data = ba4_samples)
  ba4_keep <- filterByExpr(DGEList(counts = ba4_counts), design = ba4_design)
  ba4 <- fit_model(
    ba4_counts,
    ba4_samples,
    ba4_keep,
    list(model_id = "BA4_only", formula = "~ disease_binary_w33")
  )
  ba4$feature_id <- rownames(ba4)
  ba4_small <- ba4[, c("feature_id", "logFC", "P.Value", "adj.P.Val")]
  colnames(ba4_small)[2:4] <- c("BA4_logFC", "BA4_P_Value", "BA4_adj_P_Val")
  ba4_comparison <- merge(primary_small, ba4_small, by = "feature_id", all = FALSE)
  ba4_primary <- ba4_comparison[
    ba4_comparison$feature_id %in% primary_ids,
    ,
    drop = FALSE
  ]
  ba4_primary$direction_retained <- (
    sign(ba4_primary$primary_logFC) == sign(ba4_primary$BA4_logFC)
  )
  ba4_primary$nominal_directional <- (
    ba4_primary$direction_retained & ba4_primary$BA4_P_Value < 0.05
  )
  ba4_primary$fdr05_directional <- (
    ba4_primary$direction_retained & ba4_primary$BA4_adj_P_Val < 0.05
  )
  ba4_primary$modality <- spec$modality
  ba4_feature_rows[[length(ba4_feature_rows) + 1]] <- ba4_primary
  ba4_summary_rows[[length(ba4_summary_rows) + 1]] <- data.frame(
    modality = spec$modality,
    n_libraries = nrow(ba4_samples),
    n_als_donors = sum(ba4_samples$disease_binary_w33 == "ALS"),
    n_control_donors = sum(ba4_samples$disease_binary_w33 == "Control"),
    n_features_tested = nrow(ba4),
    n_BA4_fdr05 = sum(ba4$adj.P.Val < 0.05),
    n_primary_fdr05 = length(primary_ids),
    n_primary_features_available_in_BA4 = nrow(ba4_primary),
    n_primary_direction_retained = sum(ba4_primary$direction_retained),
    n_primary_nominal_directional = sum(ba4_primary$nominal_directional),
    n_primary_fdr05_directional = sum(ba4_primary$fdr05_directional),
    primary_feature_direction_retained_fraction = mean(ba4_primary$direction_retained),
    all_feature_logFC_pearson_r_vs_primary = cor(
      ba4_comparison$primary_logFC, ba4_comparison$BA4_logFC
    ),
    all_feature_logFC_spearman_rho_vs_primary = cor(
      ba4_comparison$primary_logFC, ba4_comparison$BA4_logFC, method = "spearman"
    ),
    stringsAsFactors = FALSE
  )

  manifest <- samples[, c(
    "sample_col", "sample_id", "library", "donor", "disease_binary_w33",
    "sex", "age", "PMI", "tissue_subregion", "n_cells", "label_source_w33",
    "label_decision_w33"
  )]
  manifest$modality <- spec$modality
  manifest$used_in_nested_models <- TRUE
  manifest$used_in_BA4_model <- manifest$tissue_subregion == "Cortex_Motor_BA4"
  manifest_rows[[length(manifest_rows) + 1]] <- manifest
}

reproduction <- do.call(rbind, reproduction_rows)
nested_audit <- do.call(rbind, nested_audit_rows)
nested_features <- do.call(rbind, nested_feature_rows)
ba4_summary <- do.call(rbind, ba4_summary_rows)
ba4_features <- do.call(rbind, ba4_feature_rows)
manifest <- do.call(rbind, manifest_rows)

write.table(
  reproduction,
  file.path(out, "w38_primary_reproduction_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  nested_audit,
  file.path(out, "w38_nested_model_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  nested_features,
  gzfile(file.path(out, "tables", "w38_primary_feature_nested_sensitivity.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  ba4_summary,
  file.path(out, "w38_BA4_restriction_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  ba4_features,
  gzfile(file.path(out, "tables", "w38_BA4_primary_feature_sensitivity.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  manifest,
  file.path(out, "w38_resolved_library_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(c(
  "# W38 W34 motor-cortex astrocyte model sensitivity",
  "",
  "- W38 uses the exact W33/W34 Synapse/NYGC-resolved model libraries and disease labels.",
  "- The RNA and ATAC full models must reproduce W33 feature statistics before sensitivity outputs are accepted.",
  "- Nested models use fixed complete-case donor sets and a feature universe filtered by the W34 full model.",
  "- BA4-only analyses are post hoc directional sensitivities with three ALS and five control donors.",
  "- W38 supports model-dependence and direction-stability statements only; it does not establish subregion replication, cell-intrinsic causality, TF activity or enhancer function."
), file.path(out, "README.md"))

cat("[W38] wrote ", out, "\n", sep = "")
print(reproduction)
print(nested_audit)
print(ba4_summary)
