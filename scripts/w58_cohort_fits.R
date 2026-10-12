# Shared W58 setup: libraries, paths, module seeds, cohort loaders and the
# per-cohort limma-voom / donor-level module fit (Amendment 4). Sourced by
# run_w58_multicohort_meta.R and run_w60_age70_sensitivity.R.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(metafor)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
res_dir <- file.path(root, "results")
out <- file.path(res_dir, "W58_multicohort_meta")
seed_dir <- file.path(res_dir, "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
DELTA_AGE_ULM <- -11.1
set.seed(20261011)

seeds <- lapply(c(up = "up", down = "down"), function(d) read.delim(file.path(seed_dir, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", d, ".tsv"))))
gene_map <- read.delim(file.path(p01_source_data_root(), "Other_Datasets/Pineda_snRNAseq/synapse_download/Processed Data/191112_ALS_110_snRNA-B9/row_metadata.tsv"), stringsAsFactors = FALSE)
unique_sym <- gene_map[!duplicated(gene_map$Gene) & !gene_map$Gene %in% gene_map$Gene[duplicated(gene_map$Gene)], ]
sym_to_ensg <- setNames(unique_sym$ENSEMBL, unique_sym$Gene)
read_counts <- function(path) as.matrix(read.delim(gzfile(path), row.names = 1, check.names = FALSE))
pineda_age <- function() {
  clin <- read.delim(file.path(res_dir, "W22_Pineda_S01_clinical_pathology_context", "w22_pineda_lead_module_donor_context.tsv"), stringsAsFactors = FALSE)
  clin <- clin[!duplicated(clin$Donor), c("Donor", "Age.of.death..Y.")]
  setNames(clin$Age.of.death..Y., clin$Donor)
}

# ---- cohort loaders: counts (genes x units), samples s with unit/donor/disease ----
load_li <- function() {
  d <- file.path(res_dir, "W54_Li2023_motor_cortex_validation")
  s <- read.delim(file.path(d, "w54_li2023_astrocyte_pseudobulk_samples.tsv"), stringsAsFactors = FALSE)
  s <- s[s$disease %in% c("C9-ALS", "Control"), ]
  s$unit <- s$GSM; s$donor <- s$donor; s$disease <- ifelse(s$disease == "C9-ALS", "ALS", "Control")
  list(name = "Li2023_BA4", region = "motor_cortex", counts = read_counts(file.path(d, "w54_li2023_astrocyte_pseudobulk_counts.tsv.gz")),
       s = s, rhs = "age + sex", ids = "ensg", has_age = TRUE)
}
load_ulm <- function() {
  d <- file.path(res_dir, "W52_Ulm_motor_cortex_validation")
  s <- read.delim(file.path(d, "w52_ulm_astrocyte_pseudobulk_samples.tsv"), stringsAsFactors = FALSE)
  s <- s[s$Case_Type %in% c("ALS", "HC"), ]
  s$unit <- s$stratum; s$donor <- s$ID; s$disease <- ifelse(s$Case_Type == "ALS", "ALS", "Control")
  s$sex <- tolower(s$Sex); s$chip <- s$Chip
  main <- tapply(seq_len(nrow(s)), s$donor, function(i) s$Chip[i][which.max(s$n_cells[i])])
  s$main_chip <- main[s$donor]
  list(name = "Ulm_M1", region = "motor_cortex", counts = read_counts(file.path(d, "w52_ulm_astrocyte_pseudobulk_counts.tsv.gz")),
       s = s, rhs = "sex + chip", rhs_module = "sex + main_chip", block = TRUE, ids = "symbol", has_age = FALSE)
}
load_pineda_ba4 <- function() {
  d <- file.path(res_dir, "W51_Pineda_BA4_discovery_matched_limma")
  s <- read.delim(file.path(d, "w51_pineda_astrocyte_pseudobulk_samples.tsv"), stringsAsFactors = FALSE)
  s <- s[s$Region == "BA4" & s$Condition %in% c("ALS", "PN") & s$lib_size > 0 & substr(s$sample_folder, 1, 6) %in% c("191112", "191114") & !is.na(s$age), ]
  s$unit <- s$sample_folder; s$donor <- s$Donor; s$disease <- ifelse(s$Condition == "ALS", "ALS", "Control"); s$sex <- s$Sex
  list(name = "Pineda_BA4_2019", region = "motor_cortex", counts = read_counts(file.path(d, "w51_pineda_astrocyte_pseudobulk_counts.tsv.gz")),
       s = s, rhs = "age + sex", ids = "ensg", has_age = TRUE)
}
load_pineda_ba9 <- function() {
  d <- file.path(res_dir, "W56_cell_source_check")
  s <- read.delim(file.path(d, "w56_pineda_lineage_strata.tsv"), stringsAsFactors = FALSE)
  s <- s[s$lineage == "Astrocyte" & s$Region == "BA9" & startsWith(s$sample_folder, "210526") & s$Condition %in% c("ALS", "PN") & s$lib_size > 0, ]
  s$age <- pineda_age()[as.character(s$Donor)]
  s <- s[!is.na(s$age), ]
  s$unit <- s$stratum; s$donor <- s$Donor; s$disease <- ifelse(s$Condition == "ALS", "ALS", "Control"); s$sex <- s$Sex
  counts <- read_counts(file.path(d, "w56_pineda_lineage_counts.tsv.gz"))
  list(name = "Pineda_BA9_210526", region = "prefrontal_cortex", counts = counts[, s$unit], s = s, rhs = "age + sex", ids = "ensg", has_age = TRUE)
}
load_zelic <- function() {
  s <- read.delim(file.path(out, "w58_zelic_astrocyte_pseudobulk_samples.tsv"), stringsAsFactors = FALSE)
  s$unit <- as.character(s$donor_id); s$donor <- s$unit; s$disease <- ifelse(s$disease == "ALS", "ALS", "Control")
  counts <- read_counts(file.path(out, "w58_zelic_astrocyte_pseudobulk_counts.tsv.gz"))
  list(name = "Zelic_cervicalSC", region = "spinal_cord", counts = counts, s = s, rhs = "age + sex + PMI + RIN", ids = "symbol", has_age = TRUE)
}
load_discovery <- function() {
  a <- read.delim(file.path(res_dir, "W33_synapse_resolved_diagnosis_sensitivity", "MTC_Astrocyte_RNA.library_label_audit.tsv"), stringsAsFactors = FALSE)
  s <- a[a$used_in_w33 %in% c(TRUE, "TRUE"), ]
  s$unit <- s$library; s$disease <- ifelse(s$disease_binary_w33 == "ALS", "ALS", "Control")
  list(name = "NYGC_MTC_discovery", region = "motor_cortex", counts = read_counts(file.path(res_dir, "W57_ambient_contamination_test", "w57_discovery_astro_counts.tsv.gz")),
       s = s, rhs = "age + sex + PMI + tissue_subregion", ids = "ensg", has_age = TRUE)
}

# ---- per-cohort fit ----
fit_cohort <- function(co) {
  s <- co$s
  s$disease <- factor(s$disease, levels = c("Control", "ALS"))
  for (v in intersect(c("sex", "chip", "main_chip", "tissue_subregion"), names(s))) s[[v]] <- factor(s[[v]])
  design <- model.matrix(as.formula(paste("~ disease +", co$rhs)), data = s)
  y <- DGEList(counts = co$counts[, s$unit, drop = FALSE])
  y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(y, design)
  if (isTRUE(co$block)) {
    cons <- duplicateCorrelation(v, design, block = s$donor)$consensus.correlation
    v <- voom(y, design, block = s$donor, correlation = cons)
    cons <- duplicateCorrelation(v, design, block = s$donor)$consensus.correlation
    fit <- eBayes(lmFit(v, design, block = s$donor, correlation = cons), trend = TRUE)
  } else {
    fit <- eBayes(lmFit(v, design), trend = TRUE)
  }
  coef_table <- function(coef) {
    data.frame(gene = rownames(fit), logFC = fit$coefficients[, coef], SE = fit$stdev.unscaled[, coef] * sqrt(fit$s2.post),
               P = fit$p.value[, coef], AveExpr = fit$Amean, stringsAsFactors = FALSE)
  }
  to_ensg <- function(tab) {
    if (co$ids == "symbol") {
      tab$gene <- unname(sym_to_ensg[tab$gene])
      tab <- tab[!is.na(tab$gene) & !duplicated(tab$gene), ]
    }
    tab
  }
  disease_tab <- to_ensg(coef_table("diseaseALS"))
  age_tab <- if (co$has_age) to_ensg(coef_table("age")) else NULL

  # Module scores: mean voom log-CPM over module genes present; donor-level model.
  ids <- if (co$ids == "symbol") lapply(seeds, `[[`, "common_name") else lapply(seeds, `[[`, "gene_id")
  rhs_m <- if (!is.null(co$rhs_module)) co$rhs_module else co$rhs
  mod <- do.call(rbind, lapply(names(ids), function(m) {
    g <- intersect(ids[[m]], rownames(v))
    unit_score <- colMeans(v$E[g, , drop = FALSE])
    dd <- aggregate(data.frame(score = unit_score), by = list(donor = s$donor), FUN = mean)
    meta <- s[!duplicated(s$donor), ]
    dd <- merge(dd, meta, by = "donor")
    lmf <- lm(as.formula(paste("score ~ disease +", rhs_m)), data = dd)
    sig <- summary(lmf)$sigma
    cf <- summary(lmf)$coefficients
    data.frame(cohort = co$name, region = co$region, module = m, n_genes = length(g),
               n_als = sum(dd$disease == "ALS"), n_control = sum(dd$disease == "Control"),
               d = cf["diseaseALS", 1] / sig, se_d = cf["diseaseALS", 2] / sig,
               age_std = if ("age" %in% rownames(cf)) cf["age", 1] / sig else NA,
               se_age_std = if ("age" %in% rownames(cf)) cf["age", 2] / sig else NA)
  }))
  # Donor-level module scores (one row per donor and module), reused by W60.
  scores <- do.call(rbind, lapply(names(ids), function(m) {
    g <- intersect(ids[[m]], rownames(v))
    dd <- aggregate(data.frame(score = colMeans(v$E[g, , drop = FALSE])), by = list(donor = s$donor), FUN = mean)
    cbind(module = m, merge(dd, s[!duplicated(s$donor), ], by = "donor"))
  }))
  list(name = co$name, region = co$region, disease = disease_tab, age = age_tab, module = mod, scores = scores)
}
