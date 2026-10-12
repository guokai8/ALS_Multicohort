#!/usr/bin/env Rscript

# W57 (Amendment 3, exploratory): pseudobulk-level SoupX-style contamination
# correction of astrocyte-compartment module effects.
# rho_L = sum_M c_L / (N_L * sum_M a_L) over astrocyte-non-expressed markers M;
# corrected c* = max(0, round(c - k * rho_L * N_L * a_L)).
# Variants: (a) uncorrected, (b) k = 1, (c) k = 2, (d) microglia-only rho
# (discovery), (e) uncorrected + log(rho) covariate.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
out <- file.path(root, "results", "W57_ambient_contamination_test")
seed_dir <- file.path(root, "results", "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
seeds <- lapply(c(up = "up", down = "down"), function(d) read.delim(file.path(seed_dir, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", d, ".tsv"))))
marker_all <- c("P2RY12", "CX3CR1", "CSF1R", "CLDN5", "VWF", "MOBP", "MOG", "SLC17A7", "GAD1", "NRGN")
marker_micro <- c("P2RY12", "CX3CR1", "CSF1R")
stopifnot(!any(marker_all %in% unlist(lapply(seeds, `[[`, "common_name"))))

read_counts <- function(file) as.matrix(read.delim(gzfile(file.path(out, file)), row.names = 1, check.names = FALSE))

estimate_rho <- function(astro, ambient, marker_rows) {
  prop <- sweep(ambient, 2, colSums(ambient), "/")
  colSums(astro[marker_rows, , drop = FALSE]) / (colSums(astro) * colSums(prop[marker_rows, , drop = FALSE]))
}

correct_counts <- function(astro, ambient, rho, k = 1) {
  prop <- sweep(ambient, 2, colSums(ambient), "/")
  expected <- sweep(prop, 2, k * rho * colSums(astro), "*")
  pmax(round(astro - expected), 0)  # matrix first so pmax keeps dimensions
}

run_model <- function(counts, s, rhs, gene_index, label, block = NULL) {
  design <- model.matrix(as.formula(paste("~ disease +", rhs)), data = s)
  y <- DGEList(counts = counts[, s$unit, drop = FALSE])
  y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
  index <- lapply(gene_index, intersect, rownames(y))
  if (is.null(block)) {
    v <- voom(y, design)
    fit <- eBayes(lmFit(v, design), trend = TRUE)
    fr <- fry(v, index = index, design = design, contrast = "diseaseALS", sort = "none")
    ca <- camera(v, index = index, design = design, contrast = "diseaseALS", sort = FALSE)
  } else {
    v <- voom(y, design)
    cons <- duplicateCorrelation(v, design, block = block)$consensus.correlation
    v <- voom(y, design, block = block, correlation = cons)
    cons <- duplicateCorrelation(v, design, block = block)$consensus.correlation
    fit <- eBayes(lmFit(v, design, block = block, correlation = cons), trend = TRUE)
    fr <- fry(v, index = index, design = design, contrast = "diseaseALS", block = block, correlation = cons, sort = "none")
    ca <- NULL
  }
  tab <- topTable(fit, coef = "diseaseALS", number = Inf, sort.by = "none")
  sets <- do.call(rbind, lapply(names(index), function(m) {
    lfc <- tab[index[[m]], "logFC"]
    data.frame(analysis = label, module = m, n_genes_tested = length(index[[m]]),
               fry_direction = fr[m, "Direction"], fry_p = fr[m, "PValue"],
               camera_direction = if (is.null(ca)) NA else ca[m, "Direction"], camera_p = if (is.null(ca)) NA else ca[m, "PValue"],
               median_logFC = median(lfc), n_genes_expected_direction = sum(sign(lfc) == if (m == "up") 1 else -1))
  }))
  genes <- data.frame(analysis = label, gene = unlist(index), module = rep(names(index), lengths(index)), tab[unlist(index), c("logFC", "P.Value")])
  list(sets = sets, genes = genes, tab = tab)
}

variants <- function(astro, ambient, s, rhs, gene_index, cohort, marker_rows, micro_rows = NULL, block = NULL, aggregate_fn = identity) {
  rho <- estimate_rho(astro, ambient, marker_rows)
  runs <- list(
    a_uncorrected = astro,
    b_corrected_rho = correct_counts(astro, ambient, rho, 1),
    c_corrected_2rho = correct_counts(astro, ambient, rho, 2)
  )
  if (!is.null(micro_rows)) runs$d_corrected_micro_rho <- correct_counts(astro, ambient, estimate_rho(astro, ambient, micro_rows), 1)
  res <- lapply(names(runs), function(v) run_model(aggregate_fn(runs[[v]]), s, rhs, gene_index, paste(cohort, v, sep = ":"), block))
  names(res) <- names(runs)
  list(rho = rho, res = res)
}

collect <- function(x) list(sets = do.call(rbind, lapply(x, `[[`, "sets")), genes = do.call(rbind, lapply(x, `[[`, "genes")))

# Deviation from Amendment 3 (documented in README): log(rho) is undefined
# when rho = 0 (one Ulm unit), so every cohort uses log(rho + min positive rho / 2).
safe_log_rho <- function(rho) log(rho + min(rho[rho > 0]) / 2)

# ---- E4 discovery (GSE330130 MTC) ----
d_astro <- read_counts("w57_discovery_astro_counts.tsv.gz")
d_amb <- read_counts("w57_discovery_ambient_counts.tsv.gz")
d_ann <- read.delim(file.path(out, "w57_discovery_gene_annotation.tsv"), stringsAsFactors = FALSE)
audit <- read.delim(file.path(root, "results", "W33_synapse_resolved_diagnosis_sensitivity", "MTC_Astrocyte_RNA.library_label_audit.tsv"), stringsAsFactors = FALSE)
s <- audit[audit$used_in_w33 %in% c(TRUE, "TRUE"), ]
s$unit <- s$library
s$disease <- factor(ifelse(s$disease_binary_w33 == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
s$sex <- factor(s$sex)
s$tissue_subregion <- factor(s$tissue_subregion)
d_markers <- d_ann$id[match(marker_all, d_ann$common_name)]
d_micro <- d_ann$id[match(marker_micro, d_ann$common_name)]
stopifnot(length(d_markers) == length(marker_all), !anyNA(d_markers), length(d_micro) == length(marker_micro), !anyNA(d_micro))
d_index <- lapply(seeds, `[[`, "gene_id")
disc <- variants(d_astro[, s$unit], d_amb[, s$unit], s, "age + sex + PMI + tissue_subregion", d_index, "discovery", d_markers, d_micro)
repro <- merge(do.call(rbind, seeds)[, c("gene_id", "logFC")], data.frame(gene_id = rownames(disc$res$a_uncorrected$tab), rebuilt = disc$res$a_uncorrected$tab$logFC), by = "gene_id")
max_diff <- max(abs(repro$logFC - repro$rebuilt))
writeLines(sprintf("discovery reproduction: %d module genes, max |logFC diff| = %.3g", nrow(repro), max_diff), file.path(out, "w57_discovery_reproduction_check.txt"))
if (nrow(repro) != 29 || max_diff > 1e-6) stop("Discovery reproduction check failed: max diff ", max_diff)
s$log_rho <- safe_log_rho(disc$rho[s$unit])
disc$res$e_rho_covariate <- run_model(d_astro, s, "age + sex + PMI + tissue_subregion + log_rho", d_index, "discovery:e_rho_covariate")
rho_tab <- data.frame(cohort = "discovery", unit = s$unit, disease = s$disease, rho = disc$rho[s$unit])

# ---- E5 Ulm ----
u_all <- read_counts("w57_ulm_astro_and_ambient_counts.tsv.gz")
u_astro_cols <- grep("^ambient\\|", colnames(u_all), value = TRUE, invert = TRUE)
u_well <- sub("^.*\\|", "", u_astro_cols)
u_astro <- u_all[, u_astro_cols]
u_amb <- u_all[, paste0("ambient|", u_well)]
colnames(u_amb) <- u_astro_cols
u_s <- read.delim(file.path(root, "results", "W52_Ulm_motor_cortex_validation", "w52_ulm_astrocyte_pseudobulk_samples.tsv"), stringsAsFactors = FALSE)
u_s <- u_s[u_s$Case_Type %in% c("ALS", "HC"), ]
u_s$unit <- u_s$stratum
u_s$disease <- factor(ifelse(u_s$Case_Type == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
u_s$sex <- factor(tolower(u_s$Sex))
u_s$chip <- factor(u_s$Chip)
well_to_chip <- paste(sub("\\|.*$", "", u_astro_cols), sub("^(Chip[0-9]+).*$", "\\1", u_well), sep = "|")
to_chip <- function(m) {
  agg <- t(rowsum(t(m), group = well_to_chip))
  agg[, u_s$unit, drop = FALSE]
}
u_markers <- intersect(marker_all, rownames(u_all))
stopifnot(length(u_markers) == length(marker_all))
ulm <- variants(u_astro, u_amb, u_s, "sex + chip", lapply(seeds, `[[`, "common_name"), "Ulm", u_markers,
                block = u_s$ID, aggregate_fn = to_chip)
u_rho_chip <- tapply(ulm$rho * colSums(u_astro), well_to_chip, sum) / tapply(colSums(u_astro), well_to_chip, sum)
u_s$log_rho <- safe_log_rho(u_rho_chip[u_s$unit])
ulm$res$e_rho_covariate <- run_model(to_chip(u_astro), u_s, "sex + chip + log_rho", lapply(seeds, `[[`, "common_name"), "Ulm:e_rho_covariate", block = u_s$ID)
rho_tab <- rbind(rho_tab, data.frame(cohort = "Ulm", unit = names(u_rho_chip[u_s$unit]), disease = u_s$disease, rho = as.numeric(u_rho_chip[u_s$unit])))

# ---- E5 Pineda BA9 batch 210526 ----
w56 <- file.path(root, "results", "W56_cell_source_check")
p_line <- as.matrix(read.delim(gzfile(file.path(w56, "w56_pineda_lineage_counts.tsv.gz")), row.names = 1, check.names = FALSE))
p_amb <- read_counts("w57_pineda_ba9_210526_ambient_counts.tsv.gz")
p_strata <- read.delim(file.path(w56, "w56_pineda_lineage_strata.tsv"), stringsAsFactors = FALSE)
clin <- read.delim(file.path(root, "results", "W22_Pineda_S01_clinical_pathology_context", "w22_pineda_lead_module_donor_context.tsv"), stringsAsFactors = FALSE)
clin <- clin[!duplicated(clin$Donor), c("Donor", "Age.of.death..Y.")]
names(clin) <- c("Donor", "age")
p_s <- merge(p_strata[p_strata$lineage == "Astrocyte", ], clin, by = "Donor")
p_s <- p_s[p_s$sample_folder %in% colnames(p_amb) & p_s$Condition %in% c("ALS", "PN") & p_s$lib_size > 0 & !is.na(p_s$age), ]
p_s$unit <- p_s$sample_folder
p_s$disease <- factor(ifelse(p_s$Condition == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
p_s$sex <- factor(p_s$Sex)
p_astro <- p_line[, paste0(p_s$unit, "|Astrocyte")]
colnames(p_astro) <- p_s$unit
gene_map <- read.delim(file.path(p01_source_data_root(), "Other_Datasets/Pineda_snRNAseq/synapse_download/Processed Data/191112_ALS_110_snRNA-B9/row_metadata.tsv"), stringsAsFactors = FALSE)
p_markers <- gene_map$ENSEMBL[match(marker_all, gene_map$Gene)]
stopifnot(length(p_markers) == length(marker_all), !anyNA(p_markers))
pin <- variants(p_astro, p_amb[, p_s$unit], p_s, "age + sex", lapply(seeds, `[[`, "gene_id"), "PinedaBA9_210526", p_markers)
p_s$log_rho <- safe_log_rho(pin$rho[p_s$unit])
pin$res$e_rho_covariate <- run_model(p_astro, p_s, "age + sex + log_rho", lapply(seeds, `[[`, "gene_id"), "PinedaBA9_210526:e_rho_covariate")
rho_tab <- rbind(rho_tab, data.frame(cohort = "PinedaBA9_210526", unit = p_s$unit, disease = p_s$disease, rho = pin$rho[p_s$unit]))

# ---- outputs ----
all <- collect(c(disc$res, ulm$res, pin$res))
rho_test <- do.call(rbind, lapply(split(rho_tab, rho_tab$cohort), function(d) {
  data.frame(cohort = d$cohort[1], median_rho_ALS = median(d$rho[d$disease == "ALS"]), median_rho_control = median(d$rho[d$disease == "Control"]),
             wilcox_p = wilcox.test(rho ~ disease, data = d)$p.value)
}))
write.table(all$sets, file.path(out, "w57_module_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all$genes, file.path(out, "w57_module_gene_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rho_tab, file.path(out, "w57_rho_by_unit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rho_test, file.path(out, "w57_rho_disease_test.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(readLines(file.path(out, "w57_discovery_reproduction_check.txt")), sep = "\n")
print(rho_test, digits = 3, row.names = FALSE)
print(all$sets, digits = 3, row.names = FALSE)
