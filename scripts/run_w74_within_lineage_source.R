#!/usr/bin/env Rscript

# W74 (Amendment 16): is the ALS-up signal astrocyte-intrinsic or shared with the microglial compartment?
# A: discovery (GSE330130 MTC) microglia compartment, W33 diagnoses and covariate rules (primary).
# B: discovery astrocyte vs microglia, same libraries, ~ lineage * disease + covariates, library-blocked.
# C: Ulm, Pineda BA9 210526 and Zelic, astrocyte vs microglia in the same donors; Stouffer of signed z.
# D: discovery, ALS-up module in every lineage (descriptive).

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
res <- file.path(root, "results")
out <- file.path(res, "W74_within_lineage_source")
dir.create(out, showWarnings = FALSE)
# Per-lineage W01 pseudobulk inputs (export_w01_rna_de_inputs.py), as read by W33.
DISC <- file.path(root, "intermediate", "W01_RNA_DE", "GSE330130_MTC_RNA")

seed <- file.path(res, "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
seeds <- lapply(c(up = "up", down = "down"), function(d) read.delim(file.path(seed, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", d, ".tsv"))))

# ---- W33 diagnosis resolution and covariate rules (copied from run_w33_synapse_resolved_diagnosis_sensitivity.R) ----
bridge <- read.delim(file.path(res, "Supplemental_QC_covariate_audit", "gse330130_nygc_latest_zenodo_bridge_2026-08-12.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
bridge <- bridge[!duplicated(bridge$mo_id), , drop = FALSE]
empty_to_na <- function(x) { x <- as.character(x); x[x %in% c("", "NA", "NaN", "<NA>")] <- NA_character_; x }
resolve_diagnosis <- function(samples) {
  m <- merge(samples, bridge[, c("mo_id", "nygc_latest_direct_rna_id_match", "resolved_disease_group")], by.x = "library", by.y = "mo_id", all.x = TRUE, sort = FALSE)
  m <- m[match(samples$sample_col, m$sample_col), , drop = FALSE]
  direct <- tolower(empty_to_na(m$nygc_latest_direct_rna_id_match)) %in% c("yes", "true", "1")
  resolved <- ifelse(m$resolved_disease_group %in% c("ALS", "Control"), m$resolved_disease_group, NA_character_)
  m$disease_binary_w33 <- m$disease_binary
  m$disease_binary_w33[direct] <- resolved[direct]
  m
}
# Port of W33 build_design(): drop constant covariates, then (for power) covariates with missing values, then for rank.
covariate_rule <- function(s) {
  s$disease_binary_w33 <- factor(s$disease_binary_w33, levels = c("Control", "ALS"))
  s$sex <- factor(empty_to_na(s$sex)); s$tissue_subregion <- factor(empty_to_na(s$tissue_subregion))
  s$age <- suppressWarnings(as.numeric(s$age)); s$PMI <- suppressWarnings(as.numeric(s$PMI))
  used <- Filter(function(cv) { v <- s[[cv]]; sum(!is.na(v)) > 0 && length(unique(v[!is.na(v)])) >= 2 }, c("age", "sex", "PMI", "tissue_subregion"))
  sm <- s[complete.cases(s[, c("disease_binary_w33", used), drop = FALSE]), ]
  repeat {
    n_a <- length(unique(sm$donor[sm$disease_binary_w33 == "ALS"])); n_c <- length(unique(sm$donor[sm$disease_binary_w33 == "Control"]))
    if (n_a >= 5 && n_c >= 5) break
    removable <- intersect(c("tissue_subregion", "PMI", "sex", "age"), used[vapply(used, function(cv) any(is.na(s[[cv]])), logical(1))])
    if (!length(removable)) break
    used <- setdiff(used, removable[1]); sm <- s[complete.cases(s[, c("disease_binary_w33", used), drop = FALSE]), ]
  }
  repeat {
    d <- model.matrix(as.formula(paste("~ disease_binary_w33", if (length(used)) paste("+", paste(used, collapse = " + ")) else "")), data = sm)
    if (qr(d)$rank == ncol(d)) break
    rm1 <- intersect(c("tissue_subregion", "PMI", "sex", "age"), used)[1]; if (is.na(rm1)) break; used <- setdiff(used, rm1)
  }
  list(samples = sm, used = used)
}
read_disc <- function(lineage) {
  s <- read.delim(file.path(DISC, lineage, "samples.tsv"), stringsAsFactors = FALSE)
  tab <- read.delim(gzfile(file.path(DISC, lineage, "counts.tsv.gz")), check.names = FALSE)
  x <- as.matrix(tab[, setdiff(colnames(tab), "gene_id")]); rownames(x) <- tab$gene_id; storage.mode(x) <- "integer"
  s <- s[match(colnames(x), s$sample_col), ]
  s <- resolve_diagnosis(s)
  s <- s[s$disease_binary %in% c("ALS", "Control") & s$disease_binary_w33 %in% c("ALS", "Control"), ]
  list(x = x[, s$sample_col, drop = FALSE], s = s)
}

# ---- generic set test ----
set_test <- function(x, meta, formula, coef, index, block = NULL, label) {
  meta <- droplevels(meta)
  design <- model.matrix(as.formula(formula), data = meta); colnames(design) <- make.names(colnames(design))
  coef <- make.names(coef)
  y <- DGEList(counts = x); y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
  idx <- lapply(index, intersect, rownames(y))
  if (is.null(block)) {
    v <- voom(y, design); cons <- NULL
    fr <- fry(v, index = idx, design = design, contrast = coef, sort = "none")
    ca <- camera(v, index = idx, design = design, contrast = coef, sort = FALSE)
  } else {
    v <- voom(y, design); c1 <- duplicateCorrelation(v, design, block = block)$consensus.correlation
    v <- voom(y, design, block = block, correlation = c1); cons <- duplicateCorrelation(v, design, block = block)$consensus.correlation
    fr <- fry(v, index = idx, design = design, contrast = coef, block = block, correlation = cons, sort = "none")
    ca <- NULL  # camera has no block argument; as in W56 E2, not reported for blocked models
  }
  fit <- eBayes(lmFit(v, design, block = block, correlation = cons), trend = TRUE)
  tab <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  out_tab <- do.call(rbind, lapply(names(idx), function(m) data.frame(
    analysis = label, module = m, contrast = coef, n_units = ncol(x), n_genes_tested = length(idx[[m]]),
    fry_direction = fr[m, "Direction"], fry_p = fr[m, "PValue"], camera_p = if (is.null(ca)) NA else ca[m, "PValue"], block_correlation = if (is.null(cons)) NA else cons,
    median_logFC = if (length(idx[[m]])) median(tab[idx[[m]], "logFC"]) else NA, n_genes_up = sum(tab[idx[[m]], "logFC"] > 0))))
  genes <- do.call(rbind, lapply(names(idx), function(m) if (length(idx[[m]])) data.frame(analysis = label, module = m, gene = idx[[m]],
    logFC = tab[idx[[m]], "logFC"], P.Value = tab[idx[[m]], "P.Value"], AveExpr = tab[idx[[m]], "AveExpr"])))
  attr(out_tab, "genes") <- genes
  out_tab
}
disc_index <- lapply(seeds, `[[`, "gene_id")

# ---- A and D: discovery, every lineage ----
LINEAGES <- c("Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Excitatory_neuron", "Inhibitory_neuron")
disc <- lapply(setNames(LINEAGES, LINEAGES), read_disc)
D_list <- lapply(LINEAGES, function(l) {
  cr <- covariate_rule(disc[[l]]$s)
  f <- paste("~ disease_binary_w33", if (length(cr$used)) paste("+", paste(cr$used, collapse = " + ")) else "")
  r0 <- set_test(disc[[l]]$x[, cr$samples$sample_col], cr$samples, f, "disease_binary_w33ALS", disc_index, label = paste0("discovery_", l)); r <- r0
  r$formula <- f; r$n_als <- length(unique(cr$samples$donor[cr$samples$disease_binary_w33 == "ALS"]))
  r$n_control <- length(unique(cr$samples$donor[cr$samples$disease_binary_w33 == "Control"])); attr(r, "genes") <- attr(r0, "genes"); r
})
D <- do.call(rbind, D_list)
# descriptive, not prespecified: per-gene log2FC of module genes in each discovery lineage
D_genes <- do.call(rbind, lapply(D_list, attr, "genes"))
sym <- do.call(rbind, seeds)[, c("gene_id", "common_name")]
D_genes$symbol <- sym$common_name[match(D_genes$gene, sym$gene_id)]
A <- D[D$analysis == "discovery_Microglia", ]
A_reading <- with(A[A$module == "up", ], ifelse(fry_p < 0.05 & fry_direction == "Up" & median_logFC > 0,
  "ALS-up signal shared with the microglial compartment", "not detected in the microglial compartment"))

# ---- B: discovery astrocyte vs microglia, same libraries ----
sa <- disc$Astrocyte$s; sm <- disc$Microglia$s
common <- intersect(sa$library, sm$library)
sa <- sa[sa$library %in% common, ]; sm <- sm[sm$library %in% common, ]
genes <- intersect(rownames(disc$Astrocyte$x), rownames(disc$Microglia$x))
xb <- cbind(disc$Astrocyte$x[genes, sa$sample_col], disc$Microglia$x[genes, sm$sample_col])
colnames(xb) <- c(paste0(sa$sample_col, "_A"), paste0(sm$sample_col, "_M"))
mb <- rbind(transform(sa, lineage = "Astrocyte"), transform(sm, lineage = "Microglia"))
cr <- covariate_rule(sa)
mb$lineage <- factor(mb$lineage, levels = c("Microglia", "Astrocyte"))
mb$disease_binary_w33 <- factor(mb$disease_binary_w33, levels = c("Control", "ALS"))
mb$sex <- factor(empty_to_na(mb$sex)); mb$tissue_subregion <- factor(empty_to_na(mb$tissue_subregion))
mb$age <- as.numeric(mb$age); mb$PMI <- as.numeric(mb$PMI)
keep <- complete.cases(mb[, c("disease_binary_w33", cr$used)]); mb <- mb[keep, ]; xb <- xb[, keep]
fB <- paste("~ lineage * disease_binary_w33", if (length(cr$used)) paste("+", paste(cr$used, collapse = " + ")) else "")
B <- set_test(xb, mb, fB, "lineageAstrocyte:disease_binary_w33ALS", disc_index, block = mb$library, label = "discovery_astrocyte_minus_microglia")
ast <- D[D$analysis == "discovery_Astrocyte" & D$module == "up", ]; mic <- A[A$module == "up", ]
B$microglia_over_astrocyte_median_logFC <- ifelse(B$module == "up", mic$median_logFC / ast$median_logFC, NA)

# ---- C: independent cohorts ----
`%||%` <- function(a, b) if (is.null(a)) b else a
paired <- function(x, s, rhs, index, block, label, lin_names) {
  s$lineage <- factor(ifelse(s$lineage == lin_names[["astro"]], "Astrocyte", "Microglia"), levels = c("Microglia", "Astrocyte"))
  within <- do.call(rbind, lapply(c("Astrocyte", "Microglia"), function(l) {
    ss <- s[s$lineage == l, ]
    set_test(x[, ss$unit], ss, paste("~ disease", rhs), "diseaseALS", index, block = if (is.null(block)) NULL else ss[[block]], label = paste0(label, "_", l))
  }))
  inter <- set_test(x[, s$unit], s, paste("~ lineage * disease", rhs), "lineageAstrocyte:diseaseALS", index, block = s[[block %||% "donor"]],
                    label = paste0(label, "_astrocyte_minus_microglia"))
  rbind(within, inter)
}
# Ulm (W56)
uc <- as.matrix(read.delim(gzfile(file.path(res, "W56_cell_source_check", "w56_ulm_lineage_counts.tsv.gz")), row.names = 1, check.names = FALSE))
us <- read.delim(file.path(res, "W56_cell_source_check", "w56_ulm_lineage_strata.tsv"), stringsAsFactors = FALSE)
us <- us[us$Case_Type %in% c("ALS", "HC") & us$lib_size > 0 & us$lineage %in% c("Astrocytes", "Microglia"), ]
us$disease <- factor(ifelse(us$Case_Type == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
us$sex <- factor(tolower(us$Sex)); us$chip <- factor(us$Chip); us$donor <- us$ID; us$unit <- us$stratum
C_ulm <- paired(uc, us, "+ sex + chip", lapply(seeds, `[[`, "common_name"), "donor", "Ulm", c(astro = "Astrocytes"))
# Pineda BA9 210526 (W56)
pc <- as.matrix(read.delim(gzfile(file.path(res, "W56_cell_source_check", "w56_pineda_lineage_counts.tsv.gz")), row.names = 1, check.names = FALSE))
ps <- read.delim(file.path(res, "W56_cell_source_check", "w56_pineda_lineage_strata.tsv"), stringsAsFactors = FALSE)
clin <- read.delim(file.path(res, "W22_Pineda_S01_clinical_pathology_context", "w22_pineda_lead_module_donor_context.tsv"), stringsAsFactors = FALSE)
clin <- clin[!duplicated(clin$Donor), c("Donor", "Age.of.death..Y.")]; names(clin) <- c("Donor", "age")
ps <- merge(ps, clin, by = "Donor", all.x = TRUE)
ps <- ps[ps$Region == "BA9" & startsWith(ps$sample_folder, "210526") & ps$Condition %in% c("ALS", "PN") & ps$lib_size > 0 & !is.na(ps$age) &
           ps$lineage %in% c("Astrocyte", "Microglia"), ]
ps$disease <- factor(ifelse(ps$Condition == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
ps$sex <- factor(ps$Sex); ps$donor <- as.character(ps$Donor); ps$unit <- ps$stratum
C_pin <- paired(pc, ps, "+ age + sex", lapply(seeds, `[[`, "gene_id"), NULL, "PinedaBA9_210526", c(astro = "Astrocyte"))
# Zelic (W72), units with >= 20 nuclei
zc <- as.matrix(read.delim(gzfile(file.path(res, "W72_zelic_disease_cell_source", "inputs", "w72_lineage_pseudobulk_counts.tsv.gz")), row.names = 1, check.names = FALSE))
zs <- read.delim(file.path(res, "W72_zelic_disease_cell_source", "inputs", "w72_lineage_units.tsv"), stringsAsFactors = FALSE)
zcov <- read.delim(file.path(res, "W58_multicohort_meta", "w58_zelic_astrocyte_pseudobulk_samples.tsv"), stringsAsFactors = FALSE)[, c("donor_id", "age", "sex")]
zs <- merge(zs, zcov, by.x = "donor", by.y = "donor_id")
zs <- zs[zs$lineage %in% c("Astrocyte", "Microglia") & zs$n_nuclei >= 20, ]
zs <- zs[zs$donor %in% names(which(table(zs$donor) == 2)), ]
zs$disease <- factor(zs$disease, levels = c("Control", "ALS")); zs$sex <- factor(zs$sex); zs$donor <- as.character(zs$donor)
C_zel <- paired(zc, zs, "+ age + sex", lapply(seeds, `[[`, "common_name"), NULL, "Zelic", c(astro = "Astrocyte"))
C <- rbind(C_ulm, C_pin, C_zel)

inter <- C[grepl("astrocyte_minus_microglia", C$analysis) & C$module == "up", ]
inter$signed_z <- ifelse(inter$fry_direction == "Up", 1, -1) * qnorm(1 - inter$fry_p / 2)
inter$reading <- ifelse(inter$fry_p < 0.05 & inter$fry_direction == "Up", "larger in astrocytes",
                 ifelse(inter$fry_p < 0.05 & inter$fry_direction == "Down", "larger in microglia", "not distinguishable"))
Z <- sum(inter$signed_z) / sqrt(nrow(inter))
C_comb <- data.frame(cohorts = paste(inter$analysis, collapse = "; "), stouffer_z = Z, two_sided_p = 2 * pnorm(-abs(Z)))

summary <- data.frame(part = c("A", "B", "C"), reading = c(A_reading,
  sprintf("interaction (astrocyte minus microglia) fry P = %.3g (%s); microglia/astrocyte median logFC ratio = %.2f",
          B$fry_p[B$module == "up"], B$fry_direction[B$module == "up"], B$microglia_over_astrocyte_median_logFC[B$module == "up"]),
  sprintf("combined signed z = %.2f, two-sided P = %.3g; per cohort: %s", Z, C_comb$two_sided_p, paste(inter$analysis, inter$reading, sep = " = ", collapse = "; "))))

write.table(D, file.path(out, "w74_discovery_lineage_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(D_genes, file.path(out, "w74_discovery_lineage_module_genes_descriptive.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(B, file.path(out, "w74_discovery_astrocyte_vs_microglia.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(C, file.path(out, "w74_independent_cohort_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(inter, file.path(out, "w74_independent_cohort_interaction_reading.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out, "w74_reading.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(D[, c("analysis", "module", "n_als", "n_control", "n_genes_tested", "fry_direction", "fry_p", "camera_p", "median_logFC", "n_genes_up")], digits = 3, row.names = FALSE)
print(B, digits = 3, row.names = FALSE)
print(C[, c("analysis", "module", "n_units", "n_genes_tested", "fry_direction", "fry_p", "median_logFC")], digits = 3, row.names = FALSE)
print(inter[, c("analysis", "fry_direction", "fry_p", "signed_z", "reading")], digits = 3, row.names = FALSE)
print(summary, row.names = FALSE)
