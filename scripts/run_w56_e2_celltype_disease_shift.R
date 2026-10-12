#!/usr/bin/env Rscript

# W56 E2 (Amendment 2): ALS-vs-control module shift by cell type.
# - Ulm: Astrocytes (consistency check against W52), Microglia,
#   Oligodendrocytes, OPC; donor x chip, ~ disease + sex + chip, duplicateCorrelation.
# - Pineda BA9 batch 210526 (batch-balanced): Astrocyte, Microglia, Endothelial;
#   one profile per donor, ~ disease + age + sex.
# Blocked fry where applicable, plus camera (unblocked designs only).

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
out <- file.path(root, "results", "W56_cell_source_check")
seed_dir <- file.path(root, "results", "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
seeds <- lapply(c(up = "up", down = "down"), function(d) read.delim(file.path(seed_dir, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", d, ".tsv"))))

test_sets <- function(counts, s, rhs, gene_index, block = NULL, label) {
  design <- model.matrix(as.formula(paste("~ disease +", rhs)), data = s)
  y <- DGEList(counts = counts[, s$stratum, drop = FALSE])
  y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
  index <- lapply(gene_index, intersect, rownames(y))
  if (is.null(block)) {
    v <- voom(y, design)
    fr <- fry(v, index = index, design = design, contrast = "diseaseALS", sort = "none")
    ca <- camera(v, index = index, design = design, contrast = "diseaseALS", sort = FALSE)
    cons <- NA_real_
  } else {
    v <- voom(y, design)
    cons <- duplicateCorrelation(v, design, block = block)$consensus.correlation
    v <- voom(y, design, block = block, correlation = cons)
    cons <- duplicateCorrelation(v, design, block = block)$consensus.correlation
    fr <- fry(v, index = index, design = design, contrast = "diseaseALS", block = block, correlation = cons, sort = "none")
    ca <- NULL
  }
  fit <- eBayes(lmFit(v, design, block = block, correlation = if (is.null(block)) NULL else cons), trend = TRUE)
  tab <- topTable(fit, coef = "diseaseALS", number = Inf, sort.by = "none")
  do.call(rbind, lapply(names(index), function(m) {
    lfc <- tab[index[[m]], "logFC"]
    data.frame(analysis = label, module = m, n_als = length(unique(s$donor[s$disease == "ALS"])),
               n_control = length(unique(s$donor[s$disease == "Control"])), n_genes_tested = length(index[[m]]),
               fry_direction = fr[m, "Direction"], fry_p = fr[m, "PValue"],
               camera_direction = if (is.null(ca)) NA else ca[m, "Direction"], camera_p = if (is.null(ca)) NA else ca[m, "PValue"],
               n_genes_up = sum(lfc > 0), median_logFC = median(lfc))
  }))
}

# Ulm.
u_counts <- as.matrix(read.delim(gzfile(file.path(out, "w56_ulm_lineage_counts.tsv.gz")), row.names = 1, check.names = FALSE))
u_strata <- read.delim(file.path(out, "w56_ulm_lineage_strata.tsv"), stringsAsFactors = FALSE)
u_strata <- u_strata[u_strata$Case_Type %in% c("ALS", "HC") & u_strata$lib_size > 0, ]
u_strata$disease <- factor(ifelse(u_strata$Case_Type == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
u_strata$sex <- factor(tolower(u_strata$Sex))
u_strata$chip <- factor(u_strata$Chip)
u_strata$donor <- u_strata$ID
u_index <- lapply(seeds, function(x) x$common_name)
ulm <- do.call(rbind, lapply(c("Astrocytes", "Microglia", "Oligodendrocytes", "OPC"), function(lin) {
  s <- u_strata[u_strata$lineage == lin, ]
  s$chip <- droplevels(s$chip)
  test_sets(u_counts, s, "sex + chip", u_index, block = s$donor, label = paste0("Ulm_", lin))
}))

# Pineda BA9, batch 210526.
p_counts <- as.matrix(read.delim(gzfile(file.path(out, "w56_pineda_lineage_counts.tsv.gz")), row.names = 1, check.names = FALSE))
p_strata <- read.delim(file.path(out, "w56_pineda_lineage_strata.tsv"), stringsAsFactors = FALSE)
clin <- read.delim(file.path(root, "results", "W22_Pineda_S01_clinical_pathology_context", "w22_pineda_lead_module_donor_context.tsv"), stringsAsFactors = FALSE)
clin <- clin[!duplicated(clin$Donor), c("Donor", "Age.of.death..Y.")]
names(clin) <- c("Donor", "age")
p_strata <- merge(p_strata, clin, by = "Donor", all.x = TRUE)
p_strata <- p_strata[p_strata$Region == "BA9" & startsWith(p_strata$sample_folder, "210526") &
                       p_strata$Condition %in% c("ALS", "PN") & p_strata$lib_size > 0 & !is.na(p_strata$age), ]
p_strata$disease <- factor(ifelse(p_strata$Condition == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
p_strata$sex <- factor(p_strata$Sex)
p_strata$donor <- p_strata$Donor
p_index <- lapply(seeds, function(x) x$gene_id)
pineda <- do.call(rbind, lapply(c("Astrocyte", "Microglia", "Endothelial"), function(lin) {
  test_sets(p_counts, p_strata[p_strata$lineage == lin, ], "age + sex", p_index, label = paste0("PinedaBA9_210526_", lin))
}))

res <- rbind(ulm, pineda)
write.table(res, file.path(out, "w56_e2_celltype_disease_shift.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, digits = 3, row.names = FALSE)
