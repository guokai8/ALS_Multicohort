#!/usr/bin/env Rscript

# W68 A/B (Amendment 11): replication of the W67 motor-territory selectivity in GSE288365 Visium
# CytAssist cervical spinal cord. limma-voom ~ territory * disease (donor-blocked); interaction
# contrast tested with fry (primary) and camera. Sensitivity analyses as fixed in Amendment 11.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(Matrix)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
inp <- file.path(root, "results", "W68_GSE288365_replication", "inputs")
out <- file.path(root, "results", "W68_GSE288365_replication")

counts <- as.matrix(readMM(gzfile(file.path(inp, "w68_pseudobulk_counts.mtx.gz"))))
rownames(counts) <- readLines(gzfile(file.path(inp, "w68_genes.txt.gz")))
colnames(counts) <- readLines(file.path(inp, "w68_unit_order.txt"))
units <- read.delim(file.path(inp, "w68_units.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(units$unit, colnames(counts)))

seed <- file.path(root, "results", "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
sets <- lapply(c(up = "up", down = "down"), function(m)
  read.delim(file.path(seed, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", m, ".tsv")))$common_name)
sig <- read.delim(file.path(root, "results", "W59_NYGC_consensus_signature_holdout", "w59_frozen_consensus_signature.tsv"))
sets$up_consensus <- sig$symbol[sig$direction == "up"]
sets$down_consensus <- sig$symbol[sig$direction == "down"]

EXCL_PRIMARY <- c("NEUVM032BF2", "NEUWD570BTK")          # discovery donor; donor also in W67
W55_DONORS <- c("NEUHD481VCL", "NEUEP010ZFE", "NEUBB078VEG")

comparisons <- list(
  grey = list(unit_type = "territory", motor = "anterior_horn", sensory = "other_grey"),
  white = list(unit_type = "territory", motor = "lateral_ventral_white", sensory = "dorsal_columns"),
  mn_proximity = list(unit_type = "mn_proximity", motor = "anterior_horn_MN_adjacent", sensory = "anterior_horn_MN_distal")
)

fit_one <- function(cmp_name, analysis, keep_donor, rhs = "") {
  cmp <- comparisons[[cmp_name]]
  meta <- units[units$unit_type == cmp$unit_type & units$stratum %in% c(cmp$motor, cmp$sensory) & keep_donor(units) & units$lib_size > 0, ]
  if (nzchar(rhs)) meta <- meta[complete.cases(meta[, c("age", "sex")]), ]
  meta <- meta[meta$donor %in% names(which(table(meta$donor) == 2)), ]
  meta$disease <- factor(ifelse(meta$diagnosis == "Control", "Control", "ALS"), levels = c("Control", "ALS"))
  meta$territory <- factor(ifelse(meta$stratum == cmp$motor, "motor", "sensory"), levels = c("sensory", "motor"))
  meta$sex <- factor(meta$sex)
  design <- model.matrix(as.formula(paste("~ territory * disease", rhs)), data = meta)
  colnames(design) <- make.names(colnames(design))
  y <- DGEList(counts = counts[, meta$unit, drop = FALSE])
  y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(y, design)
  cor1 <- duplicateCorrelation(v, design, block = meta$donor)$consensus.correlation
  v <- voom(y, design, block = meta$donor, correlation = cor1)
  cons <- duplicateCorrelation(v, design, block = meta$donor)$consensus.correlation
  fit <- lmFit(v, design, block = meta$donor, correlation = cons)
  contr <- makeContrasts(interaction = territorymotor.diseaseALS, ALS_in_sensory = diseaseALS,
                         ALS_in_motor = diseaseALS + territorymotor.diseaseALS, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contr), trend = TRUE)
  index <- lapply(sets, intersect, rownames(v))
  res <- do.call(rbind, lapply(colnames(contr), function(cn) {
    fr <- fry(v, index = index, design = design, contrast = contr[, cn], block = meta$donor, correlation = cons, sort = "none")
    ca <- camera(v, index = index, design = design, contrast = contr[, cn], sort = FALSE)
    tab <- topTable(fit2, coef = cn, number = Inf, sort.by = "none")
    do.call(rbind, lapply(names(index), function(m) {
      lfc <- tab[index[[m]], "logFC"]
      data.frame(analysis = analysis, comparison = cmp_name, contrast = cn, gene_set = m,
                 n_als = length(unique(meta$donor[meta$disease == "ALS"])), n_control = length(unique(meta$donor[meta$disease == "Control"])),
                 consensus_correlation = cons, n_genes_tested = length(index[[m]]), fry_direction = fr[m, "Direction"],
                 fry_p = fr[m, "PValue"], camera_direction = ca[m, "Direction"], camera_p = ca[m, "PValue"],
                 median_logFC = median(lfc), n_genes_up = sum(lfc > 0))
    }))
  }))
  genes <- do.call(rbind, lapply(colnames(contr), function(cn) {
    tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "none")
    ids <- c(index$up, index$down)
    data.frame(analysis = analysis, comparison = cmp_name, contrast = cn, symbol = ids,
               gene_set = rep(c("up", "down"), c(length(index$up), length(index$down))), tt[ids, c("logFC", "AveExpr", "t", "P.Value")])
  }))
  list(sets = res, genes = genes)
}

not_in <- function(ids) function(u) !u$donor %in% ids
runs <- list(
  fit_one("grey", "primary", not_in(EXCL_PRIMARY)),
  fit_one("white", "primary", not_in(EXCL_PRIMARY)),
  fit_one("grey", "sens_excl_W55_donors", not_in(c(EXCL_PRIMARY, W55_DONORS))),
  fit_one("grey", "sens_all_donors", function(u) rep(TRUE, nrow(u))),
  fit_one("grey", "sens_sporadic_only", function(u) !u$donor %in% EXCL_PRIMARY & u$diagnosis != "ALS_C9"),
  fit_one("grey", "sens_sex_age", not_in(EXCL_PRIMARY), "+ sex + age"),
  fit_one("mn_proximity", "secondary_mn_proximity", not_in(EXCL_PRIMARY))
)
tab <- do.call(rbind, lapply(runs, `[[`, "sets"))
gene_tab <- do.call(rbind, lapply(runs, `[[`, "genes"))
fam <- tab$analysis == "primary" & tab$contrast == "interaction" & tab$gene_set %in% c("up", "up_consensus")
tab$fry_BH_primary_family <- NA_real_
tab$fry_BH_primary_family[fam] <- p.adjust(tab$fry_p[fam], method = "BH")
write.table(tab, file.path(out, "w68_replication_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_tab, file.path(out, "w68_module_gene_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(tab[tab$gene_set %in% c("up", "up_consensus"), c("analysis", "comparison", "contrast", "gene_set", "n_als", "n_control",
      "fry_direction", "fry_p", "camera_p", "median_logFC", "n_genes_up", "n_genes_tested", "fry_BH_primary_family")], digits = 3, row.names = FALSE)
