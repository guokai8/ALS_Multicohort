#!/usr/bin/env Rscript

# W67 (Amendment 10 B): is the ALS change of the programme larger in motor than in sensory
# territory of the same sections? Donor-level pseudobulk per territory (GSE330224 Visium),
# independent donors, limma-voom ~ territory * disease + sex + age blocked on donor; interaction
# contrast territory(motor):diseaseALS tested with fry (primary) and camera.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
inp <- file.path(root, "results", "W67_motor_selectivity", "inputs")
out <- file.path(root, "results", "W67_motor_selectivity")

seed <- file.path(root, "results", "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
sets <- lapply(c(up = "up", down = "down"), function(m)
  read.delim(file.path(seed, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", m, ".tsv")))$gene_id)
sig <- read.delim(file.path(root, "results", "W59_NYGC_consensus_signature_holdout", "w59_frozen_consensus_signature.tsv"))
sets$up_consensus <- sig$gene[sig$direction == "up"]
sets$down_consensus <- sig$gene[sig$direction == "down"]

comparisons <- list(
  list(tissue = "LSC", name = "grey_ventral_vs_dorsal_horn", motor = "ventral_horn", sensory = "dorsal_horn", family = "primary"),
  list(tissue = "LSC", name = "white_lateral_funiculus_vs_dorsal_columns", motor = "lateral_funiculus", sensory = "dorsal_columns", family = "primary"),
  list(tissue = "MTC", name = "cortex_layer5_vs_layer2_3", motor = "layer5", sensory = "layer2_3", family = "secondary")
)

run_one <- function(cmp) {
  meta <- read.delim(file.path(inp, paste0(cmp$tissue, "_territory_units.tsv")), stringsAsFactors = FALSE)
  counts <- as.matrix(read.delim(gzfile(file.path(inp, paste0(cmp$tissue, "_territory_counts.tsv.gz"))), row.names = 1, check.names = FALSE))
  meta$gse330130_donor <- meta$gse330130_donor %in% c(TRUE, "True", "TRUE")
  meta <- meta[meta$group %in% c("ALS", "Control") & !meta$gse330130_donor & meta$lib_size > 0 &
                 meta$territory %in% c(cmp$motor, cmp$sensory) & !is.na(meta$age), ]
  meta <- meta[meta$donor %in% names(which(table(meta$donor) == 2)), ]  # donors with both territories
  meta$disease <- factor(ifelse(meta$group == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
  meta$territory <- factor(ifelse(meta$territory == cmp$motor, "motor", "sensory"), levels = c("sensory", "motor"))
  meta$sex <- factor(meta$sex)
  design <- model.matrix(~ territory * disease + sex + age, data = meta)
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
    ca <- camera(v, index = index, design = design, contrast = contr[, cn], inter.gene.cor = 0.01, sort = FALSE)
    tab <- topTable(fit2, coef = cn, number = Inf, sort.by = "none")
    do.call(rbind, lapply(names(index), function(m) {
      lfc <- tab[index[[m]], "logFC"]
      data.frame(comparison = cmp$name, tissue = cmp$tissue, family = cmp$family, contrast = cn, gene_set = m,
                 n_als = length(unique(meta$donor[meta$disease == "ALS"])),
                 n_control = length(unique(meta$donor[meta$disease == "Control"])), consensus_correlation = cons,
                 n_genes_tested = length(index[[m]]), fry_direction = fr[m, "Direction"], fry_p = fr[m, "PValue"],
                 camera_direction = ca[m, "Direction"], camera_p = ca[m, "PValue"], median_logFC = median(lfc),
                 n_genes_up = sum(lfc > 0))
    }))
  }))
  genes <- do.call(rbind, lapply(colnames(contr), function(cn) {
    tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "none")
    ids <- unlist(index[c("up", "down")])
    data.frame(comparison = cmp$name, contrast = cn, gene_id = ids, gene_set = rep(c("up", "down"), lengths(index[c("up", "down")])),
               tt[ids, c("logFC", "AveExpr", "t", "P.Value")])
  }))
  list(sets = res, genes = genes)
}

runs <- lapply(comparisons, run_one)
tab <- do.call(rbind, lapply(runs, `[[`, "sets"))
gene_tab <- do.call(rbind, lapply(runs, `[[`, "genes"))
up_sym <- read.delim(file.path(seed, "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_up.tsv"))
dn_sym <- read.delim(file.path(seed, "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_down.tsv"))
gene_tab$symbol <- c(setNames(up_sym$common_name, up_sym$gene_id), setNames(dn_sym$common_name, dn_sym$gene_id))[gene_tab$gene_id]
write.table(gene_tab, file.path(out, "w67_module_gene_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
fam <- tab$family == "primary" & tab$contrast == "interaction" & tab$gene_set %in% c("up", "up_consensus")
tab$fry_BH_primary_family <- NA_real_
tab$fry_BH_primary_family[fam] <- p.adjust(tab$fry_p[fam], method = "BH")
write.table(tab, file.path(out, "w67_motor_selectivity_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(tab[, c("comparison", "contrast", "gene_set", "n_als", "n_control", "fry_direction", "fry_p", "camera_p",
              "median_logFC", "n_genes_up", "n_genes_tested", "fry_BH_primary_family")], digits = 3, row.names = FALSE)
