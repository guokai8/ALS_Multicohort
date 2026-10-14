#!/usr/bin/env Rscript

# W72 secondary (Amendment 14): within-lineage ALS vs control log2 fold change of the frozen ALS-up
# genes in GSE287257 lineage pseudobulk (limma-voom ~ disease). Lineages need >= 3 donors per group
# with >= 20 nuclei. With 4 controls the results are descriptive.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
out <- file.path(p01_project_root(), "results", "W72_zelic_disease_cell_source")
res <- file.path(p01_project_root(), "results")
up <- read.delim(file.path(res, "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists",
                           "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_up.tsv"))$common_name

counts <- as.matrix(read.delim(gzfile(file.path(out, "inputs", "w72_lineage_pseudobulk_counts.tsv.gz")), row.names = 1, check.names = FALSE))
units <- read.delim(file.path(out, "inputs", "w72_lineage_units.tsv"), stringsAsFactors = FALSE)
units <- units[units$n_nuclei >= 20, ]

rows <- list()
for (l in unique(units$lineage)) {
  u <- units[units$lineage == l, ]
  if (min(table(factor(u$disease, levels = c("ALS", "Control")))) < 3) next
  meta <- data.frame(disease = factor(u$disease, levels = c("Control", "ALS")))
  design <- model.matrix(~ disease, data = meta)
  y <- DGEList(counts = counts[, u$unit])
  y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
  fit <- eBayes(lmFit(voom(y, design), design), trend = TRUE)
  tab <- topTable(fit, coef = "diseaseALS", number = Inf, sort.by = "none")
  g <- intersect(up, rownames(tab))
  rows[[l]] <- data.frame(lineage = l, n_ALS = sum(meta$disease == "ALS"), n_control = sum(meta$disease == "Control"),
                          gene = g, logFC = tab[g, "logFC"], P = tab[g, "P.Value"])
}
de <- do.call(rbind, rows)
write.table(de, file.path(out, "w72_within_lineage_de.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
focus <- de[de$gene %in% c("GPNMB", "C3", "ANXA1", "LTF", "EMP3"), ]
focus$v <- sprintf("%.2f (%.2g)", focus$logFC, focus$P)
print(reshape(focus[, c("gene", "lineage", "v")], idvar = "gene", timevar = "lineage", direction = "wide"), row.names = FALSE)
print(unique(de[, c("lineage", "n_ALS", "n_control")]), row.names = FALSE)
