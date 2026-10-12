#!/usr/bin/env Rscript

# W54 post hoc: C9-ALS astrocytes in Li 2023 have lower per-nucleus depth
# (median UMI about 5,500 vs about 10,000 in controls), noted before the
# primary model was run but not prespecified. This refits the primary design
# with log10 astrocyte median UMI added.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
out <- file.path(p01_project_root(), "results", "W54_Li2023_motor_cortex_validation")
seed_dir <- file.path(p01_project_root(), "results", "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")

counts <- as.matrix(read.delim(gzfile(file.path(out, "w54_li2023_astrocyte_pseudobulk_counts.tsv.gz")), row.names = 1, check.names = FALSE))
s <- read.delim(file.path(out, "w54_li2023_astrocyte_pseudobulk_samples.tsv"), stringsAsFactors = FALSE)
s <- s[s$disease %in% c("C9-ALS", "Control"), ]
s$disease <- factor(ifelse(s$disease == "C9-ALS", "Case", "Control"), levels = c("Control", "Case"))
s$sex <- factor(s$sex)
s$log10_depth <- log10(s$astro_median_umi)
seeds <- lapply(c(up = "up", down = "down"), function(d) read.delim(file.path(seed_dir, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", d, ".tsv")))$gene_id)

rows <- lapply(c("~ disease + age + sex", "~ disease + age + sex + log10_depth"), function(f) {
  design <- model.matrix(as.formula(f), data = s)
  y <- DGEList(counts = counts[, s$GSM])
  y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(y, design)
  index <- lapply(seeds, intersect, rownames(v))
  fr <- fry(v, index = index, design = design, contrast = "diseaseCase", sort = "none")
  ca <- camera(v, index = index, design = design, contrast = "diseaseCase", sort = FALSE)
  data.frame(formula = f, module = names(index), fry_direction = fr$Direction, fry_p = fr$PValue,
             camera_direction = ca$Direction, camera_p = ca$PValue)
})
res <- do.call(rbind, rows)
write.table(res, file.path(out, "w54_post_hoc_depth_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Depth by group (median UMI):\n"); print(tapply(s$astro_median_umi, s$disease, median))
print(res, digits = 3, row.names = FALSE)
