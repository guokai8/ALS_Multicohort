# W51 事后诊断：区分 PMI 校正效应 与 样本组成效应（丢失 7 个缺 PMI 的对照）
suppressPackageStartupMessages({library(edgeR); library(limma)})
setwd(commandArgs(TRUE)[1])
counts <- read.delim(gzfile("w51_pineda_astrocyte_pseudobulk_counts.tsv.gz"), row.names = 1, check.names = FALSE)
s <- read.delim("w51_pineda_astrocyte_pseudobulk_samples.tsv", stringsAsFactors = FALSE)
seed <- "../W34_synapse_resolved_regulatory_reanalysis/inputs/rna_seed_lists/GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_"
up <- read.delim(paste0(seed, "up.tsv"))$gene_id; down <- read.delim(paste0(seed, "down.tsv"))$gene_id
s <- s[s$Region == "BA4" & s$Condition %in% c("ALS", "PN") & s$lib_size > 0, ]
s$disease <- factor(ifelse(s$Condition == "ALS", "ALS", "Control"), levels = c("Control", "ALS")); s$sex <- factor(s$Sex)
run <- function(id, d, f) {
  X <- model.matrix(as.formula(f), d); y <- DGEList(as.matrix(counts[, d$sample_folder]))
  y <- calcNormFactors(y[filterByExpr(y, design = X), , keep.lib.sizes = FALSE]); v <- voom(y, X)
  idx <- list(up = intersect(up, rownames(v)), down = intersect(down, rownames(v)))
  fr <- fry(v, idx, X, contrast = "diseaseALS", sort = "none"); ca <- camera(v, idx, X, contrast = "diseaseALS", sort = FALSE)
  data.frame(analysis = id, formula = f, n_als = sum(d$disease == "ALS"), n_ctrl = sum(d$disease == "Control"),
             module = names(idx), fry_dir = fr$Direction, fry_p = signif(fr$PValue, 3), camera_dir = ca$Direction, camera_p = signif(ca$PValue, 3))
}
cc <- s[!is.na(s$PMI) & !is.na(s$age), ]
res <- rbind(
  run("posthoc_same30v9_noPMI", cc, "~ disease + age + sex"),
  run("posthoc_same30v9_withPMI", cc, "~ disease + age + sex + PMI"),
  run("posthoc_all_sex_only", s, "~ disease + sex")
)
print(res, row.names = FALSE)
write.table(res, "w51_post_hoc_pmi_vs_composition.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nPMI by group (h):\n"); print(tapply(s$PMI, s$disease, function(x) round(c(n = sum(!is.na(x)), median = median(x, na.rm = TRUE)), 1)))
cat("PMI missing by group:\n"); print(table(s$disease, is.na(s$PMI)))
cat("Wilcoxon PMI ALS vs Control P =", signif(wilcox.test(PMI ~ disease, s)$p.value, 3), "\n")
miss <- s[is.na(s$PMI), c("sample_folder", "Donor", "Condition", "Sex", "age", "n_cells", "lib_size")]; print(miss, row.names = FALSE)
