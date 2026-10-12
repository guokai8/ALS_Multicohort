#!/usr/bin/env Rscript

# W62 discovery confounder audit and tissue-programme sensitivity (Amendment 7).
# Balance table by diagnosis, a CGND-only model (primary sensitivity), a
# programme-adjusted model (secondary) and a control-source comparison of
# per-library module scores. Uses the W57 rebuilt discovery astrocyte pseudobulk.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(readxl)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
res_dir <- file.path(root, "results")
out <- file.path(res_dir, "W62_discovery_confounder_audit")
dir.create(out, showWarnings = FALSE)

# ---- inputs ----
seed_dir <- file.path(res_dir, "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
seeds <- do.call(rbind, lapply(c("up", "down"), function(m) {
  x <- read.delim(file.path(seed_dir, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", m, ".tsv")))
  data.frame(module = m, gene_id = x$gene_id, symbol = x$common_name, w34_logFC = x$logFC)
}))
audit <- read.delim(file.path(res_dir, "W33_synapse_resolved_diagnosis_sensitivity", "MTC_Astrocyte_RNA.library_label_audit.tsv"), stringsAsFactors = FALSE)
s <- audit[audit$used_in_w33 %in% c(TRUE, "TRUE"), ]
bridge <- read.delim(file.path(root, "metadata", "gse330130_geo_sample_subject_bridge_2026-08-12.tsv"), stringsAsFactors = FALSE)
bridge <- bridge[bridge$geo_modality == "GEX", c("library_name", "geo_library_id")]
s <- merge(s, bridge, by.x = "library", by.y = "library_name", all.x = TRUE)
s$programme <- ifelse(grepl("-U54-", s$geo_library_id), "U54", ifelse(grepl("-CGND-", s$geo_library_id), "CGND", NA))
if (anyNA(s$programme)) stop("library without a tissue programme: ", paste(s$library[is.na(s$programme)], collapse = ", "))
st1 <- read_excel(file.path(p01_source_data_root(), "Other_Datasets", "NYGC_ALS_Consortium", "files", "NYGC_Flagship_ST_1_3_Metadata.xlsx"),
                  sheet = "ST1", skip = 1)
s$site <- st1$site_sample_collected[match(s$donor, st1$donor_id)]
s$site[is.na(s$site)] <- "not listed in NYGC ST1"
s$disease <- factor(ifelse(s$disease_binary_w33 == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
for (v in c("sex", "tissue_subregion", "programme")) s[[v]] <- factor(s[[v]])
counts <- as.matrix(read.delim(gzfile(file.path(res_dir, "W57_ambient_contamination_test", "w57_discovery_astro_counts.tsv.gz")),
                               row.names = 1, check.names = FALSE))[, s$library]
s$lib_size <- colSums(counts)

# ---- 1. balance table ----
cont <- c("age", "PMI", "n_cells", "lib_size")
categ <- c("sex", "tissue_subregion", "programme")
bal <- rbind(
  do.call(rbind, lapply(cont, function(v) {
    x <- as.numeric(s[[v]])
    data.frame(variable = v, level = "median", ALS = median(x[s$disease == "ALS"]), Control = median(x[s$disease == "Control"]),
               test = "Wilcoxon", p = wilcox.test(x ~ s$disease, exact = FALSE)$p.value)
  })),
  do.call(rbind, lapply(categ, function(v) {
    tb <- table(s[[v]], s$disease)
    data.frame(variable = v, level = rownames(tb), ALS = tb[, "ALS"], Control = tb[, "Control"], test = "Fisher",
               p = fisher.test(tb)$p.value)
  }))
)
site_tab <- as.data.frame.matrix(table(s$site, s$disease))
site_tab <- cbind(site = rownames(site_tab), site_tab)

# ---- model helper (W34 pipeline) ----
fit_model <- function(keep, rhs, label) {
  ss <- droplevels(s[keep, ])
  design <- model.matrix(as.formula(paste("~ disease +", rhs)), data = ss)
  y <- DGEList(counts = counts[, ss$library, drop = FALSE])
  y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(y, design)
  fit <- eBayes(lmFit(v, design), trend = TRUE)
  tt <- topTable(fit, coef = "diseaseALS", number = Inf, sort.by = "none")
  g <- merge(seeds, data.frame(gene_id = rownames(tt), logFC = tt$logFC, P = tt$P.Value), by = "gene_id", all.x = TRUE)
  g$direction_kept <- sign(g$logFC) == sign(g$w34_logFC)
  idx <- lapply(c(up = "up", down = "down"), function(m) intersect(seeds$gene_id[seeds$module == m], rownames(v)))
  fr <- fry(v, index = idx, design = design, contrast = "diseaseALS", sort = "none")
  ca <- camera(v, index = idx, design = design, contrast = "diseaseALS", sort = FALSE)
  sets <- data.frame(model = label, module = names(idx), n_als = sum(ss$disease == "ALS"), n_control = sum(ss$disease == "Control"),
                     residual_df = fit$df.residual[1], n_genes_tested = lengths(idx),
                     n_direction_kept = sapply(names(idx), function(m) sum(g$direction_kept[g$module == m], na.rm = TRUE)),
                     median_logFC = sapply(names(idx), function(m) median(g$logFC[g$module == m], na.rm = TRUE)),
                     fry_direction = fr[names(idx), "Direction"], fry_p = fr[names(idx), "PValue"],
                     camera_direction = ca[names(idx), "Direction"], camera_p = ca[names(idx), "PValue"])
  list(genes = cbind(model = label, g), sets = sets, v = v,
       spearman = cor(g$logFC, g$w34_logFC, method = "spearman", use = "complete.obs"))
}

primary <- fit_model(rep(TRUE, nrow(s)), "age + sex + PMI + tissue_subregion", "W34_primary_reproduction")
repro <- max(abs(primary$genes$logFC - primary$genes$w34_logFC))
if (repro > 1e-8) stop("W34 primary logFC not reproduced: ", repro)
cgnd <- fit_model(s$programme == "CGND", "age + sex", "CGND_only")
prog <- fit_model(rep(TRUE, nrow(s)), "age + sex + PMI + programme", "programme_adjusted")

# ---- 4. control-source comparison (full-cohort voom) ----
scores <- do.call(rbind, lapply(c("up", "down"), function(m) {
  g <- intersect(seeds$gene_id[seeds$module == m], rownames(primary$v))
  data.frame(library = colnames(primary$v), module = m, score = colMeans(primary$v$E[g, , drop = FALSE]))
}))
scores <- merge(scores, s[, c("library", "donor", "disease", "programme", "tissue_subregion", "site")], by = "library")
scores$group <- ifelse(scores$disease == "ALS", "ALS", paste0(scores$programme, " control"))
group_summary <- aggregate(score ~ module + group, data = scores, FUN = function(x) c(n = length(x), median = median(x), min = min(x), max = max(x)))
group_summary <- do.call(data.frame, group_summary)

# ---- reading (Amendment 7) ----
n_kept <- sum(cgnd$genes$direction_kept, na.rm = TRUE)
signs_kept <- all(sign(cgnd$sets$median_logFC) == c(up = 1, down = -1)[cgnd$sets$module])
reading <- if (n_kept >= 24 && signs_kept) "not explained by tissue programme" else if (n_kept < 15) "tissue programme is a plausible explanation" else "inconclusive"
summary_tab <- data.frame(model = c("CGND_only", "programme_adjusted"), n_genes = c(sum(!is.na(cgnd$genes$logFC)), sum(!is.na(prog$genes$logFC))),
                          n_direction_kept = c(n_kept, sum(prog$genes$direction_kept, na.rm = TRUE)),
                          spearman_vs_w34 = c(cgnd$spearman, prog$spearman),
                          reading = c(reading, "secondary; descriptive"))

write.table(s[, c("library", "geo_library_id", "donor", "disease", "programme", "tissue_subregion", "site", "age", "sex", "PMI", "n_cells", "lib_size")],
            file.path(out, "w62_discovery_libraries.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bal, file.path(out, "w62_balance_by_diagnosis.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(site_tab, file.path(out, "w62_collection_site_by_diagnosis.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rbind(cgnd$genes, prog$genes), file.path(out, "w62_module_gene_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rbind(primary$sets, cgnd$sets, prog$sets), file.path(out, "w62_module_set_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(scores, file.path(out, "w62_library_module_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(group_summary, file.path(out, "w62_module_scores_by_control_source.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_tab, file.path(out, "w62_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("W34 primary reproduction max |diff|:", repro, "\n")
print(bal, digits = 3, row.names = FALSE)
print(site_tab, row.names = FALSE)
print(rbind(primary$sets, cgnd$sets, prog$sets), digits = 3, row.names = FALSE)
print(summary_tab, digits = 3, row.names = FALSE)
print(group_summary, digits = 3, row.names = FALSE)
