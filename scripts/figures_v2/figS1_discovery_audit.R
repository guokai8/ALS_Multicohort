#!/usr/bin/env Rscript
# Supplementary Fig. S1 | Discovery-cohort confounder audit (W62, Amendment 7) and model dependence (W38).

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("FigS1_discovery_audit")
W62 <- file.path(RES, "W62_discovery_confounder_audit")
MODULE_COL <- c("ALS-up" = PAL$up, "ALS-down" = PAL$down)
module_lab <- function(m) factor(ifelse(m == "up", "ALS-up", "ALS-down"), levels = c("ALS-up", "ALS-down"))

# ---- a: libraries by tissue programme, subregion and diagnosis ----
lib <- read.delim(file.path(W62, "w62_discovery_libraries.tsv"))
lib$subregion <- sub("^Cortex_Motor_", "", lib$tissue_subregion)
a <- as.data.frame(table(programme = lib$programme, subregion = lib$subregion, disease = lib$disease))
a$stratum <- factor(paste0(a$programme, " · ", a$subregion),
                    levels = rev(c("U54 · BA4", "CGND · BA4", "CGND · Lateral", "CGND · Medial")))
a <- a[!is.na(a$stratum), ]
a$disease <- factor(a$disease, levels = c("Control", "ALS"))
write_source(a, out, "panel_a_library_composition")
pa <- ggplot(a, aes(x = disease, y = stratum, fill = Freq)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(Freq > 0, Freq, ""), colour = Freq > 3), size = 2, family = "Arial", show.legend = FALSE) +
  scale_fill_gradient(low = "#F2F2F2", high = PAL$neutral, guide = "none") +
  scale_colour_manual(values = c(`FALSE` = PAL$ink, `TRUE` = "white")) +
  labs(x = NULL, y = "Tissue programme · subregion", title = "Discovery libraries") +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

# ---- b: model dependence of the module genes (W38 nested models) ----
w38 <- read.delim(file.path(RES, "W38_W34_MTC_astrocyte_model_sensitivity", "w38_nested_model_audit.tsv"))
w38 <- w38[w38$modality == "RNA", ]
model_lab <- c(D0_disease_only = "disease", D1_disease_sex = "+ sex", D2_disease_age_sex = "+ age",
               D3_disease_age_sex_PMI = "+ PMI", D4_W34_primary_full = "+ subregion\n(primary)")
b <- rbind(data.frame(model = w38$model_id, metric = "Genes at FDR < 0.05\n(genome-wide)", n = w38$n_model_fdr05),
           data.frame(model = w38$model_id, metric = "Module genes, same direction\nand nominal P < 0.05 (of 29)", n = w38$n_primary_nominal_directional))
b$model <- factor(model_lab[b$model], levels = model_lab)
write_source(b, out, "panel_b_nested_models")
pb <- ggplot(b, aes(x = model, y = n)) +
  geom_col(fill = PAL$neutral, width = 0.65) +
  geom_text(aes(label = n), vjust = -0.4, size = 1.9, family = "Arial") +
  facet_wrap(~ metric, nrow = 1) +
  scale_y_continuous(limits = c(0, 33), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Covariates added (cumulative)", y = "Genes", title = "Dependence on the subregion term (W38)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# ---- c: CGND-only fold changes vs primary ----
g <- read.delim(file.path(W62, "w62_module_gene_results.tsv"))
g <- g[g$model == "CGND_only", ]
g$module <- module_lab(g$module)
rho <- cor(g$logFC, g$w34_logFC, method = "spearman")
write_source(g, out, "panel_c_cgnd_only_vs_primary")
lim <- range(c(g$logFC, g$w34_logFC, 0)) * 1.05
pc <- ggplot(g, aes(x = w34_logFC, y = logFC, colour = module)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.3, colour = PAL$neutral) +
  geom_point(size = 1.3) +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 1.9, family = "Arial",
           label = sprintf("%d/29 same direction\nSpearman ρ = %.2f", sum(g$direction_kept), rho)) +
  scale_colour_manual(values = MODULE_COL, name = NULL) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression("Primary model log"[2]*" fold change"), y = expression("CGND-only log"[2]*" fold change"),
       title = "Shared programme only (11 ALS, 2 controls)") +
  theme(legend.position = "bottom")

# ---- d: per-library module scores by control source ----
sc <- read.delim(file.path(W62, "w62_library_module_scores.tsv"))
sc$module <- module_lab(sc$module)
sc$group <- factor(sc$group, levels = c("U54 control", "CGND control", "ALS"))
write_source(sc, out, "panel_d_module_scores_by_source")
set.seed(20261012)
pd <- ggplot(sc, aes(x = group, y = score, colour = group)) +
  geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = PAL$ink) +
  facet_wrap(~ module, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c(`U54 control` = "#9DB4C8", `CGND control` = PAL$control, ALS = PAL$als), guide = "none") +
  labs(x = NULL, y = "Module score (mean voom log-CPM)", title = "Module scores by control source") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig <- ((pa | pb) + plot_layout(widths = c(0.7, 1.3))) / ((pc | pd) + plot_layout(widths = c(1, 1.2))) +
  plot_annotation(tag_levels = "a")
save_pub_r(fig, out, "FigureS1_discovery_audit", width_mm = 183, height_mm = 140)
cat("Fig. S1 written to", out, "\n")
