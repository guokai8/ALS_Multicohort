#!/usr/bin/env Rscript
# Fig. 2 | Discovery signal is robust to contamination correction but multi-lineage in baseline expression.

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("Fig2_contamination_cell_source")

# ---- a: discovery module effects across contamination-correction variants (W57) ----
w57 <- read.delim(file.path(RES, "W57_ambient_contamination_test", "w57_module_tests.tsv"))
a <- w57[startsWith(w57$analysis, "discovery:"), ]
variant_lab <- c(a_uncorrected = "Uncorrected", b_corrected_rho = "Corrected (rho)", c_corrected_2rho = "Corrected (2 x rho)",
                 d_corrected_micro_rho = "Corrected (microglial rho)", e_rho_covariate = "log(rho) covariate")
a$variant <- factor(variant_lab[sub("^discovery:", "", a$analysis)], levels = variant_lab)
a$module <- factor(ifelse(a$module == "up", "ALS-up", "ALS-down"), levels = c("ALS-up", "ALS-down"))
a$p_label <- formatC(a$fry_p, format = "e", digits = 1)
write_source(a[, c("variant", "module", "n_genes_tested", "median_logFC", "fry_p", "camera_p")], out, "panel_a_discovery_correction")
pa <- ggplot(a, aes(x = median_logFC, y = variant, colour = module)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_point(size = 1.8) +
  geom_text(aes(label = paste0("P = ", p_label)), nudge_y = 0.36, size = 1.7, family = "Arial", show.legend = FALSE) +
  scale_x_continuous(limits = c(-1.9, 3.8), breaks = -1:3) +
  scale_colour_manual(values = c("ALS-up" = PAL$up, "ALS-down" = PAL$down), name = NULL) +
  scale_y_discrete(limits = rev) +
  labs(x = "Median module log2 fold change", y = NULL, title = "Discovery: contamination correction") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))

# ---- b: contamination fraction by diagnosis (W57) ----
rho <- read.delim(file.path(RES, "W57_ambient_contamination_test", "w57_rho_by_unit.tsv"))
rtest <- read.delim(file.path(RES, "W57_ambient_contamination_test", "w57_rho_disease_test.tsv"))
cohort_lab <- c(discovery = "Discovery", PinedaBA9_210526 = "Pineda BA9", Ulm = "Ulm\u2020")
p_by_cohort <- setNames(formatC(rtest$wilcox_p, format = "f", digits = 2), rtest$cohort)
cohort_lab <- setNames(paste0(cohort_lab, "\nP = ", p_by_cohort[names(cohort_lab)]), names(cohort_lab))
rho$cohort_lab <- factor(cohort_lab[rho$cohort], levels = cohort_lab)
rho$disease <- factor(rho$disease, levels = c("Control", "ALS"))
write_source(rho[, c("cohort", "unit", "disease", "rho")], out, "panel_b_rho_by_unit")
set.seed(20261012)
pb <- ggplot(rho, aes(x = disease, y = rho, colour = disease)) +
  geom_jitter(width = 0.18, height = 0, size = 0.7, alpha = 0.8) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = PAL$ink) +
  facet_wrap(~ cohort_lab, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c(Control = PAL$control, ALS = PAL$als), guide = "none") +
  labs(x = NULL, y = "Contamination fraction (rho)", title = "Contamination fraction by diagnosis")

# ---- c: baseline lineage share of module genes (W56 E1, Pineda) ----
spec <- read.delim(file.path(RES, "W56_cell_source_check", "w56_e1_module_gene_specificity.tsv"), check.names = FALSE)
lin_cols <- grep("^pineda_share_", names(spec), value = TRUE)
long <- do.call(rbind, lapply(lin_cols, function(cl) data.frame(gene = spec$common_name, module = spec$module,
  lineage = sub("^pineda_share_", "", cl), share = spec[[cl]], enriched = spec$astrocyte_enriched %in% c(TRUE, "True"))))
lineage_order <- c("Astrocyte", "Microglia", "Endothelial", "Vascular/Meningeal", "Oligodendrocyte", "OPC", "Excitatory neuron", "Inhibitory neuron")
long$lineage <- factor(long$lineage, levels = lineage_order)
gene_order <- spec[order(spec$module != "up", -spec$pineda_share_Astrocyte), "common_name"]
enriched_genes <- spec$common_name[spec$astrocyte_enriched %in% c(TRUE, "True")]
long$gene <- factor(long$gene, levels = rev(gene_order))
gene_labels <- setNames(ifelse(levels(long$gene) %in% enriched_genes, paste0(levels(long$gene), " *"), levels(long$gene)), levels(long$gene))
long$module <- factor(ifelse(long$module == "up", "ALS-up module", "ALS-down module"), levels = c("ALS-up module", "ALS-down module"))
write_source(long, out, "panel_c_lineage_share")
pc <- ggplot(long, aes(x = lineage, y = gene, fill = share)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "white", high = PAL$ink, limits = c(0, 1), name = "Share of\nsummed CPM") +
  scale_y_discrete(labels = gene_labels) +
  labs(x = "Lineage (Pineda, human BA4 + BA9)", y = NULL, title = "Baseline lineage share of module genes") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0))

fig <- (pc | (pa / pb)) + plot_layout(widths = c(1.05, 1)) + plot_annotation(tag_levels = "a")
save_panels(list(pa = pa, pb = pb, pc = pc), "Figure2_contamination_cell_source")
save_pub_r(fig, out, "Figure2_contamination_cell_source", width_mm = 183, height_mm = 125)
cat("Fig. 2 written to", out, "\n")
