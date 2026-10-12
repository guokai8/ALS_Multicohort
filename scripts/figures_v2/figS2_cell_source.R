#!/usr/bin/env Rscript
# Supplementary Fig. S2 | Cell source of the module genes: Ulm lineage shares and per-lineage disease shifts (W56 E1, E2).

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("FigS2_cell_source")
W56 <- file.path(RES, "W56_cell_source_check")

# ---- a: baseline lineage share in Ulm, ordered as in Fig. 2a (Pineda astrocyte share) ----
spec <- read.delim(file.path(W56, "w56_e1_module_gene_specificity.tsv"), check.names = FALSE)
ulm_cols <- grep("^ulm_share_", names(spec), value = TRUE)
ulm_lab <- c(Astrocytes = "Astrocyte", Microglia = "Microglia", Oligodendrocytes = "Oligodendrocyte", OPC = "OPC",
             Exc_Neurons = "Excitatory neuron", Inh_Neurons = "Inhibitory neuron")
long <- do.call(rbind, lapply(ulm_cols, function(cl) data.frame(gene = spec$common_name, module = spec$module,
  lineage = unname(ulm_lab[sub("^ulm_share_", "", cl)]), share = spec[[cl]])))
long$lineage <- factor(long$lineage, levels = ulm_lab)
gene_order <- spec[order(spec$module != "up", -spec$pineda_share_Astrocyte), "common_name"]
long$gene <- factor(long$gene, levels = rev(gene_order))
enriched <- spec$common_name[spec$astrocyte_enriched %in% c(TRUE, "True")]
gene_labels <- setNames(ifelse(levels(long$gene) %in% enriched, paste0(levels(long$gene), " *"), levels(long$gene)), levels(long$gene))
long$module <- factor(ifelse(long$module == "up", "ALS-up module", "ALS-down module"), levels = c("ALS-up module", "ALS-down module"))
write_source(long, out, "panel_a_ulm_lineage_share")
pa <- ggplot(long, aes(x = lineage, y = gene, fill = share)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "white", high = PAL$ink, limits = c(0, 1), na.value = "#F2F2F2", name = "Share of\nsummed CPM") +
  scale_y_discrete(labels = gene_labels) +
  labs(x = "Lineage (Ulm/DZNE, author labels)", y = NULL, title = "Baseline lineage share in Ulm") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0))

# ---- b: Pineda vs Ulm astrocyte share per gene ----
b <- data.frame(gene = spec$common_name, module = spec$module, pineda = spec$pineda_share_Astrocyte, ulm = spec$ulm_share_Astrocytes)
b$module <- factor(ifelse(b$module == "up", "ALS-up", "ALS-down"), levels = c("ALS-up", "ALS-down"))
b <- b[!is.na(b$ulm), ]
rho_b <- cor(b$pineda, b$ulm, method = "spearman")
write_source(b, out, "panel_b_astrocyte_share_pineda_vs_ulm")
pb <- ggplot(b, aes(x = pineda, y = ulm, colour = module)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.3, colour = PAL$neutral) +
  geom_vline(xintercept = 0.5, linewidth = 0.3, colour = PAL$light) +
  geom_point(size = 1.3) +
  ggrepel::geom_text_repel(data = b[b$pineda >= 0.3 | b$ulm >= 0.5, ], aes(label = gene), size = 1.8, family = "Arial",
                           show.legend = FALSE, segment.size = 0.2, min.segment.length = 0, seed = 1) +
  annotate("text", x = 0, y = 1, hjust = 0, vjust = 1, size = 1.9, family = "Arial", label = sprintf("Spearman ρ = %.2f", rho_b)) +
  scale_colour_manual(values = c("ALS-up" = PAL$up, "ALS-down" = PAL$down), name = NULL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Astrocyte share, Pineda", y = "Astrocyte share, Ulm", title = "Astrocyte share, Pineda vs Ulm") +
  theme(legend.position = "bottom")

# ---- c: disease shift of the modules within each lineage (E2) ----
e2 <- read.delim(file.path(W56, "w56_e2_celltype_disease_shift.tsv"))
e2$cohort <- ifelse(startsWith(e2$analysis, "Ulm"), "Ulm/DZNE (motor cortex)", "Pineda BA9 (batch 210526)")
e2$lineage <- sub("^(Ulm|PinedaBA9_210526)_", "", e2$analysis)
e2$lineage <- c(Astrocytes = "Astrocyte", Astrocyte = "Astrocyte", Microglia = "Microglia", Oligodendrocytes = "Oligodendrocyte",
                OPC = "OPC", Endothelial = "Endothelial")[e2$lineage]
e2$module <- factor(ifelse(e2$module == "up", "ALS-up", "ALS-down"), levels = c("ALS-up", "ALS-down"))
e2$lineage <- factor(e2$lineage, levels = rev(c("Astrocyte", "Microglia", "Endothelial", "Oligodendrocyte", "OPC")))
e2 <- droplevels(e2)
# Up and down modules share a row per lineage (dodged); P values sit in a column to the right.
e2$label <- sprintf("P = %s (%d genes)", formatC(e2$fry_p, format = "fg", digits = 2), e2$n_genes_tested)
write_source(e2, out, "panel_c_lineage_disease_shift")
x_lab <- max(abs(e2$median_logFC)) * 1.25
dodge <- position_dodge(width = 0.55, orientation = "y")
pc <- ggplot(e2, aes(x = median_logFC, y = lineage, colour = module, group = module)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_point(size = 1.5, position = dodge) +
  geom_text(aes(x = x_lab, label = label), hjust = 0, size = 1.65, family = "Arial", position = dodge, show.legend = FALSE) +
  facet_wrap(~ cohort, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("ALS-up" = PAL$up, "ALS-down" = PAL$down), guide = "none") +
  coord_cartesian(xlim = c(-x_lab, x_lab), clip = "off") +
  labs(x = expression("Median module log"[2]*" fold change (ALS vs control)"), y = NULL,
       title = "Per-lineage disease shift (fry P)") +
  theme(strip.text = element_text(hjust = 0), plot.margin = margin(4, 62, 4, 4))

fig <- (pa | (pb / pc + plot_layout(heights = c(1, 1.15)))) + plot_layout(widths = c(0.9, 1.1)) + plot_annotation(tag_levels = "a")
save_pub_r(fig, out, "FigureS2_cell_source", width_mm = 183, height_mm = 150)
cat("Fig. S2 written to", out, "\n")
