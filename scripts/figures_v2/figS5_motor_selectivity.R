#!/usr/bin/env Rscript
# Supplementary Fig. S5 | Motor- vs sensory-territory ALS effects in GSE330224 Visium (W67, Amendment 10).

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("FigS5_motor_selectivity")
w67 <- file.path(RES, "W67_motor_selectivity")
cmp_lab <- c(grey_ventral_vs_dorsal_horn = "Spinal grey:\nventral vs\ndorsal horn",
             white_lateral_funiculus_vs_dorsal_columns = "Spinal white:\nlateral funiculus vs\ndorsal columns",
             cortex_layer5_vs_layer2_3 = "Motor cortex:\nlayer 5 vs\nlayers 2-3")
terr_lab <- c(ALS_in_motor = "Motor", ALS_in_sensory = "Sensory")
set_lab <- c(up = "Frozen ALS-up", up_consensus = "Consensus up", down_consensus = "Consensus down")

# ---- a: ALS effect in each territory, with the interaction P ----
st <- read.delim(file.path(w67, "w67_motor_selectivity_tests.tsv"))
st <- st[st$gene_set %in% names(set_lab), ]
eff <- st[st$contrast %in% names(terr_lab), ]
eff$territory <- factor(terr_lab[eff$contrast], levels = terr_lab)
eff$comparison <- factor(cmp_lab[eff$comparison], levels = cmp_lab)
eff$set <- factor(set_lab[eff$gene_set], levels = set_lab)
eff$sig <- ifelse(eff$fry_p < 0.05, "fry P < 0.05", "fry P >= 0.05")
inter <- st[st$contrast == "interaction", ]
inter$comparison <- factor(cmp_lab[inter$comparison], levels = cmp_lab)
inter$set <- factor(set_lab[inter$gene_set], levels = set_lab)
inter$lab <- paste0("interaction P = ", formatC(inter$fry_p, format = "fg", digits = 2),
                    ifelse(is.na(inter$fry_BH_primary_family), "", paste0("\nBH = ", formatC(inter$fry_BH_primary_family, format = "fg", digits = 2))))
write_source(eff[, c("comparison", "territory", "gene_set", "n_als", "n_control", "median_logFC", "fry_direction", "fry_p", "camera_p")], out, "panel_a_territory_effects")
write_source(inter[, c("comparison", "gene_set", "median_logFC", "fry_direction", "fry_p", "camera_p", "fry_BH_primary_family")], out, "panel_a_interaction")
pa <- ggplot(eff, aes(x = territory, y = median_logFC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_line(aes(group = 1), linewidth = 0.35, colour = PAL$neutral) +
  geom_point(aes(shape = sig, colour = territory), size = 1.8, stroke = 0.5) +
  geom_text(data = inter, aes(x = 1.5, y = 0.43, label = lab), size = 1.6, family = "Arial", lineheight = 0.9, inherit.aes = FALSE) +
  facet_grid(set ~ comparison) +
  scale_y_continuous(limits = c(-0.2, 0.5)) +
  scale_shape_manual(values = c(`fry P < 0.05` = 16, `fry P >= 0.05` = 1), name = NULL) +
  scale_colour_manual(values = c(Motor = PAL$als, Sensory = PAL$control), guide = "none") +
  labs(x = "Territory", y = "Median log2 fold change, ALS vs control", title = "ALS effect by territory (independent donors)") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), strip.text.y = element_text(angle = 0, hjust = 0))

# ---- b: ALS-up genes, ventral vs dorsal horn ----
g <- read.delim(file.path(w67, "w67_module_gene_results.tsv"))
g <- g[g$comparison == "grey_ventral_vs_dorsal_horn" & g$gene_set == "up" & g$contrast %in% names(terr_lab), ]
g$territory <- factor(c(ALS_in_motor = "Ventral horn", ALS_in_sensory = "Dorsal horn")[g$contrast], levels = c("Dorsal horn", "Ventral horn"))
ord <- with(g[g$territory == "Ventral horn", ], symbol[order(logFC)])
g$symbol <- factor(g$symbol, levels = ord)
g$sig <- ifelse(g$P.Value < 0.05, "P < 0.05", "P >= 0.05")
write_source(g, out, "panel_b_genes_ventral_dorsal")
pb <- ggplot(g, aes(x = logFC, y = symbol)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_line(aes(group = symbol), linewidth = 0.35, colour = PAL$light) +
  geom_point(aes(colour = territory, shape = sig), size = 1.6, stroke = 0.5) +
  scale_colour_manual(values = c(`Ventral horn` = PAL$als, `Dorsal horn` = PAL$control), name = NULL) +
  scale_shape_manual(values = c(`P < 0.05` = 16, `P >= 0.05` = 1), name = NULL) +
  labs(x = "log2 fold change, ALS vs control", y = NULL, title = "ALS-up genes in spinal grey matter") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), legend.box = "vertical")

fig <- (pa | pb) + plot_layout(widths = c(2, 1)) + plot_annotation(tag_levels = "a")
save_pub_r(fig, out, "FigureS5_motor_selectivity", width_mm = 183, height_mm = 105)
cat("Supplementary Fig. S5 written to", out, "\n")
