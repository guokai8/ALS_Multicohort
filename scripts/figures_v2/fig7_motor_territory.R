#!/usr/bin/env Rscript
# Fig. 7 | The frozen ALS-up module rises in spinal motor territory in two spatial cohorts
# (W67, GSE330224 lumbar; W68, GSE288365 cervical) and co-localises with myeloid proteins (W68 C).

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("Fig7_motor_territory")
cohort_lab <- c(W67 = "Lumbar cord\n18 ALS / 5 ctrl", W68 = "Cervical cord\n9 ALS / 7 ctrl")
set_lab <- c(up = "ALS-up\nmodule")

# ---- a: ALS effect in motor vs sensory grey matter, interaction P ----
w67 <- read.delim(file.path(RES, "W67_motor_selectivity", "w67_motor_selectivity_tests.tsv"))
w67 <- w67[w67$comparison == "grey_ventral_vs_dorsal_horn", ]
w67$cohort <- "W67"; w67$bh <- w67$fry_BH_primary_family
w68 <- read.delim(file.path(RES, "W68_GSE288365_replication", "w68_replication_tests.tsv"))
w68 <- w68[w68$analysis == "primary" & w68$comparison == "grey", ]
w68$cohort <- "W68"; w68$bh <- w68$fry_BH_primary_family
cols <- c("cohort", "contrast", "gene_set", "median_logFC", "fry_p", "camera_p", "bh", "n_als", "n_control")
st <- rbind(w67[, cols], w68[, cols])
st <- st[st$gene_set %in% names(set_lab), ]
eff <- st[st$contrast %in% c("ALS_in_motor", "ALS_in_sensory"), ]
eff$territory <- factor(ifelse(eff$contrast == "ALS_in_motor", "Ventral / anterior horn", "Dorsal horn / other grey"),
                        levels = c("Ventral / anterior horn", "Dorsal horn / other grey"))
eff$cohort_l <- factor(cohort_lab[eff$cohort], levels = cohort_lab)
eff$set <- factor(set_lab[eff$gene_set], levels = set_lab)
eff$sig <- ifelse(eff$fry_p < 0.05, "fry P < 0.05", "fry P >= 0.05")
inter <- st[st$contrast == "interaction", ]
inter$cohort_l <- factor(cohort_lab[inter$cohort], levels = cohort_lab)
inter$set <- factor(set_lab[inter$gene_set], levels = set_lab)
inter$lab <- paste0("interaction P = ", formatC(inter$fry_p, format = "fg", digits = 2), "\nBH = ", formatC(inter$bh, format = "fg", digits = 2))
write_source(eff, out, "panel_a_territory_effects")
write_source(inter, out, "panel_a_interaction")
pa <- ggplot(eff, aes(x = territory, y = median_logFC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_col(aes(fill = territory), width = 0.6) +
  geom_point(aes(shape = sig), size = 1.4, colour = PAL$ink) +
  geom_text(data = inter, aes(x = 1.5, y = 1.42, label = lab), size = 1.6, family = "Arial", lineheight = 0.9, inherit.aes = FALSE) +
  facet_grid(set ~ cohort_l) +
  scale_y_continuous(limits = c(-0.1, 1.6)) +
  scale_fill_manual(values = c(`Ventral / anterior horn` = PAL$als, `Dorsal horn / other grey` = PAL$control), guide = "none") +
  scale_shape_manual(values = c(`fry P < 0.05` = 16, `fry P >= 0.05` = 1), name = NULL) +
  scale_x_discrete(labels = c("Motor", "Sensory")) +
  labs(x = "Spinal grey matter territory", y = "Median log2 fold change, ALS vs control", title = "ALS effect by territory") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), strip.text.y = element_text(angle = 0, hjust = 0))

# ---- b: module genes in motor territory, both cohorts ----
g67 <- read.delim(file.path(RES, "W67_motor_selectivity", "w67_module_gene_results.tsv"))
g67 <- g67[g67$comparison == "grey_ventral_vs_dorsal_horn" & g67$gene_set == "up" & g67$contrast %in% c("ALS_in_motor", "ALS_in_sensory"), ]
g67$cohort <- "W67"
g68 <- read.delim(file.path(RES, "W68_GSE288365_replication", "w68_module_gene_results.tsv"))
g68 <- g68[g68$analysis == "primary" & g68$comparison == "grey" & g68$gene_set == "up" & g68$contrast %in% c("ALS_in_motor", "ALS_in_sensory"), ]
g68$cohort <- "W68"
gc <- c("cohort", "contrast", "symbol", "logFC", "P.Value")
g <- rbind(g67[, gc], g68[, gc])
shared <- intersect(g67$symbol, g68$symbol)
g <- g[g$symbol %in% shared, ]
ord <- with(g[g$cohort == "W68" & g$contrast == "ALS_in_motor", ], symbol[order(logFC)])
g$symbol <- factor(g$symbol, levels = ord)
g$territory <- factor(ifelse(g$contrast == "ALS_in_motor", "Motor", "Sensory"), levels = c("Sensory", "Motor"))
g$cohort_l <- factor(sub("\n.*", "", cohort_lab[g$cohort]), levels = sub("\n.*", "", cohort_lab))
g$sig <- ifelse(g$P.Value < 0.05, "P < 0.05", "P >= 0.05")
write_source(g, out, "panel_b_genes")
pb <- ggplot(g, aes(x = logFC, y = symbol)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_line(aes(group = symbol), linewidth = 0.35, colour = PAL$light) +
  geom_point(aes(colour = territory, shape = sig), size = 1.5, stroke = 0.5) +
  facet_grid(. ~ cohort_l) +
  scale_x_continuous(breaks = c(0, 1, 2, 3)) +
  scale_colour_manual(values = c(Motor = PAL$als, Sensory = PAL$control), name = NULL) +
  scale_shape_manual(values = c(`P < 0.05` = 16, `P >= 0.05` = 1), name = NULL) +
  labs(x = "log2 fold change, ALS vs control", y = NULL, title = "Module genes") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), legend.box = "vertical")

# ---- c: protein co-localisation (W68 C) ----
pc_d <- read.delim(file.path(RES, "W68_GSE288365_replication", "w68_protein_colocalisation_by_section.tsv"))
pc_d <- pc_d[pc_d$scope == "all_annotated", ]
prot_lab <- c(VIM.1 = "Vimentin", CD68.1 = "CD68", `HLA-DRA` = "HLA-DR", CD163.1 = "CD163", ACTA2.1 = "ACTA2 (mural)",
              PECAM1.1 = "PECAM1 (endothelial)", PTPRC.1 = "CD45", ITGAM.1 = "CD11b")
pc_d$protein <- factor(prot_lab[pc_d$protein], levels = rev(prot_lab))
write_source(pc_d, out, "panel_c_protein_colocalisation")
set.seed(20261013)
pc <- ggplot(pc_d, aes(x = rho_depth_partial, y = protein)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_jitter(aes(colour = diagnosis == "Control"), height = 0.15, width = 0, size = 0.7, alpha = 0.8) +
  stat_summary(fun = median, geom = "point", shape = 124, size = 3.2, colour = PAL$ink) +
  scale_colour_manual(values = c(`TRUE` = PAL$control, `FALSE` = PAL$als), labels = c(`TRUE` = "Control", `FALSE` = "ALS"), name = NULL) +
  labs(x = "Spearman rho with module\nscore (depth-partial; cervical)", y = NULL, title = "Protein co-localisation") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))

fig <- (pa | pb | pc) + plot_layout(widths = c(1.35, 1, 0.9)) + plot_annotation(tag_levels = "a")
save_panels(list(pa = pa, pb = pb, pc = pc), "Figure7_motor_territory")
save_pub_r(fig, out, "Figure7_motor_territory", width_mm = 183, height_mm = 95)
cat("Fig. 7 written to", out, "\n")
