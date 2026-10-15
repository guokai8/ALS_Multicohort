#!/usr/bin/env Rscript
# Supplementary Fig. S10 | SOD1-G93A spatial time course (W71 A), disease-state lineage share in human
# spinal cord (W72) and spot co-localisation in a third human spatial cohort (W71 B); Amendment 14.

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("FigS10_mouse_cell_source")
w71 <- file.path(RES, "W71_maniatis_spatial_time_course")
w72 <- file.path(RES, "W72_zelic_disease_cell_source")
ages <- c("p30", "p70", "p100", "p120")
age_lab <- c(p30 = "P30", p70 = "P70", p100 = "P100", p120 = "P120")

# ---- a: frozen module (mouse orthologs) and C3 in ventral and dorsal horn by age ----
t <- read.delim(file.path(w71, "w71_mouse_tests.tsv"))
t <- t[t$comparison == "grey" & !t$adjusted & t$gene_set == "up" & t$contrast %in% c("ventral", "dorsal"), ]
g <- read.delim(file.path(w71, "w71_mouse_gene_effects.tsv"))
g <- g[g$comparison == "grey" & !g$adjusted & g$gene == "C3" & g$contrast %in% c("ventral", "dorsal"), ]
a <- rbind(data.frame(what = "ALS-up module (median)", age = t$age, contrast = t$contrast, logFC = t$median_logFC, sig = t$fry_p < 0.05),
           data.frame(what = "C3", age = g$age, contrast = g$contrast, logFC = g$logFC, sig = g$P < 0.05))
a$age <- factor(age_lab[a$age], levels = age_lab)
a$territory <- factor(ifelse(a$contrast == "ventral", "Ventral horn", "Dorsal horn"), levels = c("Ventral horn", "Dorsal horn"))
a$what <- factor(a$what, levels = c("ALS-up module (median)", "C3"))
ti <- read.delim(file.path(w71, "w71_mouse_tests.tsv"))
ti <- ti[ti$comparison == "grey" & !ti$adjusted & ti$gene_set == "up" & ti$contrast == "interaction", ]
lab <- data.frame(what = factor("ALS-up module (median)", levels = levels(a$what)), age = factor(age_lab[ti$age], levels = age_lab),
                  label = paste0("int. P = ", formatC(ti$fry_p, format = "fg", digits = 2)), y = max(a$logFC) * 1.12)
write_source(a, out, "panel_a_mouse_time_course")
pa <- ggplot(a, aes(x = age, y = logFC, colour = territory, group = territory)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_line(linewidth = 0.5) +
  geom_point(aes(shape = sig), size = 1.8, stroke = 0.5, fill = "white") +
  geom_text(data = lab, aes(x = age, y = y, label = label), inherit.aes = FALSE, size = 1.6, family = "Arial") +
  facet_wrap(~ what, nrow = 1) +
  scale_colour_manual(values = c(`Ventral horn` = PAL$als, `Dorsal horn` = PAL$control), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "P < 0.05", `FALSE` = "P ≥ 0.05"), name = NULL) +
  labs(x = "Age (SOD1-G93A vs SOD1-WT)", y = "log2FC, G93A vs WT", title = "SOD1-G93A lumbar cord (Maniatis et al.)") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))

# ---- b: lineage share of focus genes in ALS donors (W72) ----
focus <- c("GPNMB", "C3", "LTF", "ANXA1", "EMP3")
lineages <- c("Microglia", "Astrocyte", "Vascular", "Oligodendrocyte", "OPC", "Neuron")
s <- read.delim(file.path(w72, "w72_lineage_share.tsv"))
s <- s[s$gene %in% focus, c("gene", lineages)]
b <- do.call(rbind, lapply(lineages, function(l) data.frame(gene = s$gene, lineage = l, share = s[[l]])))
b$gene <- factor(b$gene, levels = rev(focus)); b$lineage <- factor(b$lineage, levels = lineages)
b$lab <- ifelse(b$share >= 0.1, sprintf("%.2f", b$share), "")
write_source(b, out, "panel_b_lineage_share_ALS")
pb <- ggplot(b, aes(x = lineage, y = gene, fill = share)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = lab, colour = share > 0.5), size = 1.6, family = "Arial") +
  scale_fill_gradient(low = "white", high = PAL$ink, limits = c(0, 1), name = "Share") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = PAL$ink), guide = "none") +
  labs(x = NULL, y = NULL, title = "Lineage share, ALS cord") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())

# ---- c: within-lineage ALS vs control log2FC of focus genes (W72) ----
de <- read.delim(file.path(w72, "w72_within_lineage_de.tsv"))
de <- de[de$gene %in% focus & de$lineage %in% lineages, ]
de$gene <- factor(de$gene, levels = rev(focus)); de$lineage <- factor(de$lineage, levels = lineages)
de$lab <- ifelse(de$P < 0.05, "*", "")
write_source(de, out, "panel_c_within_lineage_de")
pc <- ggplot(de, aes(x = lineage, y = gene, fill = logFC)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = lab), size = 2.4, family = "Arial", vjust = 0.75) +
  scale_fill_gradient2(low = PAL$control, mid = "white", high = PAL$als, midpoint = 0, name = "log2FC") +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, title = "ALS vs control (* P < 0.05)") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())

# ---- d: human co-localisation per section (W71 B) ----
h <- read.delim(file.path(w71, "w71_human_colocalisation_sections.tsv"))
lin_lab <- c(microglia = "Microglia", astrocyte = "Astrocyte", mural_fibroblast = "Mural/fibroblast", endothelial = "Endothelial",
             oligodendrocyte = "Oligodendrocyte", neuron = "Neuron")
h$lineage <- factor(lin_lab[h$lineage], levels = lin_lab)
med <- aggregate(rho_depth_partial ~ lineage, h, median)
write_source(h, out, "panel_d_human_colocalisation")
pd <- ggplot(h, aes(x = lineage, y = rho_depth_partial)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_jitter(width = 0.18, height = 0, size = 0.6, alpha = 0.5, colour = PAL$neutral) +
  geom_crossbar(data = med, aes(ymin = rho_depth_partial, ymax = rho_depth_partial), width = 0.55, linewidth = 0.4, colour = PAL$ink) +
  labs(x = NULL, y = "Depth-partial Spearman rho\n(module vs lineage score)", title = "Human ALS cord (77 sections)") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

save_panels(list(pa = pa, pb = pb, pc = pc, pd = pd), "FigureS10_mouse_cell_source")
fig <- pa / ((pb | pc | pd) + plot_layout(widths = c(1, 1, 1.1))) + plot_layout(heights = c(1, 1.05)) + plot_annotation(tag_levels = "a")
save_pub_r(fig, out, "FigureS10_mouse_cell_source", width_mm = 183, height_mm = 140)
cat("Supplementary Fig. S10 written to", out, "\n")
