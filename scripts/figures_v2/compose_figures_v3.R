#!/usr/bin/env Rscript
# Streamlined figure set (manuscript v3): five main figures and eight supplementary figures,
# assembled from the panel objects saved by the figures_v2 scripts (run those first). The source
# data of each contributing panel are copied next to the assembled figure.

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
V3 <- file.path(RES, "Figures_v3")
P <- function(name) readRDS(file.path(FIG_OUT, "panels", paste0(name, ".rds")))
untag <- function(p) p + labs(tag = NULL)
tag <- function(p, t) p + labs(tag = t)

compose <- function(dirname, fname, fig, width, height, sources) {
  d <- file.path(V3, dirname)
  dir.create(file.path(d, "source_data"), recursive = TRUE, showWarnings = FALSE)
  for (src in sources) for (f in list.files(file.path(FIG_OUT, src, "source_data"), full.names = TRUE))
    file.copy(f, file.path(d, "source_data", paste0(src, "__", basename(f))), overwrite = TRUE)
  save_pub_r(fig, d, fname, width_mm = width, height_mm = height)
}

f1 <- P("Figure1_study_design"); f2 <- P("Figure2_contamination_cell_source"); f3 <- P("Figure3_cohort_confounding")
f4 <- P("Figure4_meta_analysis"); f5 <- P("Figure5_NYGC_bulk_holdout"); f6 <- P("Figure6_evidence_map")
f7 <- P("Figure7_motor_territory"); s1 <- P("FigureS1_discovery_audit"); s2 <- P("FigureS2_cell_source")
s3 <- P("FigureS3_holdout_genes"); s4 <- P("FigureS4_spatial"); s5 <- P("FigureS5_motor_selectivity"); s6 <- P("FigureS6_dose_response")

# ---- Fig. 1: design ----
compose("Fig1_design", "Figure1_design", (f1$pa | f1$pb) / f1$pc + plot_layout(heights = c(1.15, 1), widths = c(1.15, 1)) +
          plot_annotation(tag_levels = "a"), 183, 150, "Fig1_study_design")

# ---- Fig. 2: single-nucleus replication is shaped by cohort design ----
top <- (tag(untag(f3$pa), "a") | tag(untag(f3$pb), "b") | tag(untag(f3$pe), "c")) + plot_layout(widths = c(0.75, 1, 1.15))
mid <- (tag(untag(f3$pc1), "d") | untag(f3$pc2) | tag(untag(f3$pd1), "e") | untag(f3$pd2)) + plot_layout(widths = c(0.8, 0.8, 0.8, 1.25))
bot <- tag(untag(f4$pa), "f")
compose("Fig2_single_nucleus_replication", "Figure2_single_nucleus_replication",
        wrap_elements(full = top) / wrap_elements(full = mid) / wrap_elements(full = bot) + plot_layout(heights = c(1.15, 0.9, 1.35)),
        183, 215, c("Fig3_confounding", "Fig4_meta"))

# ---- Fig. 3: confirmation in an independent bulk hold-out ----
compose("Fig3_bulk_holdout", "Figure3_bulk_holdout",
        (tag(untag(f5$pa), "a") | tag(untag(f4$pb), "b") | tag(untag(f5$pb), "c")) + plot_layout(widths = c(1, 1.15, 1)),
        183, 85, c("Fig5_nygc_holdout", "Fig4_meta"))

# ---- Fig. 4: the signal is multicellular, not astrocyte-specific ----
spec <- read.delim(file.path(RES, "W56_cell_source_check", "w56_e1_module_gene_specificity.tsv"), check.names = FALSE)
spec <- spec[spec$module == "up", ]
lin_cols <- grep("^pineda_share_", names(spec), value = TRUE)
lineage_order <- c("Astrocyte", "Microglia", "Endothelial", "Vascular/Meningeal", "Oligodendrocyte", "OPC", "Excitatory neuron", "Inhibitory neuron")
ls_long <- do.call(rbind, lapply(lin_cols, function(cl) data.frame(gene = spec$common_name, lineage = sub("^pineda_share_", "", cl),
                                                                   share = spec[[cl]], enriched = spec$astrocyte_enriched %in% c(TRUE, "True"))))
ls_long$lineage <- factor(ls_long$lineage, levels = lineage_order)
ls_long$top <- sub("^pineda_share_", "", lin_cols[max.col(spec[, lin_cols])])[match(ls_long$gene, spec$common_name)]
gene_order <- spec$common_name[order(match(sub("^pineda_share_", "", lin_cols[max.col(spec[, lin_cols])]), lineage_order), -apply(spec[, lin_cols], 1, max))]
ls_long$gene <- factor(ls_long$gene, levels = rev(gene_order))
p4a <- ggplot(ls_long, aes(x = lineage, y = gene, fill = share)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "white", high = PAL$ink, limits = c(0, 1), name = "Share of\nsummed CPM") +
  labs(x = "Lineage (Pineda, human BA4 + BA9)", y = NULL, title = "Baseline lineage share, ALS-up genes", tag = "a") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
co <- read.delim(file.path(RES, "W66_spatial_context", "post_hoc", "w66_colocalisation_depth_partial.tsv"))
co <- co[co$lineage != "vascular", ]
lin_lab <- c(microglia = "Microglia", astrocyte = "Astrocyte", mural_fibroblast = "Mural/fibroblast", endothelial = "Endothelial",
             oligodendrocyte = "Oligodendrocyte", neuron = "Neuron")
co$lineage <- factor(lin_lab[co$lineage], levels = rev(lin_lab))
co$tissue <- factor(c(LSC = "Lumbar cord", MTC = "Motor cortex")[co$tissue], levels = c("Lumbar cord", "Motor cortex"))
p4b <- ggplot(co, aes(x = rho_depth_partial, y = lineage, colour = tissue)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_point(size = 1.8, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(`Lumbar cord` = PAL$als, `Motor cortex` = PAL$accent), name = NULL) +
  labs(x = "Median Spearman rho with module score\n(per section, depth-partial)", y = "Lineage marker score",
       title = "Spot co-localisation (GSE330224)", tag = "b") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))
d4 <- file.path(V3, "Fig4_cell_source", "source_data"); dir.create(d4, recursive = TRUE, showWarnings = FALSE)
write.table(ls_long, file.path(d4, "panel_a_lineage_share_up.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(co, file.path(d4, "panel_b_colocalisation_depth_partial.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
compose("Fig4_cell_source", "Figure4_cell_source",
        (p4a | (p4b / tag(untag(f7$pc), "c") + plot_layout(heights = c(1, 1)))) + plot_layout(widths = c(1, 1)),
        183, 110, "Fig7_motor_territory")

# ---- Fig. 5: ventral-horn response in two spatial cohorts, and working model ----
model_boxes <- data.frame(
  x = c(2, 1.05, 2.95, 1.05, 2),
  y = c(4.3, 2.95, 2.95, 1.6, 0.35),
  lab = c("ALS: lower motor neuron\nloss in the ventral horn",
          "Focal ventral-horn\nresponse (frozen module):\nGPNMB, ANXA1, LTF,\nEMP3, C3",
          "Broader astrocyte-\ncompartment remodelling\n(consensus signature)",
          "Tracks CD68+ HLA-DR+\nmyeloid cells; baseline\nmicroglial and vascular-\nmeningeal expression",
          "Detected in astrocyte-compartment\nnuclei and in bulk tissue"),
  fill = c("#EFE3DF", "#EAD7A8", "white", "#C9D6E3", "white"))
model_edges <- data.frame(x = c(1.7, 2.3, 1.05, 1.25, 2.9), y = c(3.95, 3.95, 2.3, 1.0, 2.45),
                          xend = c(1.2, 2.8, 1.05, 1.7, 2.35), yend = c(3.5, 3.5, 2.15, 0.7, 0.7))
pm <- ggplot() +
  geom_segment(data = model_edges, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.35, colour = PAL$neutral,
               arrow = arrow(length = unit(1.1, "mm"), type = "closed")) +
  geom_label(data = model_boxes, aes(x = x, y = y, label = lab, fill = I(fill)), size = 1.7, family = "Arial", lineheight = 0.92,
             linewidth = 0.25, label.r = unit(0.6, "mm"), label.padding = unit(1, "mm")) +
  coord_cartesian(xlim = c(0.3, 3.7), ylim = c(-0.05, 4.7), clip = "off") +
  labs(title = "Working model", subtitle = "Not a demonstrated mechanism", tag = "c") +
  theme_void(base_family = "Arial") + theme(plot.title = element_text(size = 6.4, face = "bold"), plot.subtitle = element_text(size = 5.6, colour = "#555555"),
                                             plot.tag = element_text(size = 8, face = "bold"))
compose("Fig5_ventral_horn", "Figure5_ventral_horn",
        (tag(untag(f7$pa), "a") | tag(untag(f7$pb), "b") | pm) + plot_layout(widths = c(1.25, 0.9, 1.3)),
        183, 95, "Fig7_motor_territory")

# ---- Supplementary figures ----
compose("FigS1_discovery_audit", "FigureS1_discovery_audit",
        ((s1$pa | s1$pb) + plot_layout(widths = c(0.7, 1.3))) / ((s1$pc | s1$pd) + plot_layout(widths = c(1, 1.2))) + plot_annotation(tag_levels = "a"),
        183, 140, "FigS1_discovery_audit")
compose("FigS2_contamination", "FigureS2_contamination", (untag(f2$pa) | untag(f2$pb)) + plot_annotation(tag_levels = "a"),
        183, 75, "Fig2_contamination_cell_source")
compose("FigS3_cell_source", "FigureS3_cell_source",
        (s2$pa | (s2$pb / s2$pc + plot_layout(heights = c(1, 1.15)))) + plot_layout(widths = c(0.9, 1.1)) + plot_annotation(tag_levels = "a"),
        183, 150, "FigS2_cell_source")
compose("FigS4_holdout_genes_clinical", "FigureS4_holdout_genes_clinical",
        untag(s3$pa) / (untag(f5$pc1) | untag(f5$pc2)) + plot_layout(heights = c(1.5, 1)) + plot_annotation(tag_levels = list(c("a", "b", ""))),
        183, 190, c("FigS3_holdout_genes", "Fig5_nygc_holdout"))
compose("FigS5_spatial_lumbar_cortex", "FigureS5_spatial_lumbar_cortex",
        (untag(s4$pa) | untag(s4$pb)) / (untag(s4$pc) | untag(s4$pd)) + plot_layout(heights = c(1.15, 1)) + plot_annotation(tag_levels = "a"),
        183, 190, "FigS4_spatial")
compose("FigS6_lumbar_territory", "FigureS6_lumbar_territory",
        (untag(s5$pa) | untag(s5$pb)) + plot_layout(widths = c(2, 1)) + plot_annotation(tag_levels = "a"), 183, 105, "FigS5_motor_selectivity")
compose("FigS7_dose_response", "FigureS7_dose_response", s6$pa, 140, 100, "FigS6_dose_response")
compose("FigS8_evidence_map", "FigureS8_evidence_map", f6$fig, 183, 130, "Fig6_evidence_map")
s9 <- P("FigureS9_microglia_adjusted")
compose("FigS9_microglia_adjusted", "FigureS9_microglia_adjusted", (s9$pa | s9$pb) + plot_layout(widths = c(1.2, 1)) + plot_annotation(tag_levels = "a"),
        183, 95, "FigS9_microglia_adjusted")
s10 <- P("FigureS10_mouse_cell_source")
compose("FigS10_mouse_cell_source", "FigureS10_mouse_cell_source",
        s10$pa / ((s10$pb | s10$pc | s10$pd) + plot_layout(widths = c(1, 1, 1.1))) + plot_layout(heights = c(1, 1.05)) + plot_annotation(tag_levels = "a"),
        183, 140, "FigS10_mouse_cell_source")
cat("Figures_v3 written to", V3, "\n")
