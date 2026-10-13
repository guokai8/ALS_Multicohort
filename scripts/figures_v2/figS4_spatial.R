#!/usr/bin/env Rscript
# Supplementary Fig. S4 | Spatial context in GSE330224 Visium (W66, Amendment 9; exploratory).

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("FigS4_spatial")
w66 <- file.path(RES, "W66_spatial_context")
tissue_lab <- c(LSC = "Lumbar spinal cord", MTC = "Motor cortex")

# ---- a: disease shift by unit (independent donors) ----
st <- read.delim(file.path(w66, "w66_set_tests.tsv"))
st <- st[st$analysis == "primary", ]
set_lab <- c(up = "Frozen ALS-up", up_consensus = "Consensus up", down_consensus = "Consensus down")
st <- st[st$gene_set %in% names(set_lab), ]
unit_order <- c("whole", "grey", "white", "edge", "L1", "L2_3", "L5_6", "WM",
                "astrocyte_only", "astrocyte_rich", "vascular_rich", "microglia_rich")
unit_lab <- c(whole = "Whole section", grey = "Grey matter", white = "White matter", edge = "Edge (pial)",
              L1 = "Layer 1", L2_3 = "Layers 2-3", L5_6 = "Layers 5-6", WM = "White matter ",
              astrocyte_only = "Astrocyte-only niche", astrocyte_rich = "Astrocyte-rich niche",
              vascular_rich = "Vascular-rich niche", microglia_rich = "Microglia-rich niche")
st$unit <- factor(unit_lab[st$stratum], levels = rev(unit_lab[unit_order]))
st$set <- factor(set_lab[st$gene_set], levels = set_lab)
st$tissue <- factor(st$tissue, levels = c("LSC", "MTC"))
st$sig <- ifelse(st$fry_p < 0.05, "fry P < 0.05", "fry P >= 0.05")
st$p_lab <- formatC(st$fry_p, format = "fg", digits = 2)
write_source(st[, c("tissue", "stratum", "gene_set", "n_als", "n_control", "median_logFC", "fry_direction", "fry_p", "camera_p")],
             out, "panel_a_disease_shift")
pa <- ggplot(st, aes(x = median_logFC, y = unit)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_point(aes(shape = sig, colour = set), size = 1.5, stroke = 0.5) +
  geom_text(aes(label = p_lab, colour = set), nudge_y = 0.38, size = 1.55, family = "Arial", show.legend = FALSE) +
  facet_grid(tissue ~ set, scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(-0.3, 0.42), breaks = c(-0.2, 0, 0.2, 0.4)) +
  scale_shape_manual(values = c(`fry P < 0.05` = 16, `fry P >= 0.05` = 1), name = NULL) +
  scale_colour_manual(values = c(`Frozen ALS-up` = PAL$up, `Consensus up` = PAL$accent, `Consensus down` = PAL$down), guide = "none") +
  labs(x = "Median log2 fold change, ALS vs control (independent donors)", y = NULL,
       title = "Disease shift in donor pseudobulk") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), strip.text.y = element_text(angle = 0, hjust = 0),
        panel.spacing.x = unit(2.5, "mm"))

# ---- b: region share of ALS-up genes ----
rs <- read.delim(file.path(w66, "w66_region_share.tsv"))
region_lab <- c(edge = "Edge", grey = "Grey", white = "White", L1 = "L1", L2_3 = "L2-3", L5_6 = "L5-6", WM = "WM")
rs$region <- factor(region_lab[rs$region], levels = region_lab)
rs$tissue <- factor(c(LSC = "Spinal cord", MTC = "Motor cortex")[rs$tissue], levels = c("Spinal cord", "Motor cortex"))
gene_order <- with(rs[rs$tissue == "Spinal cord" & rs$region == "Edge", ], symbol[order(share)])
rs$symbol <- factor(rs$symbol, levels = gene_order)
write_source(rs, out, "panel_b_region_share")
pb <- ggplot(rs, aes(x = region, y = symbol, fill = share)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
  scale_fill_gradient(low = "white", high = PAL$ink, limits = c(0, 1), name = "Share of\nsummed CPM") +
  labs(x = NULL, y = NULL, title = "Region share of ALS-up genes") +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

# ---- c: spot co-localisation with lineage markers ----
co <- read.delim(file.path(w66, "w66_colocalisation_summary.tsv"))
lin_lab <- c(astrocyte = "Astrocyte", microglia = "Microglia", endothelial = "Endothelial", mural_fibroblast = "Mural/fibroblast",
             oligodendrocyte = "Oligodendrocyte", neuron = "Neuron")
co <- co[co$lineage %in% names(lin_lab), ]
co$lineage <- factor(lin_lab[co$lineage], levels = lin_lab)
co$feature <- sub("^gene_", "", co$feature)
co$feature[co$feature == "module_up"] <- "ALS-up module"
co$feature <- factor(co$feature, levels = c(rev(gene_order), "ALS-up module"))
co$tissue <- factor(tissue_lab[co$tissue], levels = tissue_lab)
write_source(co, out, "panel_c_colocalisation")
pc <- ggplot(co, aes(x = lineage, y = feature, fill = median_rho)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(. ~ tissue) +
  scale_fill_gradient2(low = PAL$down, mid = "white", high = PAL$up, midpoint = 0, limits = c(-0.3, 0.3), oob = scales::squish,
                       name = "Median\nSpearman rho") +
  scale_y_discrete(limits = rev) +
  labs(x = "Lineage marker score", y = NULL, title = "Spot co-localisation (spots >= 500 counts)") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())

# ---- d: ligand-receptor co-detection ----
lr <- read.delim(file.path(w66, "w66_lr_codetection_by_donor.tsv"))
lr <- lr[lr$gse330130_donor %in% c(FALSE, "False"), ]
lr$tissue <- factor(tissue_lab[lr$tissue], levels = tissue_lab)
lr$group <- factor(lr$group, levels = c("Control", "ALS"))
pair_order <- with(aggregate(log2_obs_exp ~ pair, lr, median), pair[order(log2_obs_exp)])
lr$pair <- factor(lr$pair, levels = pair_order)
nul <- read.delim(file.path(w66, "post_hoc", "w66_lr_codetection_empirical_null.tsv"))
nul$tissue <- factor(tissue_lab[nul$tissue], levels = tissue_lab)
nul$pair <- factor(nul$pair, levels = pair_order)
nul$lab <- paste0("P = ", formatC(nul$empirical_p, format = "fg", digits = 2))
write_source(lr, out, "panel_d_lr_codetection")
write_source(nul, out, "panel_d_empirical_null")
set.seed(20261013)
pd_ <- ggplot(lr, aes(x = log2_obs_exp, y = pair, colour = group)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_jitter(height = 0.18, width = 0, size = 0.55, alpha = 0.75) +
  geom_point(data = nul, aes(x = random_median, y = pair), inherit.aes = FALSE, shape = 4, size = 1.3, colour = PAL$neutral) +
  stat_summary(fun = median, geom = "point", shape = 124, size = 3, colour = PAL$ink) +
  geom_text(data = nul, aes(x = 2.35, y = pair, label = lab, fontface = ifelse(empirical_p < 0.05, "bold", "plain")),
            inherit.aes = FALSE, size = 1.6, family = "Arial", hjust = 1, colour = PAL$ink) +
  facet_grid(. ~ tissue) +
  coord_cartesian(xlim = c(-0.6, 2.4)) + scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5)) +
  scale_colour_manual(values = c(Control = PAL$control, ALS = PAL$als), name = NULL) +
  labs(x = "Same-spot co-detection, log2 observed / expected", y = NULL,
       title = "Ligand-receptor co-detection (x, matched random pairs)") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))

fig <- (pa | pb) / (pc | pd_) + plot_layout(heights = c(1.15, 1)) + plot_annotation(tag_levels = "a")
fig[[1]] <- fig[[1]] + plot_layout(widths = c(1.55, 1))
save_panels(list(pa = pa, pb = pb, pc = pc, pd = pd_), "FigureS4_spatial")
save_pub_r(fig, out, "FigureS4_spatial", width_mm = 183, height_mm = 190)
cat("Supplementary Fig. S4 written to", out, "\n")
