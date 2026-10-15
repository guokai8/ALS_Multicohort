#!/usr/bin/env Rscript
# Supplementary Fig. S3 | Gene-level fold changes of the frozen modules in the NYGC bulk hold-out (W55 primary design)
# and after exploratory activated-myeloid and vascular adjustment (W56 E3a).

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("FigS3_holdout_genes")

sym <- read.delim(file.path(RES, "W56_cell_source_check", "w56_e1_module_gene_specificity.tsv"))[, c("gene_id", "common_name")]
prim <- read.delim(file.path(RES, "W55_NYGC_bulk_holdout", "w55_module_gene_results.tsv"))
prim <- prim[prim$pmi_in_model %in% c(TRUE, "TRUE"), ]
prim$design <- "Primary (confirmatory)"
e3a <- read.delim(file.path(RES, "W56_cell_source_check", "w56_e3a_module_gene_results.tsv"))
e3a$design <- "+ myeloid/vascular (exploratory)"
cols <- c("region", "gene_id", "module", "logFC", "P.Value", "design")
g <- merge(rbind(prim[, cols], e3a[, cols]), sym, by = "gene_id")
g$region <- factor(c(motor_cortex = "Motor cortex", spinal_cord = "Spinal cord")[g$region], levels = c("Motor cortex", "Spinal cord"))
g$module <- factor(ifelse(g$module == "up", "ALS-up module", "ALS-down module"), levels = c("ALS-up module", "ALS-down module"))
g$design <- factor(g$design, levels = c("Primary (confirmatory)", "+ myeloid/vascular (exploratory)"))
g$nominal <- g$P.Value < 0.05

# Order genes within each module by the primary spinal-cord fold change.
ord <- g[g$design == "Primary (confirmatory)" & g$region == "Spinal cord", ]
ord <- ord[order(ord$module, ord$logFC), "common_name"]
g$gene <- factor(g$common_name, levels = unique(c(ord, setdiff(unique(g$common_name), ord))))
write_source(g[, c("region", "module", "common_name", "gene_id", "design", "logFC", "P.Value", "nominal")], out, "gene_level_fold_changes")

wide <- reshape(g[, c("region", "module", "gene", "design", "logFC")], idvar = c("region", "module", "gene"),
                timevar = "design", direction = "wide")
names(wide) <- sub("^logFC\\.", "", names(wide))
pa <- ggplot(g, aes(x = logFC, y = gene)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_segment(data = wide, aes(x = `Primary (confirmatory)`, xend = `+ myeloid/vascular (exploratory)`, y = gene, yend = gene),
               inherit.aes = FALSE, linewidth = 0.35, colour = PAL$neutral) +
  # Exploratory (hollow) first, primary (filled) on top, so overlapping estimates stay visible.
  geom_point(data = g[g$design != "Primary (confirmatory)", ], aes(shape = design, colour = module, alpha = nominal), size = 1.6, stroke = 0.5) +
  geom_point(data = g[g$design == "Primary (confirmatory)", ], aes(shape = design, colour = module, alpha = nominal), size = 1.3) +
  facet_grid(module ~ region, scales = "free_y", space = "free_y") +
  scale_shape_manual(values = c("Primary (confirmatory)" = 16, "+ myeloid/vascular (exploratory)" = 1), name = NULL) +
  scale_colour_manual(values = c("ALS-up module" = PAL$up, "ALS-down module" = PAL$down), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), labels = c(`TRUE` = "P < 0.05", `FALSE` = "P ≥ 0.05"), name = NULL) +
  labs(x = expression("ALS vs control log"[2]*" fold change (NYGC bulk)"), y = NULL,
       title = "Module genes in the bulk hold-out") +
  theme(legend.position = "bottom", strip.text.y = element_text(angle = 0, hjust = 0))

save_panels(list(pa = pa), "FigureS3_holdout_genes")
save_pub_r(pa, out, "FigureS3_holdout_genes", width_mm = 183, height_mm = 130)
cat("Fig. S3 written to", out, "\n")
