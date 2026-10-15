#!/usr/bin/env Rscript
# Fig. 5 | NYGC bulk hold-out: frozen modules, frozen consensus signature and clinical association.

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("Fig5_NYGC_bulk_holdout")
W55 <- file.path(RES, "W55_NYGC_bulk_holdout")
region_lab <- c(motor_cortex = "Motor cortex", spinal_cord = "Spinal cord")
cols <- c("region", "module", "n_samples", "n_als_donors", "n_control_donors", "n_module_genes_tested", "n_genes_concordant",
          "median_logFC", "fry_direction", "fry_p", "camera_p")
fmt_p <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), formatC(p, format = "fg", digits = 2))

# Shared encoding: x = analysis, y = -log10 fry P, filled = confirmatory, open = exploratory/sensitivity.
p_dotplot <- function(df, title, colour_values) {
  ggplot(df, aes(x = analysis, y = neglog10(fry_p), colour = arm, group = arm)) +
    geom_hline(yintercept = -log10(0.05), linetype = "22", linewidth = 0.3, colour = PAL$neutral) +
    geom_line(linewidth = 0.3) +
    geom_point(aes(shape = role), size = 1.7, fill = "white", stroke = 0.5) +
    geom_label(data = df[!is.na(df$label), ], aes(label = label), vjust = -0.55, hjust = 0.3, size = 1.7, family = "Arial",
               fill = "white", linewidth = 0, label.padding = unit(0.4, "mm"), show.legend = FALSE) +
    facet_wrap(~ region, nrow = 1) +
    scale_colour_manual(values = colour_values, name = NULL) +
    scale_shape_manual(values = c(Confirmatory = 16, Exploratory = 21, Sensitivity = 24), name = NULL, drop = FALSE) +
    scale_x_discrete(expand = expansion(add = 0.45)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
    labs(x = NULL, y = expression(-log[10]~italic(P)[fry]), title = title) +
    guides(colour = guide_legend(order = 1, nrow = 2), shape = guide_legend(order = 2, nrow = 3)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom", legend.box = "horizontal",
          legend.key.size = unit(2.4, "mm"), legend.spacing.y = unit(0.3, "mm"), legend.margin = margin(0, 6, 0, 0))
}

# ---- a: frozen modules (W55 primary, W55 no-PMI sensitivity, W56 E3a myeloid + vascular adjustment) ----
prim <- read.delim(file.path(W55, "w55_primary_module_tests.tsv"))
nopmi <- read.delim(file.path(W55, "w55_sensitivity_no_pmi_module_tests.tsv"))
e3 <- read.delim(file.path(RES, "W56_cell_source_check", "w56_e3_nygc_composition_retest.tsv"))
e3a <- e3[e3$analysis == "E3a_myeloid_vascular_adjusted", ]
a <- rbind(cbind(analysis = "Primary", role = "Confirmatory", prim[, cols], fry_BH = prim$fry_BH),
           cbind(analysis = "No PMI", role = "Sensitivity", nopmi[, cols], fry_BH = nopmi$fry_BH),
           cbind(analysis = "+ myeloid/vascular", role = "Exploratory", e3a[, cols], fry_BH = NA))
a$analysis <- factor(a$analysis, levels = c("Primary", "No PMI", "+ myeloid/vascular"))
a$role <- factor(a$role, levels = c("Confirmatory", "Exploratory", "Sensitivity"))
a$region <- factor(region_lab[a$region], levels = region_lab)
a$arm <- factor(ifelse(a$module == "up", "ALS-up module", "ALS-down module"), levels = c("ALS-up module", "ALS-down module"))
a$label <- ifelse(a$analysis == "Primary", paste0("FDR ", fmt_p(a$fry_BH)), NA)
write_source(a, out, "panel_a_frozen_module_tests")
pa <- p_dotplot(a, "ALS-up and ALS-down modules", c("ALS-up module" = PAL$up, "ALS-down module" = PAL$down))

# ---- b: frozen consensus signature (W59) ----
w59 <- read.delim(list.files(RES, pattern = "^w59_signature_tests.tsv$", recursive = TRUE, full.names = TRUE)[1])
an_lab <- c("1_primary_W55_design" = "Primary", "2_secondary_myeloid_vascular" = "+ myeloid/vascular",
            "3_sensitivity_signature_cohort_overlap_excluded" = "Overlap excluded")
b <- cbind(analysis = an_lab[w59$analysis], w59[, cols], fry_BH = w59$fry_BH)
b$role <- factor(ifelse(b$analysis == "Primary", "Confirmatory", ifelse(b$analysis == "Overlap excluded", "Sensitivity", "Exploratory")),
                 levels = c("Confirmatory", "Exploratory", "Sensitivity"))
b$analysis <- factor(b$analysis, levels = an_lab)
b$region <- factor(region_lab[b$region], levels = region_lab)
b$arm <- factor(ifelse(b$module == "up_consensus", "Consensus up (85)", "Consensus down (73)"),
                levels = c("Consensus up (85)", "Consensus down (73)"))
b$label <- ifelse(b$analysis == "Primary", paste0("FDR ", fmt_p(b$fry_BH)), NA)
write_source(b, out, "panel_b_consensus_signature_tests")
pb <- p_dotplot(b, "Consensus signature", c("Consensus up (85)" = PAL$up, "Consensus down (73)" = PAL$down))

# ---- c: clinical association within ALS donors (W55 clinical models) ----
coef <- read.delim(file.path(W55, "w55_clinical_coefficients_with_ci.tsv"))
assoc <- read.delim(file.path(W55, "w55_clinical_association.tsv"))
coef <- merge(coef, assoc[, c("region", "module", "lm_BH")], by = c("region", "module"))
coef$row <- factor(paste0(region_lab[coef$region], ", ", ifelse(coef$module == "up", "ALS-up", "ALS-down")),
                   levels = rev(c("Motor cortex, ALS-up", "Motor cortex, ALS-down", "Spinal cord, ALS-up", "Spinal cord, ALS-down")))
coef$arm <- ifelse(coef$module == "up", "up", "down")
coef$label <- paste0("FDR ", fmt_p(coef$lm_BH))
write_source(coef, out, "panel_c1_clinical_coefficients")
pc1 <- ggplot(coef, aes(x = beta_per_sd, y = row, colour = arm)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = PAL$neutral) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y", linewidth = 0.4) +
  geom_point(size = 1.6) +
  geom_text(aes(x = 0.37, label = label), hjust = 0, size = 1.7, family = "Arial", colour = PAL$ink) +
  scale_colour_manual(values = c(up = PAL$up, down = PAL$down), guide = "none") +
  scale_x_continuous(breaks = c(-0.3, 0, 0.3)) +
  coord_cartesian(xlim = c(-0.4, 0.36), clip = "off") +
  labs(x = "log disease duration per SD\nof module score (95% CI)", y = NULL, title = "Score vs disease duration") +
  theme(plot.margin = margin(4, 30, 4, 4))

don <- read.delim(file.path(W55, "w55_clinical_donor_scores.tsv"))
sc <- don[don$region == "spinal_cord" & don$module == "down", ]
write_source(sc, out, "panel_c2_spinal_down_donor_scores")
pc2 <- ggplot(sc, aes(x = score_sd, y = duration_months)) +
  geom_point(size = 0.8, colour = PAL$down, alpha = 0.75) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = PAL$ink, fill = PAL$light, linewidth = 0.4) +
  scale_y_log10(breaks = c(6, 12, 24, 48, 96, 192)) +
  labs(x = "Spinal ALS-down score (SD units)", y = "Disease duration (months, log scale)",
       title = sprintf("Spinal cord, ALS-down (n = %d)", nrow(sc)))

fig <- ((pa | pb) / (pc1 | pc2)) + plot_layout(heights = c(1.25, 1)) + plot_annotation(tag_levels = list(c("a", "b", "c", "")))
save_panels(list(pa = pa, pb = pb, pc1 = pc1, pc2 = pc2), "Figure5_NYGC_bulk_holdout")
save_pub_r(fig, out, "Figure5_NYGC_bulk_holdout", width_mm = 183, height_mm = 140)
cat("Fig. 5 written to", out, "\n")
