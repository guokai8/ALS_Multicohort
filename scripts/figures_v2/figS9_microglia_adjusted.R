#!/usr/bin/env Rscript
# Supplementary Fig. S9 | Territory x disease interaction with and without adjustment for microglial content (W70, Amendment 13).

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("FigS9_microglia_adjusted")
w70 <- file.path(RES, "W70_microglia_adjusted_territory")
cohort_lab <- c(lumbar_GSE330224 = "Lumbar cord", cervical_GSE288365 = "Cervical cord")
set_lab <- c(up = "ALS-up module", up_consensus = "Consensus up genes")
adj_lab <- c(none = "Unadjusted", microglia = "All (primary)", homeostatic = "Homeostatic", activated = "Activated")

t <- read.delim(file.path(w70, "w70_adjusted_interaction_tests.tsv"))
un <- unique(t[, c("cohort", "gene_set", "unadjusted_fry_p", "unadjusted_median_logFC")])
names(un)[3:4] <- c("fry_p", "median_logFC_interaction"); un$score <- "none"
d <- rbind(un, t[, c("cohort", "gene_set", "fry_p", "median_logFC_interaction", "score")])
d$adjust <- factor(adj_lab[d$score], levels = adj_lab)
d$cohort_l <- factor(cohort_lab[d$cohort], levels = cohort_lab); d$set <- factor(set_lab[d$gene_set], levels = set_lab)
d$lab <- formatC(d$fry_p, format = "fg", digits = 2)
write_source(d, out, "panel_a_interaction")
pa <- ggplot(d, aes(x = adjust, y = median_logFC_interaction, fill = adjust == "Unadjusted")) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_col(width = 0.65) +
  geom_text(aes(label = lab), vjust = -0.4, size = 1.6, family = "Arial") +
  facet_grid(set ~ cohort_l, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  scale_fill_manual(values = c(`TRUE` = PAL$neutral, `FALSE` = PAL$als), guide = "none") +
  labs(x = "Adjustment for microglial marker score", y = "Interaction median log2FC\n(motor minus sensory, ALS vs control)",
       title = "Territory x disease interaction (labels: fry P)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), strip.text.y = element_text(angle = 0, hjust = 0))

m <- read.delim(file.path(w70, "w70_microglial_score_by_territory.tsv"))
m$cohort_l <- factor(cohort_lab[m$cohort], levels = cohort_lab)
m$score_l <- factor(adj_lab[m$score], levels = adj_lab[-1])
m$territory <- factor(ifelse(m$territory == "motor", "Motor", "Sensory"), levels = c("Motor", "Sensory"))
m$lab <- formatC(m$wilcoxon_p, format = "fg", digits = 2)
write_source(m, out, "panel_b_microglial_score")
pb <- ggplot(m, aes(x = score_l, y = mean_ALS_minus_control, fill = territory)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = lab, vjust = ifelse(mean_ALS_minus_control >= 0, -0.4, 1.3)), position = position_dodge(width = 0.75), size = 1.5, family = "Arial") +
  facet_grid(. ~ cohort_l) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_fill_manual(values = c(Motor = PAL$als, Sensory = PAL$control), name = NULL) +
  labs(x = "Marker score", y = "Mean score, ALS minus control\n(log2-CPM)", title = "Microglial score by territory (labels: Wilcoxon P)") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), axis.text.x = element_text(angle = 30, hjust = 1))

save_panels(list(pa = pa, pb = pb), "FigureS9_microglia_adjusted")
fig <- (pa | pb) + plot_layout(widths = c(1.2, 1)) + plot_annotation(tag_levels = "a")
save_pub_r(fig, out, "FigureS9_microglia_adjusted", width_mm = 183, height_mm = 95)
cat("Supplementary Fig. S9 written to", out, "\n")
