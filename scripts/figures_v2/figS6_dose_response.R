#!/usr/bin/env Rscript
# Supplementary Fig. S6 | Motor-territory module response against a motor-neuron index (W69 A, Amendment 12).

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("FigS6_dose_response")
w69 <- file.path(RES, "W69_dose_response_and_combination")
d <- read.delim(file.path(w69, "w69a_donor_values.tsv"))
s <- read.delim(file.path(w69, "w69a_dose_response_by_cohort.tsv"))
cohort_lab <- c(lumbar_GSE330224 = "Lumbar cord (GSE330224)", cervical_GSE288365 = "Cervical cord (GSE288365)")
set_lab <- c(up = "Frozen ALS-up module", up_consensus = "Consensus up arm")
d$cohort_l <- factor(cohort_lab[d$cohort], levels = cohort_lab); d$set <- factor(set_lab[d$gene_set], levels = set_lab)
d$group <- factor(ifelse(d$als, "ALS", "Control"), levels = c("Control", "ALS"))
s$cohort_l <- factor(cohort_lab[s$cohort], levels = cohort_lab); s$set <- factor(set_lab[s$gene_set], levels = set_lab)
s$lab <- sprintf("ALS: rho = %.2f, one-sided P = %.2f\nrandom-set median rho = %.2f", s$rho, s$p_one_sided, s$random_median_rho)
write_source(d, out, "donor_values"); write_source(s, out, "summary")
p <- ggplot(d, aes(x = mn_index, y = module_response)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_point(aes(colour = group), size = 1.4) +
  geom_text(data = s, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.2, size = 1.6, family = "Arial", lineheight = 0.9) +
  facet_grid(set ~ cohort_l, scales = "free") +
  scale_colour_manual(values = c(Control = PAL$control, ALS = PAL$als), name = NULL) +
  labs(x = "Motor-neuron index (motor-territory mean log2-CPM; lower = more loss)",
       y = "Module response (motor minus sensory territory, log2-CPM)") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), strip.text.y = element_text(angle = 0, hjust = 0))
save_pub_r(p, out, "FigureS6_dose_response", width_mm = 140, height_mm = 100)
cat("Supplementary Fig. S6 written to", out, "\n")
