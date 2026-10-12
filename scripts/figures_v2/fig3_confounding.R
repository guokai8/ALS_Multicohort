#!/usr/bin/env Rscript
# Fig. 3 | Each discordant single-nucleus cohort carries an identifiable design confounder.

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("Fig3_cohort_confounding")
W51 <- file.path(RES, "W51_Pineda_BA4_discovery_matched_limma")
MODULE_COL <- c("ALS-up" = PAL$up, "ALS-down" = PAL$down)
module_lab <- function(m) factor(ifelse(m == "up", "ALS-up", "ALS-down"), levels = c("ALS-up", "ALS-down"))

# ---- a: Pineda processing date (sample-folder prefix) by diagnosis ----
pref <- read.delim(file.path(W51, "w51_pineda_prefix_by_condition.tsv"), colClasses = c(prefix = "character"))
a <- do.call(rbind, lapply(c("ALS", "FTLD", "PN"), function(cond) data.frame(region = pref$Region, prefix = pref$prefix,
  condition = cond, n = pref[[cond]])))
a$condition <- factor(a$condition, levels = c("ALS", "FTLD", "PN"), labels = c("ALS", "FTLD", "Control"))
a$prefix <- factor(a$prefix, levels = rev(sort(unique(a$prefix))))
a$region <- factor(a$region, levels = c("BA4", "BA9"), labels = c("BA4 (motor)", "BA9 (prefrontal)"))
# Cells that form the analysed contrasts: BA4 2019 campaign (191112 ALS vs 191114 control) and BA9 batch 210526.
a$analysed <- (a$region == "BA4 (motor)" & ((a$prefix == "191112" & a$condition == "ALS") | (a$prefix == "191114" & a$condition == "Control"))) |
  (a$region == "BA9 (prefrontal)" & a$prefix == "210526" & a$condition %in% c("ALS", "Control"))
write_source(a, out, "panel_a_pineda_prefix_by_condition")
pa <- ggplot(a, aes(x = condition, y = prefix)) +
  geom_tile(aes(fill = n), colour = "white", linewidth = 0.4) +
  geom_tile(data = a[a$analysed, ], fill = NA, colour = PAL$ink, linewidth = 0.45, width = 0.92, height = 0.88) +
  geom_text(aes(label = ifelse(n > 0, n, ""), colour = n > 12), size = 1.9, family = "Arial", show.legend = FALSE) +
  facet_grid(region ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "#F2F2F2", high = PAL$neutral, guide = "none") +
  scale_colour_manual(values = c(`FALSE` = PAL$ink, `TRUE` = "white")) +
  labs(x = NULL, y = "Processing date (yymmdd)", title = "Pineda: processing date by diagnosis") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), strip.text.y = element_text(angle = 0, hjust = 0),
        legend.position = "bottom", legend.key.height = unit(1.6, "mm"), legend.key.width = unit(4, "mm"))

# ---- b: variance explained by processing date within diagnosis ----
bv <- read.delim(file.path(W51, "w51_post_hoc_batch_verification.tsv"))
b_lab <- c("BA4 ALS, all dates", "BA4 C9ALS,\n191112 vs 210910", "BA4 sALS,\n191112 vs 201019", "BA4 FTLD, all dates", "BA9 ALS,\n210526 vs 210910")
bv$label <- factor(b_lab, levels = rev(b_lab))
bv$p_label <- paste0("P = ", formatC(bv$perm_p, format = "fg", digits = 2))
write_source(bv, out, "panel_b_within_diagnosis_date_R2")
pb <- ggplot(bv, aes(x = prefix_R2_top5PC, y = label)) +
  geom_col(aes(fill = perm_p < 0.05), width = 0.62) +
  geom_text(aes(label = p_label), hjust = -0.1, size = 1.8, family = "Arial") +
  scale_fill_manual(values = c(`TRUE` = PAL$ink, `FALSE` = PAL$light), guide = "none") +
  scale_x_continuous(limits = c(0, 0.8), breaks = seq(0, 0.6, 0.2), expand = expansion(mult = c(0, 0))) +
  labs(x = expression(R^2~"of processing date (top 5 PCs)"), y = NULL, title = "Date effect within diagnosis")

# ---- c: Ulm age imbalance and age-bias correction ----
ages <- read.delim(file.path(dirname(FIG_ROOT), "Other_Datasets", "Multiome_Dataset", "files", "Multiome_Dataset_Samples_Summary.txt"))
ages <- ages[ages$Case %in% c("ALS", "Control"), ]
ages$Age <- suppressWarnings(as.numeric(ages$Age))
ages <- ages[!is.na(ages$Age), ]
ages$group <- factor(ages$Case, levels = c("Control", "ALS"))
age_med <- aggregate(Age ~ group, data = ages, FUN = mean)
write_source(ages[, c("Case", "Sex", "Age")], out, "panel_c1_ulm_age_by_group")
set.seed(20261012)
pc1 <- ggplot(ages, aes(x = group, y = Age, colour = group)) +
  geom_jitter(width = 0.18, height = 0, size = 0.7, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = PAL$ink) +
  annotate("text", x = 1.5, y = max(ages$Age, na.rm = TRUE) + 4,
           label = sprintf("Mean %.1f vs %.1f y", age_med$Age[age_med$group == "Control"], age_med$Age[age_med$group == "ALS"]),
           size = 1.8, family = "Arial") +
  scale_colour_manual(values = c(Control = PAL$control, ALS = PAL$als), guide = "none") +
  labs(x = NULL, y = "Age at death (years)", title = "Ulm: age")

eff <- read.delim(file.path(RES, "W58_multicohort_meta", "w58_module_effects_by_cohort.tsv"))
ulm <- eff[eff$cohort %in% c("Ulm_M1", "Ulm_M1_agecorrected"), ]
ulm$version <- factor(ifelse(ulm$cohort == "Ulm_M1", "Raw", "Age-\ncorrected"), levels = c("Raw", "Age-\ncorrected"))
ulm$module <- module_lab(ulm$module)
write_source(ulm[, c("cohort", "module", "d", "se_d")], out, "panel_c2_ulm_age_correction")
dodge <- position_dodge(width = 0.45)
pc2 <- ggplot(ulm, aes(x = version, y = d, colour = module, group = module)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = PAL$light) +
  geom_line(position = dodge, linewidth = 0.3) +
  geom_pointrange(aes(ymin = d - 1.96 * se_d, ymax = d + 1.96 * se_d), position = dodge, size = 0.18, linewidth = 0.35) +
  scale_colour_manual(values = MODULE_COL, name = NULL) +
  labs(x = NULL, y = "Standardised module effect (d)", title = "Ulm: module effect") +
  theme(legend.position = "none")

# ---- d: Li 2023 depth imbalance and depth adjustment ----
li <- read.delim(file.path(RES, "W54_Li2023_motor_cortex_validation", "w54_li2023_astrocyte_pseudobulk_samples.tsv"))
li <- li[li$disease %in% c("C9-ALS", "Control"), ]
li$group <- factor(li$disease, levels = c("Control", "C9-ALS"), labels = c("Control", "C9-ALS"))
write_source(li[, c("GSM", "donor", "disease", "n_astrocytes", "astro_median_umi")], out, "panel_d1_li_astrocyte_depth")
li_med <- aggregate(astro_median_umi ~ group, data = li, FUN = median)
pd1 <- ggplot(li, aes(x = group, y = astro_median_umi / 1000, colour = group)) +
  geom_jitter(width = 0.15, height = 0, size = 0.9, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = PAL$ink) +
  scale_colour_manual(values = c(Control = PAL$control, `C9-ALS` = PAL$als), guide = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Median UMI per astrocyte (x1,000)", title = "Li: depth")

dep <- read.delim(file.path(RES, "W54_Li2023_motor_cortex_validation", "w54_post_hoc_depth_sensitivity.tsv"))
dep <- rbind(data.frame(model = dep$formula, module = dep$module, test = "fry", p = dep$fry_p),
             data.frame(model = dep$formula, module = dep$module, test = "camera", p = dep$camera_p))
dep$model <- factor(ifelse(grepl("depth", dep$model), "+ depth", "Primary"), levels = c("Primary", "+ depth"))
dep$module <- module_lab(dep$module)
write_source(dep, out, "panel_d2_li_depth_adjustment")
pd2 <- ggplot(dep, aes(x = model, y = neglog10(p), colour = module, shape = test, group = interaction(module, test))) +
  geom_hline(yintercept = -log10(0.05), linetype = "22", linewidth = 0.3, colour = PAL$neutral) +
  geom_line(linewidth = 0.3) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = MODULE_COL, name = NULL) +
  scale_shape_manual(values = c(fry = 16, camera = 1), name = NULL) +
  labs(x = NULL, y = expression(-log[10]~italic(P)), title = "Li: depth adjustment") +
  theme(legend.position = "right", legend.key.size = unit(2.2, "mm"), legend.spacing.y = unit(0.5, "mm"))

# ---- e: cross-cohort reproducibility of top-200 genes (W58) ----
rep <- read.delim(file.path(RES, "W58_multicohort_meta", "w58_reproducibility_matrix.tsv"))
cohort_short <- c(NYGC_MTC_discovery = "Discovery", Li2023_BA4 = "Li", Ulm_M1 = "Ulm", Zelic_cervicalSC = "Zelic",
                  Pineda_BA4_2019 = "Pineda BA4", Pineda_BA9_210526 = "Pineda BA9")
rep$source <- factor(cohort_short[rep$source], levels = cohort_short)
rep$target <- factor(cohort_short[rep$target], levels = cohort_short)
rep$star <- ifelse(rep$empirical_p < 0.05, "*", "")
write_source(rep, out, "panel_e_reproducibility_matrix")
pe <- ggplot(rep, aes(x = target, y = source, fill = replicated_fraction)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = paste0(round(100 * replicated_fraction), star), colour = replicated_fraction > 0.3),
            size = 1.8, family = "Arial", show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = PAL$up, limits = c(0, 0.6), name = "Replicated\nfraction") +
  scale_colour_manual(values = c(`FALSE` = PAL$ink, `TRUE` = "white")) +
  scale_y_discrete(limits = rev) +
  labs(x = "Replication cohort", y = "Top-200 source cohort", title = "Cross-cohort replication of top genes (%)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key.width = unit(2, "mm"))

# Tags are added per panel so that c and d each span two sub-plots.
pa <- pa + labs(tag = "a"); pb <- pb + labs(tag = "b"); pe <- pe + labs(tag = "e")
pc1 <- pc1 + labs(tag = "c"); pd1 <- pd1 + labs(tag = "d")
top <- (pa | pb | pe) + plot_layout(widths = c(0.75, 1, 1.15))
bottom <- (pc1 | pc2 | pd1 | pd2) + plot_layout(widths = c(0.8, 0.8, 0.8, 1.25))
fig <- top / bottom + plot_layout(heights = c(1.25, 1))
save_pub_r(fig, out, "Figure3_cohort_confounding", width_mm = 183, height_mm = 125)
cat("Fig. 3 written to", out, "\n")
