#!/usr/bin/env Rscript
# Fig. 4 | Random-effects meta-analysis of module effects and discovery-agnostic consensus genes.

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("Fig4_meta_analysis")
W58 <- file.path(RES, "W58_multicohort_meta")
W60 <- file.path(RES, "W60_age70_sensitivity")

# ---- a: forest plot of standardised module effects ----
eff <- read.delim(file.path(W58, "w58_module_effects_by_cohort.tsv"))
meta <- read.delim(file.path(W58, "w58_module_meta.tsv"))
cohort_lab <- c(Li2023_BA4 = "Li 2023 (C9, motor)", Ulm_M1 = "Ulm (raw)", Ulm_M1_agecorrected = "Ulm (age-corrected)",
                Ulm_M1_age70 = "Ulm (age-corrected, \u226570 slope)",
                Pineda_BA4_2019 = "Pineda BA4 (2019 campaign)", Pineda_BA9_210526 = "Pineda BA9 (prefrontal)",
                Zelic_cervicalSC = "Zelic (cervical cord)")
pool_lab <- c(primary_motor_independent = "Primary pool: motor cortex (k = 3)",
              sens_motor_ulm_uncorrected = "Sensitivity: Ulm raw (k = 3)",
              sens_motor_ulm_age70 = "Sensitivity: Ulm \u226570 slope (k = 3)",
              secondary_all_independent = "Secondary pool: all contrasts (k = 5)")
rows <- c(cohort_lab, pool_lab)
# W60 (Amendment 6): Ulm corrected with beta_age from donors aged >= 70, and the re-run primary pool.
u70 <- read.delim(file.path(W60, "w60_ulm_corrected_effects.tsv"))
u70 <- u70[u70$estimate == "primary_age70_BA4_BA9", ]
m70 <- read.delim(file.path(W60, "w60_module_meta.tsv"))
m70 <- m70[m70$estimate == "primary_age70_BA4_BA9" & m70$pool == "primary_motor_independent", ]
eff <- rbind(eff[, c("cohort", "module", "d", "se_d")],
             data.frame(cohort = "Ulm_M1_age70", module = u70$module, d = u70$d_corrected, se_d = u70$se_corrected))
meta <- rbind(meta[, c("analysis", "module", "pooled_d", "ci_low", "ci_high", "p", "I2")],
              data.frame(analysis = "sens_motor_ulm_age70", module = m70$module, pooled_d = m70$pooled_d,
                         ci_low = m70$ci_low, ci_high = m70$ci_high, p = m70$p, I2 = m70$I2))
coh <- eff[eff$cohort %in% names(cohort_lab), ]
f_coh <- data.frame(key = coh$cohort, module = coh$module, d = coh$d, lo = coh$d - 1.96 * coh$se_d, hi = coh$d + 1.96 * coh$se_d,
                    kind = ifelse(coh$cohort %in% c("Ulm_M1", "Ulm_M1_age70"), "cohort_excluded", "cohort"), note = "")
pm <- meta[meta$analysis %in% names(pool_lab), ]
# W61 (Amendment 7): Knapp-Hartung P for each pool, shown next to the prespecified z-test P.
w61 <- read.delim(file.path(RES, "W61_small_k_meta_inference", "w61_small_k_meta_inference.tsv"))
w61$analysis <- ifelse(w61$pool == "w60_primary_motor_age70", "sens_motor_ulm_age70", w61$pool)
pm$hksj_p <- w61$hksj_p[match(paste(pm$analysis, pm$module), paste(w61$analysis, w61$module))]
if (anyNA(pm$hksj_p)) stop("HKSJ P missing for a plotted pool")
f_pool <- data.frame(key = pm$analysis, module = pm$module, d = pm$pooled_d, lo = pm$ci_low, hi = pm$ci_high,
                     kind = ifelse(startsWith(pm$analysis, "sens_"), "pool_sens", "pool"),
                     note = sprintf("P = %s (HKSJ %s), I² = %.0f%%", formatC(pm$p, format = "fg", digits = 2),
                                    formatC(pm$hksj_p, format = "fg", digits = 2), pm$I2))
fa <- rbind(f_coh, f_pool)
fa$row <- factor(rows[fa$key], levels = rev(rows))
fa$module <- factor(ifelse(fa$module == "up", "ALS-up module", "ALS-down module"), levels = c("ALS-up module", "ALS-down module"))
disc <- eff[eff$cohort == "NYGC_MTC_discovery", ]
disc_note <- setNames(sprintf("Discovery d = %.1f (not pooled)", disc$d), ifelse(disc$module == "up", "ALS-up module", "ALS-down module"))
write_source(fa[, c("row", "module", "kind", "d", "lo", "hi", "note")], out, "panel_a_forest")  # note carries z-test and HKSJ P
diamond <- function(df) {
  y <- as.numeric(df$row)
  do.call(rbind, lapply(seq_len(nrow(df)), function(i) data.frame(id = i, module = df$module[i], kind = df$kind[i],
    x = c(df$lo[i], df$d[i], df$hi[i], df$d[i]), y = y[i] + c(0, 0.28, 0, -0.28))))
}
pools <- fa[fa$kind %in% c("pool", "pool_sens"), ]
dia <- diamond(pools)
cohorts <- fa[!fa$kind %in% c("pool", "pool_sens"), ]
# Per-module x range (all CIs fit) and the x position of the right-hand annotations.
x_range <- data.frame(module = factor(rep(levels(fa$module), each = 2), levels = levels(fa$module)), x = c(-1.2, 4.6, -7, 1.5))
note_x <- c("ALS-up module" = 1.5, "ALS-down module" = -6.9)  # inside the empty side of each panel
mod_col <- c("ALS-up module" = PAL$up, "ALS-down module" = PAL$down)
pa <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = PAL$neutral) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = length(pool_lab) + 0.5, fill = "#F3F3F3") +
  geom_errorbar(data = cohorts, aes(y = row, xmin = lo, xmax = hi, colour = module), width = 0, orientation = "y", linewidth = 0.35) +
  geom_point(data = cohorts, aes(y = row, x = d, colour = module, shape = kind), size = 1.5, fill = "white", stroke = 0.45) +
  geom_polygon(data = dia, aes(x = x, y = y, group = id, fill = module, alpha = kind), colour = PAL$ink, linewidth = 0.25) +
  geom_blank(data = x_range, aes(x = x)) +
  geom_text(data = pools, aes(y = row, x = note_x[as.character(module)], label = note), hjust = 0, size = 1.75, family = "Arial") +
  geom_text(data = data.frame(module = factor(names(disc_note), levels = levels(fa$module)), label = disc_note),
            aes(x = note_x[as.character(module)], y = length(rows) + 0.45, label = label), hjust = 0, size = 1.75, family = "Arial", colour = PAL$neutral) +
  facet_wrap(~ module, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = mod_col, guide = "none") +
  scale_fill_manual(values = mod_col, guide = "none") +
  scale_alpha_manual(values = c(pool = 1, pool_sens = 0.25), guide = "none") +
  scale_shape_manual(values = c(cohort = 16, cohort_excluded = 21), guide = "none") +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(breaks = seq(-6, 4, 2), expand = expansion(mult = 0)) +
  coord_cartesian(ylim = c(0.5, length(rows) + 0.7), clip = "off") +
  labs(x = "Standardised ALS - control module effect (d, 95% CI)", y = NULL, title = "Module effects in independent snRNA-seq contrasts") +
  theme(panel.spacing = unit(8, "mm"))

# ---- b: genome-wide meta-analysis without discovery (run ii) ----
gw <- read.delim(gzfile(file.path(W58, "w58_genome_wide_meta.tsv.gz")), na.strings = "NA")
gw <- gw[gw$analysis == "ii_without_discovery" & !is.na(gw$meta_p), ]
gw$module[is.na(gw$module)] <- ""
gw$class <- ifelse(gw$consensus %in% c(TRUE, "TRUE"), ifelse(gw$meta_logFC > 0, "Consensus up", "Consensus down"), "Other")
gw$class <- factor(gw$class, levels = c("Other", "Consensus up", "Consensus down"))
gw <- gw[order(gw$class), ]
named <- c("PXDN", "TNS1", "SPP1", "MYOF", "MYLK", "CORO1C", "SNED1", "BCL6", "SLC1A4", "ID4", "SLITRK2", "FKBP8", "CEND1", "CDH2", "ZNF462", "C3")
lab <- gw[gw$symbol %in% named, ]
lab$face <- ifelse(lab$module != "", "bold", "plain")
write_source(gw[, c("gene", "symbol", "k", "meta_logFC", "meta_se", "meta_p", "meta_BH", "I2", "frac_same_sign", "loo_max_BH", "consensus", "module", "class")],
             out, "panel_b_genome_wide_meta_run_ii")
n_up <- sum(gw$class == "Consensus up"); n_down <- sum(gw$class == "Consensus down")
pb <- ggplot(gw, aes(x = meta_logFC, y = neglog10(meta_p))) +
  geom_point(aes(colour = class, size = class), alpha = 0.8, stroke = 0) +
  ggrepel::geom_text_repel(data = lab, aes(label = symbol, fontface = face), size = 1.85, family = "Arial", min.segment.length = 0,
                           segment.size = 0.2, box.padding = 0.25, max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = c(Other = PAL$light, `Consensus up` = PAL$up, `Consensus down` = PAL$down), name = NULL,
                      labels = c(Other = "Other genes", `Consensus up` = paste0("Consensus up (", n_up, ")"),
                                 `Consensus down` = paste0("Consensus down (", n_down, ")"))) +
  scale_size_manual(values = c(Other = 0.35, `Consensus up` = 0.8, `Consensus down` = 0.8), guide = "none") +
  guides(colour = guide_legend(override.aes = list(size = 1.6))) +
  scale_x_continuous(limits = c(-1.75, 1.75)) +
  labs(x = "Meta-analysis log2 fold change (five contrasts, discovery excluded)", y = expression(-log[10]~italic(P)[meta]),
       title = "Discovery-agnostic genome-wide meta-analysis") +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1), legend.key.size = unit(2.4, "mm"))

fig <- (pa / free(pb)) + plot_layout(heights = c(1.15, 1.1)) + plot_annotation(tag_levels = "a")
save_pub_r(fig, out, "Figure4_meta_analysis", width_mm = 183, height_mm = 160)
cat("Fig. 4 written to", out, "\n")
