#!/usr/bin/env Rscript
# Fig. 1 | Study design: workflow schematic, plan/amendment commit order, cohort table.

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("Fig1_study_design")

# ---- a: workflow schematic ----
boxes <- data.frame(
  id = c("disc", "mod", "plan", "sn", "meta", "bulk1", "bulk2", "expl"),
  x = c(1, 1, 1, 0.35, 0.35, 1.65, 1.0, 1.65),
  y = c(6, 5, 4, 3, 2, 3, 1, 2),
  label = c("Discovery: GSE330130 motor cortex\nastrocyte-compartment pseudobulk (11 ALS / 7 ctrl)",
            "Frozen modules: 15 ALS-up / 14 ALS-down genes",
            "Frozen analysis plan v1.0 + Amendments 1-9\n(committed before each governed analysis)",
            "Independent snRNA-seq:\n5 contrasts from 4 studies",
            "Meta-analysis +\n159-gene consensus signature",
            "NYGC bulk hold-out:\nfrozen modules (test 1)",
            "NYGC bulk hold-out:\nconsensus signature (test 2)",
            "Exploratory: cell source, ambient RNA, composition,\npathways, ligand-receptor, Visium spatial context"),
  role = c("discovery", "discovery", "plan", "consolidation", "consolidation", "confirmatory", "confirmatory", "exploratory")
)
edges <- data.frame(from = c("disc", "mod", "plan", "plan", "sn", "meta", "bulk1"),
                    to = c("mod", "plan", "sn", "bulk1", "meta", "bulk2", "expl"))
edges <- merge(merge(edges, boxes[, c("id", "x", "y")], by.x = "from", by.y = "id"),
               boxes[, c("id", "x", "y")], by.x = "to", by.y = "id", suffixes = c("", "_end"))
role_fill <- c(discovery = "#EFE3DF", plan = "#C8C8C8", consolidation = "#C9D6E3", confirmatory = "#EAD7A8", exploratory = "white")
write_source(boxes, out, "panel_a_workflow_boxes")
pa <- ggplot() +
  geom_segment(data = edges, aes(x = x, y = y - 0.3, xend = x_end, yend = y_end + 0.3),
               linewidth = 0.35, colour = PAL$neutral, arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  geom_label(data = boxes, aes(x = x, y = y, label = label, fill = role), size = 1.95, family = "Arial", key_glyph = "rect",
             linewidth = 0.25, label.r = unit(0.6, "mm"), lineheight = 0.95, label.padding = unit(1.1, "mm")) +
  scale_fill_manual(values = role_fill, name = NULL,
                    labels = c(confirmatory = "Confirmatory", consolidation = "Consolidation", discovery = "Discovery",
                               exploratory = "Exploratory", plan = "Prespecification")) +
  guides(fill = guide_legend(override.aes = list(colour = PAL$ink, linewidth = 0.25))) +
  coord_cartesian(xlim = c(-0.25, 2.25), ylim = c(0.5, 6.5), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.6), legend.key.size = unit(2.5, "mm"),
        plot.tag = element_text(size = 8, face = "bold"))

# ---- b: commit order of plan, amendments and governed analyses ----
log_lines <- system2("git", c("-C", shQuote(FIG_ROOT), "log", "--reverse", shQuote("--format=%h|%ad|%s"), shQuote("--date=format:%Y-%m-%d %H:%M"), "841cd1a..HEAD"), stdout = TRUE)
log <- do.call(rbind, lapply(strsplit(log_lines, "\\|"), function(x) data.frame(hash = x[1], time = x[2], subject = x[3])))
events <- data.frame(
  hash = c("267949d", "b7266f0", "02f1136", "2861aef", "dfa6954", "9b1fd37", "bbd5e98", "8ce592d", "4624cc1", "957a06a", "5605120", "fedddd5", "896ce41", "7f16fd3", "284d22c", "9342a49", "b6b8c94", "73864d2", "0576362", "89c7c6a"),
  label = c("W50-W54 external checks (observed before plan)", "Plan v1.0 frozen", "Amendment 1 (hold-out implementation)",
            "W55 bulk hold-out: frozen modules", "Amendment 2 (cell source)", "W56 cell-source checks",
            "Amendment 3 (ambient RNA)", "W57 ambient-RNA test", "Amendment 4 (meta-analysis)", "W58 meta-analysis",
            "Amendment 5 (frozen consensus signature)", "W59 bulk hold-out: consensus signature",
            "Amendment 6 (age >= 70 sensitivity)", "W60 age >= 70 sensitivity of W58",
            "Amendment 7 (review-prompted checks)", "W61 HKSJ inference; W62 discovery audit",
            "Amendment 8 (exploratory biological context)", "W63-W65 pathways, composition, ligand-receptor",
            "Amendment 9 (spatial context)", "W66 Visium spatial context"),
  kind = c("prior", "plan", "plan", "confirmatory", "plan", "exploratory", "plan", "exploratory", "plan", "consolidation", "plan", "confirmatory", "plan", "consolidation", "plan", "consolidation", "plan", "exploratory", "plan", "exploratory")
)
log$commit_index <- seq_len(nrow(log))
events <- merge(events, log, by = "hash")
events <- events[order(events$commit_index), ]  # git commit order; several commits share a minute
events$order <- rev(seq_len(nrow(events)))
write_source(events, out, "panel_b_commit_order")
kind_col <- c(prior = "#8C8C8C", plan = role_fill[["plan"]], confirmatory = role_fill[["confirmatory"]],
              exploratory = role_fill[["exploratory"]], consolidation = role_fill[["consolidation"]])
pb <- ggplot(events, aes(x = 0, y = order)) +
  annotate("segment", x = 0, xend = 0, y = min(events$order), yend = max(events$order), colour = PAL$light, linewidth = 0.5) +
  geom_point(aes(fill = kind), shape = 21, size = 2.1, stroke = 0.35, colour = PAL$ink) +
  geom_text(aes(x = 0.08, label = paste0(label, "  [", hash, "]")), hjust = 0, size = 1.95, family = "Arial") +
  scale_fill_manual(values = kind_col, guide = "none") +
  coord_cartesian(xlim = c(-0.05, 1.6), clip = "off") +
  labs(title = "Version-control commit order") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 6.4, face = "bold"), plot.tag = element_text(size = 8, face = "bold"))

# ---- c: cohort table (counts computed from the analysis sample tables) ----
count_groups <- function(df, donor, group, als_value) {
  d <- df[!duplicated(df[[donor]]), ]
  c(als = sum(d[[group]] %in% als_value), ctrl = sum(!d[[group]] %in% als_value))
}
w33 <- read.delim(file.path(RES, "W33_synapse_resolved_diagnosis_sensitivity", "MTC_Astrocyte_RNA.library_label_audit.tsv"))
w33 <- w33[w33$used_in_w33 %in% c(TRUE, "TRUE"), ]
zel <- read.delim(file.path(RES, "W58_multicohort_meta", "w58_zelic_astrocyte_pseudobulk_samples.tsv"))
li <- read.delim(file.path(RES, "W54_Li2023_motor_cortex_validation", "w54_li2023_astrocyte_pseudobulk_samples.tsv"))
li <- li[li$disease %in% c("C9-ALS", "Control"), ]
ulm <- read.delim(file.path(RES, "W52_Ulm_motor_cortex_validation", "w52_ulm_astrocyte_pseudobulk_samples.tsv"))
ulm <- ulm[ulm$Case_Type %in% c("ALS", "HC"), ]
pin <- read.delim(file.path(RES, "W51_Pineda_BA4_discovery_matched_limma", "w51_pineda_astrocyte_pseudobulk_samples.tsv"))
pin$prefix <- substr(pin$sample_folder, 1, 6)
ba4 <- pin[pin$Region == "BA4" & pin$prefix %in% c("191112", "191114") & pin$Condition %in% c("ALS", "PN") & pin$lib_size > 0 & !is.na(pin$age), ]
ba9 <- pin[pin$Region == "BA9" & pin$prefix == "210526" & pin$Condition %in% c("ALS", "PN") & pin$lib_size > 0 & !is.na(pin$age), ]
nygc <- read.delim(file.path(RES, "W55_NYGC_bulk_holdout", "w55_primary_module_tests.tsv"))
nygc <- nygc[nygc$module == "up", ]
visium <- lapply(c(LSC = "LSC", MTC = "MTC"), function(t) {
  u <- read.delim(file.path(RES, "W66_spatial_context", "inputs", paste0(t, "_units.tsv")))
  u <- u[u$unit_type == "whole" & u$group %in% c("ALS", "Control") & !(u$gse330130_donor %in% c(TRUE, "True")), ]
  c(als = sum(u$group == "ALS"), ctrl = sum(u$group == "Control"))
})
counts <- rbind(
  count_groups(w33, "donor", "disease_binary_w33", "ALS"),
  count_groups(zel, "donor_id", "disease", "ALS"),
  count_groups(li, "donor", "disease", "C9-ALS"),
  count_groups(ulm, "ID", "Case_Type", "ALS"),
  count_groups(ba4, "Donor", "Condition", "ALS"),
  count_groups(ba9, "Donor", "Condition", "ALS"),
  c(als = 5, ctrl = 5),
  c(als = nygc$n_als_donors[nygc$region == "motor_cortex"], ctrl = nygc$n_control_donors[nygc$region == "motor_cortex"]),
  c(als = nygc$n_als_donors[nygc$region == "spinal_cord"], ctrl = nygc$n_control_donors[nygc$region == "spinal_cord"]),
  visium$LSC, visium$MTC
)
cohorts <- data.frame(
  cohort = c("NYGC GSE330130 (discovery)", "Zelic 2025", "Li 2023", "Ulm/DZNE (Ruf 2026)", "Pineda 2024 BA4 (2019 campaign)",
             "Pineda 2024 BA9 (batch 210526)", "Masrori 2025", "NYGC bulk hold-out", "NYGC bulk hold-out", "GSE330224 Visium", "GSE330224 Visium"),
  region = c("Motor cortex", "Cervical spinal cord", "Motor cortex (C9)", "Motor cortex", "Motor cortex", "Prefrontal cortex",
             "Motor cortex + spinal cord", "Motor cortex", "Spinal cord", "Lumbar spinal cord", "Motor cortex"),
  astrocytes = c("Author labels", "Marker rule", "Marker rule", "Author WNN labels", "Author labels", "Author labels",
                 "Author labels (summary DE)", "Bulk tissue", "Bulk tissue", "Spatial spots", "Spatial spots"),
  role = c("Discovery", "Independent snRNA", "Independent snRNA", "Independent snRNA", "Independent snRNA",
           "Independent snRNA (shares donors with BA4)", "Direction only", "Confirmatory hold-out", "Confirmatory hold-out",
           "Exploratory spatial", "Exploratory spatial"),
  als = counts[, "als"], ctrl = counts[, "ctrl"]
)
cohorts$donors <- paste0(cohorts$als, " / ", cohorts$ctrl)
write_source(cohorts, out, "panel_c_cohorts")
long <- rbind(
  data.frame(row = seq_len(nrow(cohorts)), col = 1, text = cohorts$cohort),
  data.frame(row = seq_len(nrow(cohorts)), col = 2, text = cohorts$region),
  data.frame(row = seq_len(nrow(cohorts)), col = 3, text = cohorts$donors),
  data.frame(row = seq_len(nrow(cohorts)), col = 4, text = cohorts$astrocytes),
  data.frame(row = seq_len(nrow(cohorts)), col = 5, text = cohorts$role)
)
col_x <- c(0, 2.35, 3.95, 4.75, 6.25)
long$x <- col_x[long$col]
long$y <- -long$row
header <- data.frame(x = col_x, y = 0, text = c("Cohort", "Region", "ALS / ctrl", "Astrocyte definition", "Role"))
role_rows <- data.frame(y = -seq_len(nrow(cohorts)), role = cohorts$role)
role_rows$band <- ifelse(grepl("Confirmatory", role_rows$role), "confirmatory", ifelse(role_rows$role == "Discovery", "discovery", "other"))
pc <- ggplot() +
  geom_tile(data = role_rows, aes(x = 4.1, y = y, fill = band), width = 8.6, height = 0.9, alpha = 0.6) +
  geom_text(data = header, aes(x = x, y = y, label = text), hjust = 0, fontface = "bold", size = 2.0, family = "Arial") +
  geom_text(data = long, aes(x = x, y = y, label = text), hjust = 0, size = 1.9, family = "Arial") +
  scale_fill_manual(values = c(confirmatory = "#F3E9D2", discovery = "#EFE3DF", other = "white"), guide = "none") +
  coord_cartesian(xlim = c(-0.1, 8.4), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.tag = element_text(size = 8, face = "bold"))

fig <- (pa | pb) / pc + plot_layout(heights = c(1.15, 1), widths = c(1.15, 1)) + plot_annotation(tag_levels = "a")
save_pub_r(fig, out, "Figure1_study_design", width_mm = 183, height_mm = 150)
cat("Fig. 1 written to", out, "\n")
