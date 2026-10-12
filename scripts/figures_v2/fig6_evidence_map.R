#!/usr/bin/env Rscript
# Fig. 6 | Evidence map of the reproducible ALS programme. Each box is labelled with the tier of
# evidence behind it; arrows show the order of the analyses, not signalling. All numbers are read
# from committed result tables (W55, W56, W58, W59, W61, W63, W64, W65).

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "fig_common.R"))
out <- new_figure_dir("Fig6_evidence_map")
rd <- function(...) read.delim(file.path(RES, ...), check.names = FALSE)
fp <- function(p) if (p < 1e-3) formatC(p, format = "e", digits = 1) else formatC(signif(p, 2), format = "fg")

# ---- numbers from committed tables ----
w55 <- rd("W55_NYGC_bulk_holdout", "w55_primary_module_tests.tsv")
w55_up <- setNames(w55$fry_p[w55$module == "up"], w55$region[w55$module == "up"])
w59 <- rd("W59_NYGC_consensus_signature_holdout", "w59_signature_tests.tsv")
w59 <- w59[w59$analysis == "1_primary_W55_design", ]
w59_p <- setNames(w59$fry_p, paste(w59$region, w59$module))
w61 <- rd("W61_small_k_meta_inference", "w61_small_k_meta_inference.tsv")
w61 <- w61[w61$pool == "primary_motor_independent" & w61$module == "up", ]
w64 <- rd("W64_bulk_published_composition", "w64_set_tests.tsv")
w64 <- w64[w64$design == "a_st7_replaces_markers" & w64$gene_set == "up", ]
w64_p <- setNames(w64$fry_p, w64$region)
w63 <- rbind(rd("W63_pathway_context", "w63_themed_gsea.tsv"), rd("W63_pathway_context", "w63_hallmark_gsea.tsv"))
g <- function(set) w63[w63$pathway == set, ][1, ]
spec <- rd("W56_cell_source_check", "w56_e1_module_gene_specificity.tsv")

# ---- column 1: lineage with the highest baseline share of each ALS-up gene (W56 E1, Pineda) ----
share_cols <- grep("^pineda_share_", names(spec), value = TRUE)
up <- spec[spec$module == "up", ]
up$top <- sub("^pineda_share_", "", share_cols[max.col(up[, share_cols], ties.method = "first")])
up$top_share <- apply(up[, share_cols], 1, max)
up$group <- ifelse(up$top %in% c("Excitatory neuron", "Inhibitory neuron", "Oligodendrocyte"), "Neuron / oligodendrocyte", up$top)
lineage_order <- c("Endothelial", "Vascular/Meningeal", "Microglia", "Astrocyte", "OPC", "Neuron / oligodendrocyte")
up <- up[order(match(up$group, lineage_order), -up$top_share), ]
col1 <- do.call(rbind, lapply(lineage_order, function(l) {
  d <- up[up$group == l, ]
  genes <- paste0(d$common_name, " (", sprintf("%.2f", d$top_share), ")")
  body <- paste(vapply(split(genes, ceiling(seq_along(genes) / 2)), paste, "", collapse = ", "), collapse = "\n")
  data.frame(title = paste0(l, ": ", nrow(d), if (nrow(d) == 1) " gene" else " genes"), body = body, tier = "exploratory")
}))

# ---- column 2: what is reproducible ----
col2 <- data.frame(
  title = c("Frozen ALS-up module (15 genes)",
            "Bulk hold-out, frozen module (W55)",
            "Bulk hold-out, 158-gene consensus (W59)",
            "snRNA-seq meta-analysis (W58)",
            "Published deconvolution (W64)"),
  body = c("Discovery cohort defines candidates only",
           sprintf("Motor cortex P = %s; spinal cord P = %s", fp(w55_up[["motor_cortex"]]), fp(w55_up[["spinal_cord"]])),
           sprintf("Up arm: SC P = %s; MC P = %s\nDown arm: MC P = %s; SC P = %s",
                   fp(w59_p[["spinal_cord up_consensus"]]), fp(w59_p[["motor_cortex up_consensus"]]),
                   fp(w59_p[["motor_cortex down_consensus"]]), fp(w59_p[["spinal_cord down_consensus"]])),
           sprintf("Independent cohorts, primary pool\nd = %.2f, z-test P = %s; HKSJ P = %s",
                   w61$pooled_d, fp(w61$z_p), fp(w61$hksj_p)),
           sprintf("Motor cortex P = %s (median log2FC %.2f)\nSpinal cord P = %s (median log2FC %.2f)",
                   fp(w64_p[["motor_cortex"]]), w64$median_logFC[w64$region == "motor_cortex"],
                   fp(w64_p[["spinal_cord"]]), w64$median_logFC[w64$region == "spinal_cord"])),
  tier = c("candidate", "confirmed", "confirmed", "model", "model")
)

# ---- column 3: pathway context (W63, preranked GSEA on the discovery-agnostic meta statistic) ----
emt <- g("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"); ecm <- g("GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION")
ifn <- g("HALLMARK_INTERFERON_GAMMA_RESPONSE"); tnf <- g("HALLMARK_TNFA_SIGNALING_VIA_NFKB")
cmp <- g("GOBP_COMPLEMENT_ACTIVATION"); ang <- g("HALLMARK_ANGIOGENESIS")
oxp <- g("HALLMARK_OXIDATIVE_PHOSPHORYLATION"); dna <- g("HALLMARK_DNA_REPAIR")
not_sets <- c("GOBP_ASTROCYTE_ACTIVATION", "GOBP_GLIAL_CELL_ACTIVATION", "GOBP_MAINTENANCE_OF_BLOOD_BRAIN_BARRIER")
col3 <- data.frame(
  title = c("ECM remodelling (up)", "Interferon and NF-κB (up)", "Complement activation (up)", "Angiogenesis (up)",
            "OXPHOS, DNA repair (down)", "Astrocyte/glial activation, BBB upkeep"),
  body = c(sprintf("EMT NES %.2f, FDR %s\nECM organisation FDR %s", emt$NES, fp(emt$padj), fp(ecm$padj)),
           sprintf("IFN-γ FDR %s; TNFα–NF-κB FDR %s", fp(ifn$padj), fp(tnf$padj)),
           sprintf("NES %.2f, FDR %s (C1R, C1S, C3, C1QB, C2)", cmp$NES, fp(cmp$padj)),
           sprintf("NES %.2f, FDR %s", ang$NES, fp(ang$padj)),
           sprintf("NES %.2f and %.2f, FDR %s and %s", oxp$NES, dna$NES, fp(oxp$padj), fp(dna$padj)),
           sprintf("Not enriched (FDR %s–%s)", fp(min(w63$padj[w63$pathway %in% not_sets])),
                   fp(max(w63$padj[w63$pathway %in% not_sets])))),
  tier = c(rep("exploratory", 5), "not_supported")
)

# ---- column 4: candidate communication axes (W65; liana Consensus; no tests) ----
col4 <- data.frame(
  title = c("Complement", "Angiopoietin", "Stem-cell factor", "Ephrin and EGF", "Annexin A1"),
  body = c("C3→CR1 higher; C3→C3AR1, LRP1 lower", "ANGPT1→TEK, ITGA5/ITGB1 lower",
           "KITLG→KIT lower; KITLG→EPOR higher", "EFNA5→EPH, EGF→ERBB3/4: mixed",
           "ANXA1→EGFR, DYSF higher"),
  tier = "hypothesis"
)

# ---- layout (1 unit = 1 mm) ----
W <- 183; COL_W <- 41.5; GAP <- (W - 4 * COL_W - 4) / 3
col_x0 <- 2 + (0:3) * (COL_W + GAP)
LINE <- 2.35; TOP <- 78
stack <- function(d, col) {
  n_lines <- lengths(strsplit(d$body, "\n", fixed = TRUE))
  h <- 3.6 + n_lines * LINE + 1.2
  gap <- 1.8
  y1 <- TOP - c(0, cumsum(h + gap))[seq_len(nrow(d))]
  cbind(d, col = col, x0 = col_x0[col], x1 = col_x0[col] + COL_W, y1 = y1, y0 = y1 - h)
}
boxes <- rbind(stack(col1, 1), stack(col2, 2), stack(col3, 3), stack(col4, 4))
boxes$xm <- (boxes$x0 + boxes$x1) / 2
BOTTOM <- min(boxes$y0) - 6
H <- TOP + 10
write_source(boxes, out, "evidence_map_boxes")

headers <- data.frame(
  x = col_x0 + COL_W / 2, y = TOP + 6.5,
  title = c("Cell sources", "Reproducible ALS programme", "Pathway context", "Candidate communication axes"),
  sub = c("Lineage with highest baseline share (W56)", "Prespecified tests", "Astrocyte-compartment meta-analysis (W63)", "Ligand–receptor pairs (W65)")
)
arrows <- data.frame(x = col_x0[1:3] + COL_W + 0.6, xend = col_x0[2:4] - 0.6, y = (TOP + BOTTOM) / 2 + 3, yend = (TOP + BOTTOM) / 2 + 3)

tier_fill <- c(candidate = "#EFE3DF", confirmed = "#EAD7A8", model = "#C9D6E3", exploratory = "white",
               hypothesis = "white", not_supported = "#EDEDED")
tier_lab <- c(candidate = "Candidate definition (discovery)", confirmed = "Confirmed (prespecified hold-out)",
              model = "Supportive, model-dependent", exploratory = "Exploratory association",
              hypothesis = "Hypothesis only (no test)", not_supported = "Not supported")
boxes$tier <- factor(boxes$tier, levels = names(tier_fill))
boxes$lty <- ifelse(boxes$tier == "hypothesis", "dashed", "solid")
boxes$text_col <- ifelse(boxes$tier == "not_supported", "#6B6B6B", PAL$ink)

fig <- ggplot() +
  geom_segment(data = arrows, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4, colour = PAL$neutral,
               arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  geom_rect(data = boxes, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = tier, linetype = lty),
            colour = PAL$ink, linewidth = 0.25) +
  geom_text(data = boxes, aes(x = xm, y = y1 - 2.1, label = title, colour = text_col), fontface = "bold", size = 1.95,
            family = "Arial", vjust = 0.5) +
  geom_text(data = boxes, aes(x = xm, y = y1 - 3.9, label = body, colour = text_col), size = 1.8, family = "Arial",
            vjust = 1, lineheight = 0.92) +
  geom_text(data = headers, aes(x = x, y = y, label = title), fontface = "bold", size = 2.3, family = "Arial") +
  geom_text(data = headers, aes(x = x, y = y - 3, label = sub), size = 1.85, family = "Arial", colour = "#555555") +
  annotate("text", x = W / 2, y = BOTTOM + 2.5, size = 1.8, family = "Arial", colour = "#555555",
           label = "Evidence map, not a mechanism: arrows give the order of analyses, not signalling. Numbers in brackets are baseline lineage shares. MC, motor cortex; SC, spinal cord.") +
  scale_fill_manual(values = tier_fill, labels = tier_lab, name = NULL, drop = FALSE) +
  scale_linetype_identity() +
  scale_colour_identity() +
  guides(fill = guide_legend(nrow = 1, override.aes = list(colour = PAL$ink, linewidth = 0.25,
                                                           linetype = c(rep("solid", 4), "dashed", "solid")))) +
  coord_fixed(xlim = c(0, W), ylim = c(BOTTOM, H), expand = FALSE, clip = "off") +
  theme_void(base_family = "Arial") +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.6), legend.key.size = unit(2.6, "mm"),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(1, 1, 1, 1))

save_pub_r(fig, out, "Figure6_evidence_map", width_mm = 183, height_mm = H - BOTTOM + 8)
cat("Fig. 6 written to", out, "\n")
