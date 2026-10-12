#!/usr/bin/env Rscript

# Figure source data for Supplementary Fig. S3: gene-level results of the W56 E3a
# design (W55 primary design + activated-myeloid and vascular scores), which W56
# saved only at gene-set level. Re-runs the E3a fits unchanged and stops unless the
# set-level results equal the committed w56_e3_nygc_composition_retest.tsv.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "run_w55_nygc_bulk_holdout.R"))

w56 <- file.path(root, "results", "W56_cell_source_check")
extra <- list(
  myeloid_act = symbols_to_ensg(c("CD68", "TREM2", "TYROBP", "SPP1", "CD163", "ITGAX")),
  vascular = symbols_to_ensg(c("CLDN5", "PECAM1", "VWF", "PDGFRB"))
)
fits <- lapply(c("motor_cortex", "spinal_cord"), fit_region, extra_markers = extra)
sets <- do.call(rbind, lapply(fits, `[[`, "sets"))
genes <- do.call(rbind, lapply(fits, `[[`, "genes"))

committed <- read.delim(file.path(w56, "w56_e3_nygc_composition_retest.tsv"))
committed <- committed[committed$analysis == "E3a_myeloid_vascular_adjusted", ]
chk <- merge(sets, committed, by = c("region", "module"), suffixes = c("", "_w56"))
if (nrow(chk) != 4) stop("E3a set rows not matched")
max_diff <- max(abs(chk$fry_p - chk$fry_p_w56), abs(chk$median_logFC - chk$median_logFC_w56))
if (max_diff > 1e-10) stop("E3a gene source data do not reproduce committed W56 E3a values: ", max_diff)

genes$design <- "E3a_myeloid_vascular_adjusted"
write.table(genes, file.path(w56, "w56_e3a_module_gene_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("E3a gene source data written; max diff vs committed set results =", max_diff, "\n")
