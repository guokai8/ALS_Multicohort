#!/usr/bin/env Rscript

# W56 E3 (Amendment 2): NYGC bulk re-test of the W55 ALS-up signal.
# E3a: W55 primary design + activated-myeloid and vascular scores.
# E3b: astrocyte-enriched up subset (from E1) with the W55 primary design and
#      with the E3a design.
# Blocked fry is primary; exploratory, so no new multiplicity family is formed.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "run_w55_nygc_bulk_holdout.R"))

w56 <- file.path(root, "results", "W56_cell_source_check")
extra <- list(
  myeloid_act = symbols_to_ensg(c("CD68", "TREM2", "TYROBP", "SPP1", "CD163", "ITGAX")),
  vascular = symbols_to_ensg(c("CLDN5", "PECAM1", "VWF", "PDGFRB"))
)
stopifnot(!any(unlist(extra) %in% unlist(modules)))
subset_ids <- read.delim(file.path(w56, "w56_e1_up_astrocyte_enriched_subset.tsv"))$gene_id
subset_sets <- list(up_astro_enriched = subset_ids)

runs <- list(
  list(label = "E3a_myeloid_vascular_adjusted", extra = extra, sets = modules),
  list(label = "E3b_astro_subset_primary_design", extra = list(), sets = subset_sets),
  list(label = "E3b_astro_subset_E3a_design", extra = extra, sets = subset_sets)
)
res <- do.call(rbind, lapply(runs, function(r) {
  do.call(rbind, lapply(c("motor_cortex", "spinal_cord"), function(region) {
    x <- fit_region(region, extra_markers = r$extra, gene_sets = r$sets)$sets
    cbind(analysis = r$label, x)
  }))
}))
write.table(res, file.path(w56, "w56_e3_nygc_composition_retest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res[, c("analysis", "region", "module", "n_samples", "n_module_genes_tested", "fry_direction", "fry_p", "camera_p", "n_genes_concordant", "median_logFC")], digits = 3, row.names = FALSE)
