#!/usr/bin/env Rscript
# W68 preparation: spot metadata (annotations, condition, motor-neuron distance) from the GSE288365
# Seurat objects. Metadata only; no expression is read here.
suppressPackageStartupMessages(library(SeuratObject))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
src <- file.path(p01_source_data_root(), "GSE288365", "seurat")
out <- file.path(p01_project_root(), "results", "W68_GSE288365_replication", "inputs")
tmp <- tempfile(fileext = ".rds")
cols <- c("sample", "barcode", "condition", "condition_general", "manual", "MN_dist_group", "MN_adjacent", "MN_dist_map",
          "IBA1_positive", "nCount_Spatial", "nFeature_Spatial")
md <- do.call(rbind, lapply(list.files(src, pattern = "SeuratObject\\.rds\\.gz$", full.names = TRUE), function(f) {
  system2("gzip", c("-dc", shQuote(f)), stdout = tmp)
  o <- readRDS(tmp)
  m <- o@meta.data
  m <- m[, intersect(cols, names(m))]
  m$file <- basename(f); m$spot <- rownames(o@meta.data)
  m
}))
unlink(tmp)
write.table(md, gzfile(file.path(out, "w68_spot_metadata.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
print(table(md$file, md$condition))
print(table(md$file, md$manual))
print(table(md$MN_dist_group, useNA = "ifany")); print(table(md$MN_adjacent, useNA = "ifany"))
