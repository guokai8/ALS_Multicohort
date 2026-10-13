#!/usr/bin/env Rscript
# W68 inputs (Amendment 11): donor x territory pseudobulk of raw Spatial counts from the GSE288365
# Seurat objects, and a spot table (frozen ALS-up module score, spot depth, protein CLR values) for
# the protein co-localisation analysis. No disease comparison is made here.

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
src <- file.path(p01_source_data_root(), "GSE288365", "seurat")
out <- file.path(root, "results", "W68_GSE288365_replication", "inputs")

# Sample-to-donor map from Supplementary Table 1 of Zhang et al., Nat Neurosci 2026 (sample prefix -> donor, diagnosis, age, sex)
donor_map <- read.delim(text = "prefix\tdonor\tdiagnosis\tage\tsex
GWF.17.23_\tNEUXW599BNC\tControl\t74\tf
GBB20.16.A13\tNEUKJ976GGT\tControl\t71\tm
NEUDW867LF9\tNEUDW867LF9\tControl\t61\tm
NEUHD481VCL\tNEUHD481VCL\tControl\t54\tf
GBB21.03.A12\tNEUJV725KH1\tControl\t51\tm
GBB.18.01.A9\tNEUWD570BTK\tControl\t68\tm
GBB19.13.A13\tNEUVM032BF2\tControl\t75\tf
GBB.17.21.B10\tNEUEP010ZFE\tControl\t51\tm
NEUUJ189PPH\tNEUUJ189PPH\tControl\tNA\tNA
GWF.17.24.10B\tNEUBB078VEG\tALS_C9\t51\tf
MILLERMINEUWM\tNEUWM344ZLM\tALS_C9\t55\tf
MILLERMI.AU.87\tNEUBF291JXB\tALS_C9\tNA\tNA
MILLERMI.N14.280\tNEUEN501PPU\tALS_C9\tNA\tNA
GWF.20.50.10B\tNEUBD692FG6\tALS_sporadic\t67\tm
GWF.19.37.10B\tNEUPG876HJT\tALS_sporadic\t80\tm
GWF.18.33.10B\tNEUDK592BDZ\tALS_sporadic\t62\tf
GWF_19.35_10C\tNEUYT145HVC\tALS_sporadic\t40\tm
GWF_19.47_10B\tNEUCP856ZH7\tALS_sporadic\tNA\tm", stringsAsFactors = FALSE)

seed <- file.path(root, "results", "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
up_sym <- read.delim(file.path(seed, "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_up.tsv"))$common_name
territory_of <- c(Anterior_horns = "anterior_horn", Grey_matter = "other_grey", White_matter = "lateral_ventral_white",
                  Posterior_white_matter = "dorsal_columns")
proteins <- c("CD68.1", "CD163.1", "ITGAM.1", "HLA-DRA", "PTPRC.1", "PECAM1.1", "ACTA2.1", "VIM.1")
isotypes <- c("mouse-IgG2a", "mouse-IgG1k", "mouse-IgG2bk", "rat-IgG2a")

tmp <- tempfile(fileext = ".rds")
blocks <- list(); metas <- list(); spots <- list(); genes <- NULL
for (f in list.files(src, pattern = "SeuratObject\\.rds\\.gz$", full.names = TRUE)) {
  sample <- sub("^GSE288365_", "", sub("_SeuratObject\\.rds\\.gz$", "", basename(f)))
  hit <- donor_map[startsWith(paste0(sample, "_"), donor_map$prefix) | startsWith(sample, donor_map$prefix), ]
  stopifnot(nrow(hit) == 1)
  system2("gzip", c("-dc", shQuote(f)), stdout = tmp)
  o <- readRDS(tmp)
  for (a in c("Spatial", "Protein")) if (length(Layers(o[[a]], search = "counts")) > 1) o[[a]] <- JoinLayers(o[[a]])  # two capture areas
  counts <- LayerData(o, assay = "Spatial", layer = "counts")[, rownames(o@meta.data)]
  if (is.null(genes)) genes <- rownames(counts)
  stopifnot(identical(rownames(counts), genes))
  md <- o@meta.data
  terr <- unname(territory_of[as.character(md$manual)])
  adj <- ifelse(md$manual == "Anterior_horns" & md$MN_adjacent %in% "MN_adjacent", "anterior_horn_MN_adjacent",
                ifelse(md$manual == "Anterior_horns" & md$MN_adjacent %in% "not_MN_adjacent", "anterior_horn_MN_distal", NA))
  for (lab in list(list(type = "territory", v = terr), list(type = "mn_proximity", v = adj))) {
    keep <- !is.na(lab$v)
    f_ <- factor(lab$v[keep])
    ind <- sparseMatrix(i = seq_len(sum(keep)), j = as.integer(f_), x = 1, dims = c(sum(keep), nlevels(f_)))
    pb <- counts[, keep, drop = FALSE] %*% ind
    colnames(pb) <- paste(lab$type, hit$donor, levels(f_), sep = ":")
    blocks[[length(blocks) + 1]] <- pb
    metas[[length(metas) + 1]] <- data.frame(unit = colnames(pb), unit_type = lab$type, stratum = levels(f_), donor = hit$donor,
                                             diagnosis = hit$diagnosis, age = hit$age, sex = hit$sex, sample = sample,
                                             n_spots = as.vector(table(f_)))
  }
  # spot table for W68 C
  total <- colSums(counts)
  up_rows <- intersect(up_sym, rownames(counts))
  prot <- LayerData(o, assay = "Protein", layer = "counts")[, rownames(o@meta.data)]
  prot <- prot[setdiff(rownames(prot), isotypes), , drop = FALSE]
  clr <- apply(as.matrix(prot), 2, function(x) log1p(x) - mean(log1p(x)))
  spots[[length(spots) + 1]] <- data.frame(sample = sample, donor = hit$donor, diagnosis = hit$diagnosis, manual = md$manual,
                                           total_counts = total, module_up = log1p(1e4 * colSums(counts[up_rows, , drop = FALSE]) / pmax(total, 1)),
                                           t(clr[proteins, , drop = FALSE]), check.names = FALSE)
}
unlink(tmp)
pb <- do.call(cbind, blocks)
meta <- do.call(rbind, metas)
# donors with two Seurat objects are summed per unit (none expected; checked here)
stopifnot(!anyDuplicated(meta$unit))
meta$lib_size <- colSums(pb)
write.table(meta, file.path(out, "w68_units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeMM(pb, file.path(out, "w68_pseudobulk_counts.mtx")); system2("gzip", c("-f", file.path(out, "w68_pseudobulk_counts.mtx")))
writeLines(genes, gzfile(file.path(out, "w68_genes.txt.gz")))
writeLines(colnames(pb), file.path(out, "w68_unit_order.txt"))
write.table(do.call(rbind, spots), gzfile(file.path(out, "w68_spot_table.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
print(table(meta$unit_type, meta$stratum)); cat("genes:", length(genes), "\n")
