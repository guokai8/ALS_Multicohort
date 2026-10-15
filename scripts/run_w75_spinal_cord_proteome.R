#!/usr/bin/env Rscript

# W75 (Amendment 17): ALS-up module proteins in an independent human spinal-cord proteome
# (PRIDE PXD067060; Feneberg et al. 2025). A: ALS vs control, spinal cord, detergent-soluble (sarkosyl)
# fraction (primary). B: ALS vs PD and AD, PD and AD vs control (disease specificity). C: urea fraction,
# motor cortex and the ALS-down module (descriptive).

suppressPackageStartupMessages(library(limma))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
res <- file.path(root, "results")
out <- file.path(res, "W75_spinal_cord_proteome")
dir.create(out, showWarnings = FALSE)
SRC <- file.path(p01_source_data_root(), "Other_Datasets", "Feneberg2025_PXD067060")

seed <- file.path(res, "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
modules <- lapply(c(up = "up", down = "down"), function(d)
  read.delim(file.path(seed, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", d, ".tsv")))$common_name)
LEADING <- c("GPNMB", "C3", "LTF", "ANXA1", "EMP3")
# Experiment suffix -> diagnosis. Urea: 2/4/6/8 = CTL/PD/AD/ALS (sample sheet). Sarkosyl: 1/3/5/7, verified by
# one-to-one pairing of tissue x donor index with the urea fraction (metadata only).
DX <- list(Sark = c(`1` = "Control", `3` = "PD", `5` = "AD", `7` = "ALS"),
           Urea = c(`2` = "Control", `4` = "PD", `6` = "AD", `8` = "ALS"))

load_fraction <- function(fraction, tissue) {
  pg <- read.delim(file.path(SRC, fraction, "proteinGroups.txt"), check.names = FALSE, quote = "", stringsAsFactors = FALSE)
  flag <- function(col) if (col %in% names(pg)) pg[[col]] %in% "+" else rep(FALSE, nrow(pg))
  pg <- pg[!flag("Reverse") & !flag("Potential contaminant") & !flag("Only identified by site") & pg[["Unique peptides"]] >= 2, ]
  qcol <- grep("^LFQ intensity ", names(pg), value = TRUE)
  if (!length(qcol)) qcol <- grep("^Intensity ", names(pg), value = TRUE)
  ex <- sub("^(LFQ intensity|Intensity) ", "", qcol)
  m <- regmatches(ex, regexec("^(CTX|SPC)_P(\\d+)_(\\d+)$", ex))
  keep <- vapply(m, function(z) length(z) == 4 && z[2] == tissue && z[4] %in% names(DX[[fraction]]), logical(1))
  x <- as.matrix(pg[, qcol[keep]]); x[x <= 0] <- NA; x <- log2(x)
  colnames(x) <- ex[keep]
  samples <- data.frame(sample = ex[keep], diagnosis = unname(DX[[fraction]][vapply(m[keep], `[`, "", 4)]))
  # normalisation as in the source article: median of the 90% lowest-variance proteins
  complete <- rowSums(is.na(x)) == 0
  v <- apply(x[complete, , drop = FALSE], 1, var)
  ref <- names(v)[v <= quantile(v, 0.9)]
  x <- sweep(x, 2, apply(x[ref, , drop = FALSE], 2, median))
  present <- rowMeans(!is.na(x)) >= 0.5
  rownames(x) <- paste0("pg", seq_len(nrow(x)))
  x <- x[present, , drop = FALSE]
  genes <- pg[["Gene names"]][present]; upep <- pg[["Unique peptides"]][present]
  list(x = x, samples = samples, genes = genes, upep = upep, n_proteins = nrow(x))
}

gene_rows <- function(d, symbols) {
  idx <- vapply(symbols, function(g) {
    hit <- which(vapply(strsplit(d$genes, ";"), function(s) g %in% s, logical(1)))
    if (!length(hit)) NA_integer_ else hit[which.max(d$upep[hit])]
  }, integer(1))
  idx <- idx[!is.na(idx)]
  setNames(rownames(d$x)[idx], names(idx))
}

run_contrasts <- function(d, label, groups) {
  # A: two-group model (ALS, Control); B: four-group model (Control, PD, AD, ALS)
  keep <- d$samples$diagnosis %in% groups
  d$x <- d$x[, keep, drop = FALSE]; s <- d$samples[keep, ]
  s$diagnosis <- factor(s$diagnosis, levels = intersect(c("Control", "PD", "AD", "ALS"), groups))
  design <- model.matrix(~ 0 + diagnosis, data = s); colnames(design) <- levels(s$diagnosis)
  cns <- if (length(groups) == 2) c(ALS_vs_Control = "ALS - Control") else
    c(ALS_vs_PD = "ALS - PD", ALS_vs_AD = "ALS - AD", PD_vs_Control = "PD - Control", AD_vs_Control = "AD - Control")
  contr <- makeContrasts(contrasts = cns, levels = design); colnames(contr) <- names(cns)
  fit <- eBayes(contrasts.fit(lmFit(d$x, design), contr), trend = TRUE)
  idx <- lapply(modules, function(g) gene_rows(d, g))
  sets <- do.call(rbind, lapply(colnames(contr), function(cn) do.call(rbind, lapply(names(idx), function(mod) {
    rows <- idx[[mod]]
    # fry and camera need complete data: the set test uses module proteins quantified in every analysed sample
    complete <- rowSums(is.na(d$x)) == 0
    rows_cc <- rows[complete[unname(rows)]]
    lfc <- fit$coefficients[unname(rows), cn]
    if (length(rows_cc) < 3) return(data.frame(analysis = label, contrast = cn, module = mod, n_proteins = length(rows),
      n_proteins_set_test = length(rows_cc), fry_direction = NA, fry_p = NA, camera_p = NA,
      median_logFC = if (length(rows)) median(lfc) else NA, n_up = sum(lfc > 0)))
    xc <- d$x[complete, , drop = FALSE]
    set_idx <- match(unname(rows_cc), rownames(xc))
    fr <- fry(xc, index = list(set = set_idx), design = design, contrast = contr[, cn])
    ca <- camera(xc, index = list(set = set_idx), design = design, contrast = contr[, cn])
    data.frame(analysis = label, contrast = cn, module = mod, n_proteins = length(rows), n_proteins_set_test = length(rows_cc),
               fry_direction = fr$Direction, fry_p = fr$PValue, camera_p = ca$PValue, median_logFC = median(lfc), n_up = sum(lfc > 0))
  }))))
  genes <- do.call(rbind, lapply(names(idx), function(mod) {
    rows <- idx[[mod]]; if (!length(rows)) return(NULL)
    do.call(rbind, lapply(colnames(contr), function(cn) {
      tt <- topTable(fit, coef = cn, number = Inf, sort.by = "none")[unname(rows), ]
      data.frame(analysis = label, contrast = cn, module = mod, gene = names(rows), logFC = tt$logFC, P.Value = tt$P.Value,
                 n_quantified = rowSums(!is.na(d$x[unname(rows), , drop = FALSE])), n_samples = ncol(d$x))
    }))
  }))
  genes$leading <- genes$gene %in% LEADING
  genes$BH_leading <- NA; genes$BH_module <- NA
  for (k in unique(paste(genes$contrast, genes$module))) {
    w <- paste(genes$contrast, genes$module) == k
    genes$BH_module[w] <- p.adjust(genes$P.Value[w], "BH")
    wl <- w & genes$leading
    genes$BH_leading[wl] <- p.adjust(genes$P.Value[wl], "BH")
  }
  counts <- as.data.frame(table(s$diagnosis)); names(counts) <- c("diagnosis", "n"); counts$analysis <- label
  counts$n_proteins_after_filter <- d$n_proteins
  list(sets = sets, genes = genes, counts = counts)
}

runs <- list(c("Sark", "SPC", "spinal_cord_soluble"), c("Urea", "SPC", "spinal_cord_insoluble"),
             c("Sark", "CTX", "cortex_soluble"), c("Urea", "CTX", "cortex_insoluble"))
allr <- unlist(lapply(runs, function(r) {
  d <- load_fraction(r[1], r[2])
  list(run_contrasts(d, r[3], c("Control", "ALS")), run_contrasts(d, r[3], c("Control", "PD", "AD", "ALS")))
}), recursive = FALSE)
sets <- do.call(rbind, lapply(allr, `[[`, "sets")); genes <- do.call(rbind, lapply(allr, `[[`, "genes"))
counts <- unique(do.call(rbind, lapply(allr[seq(2, length(allr), 2)], `[[`, "counts")))

# ---- prespecified readings ----
pa <- sets[sets$analysis == "spinal_cord_soluble" & sets$contrast == "ALS_vs_Control" & sets$module == "up", ]
A_reading <- if (is.na(pa$fry_p)) "fewer than 3 module proteins measured; gene level only" else
  if (pa$fry_p < 0.05 && pa$fry_direction == "Up" && pa$median_logFC > 0) "protein-level support" else
  if (pa$median_logFC > 0) "directionally consistent" else "not supported"
b <- sets[sets$analysis == "spinal_cord_soluble" & sets$module == "up", ]; g <- function(cn) b[b$contrast == cn, ]
up_sig <- function(r) !is.na(r$fry_p) && r$fry_p < 0.05 && r$fry_direction == "Up"
B_reading <- if (up_sig(g("ALS_vs_PD")) && up_sig(g("ALS_vs_AD"))) "ALS-selective among these diseases" else
  if ((up_sig(g("PD_vs_Control")) && !up_sig(g("ALS_vs_PD"))) || (up_sig(g("AD_vs_Control")) && !up_sig(g("ALS_vs_AD"))))
    "shared with other neurodegeneration" else "undetermined"
reading <- data.frame(part = c("A", "B"), reading = c(A_reading, B_reading))

write.table(counts, file.path(out, "w75_sample_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sets, file.path(out, "w75_module_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes, file.path(out, "w75_module_protein_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(reading, file.path(out, "w75_reading.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(counts, row.names = FALSE)
print(sets[sets$module == "up", ], digits = 3, row.names = FALSE)
print(genes[genes$analysis == "spinal_cord_soluble" & genes$module == "up", ], digits = 3, row.names = FALSE)
print(reading, row.names = FALSE)
