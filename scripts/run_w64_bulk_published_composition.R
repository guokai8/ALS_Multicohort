#!/usr/bin/env Rscript

# W64 (Amendment 8; exploratory): NYGC bulk hold-out re-tested against the consortium's
# published cell-type deconvolution (Flagship ST7). Design (a): W55 primary design with
# the marker-based composition scores replaced by ST7 proportions; design (b): W55
# primary design plus ST7 proportions. Gene sets: frozen modules and W59 consensus arms.
# Also: overlap of the frozen module genes with the consortium's own DE (Flagship ST4).

suppressPackageStartupMessages(library(readxl))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "run_w55_nygc_bulk_holdout.R"))

w64 <- file.path(root, "results", "W64_bulk_published_composition")
dir.create(w64, showWarnings = FALSE)
xlsx <- file.path(p01_source_data_root(), "Other_Datasets", "NYGC_ALS_Consortium", "files", "NYGC_Flagship_ST_4_16_expression.xlsx")

# ---- ST7 proportions, one column per cell type (reference type dropped) ----
st7 <- as.data.frame(read_excel(xlsx, sheet = "ST7"))
names(st7) <- c("rna_id", "cell_type", "estimate", "tissue")
# Some estimates are stored as text in the workbook (e.g. "1.1E-4"); coerce to numeric.
st7$estimate <- as.numeric(st7$estimate)
stopifnot(!anyNA(st7$estimate))
ref_type <- c(Cortex = "Excitatory", `Spinal Cord` = "Neurons")
tissue_of <- c(motor_cortex = "Cortex", spinal_cord = "Spinal Cord")
st7_wide <- function(region) {
  d <- st7[st7$tissue == tissue_of[[region]], ]
  w <- reshape(d[, c("rna_id", "cell_type", "estimate")], idvar = "rna_id", timevar = "cell_type", direction = "wide")
  names(w) <- sub("^estimate\\.", "", names(w))
  w <- w[, setdiff(names(w), ref_type[[tissue_of[[region]]]])]
  names(w)[-1] <- paste0("st7_", make.names(names(w)[-1]))
  w
}

sig <- read.delim(file.path(root, "results", "W59_NYGC_consensus_signature_holdout", "w59_frozen_consensus_signature.tsv"))
gene_sets <- c(modules, list(up_consensus = sig$gene[sig$direction == "up"], down_consensus = sig$gene[sig$direction == "down"]))

fit_w64 <- function(region, design_id) {
  p <- prepare(region)
  w <- st7_wide(region)
  s <- merge(p$s, w, by = "rna_id")
  st7_cols <- grep("^st7_", names(s), value = TRUE)
  s <- droplevels(s[complete.cases(s[, st7_cols]), ])  # drop subregion levels left empty by the ST7 restriction
  y <- DGEList(counts = counts_all[, s$rna_id])
  y <- calcNormFactors(y[filterByExpr(y, group = s$disease), , keep.lib.sizes = FALSE], method = "TMM")
  marker_cols <- paste0(names(marker_ensg), "_score")
  comp <- if (design_id == "a_st7_replaces_markers") st7_cols else c(marker_cols, st7_cols)
  rhs <- c(p$needed, comp)
  design <- model.matrix(as.formula(paste("~ disease +", paste(rhs, collapse = " + "))), data = s)
  v <- voom(y, design)
  cor1 <- duplicateCorrelation(v, design, block = s$donor_id)$consensus.correlation
  v <- voom(y, design, block = s$donor_id, correlation = cor1)
  cons <- duplicateCorrelation(v, design, block = s$donor_id)$consensus.correlation
  fit <- eBayes(lmFit(v, design, block = s$donor_id, correlation = cons), trend = TRUE)
  index <- lapply(gene_sets, intersect, rownames(v))
  fr <- fry(v, index = index, design = design, contrast = "diseaseALS", block = s$donor_id, correlation = cons, sort = "none")
  ca <- camera(v, index = index, design = design, contrast = "diseaseALS", sort = FALSE)
  tab <- topTable(fit, coef = "diseaseALS", number = Inf, sort.by = "none")
  sets <- do.call(rbind, lapply(names(index), function(m) {
    lfc <- tab[index[[m]], "logFC"]
    expected <- if (startsWith(m, "up")) 1 else -1
    data.frame(design = design_id, region = region, gene_set = m, n_samples = ncol(v),
               n_als_donors = length(unique(s$donor_id[s$disease == "ALS"])),
               n_control_donors = length(unique(s$donor_id[s$disease == "Control"])),
               n_composition_covariates = length(comp), n_genes_tested = length(index[[m]]),
               fry_direction = fr[m, "Direction"], fry_p = fr[m, "PValue"], camera_direction = ca[m, "Direction"],
               camera_p = ca[m, "PValue"], n_genes_concordant = sum(sign(lfc) == expected), median_logFC = median(lfc))
  }))
  genes <- data.frame(design = design_id, region = region, gene_id = unlist(index), gene_set = rep(names(index), lengths(index)),
                      tab[unlist(index), c("logFC", "AveExpr", "t", "P.Value")])
  list(sets = sets, genes = genes)
}

runs <- expand.grid(region = c("motor_cortex", "spinal_cord"), design = c("a_st7_replaces_markers", "b_markers_plus_st7"),
                    stringsAsFactors = FALSE)
fits <- Map(fit_w64, runs$region, runs$design)
sets <- do.call(rbind, lapply(fits, `[[`, "sets"))
genes <- do.call(rbind, lapply(fits, `[[`, "genes"))

# Reading (Amendment 8) for the ALS-up module, design (a)
up_a <- sets[sets$design == "a_st7_replaces_markers" & sets$gene_set == "up", ]
hit <- up_a$fry_p < 0.05 & up_a$fry_direction == "Up"
reading <- if (all(hit)) "not explained by published composition" else if (any(hit)) "partly explained" else "not separable from published composition"

# ---- ST4 overlap (context; ST4 includes discovery donors) ----
st4 <- as.data.frame(read_excel(xlsx, sheet = "ST4"))
names(st4)[1] <- "gene_versioned"
st4$gene_id <- sub("\\..*$", "", st4$gene_versioned)
num_cols <- grep("log2FoldChange|pvalue|padj", names(st4), value = TRUE)
st4[num_cols] <- lapply(st4[num_cols], as.numeric)  # some values are stored as text in the workbook
mod_tab <- data.frame(module = rep(names(modules), lengths(modules)), gene_id = unlist(modules))
mod_tab <- merge(mod_tab, st4[, c("gene_id", "gene_symbol", "MCX.log2FoldChange", "MCX.padj", "CSC.log2FoldChange", "CSC.padj")], by = "gene_id", all.x = TRUE)
dir_ok <- function(lfc, m) ifelse(m == "up", lfc > 0, lfc < 0)
mod_tab$MCX_DE_hypothesised <- !is.na(mod_tab$MCX.padj) & mod_tab$MCX.padj < 0.05 & dir_ok(mod_tab$MCX.log2FoldChange, mod_tab$module)
mod_tab$CSC_DE_hypothesised <- !is.na(mod_tab$CSC.padj) & mod_tab$CSC.padj < 0.05 & dir_ok(mod_tab$CSC.log2FoldChange, mod_tab$module)
overlap <- aggregate(cbind(MCX_DE_hypothesised, CSC_DE_hypothesised) ~ module, data = mod_tab, FUN = sum)
overlap$n_genes <- as.vector(table(mod_tab$module)[overlap$module])

write.table(sets, file.path(w64, "w64_set_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes, file.path(w64, "w64_gene_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mod_tab, file.path(w64, "w64_module_genes_in_nygc_st4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(overlap, file.path(w64, "w64_st4_overlap_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("reading (ALS-up, design a):", reading)), file.path(w64, "w64_reading.txt"))
print(sets[, c("design", "region", "gene_set", "n_samples", "n_als_donors", "n_control_donors", "fry_direction", "fry_p", "camera_p", "n_genes_concordant", "n_genes_tested", "median_logFC")], digits = 3, row.names = FALSE)
print(overlap, row.names = FALSE)
cat("reading (ALS-up, design a):", reading, "\n")
