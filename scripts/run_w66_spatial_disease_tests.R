#!/usr/bin/env Rscript

# W66 B/C (Amendment 9; exploratory): frozen modules and consensus arms in GSE330224 Visium
# donor-level pseudobulk (LSC lumbar spinal cord, MTC motor cortex). One profile per donor and
# unit (whole section, anatomical region group, lineage niche). Primary: donors not in GSE330130,
# ALS (Level 1) vs non-neurological control, limma-voom ~ disease + sex + age; blocked by nothing
# (one profile per donor per unit); fry primary, camera reported.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
inp <- file.path(root, "results", "W66_spatial_context", "inputs")
out <- file.path(root, "results", "W66_spatial_context")

seed <- file.path(root, "results", "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
sets <- lapply(c(up = "up", down = "down"), function(m)
  read.delim(file.path(seed, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", m, ".tsv")))$gene_id)
sig <- read.delim(file.path(root, "results", "W59_NYGC_consensus_signature_holdout", "w59_frozen_consensus_signature.tsv"))
sets$up_consensus <- sig$gene[sig$direction == "up"]
sets$down_consensus <- sig$gene[sig$direction == "down"]

fit_unit <- function(counts, meta, analysis, tissue, unit_type, stratum, formula_rhs) {
  meta$disease <- factor(ifelse(meta$group == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
  meta$sex <- factor(meta$sex)
  meta <- meta[complete.cases(meta[, all.vars(as.formula(paste("~", formula_rhs)))]), ]  # complete cases (Amendment 9)
  y <- DGEList(counts = counts[, meta$unit, drop = FALSE])
  design <- model.matrix(as.formula(paste("~ disease +", formula_rhs)), data = meta)
  y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(y, design)
  fit <- eBayes(lmFit(v, design), trend = TRUE)
  tab <- topTable(fit, coef = "diseaseALS", number = Inf, sort.by = "none")
  index <- lapply(sets, intersect, rownames(v))
  fr <- fry(v, index = index, design = design, contrast = "diseaseALS", sort = "none")
  ca <- camera(v, index = index, design = design, contrast = "diseaseALS", sort = FALSE)
  res <- do.call(rbind, lapply(names(index), function(m) {
    lfc <- tab[index[[m]], "logFC"]
    expected <- if (startsWith(m, "up")) 1 else -1
    data.frame(analysis = analysis, tissue = tissue, unit_type = unit_type, stratum = stratum, gene_set = m,
               design = paste("~ disease +", formula_rhs), n_als = sum(meta$disease == "ALS"),
               n_control = sum(meta$disease == "Control"), median_spots = median(meta$n_spots),
               n_genes_tested = length(index[[m]]), fry_direction = fr[m, "Direction"], fry_p = fr[m, "PValue"],
               camera_direction = ca[m, "Direction"], camera_p = ca[m, "PValue"],
               n_genes_concordant = sum(sign(lfc) == expected), median_logFC = median(lfc))
  }))
  genes <- data.frame(analysis = analysis, tissue = tissue, unit_type = unit_type, stratum = stratum,
                      gene_id = unlist(index[c("up", "down")]), gene_set = rep(c("up", "down"), lengths(index[c("up", "down")])),
                      tab[unlist(index[c("up", "down")]), c("logFC", "AveExpr", "t", "P.Value")])
  list(sets = res, genes = genes)
}

runs <- list()
for (tissue in c("LSC", "MTC")) {
  meta_all <- read.delim(file.path(inp, paste0(tissue, "_units.tsv")), stringsAsFactors = FALSE)
  counts <- read.delim(gzfile(file.path(inp, paste0(tissue, "_pseudobulk_counts.tsv.gz"))), row.names = 1, check.names = FALSE)
  counts <- as.matrix(counts[, -1])
  meta_all <- meta_all[meta_all$group %in% c("ALS", "Control") & meta_all$lib_size > 0, ]
  meta_all$gse330130_donor <- meta_all$gse330130_donor %in% c(TRUE, "True", "TRUE")  # written by pandas
  indep <- !meta_all$gse330130_donor
  strata <- unique(meta_all[, c("unit_type", "stratum")])
  for (i in seq_len(nrow(strata))) {
    ut <- strata$unit_type[i]; st <- strata$stratum[i]
    m <- meta_all[meta_all$unit_type == ut & meta_all$stratum == st, ]
    add <- function(analysis, keep, rhs) {
      mm <- m[keep, ]
      if (sum(mm$group == "ALS") < 3 || sum(mm$group == "Control") < 3) return(NULL)
      runs[[length(runs) + 1]] <<- fit_unit(counts, mm, analysis, tissue, ut, st, rhs)
    }
    add("primary", indep[meta_all$unit_type == ut & meta_all$stratum == st], "sex + age")
    if (ut == "whole") {
      ind <- !m$gse330130_donor
      add("sens_ALS_only", ind & (m$diagnosis == "ALS" | m$group == "Control"), "sex + age")
      add("sens_PMI", ind & !is.na(m$pmi), "sex + age + pmi")
      add("sens_all_donors", rep(TRUE, nrow(m)), "sex + age")
    }
  }
}
sets_tab <- do.call(rbind, lapply(runs, `[[`, "sets"))
genes_tab <- do.call(rbind, lapply(runs, `[[`, "genes"))

# Primary family (Amendment 9): frozen up and down modules, whole-section units, primary analysis, both tissues.
fam <- sets_tab$analysis == "primary" & sets_tab$unit_type == "whole" & sets_tab$gene_set %in% c("up", "down")
sets_tab$fry_BH_primary_family <- NA_real_
sets_tab$fry_BH_primary_family[fam] <- p.adjust(sets_tab$fry_p[fam], method = "BH")

write.table(sets_tab, file.path(out, "w66_set_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes_tab, file.path(out, "w66_module_gene_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(sets_tab[sets_tab$gene_set %in% c("up", "down"),
               c("analysis", "tissue", "unit_type", "stratum", "gene_set", "n_als", "n_control", "fry_direction", "fry_p",
                 "camera_p", "n_genes_concordant", "n_genes_tested", "median_logFC", "fry_BH_primary_family")],
      digits = 3, row.names = FALSE)
