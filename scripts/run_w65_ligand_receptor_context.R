#!/usr/bin/env Rscript

# W65 (Amendment 8; exploratory, hypothesis-generating): ligand-receptor context of the
# frozen module genes in the W56 lineage pseudobulk of Pineda BA9 batch 210526 and Ulm.
# Pairs from the liana Consensus resource in which a module gene is the ligand or the
# receptor (or a subunit). For every sender-receiver lineage combination: expression
# (mean CPM >= 1 for ligand in sender and receptor in receiver; complexes use the
# minimum over subunits) and the disease change (ligand logFC in sender + receptor
# logFC in receiver, from per-lineage limma-voom ~ disease + sex). No P values are used.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
w56 <- file.path(root, "results", "W56_cell_source_check")
out <- file.path(root, "results", "W65_ligand_receptor_context")
dir.create(out, showWarnings = FALSE)
CPM_MIN <- 1

seed_dir <- file.path(root, "results", "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
mods <- do.call(rbind, lapply(c("up", "down"), function(m) {
  x <- read.delim(file.path(seed_dir, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", m, ".tsv")))
  data.frame(module = m, symbol = x$common_name)
}))
module_of <- setNames(mods$module, mods$symbol)

lr <- liana::select_resource("Consensus")[[1]]
lr <- unique(as.data.frame(lr)[, c("source_genesymbol", "target_genesymbol")])
names(lr) <- c("ligand", "receptor")
subunits <- function(x) strsplit(x, "_", fixed = TRUE)[[1]]
lr$involves <- vapply(seq_len(nrow(lr)), function(i) paste(intersect(c(subunits(lr$ligand[i]), subunits(lr$receptor[i])), mods$symbol), collapse = ";"), "")
lr <- lr[lr$involves != "", ]

lineage_key <- c(Astrocyte = "Astrocyte", Astrocytes = "Astrocyte", Microglia = "Microglia", Endothelial = "Endothelial",
                 `Vascular/Meningeal` = "Vascular/Meningeal", Oligodendrocyte = "Oligodendrocyte", Oligodendrocytes = "Oligodendrocyte",
                 OPC = "OPC", `Excitatory neuron` = "Excitatory neuron", Exc_Neurons = "Excitatory neuron",
                 `Inhibitory neuron` = "Inhibitory neuron", Inh_Neurons = "Inhibitory neuron")

# Per-lineage mean CPM and disease logFC (gene symbols as row names).
lineage_stats <- function(counts, strata) {
  do.call(rbind, lapply(split(strata, strata$lineage), function(s) {
    y <- DGEList(counts = counts[, s$unit, drop = FALSE])
    mean_cpm <- rowMeans(cpm(y))
    design <- model.matrix(~ disease + sex, data = s)
    keep <- filterByExpr(y, design = design)
    y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
    fit <- eBayes(lmFit(voom(y, design), design), trend = TRUE)
    lfc <- setNames(rep(NA_real_, nrow(counts)), rownames(counts))
    lfc[rownames(fit)] <- fit$coefficients[, "diseaseALS"]
    data.frame(lineage = lineage_key[[s$lineage[1]]], gene = rownames(counts), mean_cpm = mean_cpm, logFC = lfc,
               n_als = sum(s$disease == "ALS"), n_control = sum(s$disease == "Control"))
  }))
}

# ---- Pineda BA9 batch 210526: one profile per donor and lineage ----
p_counts <- as.matrix(read.delim(gzfile(file.path(w56, "w56_pineda_lineage_counts.tsv.gz")), row.names = 1, check.names = FALSE))
gene_map <- read.delim(file.path(p01_source_data_root(), "Other_Datasets/Pineda_snRNAseq/synapse_download/Processed Data/191112_ALS_110_snRNA-B9/row_metadata.tsv"), stringsAsFactors = FALSE)
sym <- gene_map$Gene[match(rownames(p_counts), gene_map$ENSEMBL)]
keep_rows <- !is.na(sym) & !duplicated(sym)
p_counts <- p_counts[keep_rows, ]
rownames(p_counts) <- sym[keep_rows]
p_strata <- read.delim(file.path(w56, "w56_pineda_lineage_strata.tsv"), stringsAsFactors = FALSE)
p_strata <- p_strata[p_strata$Region == "BA9" & startsWith(p_strata$sample_folder, "210526") &
                       p_strata$Condition %in% c("ALS", "PN") & p_strata$lib_size > 0, ]
p_strata$unit <- p_strata$stratum
p_strata$disease <- factor(ifelse(p_strata$Condition == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
p_strata$sex <- factor(p_strata$Sex)
pineda <- lineage_stats(p_counts, p_strata)

# ---- Ulm: sum donor x chip profiles to one profile per donor and lineage ----
u_counts <- as.matrix(read.delim(gzfile(file.path(w56, "w56_ulm_lineage_counts.tsv.gz")), row.names = 1, check.names = FALSE))
u_strata <- read.delim(file.path(w56, "w56_ulm_lineage_strata.tsv"), stringsAsFactors = FALSE)
u_strata <- u_strata[u_strata$Case_Type %in% c("ALS", "HC") & u_strata$lib_size > 0, ]
u_strata$unit <- paste(u_strata$ID, u_strata$lineage, sep = "|")
grp <- split(u_strata$stratum, u_strata$unit)
u_donor <- sapply(grp, function(cols) rowSums(u_counts[, cols, drop = FALSE]))
u_meta <- u_strata[!duplicated(u_strata$unit), c("unit", "ID", "lineage", "Sex", "Case_Type")]
u_meta <- u_meta[match(colnames(u_donor), u_meta$unit), ]
u_meta$disease <- factor(ifelse(u_meta$Case_Type == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))
u_meta$sex <- factor(tolower(u_meta$Sex))
ulm <- lineage_stats(u_donor, u_meta)

# ---- pair scores per cohort ----
score_pairs <- function(st, cohort) {
  lineages <- unique(st$lineage)
  get <- function(genes, lin, col) {
    v <- st[st$lineage == lin & st$gene %in% genes, col]
    if (length(v) < length(genes) || anyNA(v)) return(NA_real_)
    if (col == "mean_cpm") min(v) else mean(v)
  }
  rows <- list()
  for (i in seq_len(nrow(lr))) {
    lg <- subunits(lr$ligand[i]); rc <- subunits(lr$receptor[i])
    for (snd in lineages) for (rcv in lineages) {
      lc <- get(lg, snd, "mean_cpm"); rcp <- get(rc, rcv, "mean_cpm")
      if (is.na(lc) || is.na(rcp) || lc < CPM_MIN || rcp < CPM_MIN) next
      rows[[length(rows) + 1]] <- data.frame(cohort = cohort, ligand = lr$ligand[i], receptor = lr$receptor[i],
        module_genes = lr$involves[i], sender = snd, receiver = rcv, ligand_cpm = lc, receptor_cpm = rcp,
        ligand_logFC = get(lg, snd, "logFC"), receptor_logFC = get(rc, rcv, "logFC"))
    }
  }
  r <- do.call(rbind, rows)
  r$score_change <- r$ligand_logFC + r$receptor_logFC
  r
}
pairs <- rbind(score_pairs(pineda, "PinedaBA9_210526"), score_pairs(ulm, "Ulm"))
pairs$module <- vapply(pairs$module_genes, function(g) paste(unique(module_of[strsplit(g, ";")[[1]]]), collapse = ";"), "")

key <- c("ligand", "receptor", "sender", "receiver")
both <- merge(pairs[pairs$cohort == "PinedaBA9_210526", c(key, "module_genes", "module", "score_change", "ligand_logFC", "receptor_logFC")],
              pairs[pairs$cohort == "Ulm", c(key, "score_change", "ligand_logFC", "receptor_logFC")], by = key, suffixes = c("_pineda", "_ulm"))
both$same_sign <- !is.na(both$score_change_pineda) & !is.na(both$score_change_ulm) & sign(both$score_change_pineda) == sign(both$score_change_ulm)
both <- both[order(-both$same_sign, -abs(both$score_change_pineda + both$score_change_ulm)), ]

write.table(pairs, file.path(out, "w65_expressed_pairs_by_cohort.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(both, file.path(out, "w65_pairs_in_both_cohorts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lr, file.path(out, "w65_module_gene_lr_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("LR pairs involving module genes:", nrow(lr), "\n")
print(table(pairs$cohort))
cat("pairs x lineage combos in both cohorts:", nrow(both), "| same sign:", sum(both$same_sign), "\n")
print(head(both[both$same_sign, c(key, "module_genes", "score_change_pineda", "score_change_ulm")], 25), digits = 2, row.names = FALSE)
