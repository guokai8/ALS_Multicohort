#!/usr/bin/env Rscript

# W63 pathway context (Amendment 8; exploratory). Preranked GSEA of the W58 run (ii)
# meta-analysis statistic and over-representation of the consensus arms, on the
# prespecified themed gene sets; Hallmark GSEA as a broader view; theme membership of
# the frozen module genes; NYGC Flagship ST5 GSEA for the themed GO terms as context.

suppressPackageStartupMessages({
  library(msigdbr)
  library(fgsea)
  library(readxl)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
res_dir <- file.path(p01_project_root(), "results")
out <- file.path(res_dir, "W63_pathway_context")
dir.create(out, showWarnings = FALSE)
set.seed(20261012)

themes <- list(
  Complement = c("HALLMARK_COMPLEMENT", "GOBP_COMPLEMENT_ACTIVATION"),
  `Extracellular matrix` = c("GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  `Vascular and BBB` = c("HALLMARK_ANGIOGENESIS", "GOBP_BLOOD_VESSEL_MORPHOGENESIS", "GOBP_MAINTENANCE_OF_BLOOD_BRAIN_BARRIER"),
  `Inflammation and glial activation` = c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                                          "HALLMARK_INTERFERON_GAMMA_RESPONSE", "GOBP_GLIAL_CELL_ACTIVATION", "GOBP_ASTROCYTE_ACTIVATION"),
  `Homeostatic and neuronal support` = c("GOBP_SYNAPSE_ORGANIZATION", "GOBP_AMINO_ACID_TRANSPORT")
)
theme_of <- setNames(rep(names(themes), lengths(themes)), unlist(themes))
msig <- rbind(msigdbr(species = "Homo sapiens", collection = "H")[, c("gs_name", "gene_symbol")],
              msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")[, c("gs_name", "gene_symbol")])
sets_all <- lapply(split(msig$gene_symbol, msig$gs_name), unique)
stopifnot(all(unlist(themes) %in% names(sets_all)))

# ---- W58 run (ii) statistic ----
gw <- read.delim(gzfile(file.path(res_dir, "W58_multicohort_meta", "w58_genome_wide_meta.tsv.gz")), na.strings = "NA")
gw <- gw[gw$analysis == "ii_without_discovery" & !is.na(gw$meta_p) & !is.na(gw$symbol) & gw$symbol != "", ]
gw <- gw[!duplicated(gw$symbol), ]
stat <- setNames(gw$meta_logFC / gw$meta_se, gw$symbol)
universe <- names(stat)

# ---- 1. preranked GSEA: themed family (primary) and Hallmark ----
themed <- fgsea(pathways = sets_all[unlist(themes)], stats = stat, minSize = 5, maxSize = 2000, nPermSimple = 10000)
themed_tab <- data.frame(family = "themed", theme = theme_of[themed$pathway], pathway = themed$pathway, size = themed$size,
                         NES = themed$NES, pval = themed$pval, padj = p.adjust(themed$pval, method = "BH"),
                         leading_edge = vapply(themed$leadingEdge, function(x) paste(head(x, 15), collapse = ";"), ""))
hall_names <- grep("^HALLMARK_", names(sets_all), value = TRUE)
hall <- fgsea(pathways = sets_all[hall_names], stats = stat, minSize = 5, maxSize = 2000, nPermSimple = 10000)
hall_tab <- data.frame(family = "hallmark", theme = ifelse(hall$pathway %in% names(theme_of), theme_of[hall$pathway], ""),
                       pathway = hall$pathway, size = hall$size, NES = hall$NES, pval = hall$pval,
                       padj = p.adjust(hall$pval, method = "BH"),
                       leading_edge = vapply(hall$leadingEdge, function(x) paste(head(x, 15), collapse = ";"), ""))

# ---- 2. over-representation of the consensus arms ----
cons <- gw[gw$consensus %in% c(TRUE, "TRUE"), ]
arms <- list(consensus_up = cons$symbol[cons$meta_logFC > 0], consensus_down = cons$symbol[cons$meta_logFC < 0])
ora <- do.call(rbind, lapply(names(arms), function(a) {
  hits <- arms[[a]]
  rows <- lapply(unlist(themes), function(s) {
    set <- intersect(sets_all[[s]], universe)
    k <- sum(hits %in% set)
    tab <- matrix(c(k, length(hits) - k, length(set) - k, length(universe) - length(set) - length(hits) + k), 2)
    data.frame(arm = a, theme = theme_of[[s]], pathway = s, set_size_in_universe = length(set), n_arm = length(hits), overlap = k,
               fold = (k / length(hits)) / (length(set) / length(universe)),
               p = fisher.test(tab, alternative = "greater")$p.value, genes = paste(sort(hits[hits %in% set]), collapse = ";"))
  })
  r <- do.call(rbind, rows)
  r$padj <- p.adjust(r$p, method = "BH")
  r
}))

# ---- 3. theme membership of the frozen module genes ----
seed_dir <- file.path(res_dir, "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
mods <- do.call(rbind, lapply(c("up", "down"), function(m) {
  x <- read.delim(file.path(seed_dir, paste0("GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_", m, ".tsv")))
  data.frame(module = m, symbol = x$common_name)
}))
mods$themes <- vapply(mods$symbol, function(g) paste(unique(theme_of[names(theme_of)[vapply(names(theme_of),
  function(s) g %in% sets_all[[s]], logical(1))]]), collapse = ";"), "")

# ---- 4. NYGC Flagship ST5 GSEA for the themed GO terms (context) ----
st5_path <- file.path(p01_source_data_root(), "Other_Datasets", "NYGC_ALS_Consortium", "files", "NYGC_Flagship_ST_4_16_expression.xlsx")
st5 <- as.data.frame(read_excel(st5_path, sheet = "ST5"))
go_terms <- tolower(gsub("_", " ", sub("^GOBP_", "", grep("^GOBP_", unlist(themes), value = TRUE))))
st5_ctx <- st5[tolower(st5$Description) %in% go_terms, , drop = FALSE]

write.table(themed_tab, file.path(out, "w63_themed_gsea.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(hall_tab[order(hall_tab$pval), ], file.path(out, "w63_hallmark_gsea.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ora, file.path(out, "w63_consensus_ora.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mods, file.path(out, "w63_module_gene_themes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(st5_ctx, file.path(out, "w63_nygc_st5_themed_go_context.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("universe:", length(universe), "| consensus up/down:", lengths(arms), "\n")
print(themed_tab[order(themed_tab$pval), c("theme", "pathway", "size", "NES", "pval", "padj")], digits = 3, row.names = FALSE)
print(head(hall_tab[order(hall_tab$pval), c("pathway", "NES", "pval", "padj")], 12), digits = 3, row.names = FALSE)
print(ora[ora$p < 0.05, c("arm", "pathway", "overlap", "fold", "p", "padj", "genes")], digits = 3, row.names = FALSE)
print(mods, row.names = FALSE)
cat("ST5 context rows:", nrow(st5_ctx), "\n")
