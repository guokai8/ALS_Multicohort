#!/usr/bin/env Rscript

# W71 A (Amendment 14): territory x genotype interaction of the frozen ALS-up module (mouse orthologs)
# in SOD1-G93A vs SOD1-WT lumbar spinal cord (Maniatis et al. 2019) at P30, P70, P100 and P120.
# Primary family: frozen-module interaction at the four ages (BH over four). Secondary: consensus up
# arm, microglia-adjusted interaction (W70 method), gene-level effects, white-matter comparison.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(babelgene)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
root <- p01_project_root()
res <- file.path(root, "results")
out <- file.path(res, "W71_maniatis_spatial_time_course")
inp <- file.path(out, "inputs")

# ---- gene sets: babelgene mouse orthologs with support from >= 3 sources ----
to_mouse <- function(h) {
  o <- orthologs(h, species = "mouse", human = TRUE)
  unique(o$symbol[o$support_n >= 3])
}
seed <- file.path(res, "W34_synapse_resolved_regulatory_reanalysis", "inputs", "rna_seed_lists")
up <- read.delim(file.path(seed, "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_up.tsv"))$common_name
sig <- read.delim(file.path(res, "W59_NYGC_consensus_signature_holdout", "w59_frozen_consensus_signature.tsv"))
sets <- list(up = to_mouse(up), up_consensus = to_mouse(sig$symbol[sig$direction == "up"]))
MICROGLIA <- to_mouse(c("P2RY12", "CX3CR1", "CSF1R", "PTPRC", "CD68", "TREM2", "TYROBP", "CD163", "ITGAX"))
GENES <- c("Gpnmb", "Anxa1", "Ltf", "Emp3", "C3")
stopifnot(length(sets$up) == 11, length(sets$up_consensus) == 74)

counts <- as.matrix(read.delim(gzfile(file.path(inp, "w71_mouse_unit_counts.tsv.gz")), row.names = 1, check.names = FALSE))
units <- read.delim(file.path(inp, "w71_mouse_units.tsv"), stringsAsFactors = FALSE)
units$animal <- as.character(units$animal)
COMPARISONS <- list(grey = c(motor = "Vent_Horn", sensory = "Dors_Horn"),
                    white = c(motor = "Vent_Lat_White", sensory = "Dors_Med_White"))

fit_age <- function(age, comparison, adjust = FALSE) {
  reg <- COMPARISONS[[comparison]]
  u <- units[units$age == age & units$region %in% reg & units$n_spots >= 5 & units$lib_size > 0, ]
  u <- u[u$animal %in% names(which(table(u$animal) == 2)), ]
  x <- counts[, u$unit]
  meta <- data.frame(animal = u$animal,
                     territory = factor(ifelse(u$region == reg[["motor"]], "motor", "sensory"), levels = c("sensory", "motor")),
                     genotype = factor(u$genotype, levels = c("WT", "G93A")), sex = factor(u$sex))
  rhs <- if (nlevels(droplevels(meta$sex)) > 1) "+ sex" else ""
  if (adjust) {
    lc <- cpm(calcNormFactors(DGEList(counts = x), method = "TMM"), log = TRUE, prior.count = 1)
    meta$mg <- colMeans(lc[intersect(MICROGLIA, rownames(lc)), , drop = FALSE]); meta$mg <- meta$mg - mean(meta$mg)
    rhs <- paste(rhs, "+ mg")
  }
  design <- model.matrix(as.formula(paste("~ territory * genotype", rhs)), data = meta)
  colnames(design) <- make.names(colnames(design))
  y <- DGEList(counts = x)
  y <- calcNormFactors(y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(y, design)
  cor1 <- duplicateCorrelation(v, design, block = meta$animal)$consensus.correlation
  v <- voom(y, design, block = meta$animal, correlation = cor1)
  cons <- duplicateCorrelation(v, design, block = meta$animal)$consensus.correlation
  fit <- lmFit(v, design, block = meta$animal, correlation = cons)
  cm <- cbind(interaction = setNames(as.numeric(colnames(design) == "territorymotor.genotypeG93A"), colnames(design)),
              ventral = setNames(as.numeric(colnames(design) %in% c("genotypeG93A", "territorymotor.genotypeG93A")), colnames(design)),
              dorsal = setNames(as.numeric(colnames(design) == "genotypeG93A"), colnames(design)))
  fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE)
  index <- lapply(sets, intersect, rownames(v))
  tests <- do.call(rbind, lapply(colnames(cm), function(k) {
    fr <- fry(v, index = index, design = design, contrast = cm[, k], block = meta$animal, correlation = cons, sort = "none")
    ca <- camera(v, index = index, design = design, contrast = cm[, k], sort = FALSE)
    tab <- topTable(fit2, coef = k, number = Inf, sort.by = "none")
    do.call(rbind, lapply(names(index), function(m) data.frame(
      age = age, comparison = comparison, adjusted = adjust, contrast = k, gene_set = m, n_genes_tested = length(index[[m]]),
      n_G93A = length(unique(meta$animal[meta$genotype == "G93A"])), n_WT = length(unique(meta$animal[meta$genotype == "WT"])),
      fry_direction = fr[m, "Direction"], fry_p = fr[m, "PValue"], camera_p = ca[m, "PValue"],
      median_logFC = median(tab[index[[m]], "logFC"]))))
  }))
  genes <- do.call(rbind, lapply(colnames(cm), function(k) {
    tab <- topTable(fit2, coef = k, number = Inf, sort.by = "none")
    g <- intersect(sets$up, rownames(tab))
    data.frame(age = age, comparison = comparison, adjusted = adjust, contrast = k, gene = g,
               logFC = tab[g, "logFC"], P = tab[g, "P.Value"])
  }))
  list(tests = tests, genes = genes)
}

ages <- c("p30", "p70", "p100", "p120")
runs <- unlist(lapply(ages, function(a) lapply(names(COMPARISONS), function(cmp) fit_age(a, cmp))), recursive = FALSE)
tests <- do.call(rbind, lapply(runs, `[[`, "tests"))
genes <- do.call(rbind, lapply(runs, `[[`, "genes"))

# primary family: frozen-module grey-matter interaction at four ages; consensus arm BH separately
for (gs in names(sets)) {
  i <- tests$comparison == "grey" & tests$contrast == "interaction" & tests$gene_set == gs & !tests$adjusted
  tests$fry_BH_family[i] <- p.adjust(tests$fry_p[i], "BH")
}

# microglia adjustment where the unadjusted frozen-module interaction has fry P < 0.05
prim <- tests[tests$comparison == "grey" & tests$contrast == "interaction" & tests$gene_set == "up" & !tests$adjusted, ]
adj_ages <- prim$age[prim$fry_p < 0.05]
adj <- do.call(rbind, lapply(adj_ages, function(a) fit_age(a, "grey", adjust = TRUE)$tests))
if (!is.null(adj)) {
  adj <- adj[adj$contrast == "interaction", ]
  ref <- tests[tests$comparison == "grey" & tests$contrast == "interaction" & !tests$adjusted, c("age", "gene_set", "median_logFC")]
  adj <- merge(adj, setNames(ref, c("age", "gene_set", "unadjusted_median_logFC")), by = c("age", "gene_set"))
  adj$attenuation <- 1 - adj$median_logFC / adj$unadjusted_median_logFC
}

# ---- prespecified readings ----
late <- prim[prim$age %in% c("p100", "p120"), ]
rep_reading <- if (any(late$fry_direction == "Up" & late$fry_BH_family < 0.05)) "replicated in SOD1-G93A" else
  if (any(late$fry_direction == "Up" & late$fry_p < 0.05)) "nominal" else "not replicated"
vh <- tests[tests$comparison == "grey" & tests$contrast == "ventral" & tests$gene_set == "up" & !tests$adjusted, ]
vh <- vh[match(ages, vh$age), ]
first <- vh$age[which(vh$fry_direction == "Up" & vh$fry_p < 0.05)[1]]
timing <- if (is.na(first)) "no ventral-horn response" else
  c(p30 = "presymptomatic", p70 = "at onset", p100 = "symptomatic or late", p120 = "symptomatic or late")[[first]]
reading <- data.frame(question = c("cross-species replication", "timing"), reading = c(rep_reading, timing),
                      basis = c("frozen ALS-up grey interaction, P100/P120", paste("earliest Up ventral-horn effect:", ifelse(is.na(first), "none", first))))

dir.create(out, showWarnings = FALSE)
write.table(tests, file.path(out, "w71_mouse_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes, file.path(out, "w71_mouse_gene_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(adj)) write.table(adj, file.path(out, "w71_mouse_microglia_adjusted.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(reading, file.path(out, "w71_mouse_reading.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(tests[!tests$adjusted, c("age", "comparison", "contrast", "gene_set", "n_G93A", "n_WT", "fry_direction", "fry_p", "fry_BH_family", "camera_p", "median_logFC")], digits = 3, row.names = FALSE)
if (!is.null(adj)) print(adj[, c("age", "gene_set", "fry_direction", "fry_p", "median_logFC", "unadjusted_median_logFC", "attenuation")], digits = 3, row.names = FALSE)
g <- genes[genes$gene %in% GENES & genes$comparison == "grey" & genes$contrast == "ventral", ]
print(reshape(g[, c("age", "gene", "logFC")], idvar = "gene", timevar = "age", direction = "wide"), digits = 2, row.names = FALSE)
print(reading, row.names = FALSE)
