#!/usr/bin/env Rscript

# W68 C (Amendment 11; descriptive): spot-level Spearman correlation between the frozen ALS-up module
# score and CLR-transformed protein counts (CD68, CD163, ITGAM, HLA-DRA, PTPRC, PECAM1, ACTA2, VIM)
# within each GSE288365 section, over all annotated spots and over anterior-horn spots; raw and
# with ranks residualised on log Spatial counts. Median over sections.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "p01_paths.R"))
out <- file.path(p01_project_root(), "results", "W68_GSE288365_replication")
s <- read.delim(gzfile(file.path(out, "inputs", "w68_spot_table.tsv.gz")), check.names = FALSE)
proteins <- c("CD68.1", "CD163.1", "ITGAM.1", "HLA-DRA", "PTPRC.1", "PECAM1.1", "ACTA2.1", "VIM.1")
annotated <- c("Anterior_horns", "Grey_matter", "White_matter", "Posterior_white_matter")
s <- s[s$manual %in% annotated & s$total_counts > 0, ]
resid <- function(v, d) residuals(lm(v ~ d))
rows <- list()
for (scope in c("all_annotated", "anterior_horn")) {
  ss <- if (scope == "anterior_horn") s[s$manual == "Anterior_horns", ] else s
  for (smp in unique(ss$sample)) {
    d <- ss[ss$sample == smp, ]
    rd <- rank(log(d$total_counts)); ry <- rank(d$module_up)
    for (p in proteins) {
      rx <- rank(d[[p]])
      rows[[length(rows) + 1]] <- data.frame(scope = scope, sample = smp, donor = d$donor[1], diagnosis = d$diagnosis[1], protein = p,
                                             n_spots = nrow(d), rho = cor(ry, rx), rho_depth_partial = cor(resid(ry, rd), resid(rx, rd)))
    }
  }
}
per <- do.call(rbind, rows)
summ <- aggregate(cbind(rho, rho_depth_partial) ~ scope + protein, data = per, FUN = median)
summ$frac_positive <- aggregate(rho ~ scope + protein, data = per, FUN = function(r) mean(r > 0))$rho
write.table(per, file.path(out, "w68_protein_colocalisation_by_section.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summ, file.path(out, "w68_protein_colocalisation_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ[order(summ$scope, -summ$rho), ], digits = 2, row.names = FALSE)
