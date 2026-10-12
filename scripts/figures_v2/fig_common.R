# Shared theme, palette and export helpers for the v2 manuscript figures.
# All figures are drawn and exported with R (ggplot2 + patchwork).

invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))  # Rscript may start in the C locale; labels use Unicode

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

fig_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
}
FIG_ROOT <- normalizePath(file.path(fig_script_dir(), "..", ".."), mustWork = FALSE)
RES <- file.path(FIG_ROOT, "results")
FIG_OUT <- file.path(RES, "Figures_v2")

theme_set(
  theme_classic(base_size = 6.5, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(size = 5.8, colour = "black"),
      axis.title = element_text(size = 6.4),
      legend.title = element_text(size = 6.2),
      legend.text = element_text(size = 5.8),
      legend.key.size = unit(3, "mm"),
      strip.text = element_text(size = 6.2, face = "bold"),
      strip.background = element_blank(),
      plot.title = element_text(size = 7, face = "bold"),
      plot.tag = element_text(size = 8, face = "bold"),
      panel.grid = element_blank()
    )
)

# One neutral, one signal and one accent family, reused in every figure.
PAL <- list(
  als = "#B5523B", control = "#4F6D8A", neutral = "#8C8C8C", light = "#D9D9D9",
  up = "#B5523B", down = "#4F6D8A", accent = "#D9A441", ink = "#222222"
)
REGION_SHAPE <- c(motor_cortex = 16, spinal_cord = 17, prefrontal_cortex = 15)

new_figure_dir <- function(name) {
  d <- file.path(FIG_OUT, name)
  dir.create(file.path(d, "source_data"), recursive = TRUE, showWarnings = FALSE)
  d
}

write_source <- function(df, dir, panel) {
  write.table(df, file.path(dir, "source_data", paste0(panel, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(df)
}

save_pub_r <- function(plot, dir, name, width_mm = 183, height_mm = 120, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  base <- file.path(dir, name)
  svglite::svglite(paste0(base, ".svg"), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(base, "_preview.png"), width = w, height = h, units = "in", res = 300)
  print(plot)
  dev.off()
}

neglog10 <- function(p) -log10(pmax(p, 1e-300))
