#!/usr/bin/env Rscript
#
# visualize.R — Generate Manhattan plot (genome-wide) and per-candidate
# multi-line plots from normalized coverage data.
#
# Extracted from ECCsplorer/lib/eccDNA_Rcodes.py:
#   rmanhattan_plot     — genome-wide chromosome Manhattan plot
#   rline_multiplot     — per-candidate 6-panel multi-line plots (TR + CO)
#   rline_multiplot_noco — per-candidate 3-panel plots (TR only)
#
# Usage:
#   visualize.R --data <normalized_coverage.csv> \
#               --outdir <output_directory> \
#               --prefix <output_prefix> \
#               [--window_size <bp>] \
#               [--hiconf <hiconf_win.bed>] \
#               [--mode manhattan|line|both]

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(dplyr))

# ---- Argument parsing ----
args <- commandArgs(trailingOnly = TRUE)

data_path <- NULL
outdir <- NULL
prefix <- "eccsplorer"
window_size <- 100
hiconf_path <- NULL
mode <- "both"
image_res <- 300
image_width <- 1625
image_height <- 925
image_type <- "png"
image_points <- 10

i <- 1
while (i <= length(args)) {
  if (args[i] == "--data") {
    data_path <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--outdir") {
    outdir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--prefix") {
    prefix <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--window_size") {
    window_size <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--hiconf") {
    hiconf_path <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--mode") {
    mode <- args[i + 1]
    i <- i + 2
  } else {
    cat(sprintf("WARNING: Unknown argument: %s\n", args[i]), file = stderr())
    i <- i + 1
  }
}

if (is.null(data_path) || is.null(outdir)) {
  cat("USAGE: visualize.R --data <normalized.csv> --outdir <DIR> [--prefix <STR>] [--window_size <INT>] [--hiconf <BED>] [--mode manhattan|line|both]\n",
      file = stderr())
  quit(status = 1)
}

if (!file.exists(data_path)) {
  cat(sprintf("ERROR: Data file '%s' not found.\n", data_path), file = stderr())
  quit(status = 1)
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Read data ----
coverage.data <- read.table(data_path, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, check.names = FALSE)

# ---- Manhattan plot ----
gen_manhattan <- function() {
  cat("Generating Manhattan plot...\n", file = stderr())

  # Read hiconf data if provided
  if (!is.null(hiconf_path) && file.exists(hiconf_path)) {
    hiconf_data <- read.table(hiconf_path, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE)
    coverage.data$hiconf <- hiconf_data$hiconf
  } else {
    coverage.data$hiconf <- 0
  }

  coverage.data$hiconf_map.all <- ifelse(coverage.data$hiconf >= 1,
                                         coverage.data[, 4], NA)

  # Determine y-limits
  y_up_lim <- max(coverage.data[, 4], na.rm = TRUE) * 1.1
  y_dw_lim <- 0
  if (ncol(coverage.data) >= 5 && !all(is.na(coverage.data[, 5]))) {
    y_dw_lim <- max(coverage.data[, 5], na.rm = TRUE) * 1.1
  }

  chr_colors <- c("#555555", "#111111")

  # Calculate cumulative position
  coverage.data.tmp <- coverage.data %>%
    group_by(chr) %>%
    summarise(chr_len = max(start), .groups = "drop") %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(coverage.data, ., by = c("chr" = "chr")) %>%
    arrange(chr, start) %>%
    mutate(poscum = start + tot)

  axis.label <- coverage.data.tmp %>%
    group_by(chr) %>%
    summarize(center = (max(poscum) + min(poscum)) / 2, .groups = "drop")

  # Output file
  png_path <- file.path(outdir, paste0(prefix, "_manhattan.png"))

  png(filename = png_path, width = image_width, height = image_height,
      units = "px", pointsize = image_points, res = image_res, bg = "white")

  p <- ggplot(coverage.data.tmp, aes(x = poscum)) +
    geom_point(y = -coverage.data.tmp[, 5], alpha = 0.8, size = 0.5, color = "grey") +
    geom_point(aes(y = coverage.data.tmp[, 4], color = as.factor(chr)),
               alpha = 1.0, size = 0.5) +
    geom_point(y = coverage.data.tmp$hiconf_map.all, alpha = 0.8,
               size = 0.5, color = "red") +
    scale_color_manual(values = rep(chr_colors, length(unique(coverage.data$chr)))) +
    scale_x_continuous(label = axis.label$chr, breaks = axis.label$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(-y_dw_lim, y_up_lim)) +
    labs(x = paste0("Chromosome [", window_size, "bp windows]"),
         y = "Mean depth [RPM]") +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  print(p)
  graphics.off()
  cat(sprintf("Manhattan plot saved: %s\n", png_path), file = stderr())
}

# ---- Line multiplot (with control) ----
gen_line_multiplot <- function(cand_data, candidate_name) {
  cat(sprintf("Generating line plot for %s...\n", candidate_name), file = stderr())

  out_png <- file.path(outdir, paste0(candidate_name, "_multiplot.png"))

  ncols <- ncol(cand_data)

  # Determine which columns to use
  has_control <- ncols >= 10
  y_up_lim_a <- max(c(cand_data[, 5], cand_data[, 8]), na.rm = TRUE) * 1.1
  y_up_lim_c <- max(c(cand_data[, 6], cand_data[, 9]), na.rm = TRUE) * 1.1
  y_up_lim_e <- max(c(cand_data[, 7], cand_data[, 10]), na.rm = TRUE) * 1.1

  theme_set(
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  )

  if (has_control) {
    # 6-panel layout: TR.all, CO.all, TR.SR, CO.SR, TR.DR, CO.DR
    plot.TR.all <- ggplot(cand_data, aes(x = pos)) +
      geom_line(aes(y = cand_data[, 5]), size = 0.25) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, y_up_lim_a)) +
      labs(x = colnames(cand_data)[5], y = NULL)

    plot.TR.SR <- ggplot(cand_data, aes(x = pos)) +
      geom_line(aes(y = cand_data[, 6]), size = 0.25) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, y_up_lim_c)) +
      labs(x = colnames(cand_data)[6], y = NULL)

    plot.TR.DR <- ggplot(cand_data, aes(x = pos)) +
      geom_line(aes(y = cand_data[, 7]), size = 0.25) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, y_up_lim_e)) +
      labs(x = colnames(cand_data)[7], y = NULL)

    plot.CO.all <- ggplot(cand_data, aes(x = pos)) +
      geom_line(aes(y = cand_data[, 8]), size = 0.25) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, y_up_lim_a)) +
      labs(x = colnames(cand_data)[8], y = NULL)

    plot.CO.SR <- ggplot(cand_data, aes(x = pos)) +
      geom_line(aes(y = cand_data[, 9]), size = 0.25) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, y_up_lim_c)) +
      labs(x = colnames(cand_data)[9], y = NULL)

    plot.CO.DR <- ggplot(cand_data, aes(x = pos)) +
      geom_line(aes(y = cand_data[, 10]), size = 0.25) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, y_up_lim_e)) +
      labs(x = colnames(cand_data)[10], y = NULL)

    png(filename = out_png, width = image_width, height = image_height,
        units = "px", pointsize = image_points, res = image_res, bg = "white")

    grid.arrange(plot.TR.all, plot.CO.all, plot.TR.SR, plot.CO.SR,
                 plot.TR.DR, plot.CO.DR,
                 ncol = 2, left = "Coverage depth [RPM]")

  } else {
    # 3-panel layout: TR only
    y_up_lim_a <- max(cand_data[, 5], na.rm = TRUE) * 1.1
    y_up_lim_c <- max(cand_data[, 6], na.rm = TRUE) * 1.1
    y_up_lim_e <- max(cand_data[, 7], na.rm = TRUE) * 1.1

    plot.TR.all <- ggplot(cand_data, aes(x = pos)) +
      geom_line(aes(y = cand_data[, 5]), size = 0.25) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, y_up_lim_a)) +
      labs(x = colnames(cand_data)[5], y = NULL)

    plot.TR.SR <- ggplot(cand_data, aes(x = pos)) +
      geom_line(aes(y = cand_data[, 6]), size = 0.25) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, y_up_lim_c)) +
      labs(x = colnames(cand_data)[6], y = NULL)

    plot.TR.DR <- ggplot(cand_data, aes(x = pos)) +
      geom_line(aes(y = cand_data[, 7]), size = 0.25) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, y_up_lim_e)) +
      labs(x = colnames(cand_data)[7], y = NULL)

    png(filename = out_png, width = image_width, height = image_height,
        units = "px", pointsize = image_points, res = image_res, bg = "white")

    grid.arrange(plot.TR.all, plot.TR.SR, plot.TR.DR,
                 ncol = 1, left = "Coverage depth [RPM]")
  }

  graphics.off()
  cat(sprintf("Line plot saved: %s\n", out_png), file = stderr())
}

# ---- Main execution ----
if (mode == "manhattan" || mode == "both") {
  gen_manhattan()
}

if (mode == "line" || mode == "both") {
  # Check if data has pos column (per-candidate data) or chr/start/end (genome-wide)
  # For per-candidate data, each file should have a pos column
  cat("NOTE: For per-candidate line plots, use individual candidate coverage files.\n", file = stderr())
}

cat("Visualization complete.\n", file = stderr())
