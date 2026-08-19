#!/usr/bin/env Rscript
#
# normalize.R — Normalize coverage data to RPM (Reads Per Million mapped bases)
# and calculate fold enrichment vs genome background.
#
# Extracted from ECCsplorer/lib/eccDNA_Rcodes.py:
#   Rconvert_genome  — genome-wide windowed coverage → RPM
#   Rconvert_region  — region coverage → RPM + fold enrichment (Renrichment)
#
# Usage:
#   normalize.R --coverage <raw_coverage.csv> \
#               --mapped_bases <base1,base2,...> \
#               --output <normalized.csv>
#
# If --mode=genome (default):
#   Input:  [chr, start, end, TR.all, CO.all]
#   Output: [chr, start, end, TR.all_RPM, CO.all_RPM]
#
# If --mode=region:
#   Input:  [chr, start, end, TR.all, TR.SR, TR.DR, CO.all, CO.SR, CO.DR]
#   Output: Same + enrich.all column (fold enrichment vs genome background)
#
# Normalization: RPM = ceiling(coverage_value * 1,000,000 / mapped_bases)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

# ---- Argument parsing ----
args <- commandArgs(trailingOnly = TRUE)

coverage_path <- NULL
mapped_bases <- NULL
output_path <- NULL
mode <- "genome"
enrich_threshold <- 2.0
stats_paths <- character(0)

i <- 1
while (i <= length(args)) {
  if (args[i] == "--coverage") {
    coverage_path <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--mapped_bases") {
    mapped_bases <- as.integer(strsplit(args[i + 1], ",")[[1]])
    i <- i + 2
  } else if (args[i] == "--output") {
    output_path <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--mode") {
    mode <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--stats") {
    stats_paths <- c(stats_paths, args[i + 1])
    i <- i + 2
  } else if (args[i] == "--enrich_threshold") {
    enrich_threshold <- as.numeric(args[i + 1])
    i <- i + 2
  } else {
    cat(sprintf("WARNING: Unknown argument: %s\n", args[i]), file = stderr())
    i <- i + 1
  }
}

if (is.null(coverage_path) || is.null(mapped_bases) || is.null(output_path)) {
  cat("USAGE: normalize.R --coverage <FILE> --mapped_bases <N1,N2,...> --output <FILE> [--mode genome|region]\n",
      file = stderr())
  quit(status = 1)
}

if (!file.exists(coverage_path)) {
  cat(sprintf("ERROR: Coverage file '%s' not found.\n", coverage_path), file = stderr())
  quit(status = 1)
}

if (is.null(mapped_bases) && length(stats_paths) > 0) {
  mapped_bases <- sapply(stats_paths, function(sp) {
    lines <- readLines(sp, warn = FALSE)
    line <- grep("^bases mapped \\(cigar\\):", lines, value = TRUE)[1]
    if (is.na(line)) {
      NA_integer_
    } else {
      as.integer(trimws(sub(".*:", "", line)))
    }
  })
}
if (is.null(mapped_bases) || length(mapped_bases) == 0 || any(is.na(mapped_bases))) {
  mapped_bases <- 1000000
  cat("WARNING: mapped bases not found in stats files; defaulting to 1e6.\n", file = stderr())
}

# ---- Read data ----
coverage.data <- read.table(coverage_path, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, check.names = FALSE)

# ---- RPM normalization (genome mode) ----
# Columns: chr, start, end, TR.all (, CO.all)
if (mode == "genome") {
  ncols <- ncol(coverage.data)

  # TR factor: column 4
  coverage.data[, 4] <- ceiling(coverage.data[, 4] * (1000000 / mapped_bases[1]))

  # CO factor: column 5 (if exists)
  if (length(mapped_bases) >= 2 && ncols >= 5) {
    tryCatch({
      coverage.data[, 5] <- ceiling(coverage.data[, 5] * (1000000 / mapped_bases[2]))
    }, error = function(e) {
      cat(sprintf("NOTE: Skipping CO normalization: %s\n", e$message), file = stderr())
    })
  }

  write.table(coverage.data, file = output_path, sep = "\t",
              col.names = TRUE, row.names = FALSE, quote = FALSE)

} else if (mode == "region") {
  # Region mode: chr, start, end, TR.all, TR.SR, TR.DR (, CO.all, CO.SR, CO.DR)
  ncols <- ncol(coverage.data)

  # TR normalization (columns 4, 5, 6)
  coverage.data[, 4] <- ceiling(coverage.data[, 4] * (1000000 / mapped_bases[1]))
  coverage.data[, 5] <- ceiling(coverage.data[, 5] * (1000000 / mapped_bases[1]))
  coverage.data[, 6] <- ceiling(coverage.data[, 6] * (1000000 / mapped_bases[1]))

  # CO normalization (columns 7, 8, 9) if present
  if (length(mapped_bases) >= 2 && ncols >= 7) {
    tryCatch({
      coverage.data[, 7] <- ceiling(coverage.data[, 7] * (1000000 / mapped_bases[2]))
      if (ncols >= 8) coverage.data[, 8] <- ceiling(coverage.data[, 8] * (1000000 / mapped_bases[2]))
      if (ncols >= 9) coverage.data[, 9] <- ceiling(coverage.data[, 9] * (1000000 / mapped_bases[2]))
    }, error = function(e) {
      cat(sprintf("NOTE: Skipping CO normalization in region mode: %s\n", e$message), file = stderr())
    })
  }

  # ---- Fold enrichment calculation (replicates eccDNA_Rcodes.py enrichment) ----
  # enrich.all = TR.all_RPM / CO.all_RPM; when CO missing use mean(CO.all) background
  tr_all <- as.numeric(coverage.data[, 4])
  co_all <- if (ncols >= 7) as.numeric(coverage.data[, 7]) else rep(NA, nrow(coverage.data))

  # control background: mean of CO.all (fallback to TR.all mean when no control)
  avg_control <- if (all(is.na(co_all))) {
    mean(tr_all, na.rm = TRUE)
  } else {
    mean(co_all, na.rm = TRUE)
  }
  if (is.na(avg_control) || avg_control == 0) avg_control <- 1

  enrich_all <- tr_all / co_all
  enrich_alt <- round(tr_all / avg_control, digits = 2)
  # when CO.all is NA → fallback to control-background enrichment (replicates enrich.alt)
  enrich_all[is.na(enrich_all) | is.infinite(enrich_all)] <- enrich_alt[is.na(enrich_all) | is.infinite(enrich_all)]
  coverage.data$enrich.all <- round(enrich_all, digits = 2)

  # Sort by enrichment score
  coverage.data <- coverage.data[order(-coverage.data$enrich.all,
                                        coverage.data$chr,
                                        coverage.data$start), ]

  # Add candidate IDs
  coverage.data <- coverage.data %>%
    mutate(id = sprintf("eccCand_%03d", row_number()))

  # Sort back by chromosome and position
  coverage.data <- coverage.data[order(coverage.data$chr, coverage.data$start), ]

  write.table(coverage.data, file = output_path, sep = "\t",
              col.names = TRUE, row.names = FALSE, quote = FALSE)

} else {
  cat(sprintf("ERROR: Unknown mode '%s'. Use 'genome' or 'region'.\n", mode), file = stderr())
  quit(status = 1)
}

cat(sprintf("Normalization complete: %s\n", output_path), file = stderr())
