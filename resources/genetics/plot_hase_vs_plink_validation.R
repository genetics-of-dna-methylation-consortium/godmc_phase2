#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "Rscript resources/genetics/plot_hase_vs_plink_validation.R <hase_csv_gz> <plink_glm_gz> <out_dir> [prefix]",
  sep = "\n"
)

if (length(args) < 3 || length(args) > 4) {
  stop(usage, call. = FALSE)
}

hase_file <- args[[1]]
plink_file <- args[[2]]
out_dir <- args[[3]]
prefix <- if (length(args) >= 4) args[[4]] else "hase_vs_plink"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

require_columns <- function(dt, cols, label) {
  missing_cols <- setdiff(cols, names(dt))
  if (length(missing_cols) > 0) {
    fail(label, " is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
}

clean_allele <- function(x) {
  toupper(trimws(as.character(x)))
}

normalize_plink_chromosome <- function(plink) {
  plink[, plink_ID_original := ID]

  chrom_col <- intersect(c("#CHROM", "CHROM"), names(plink))
  if (length(chrom_col) > 0) {
    chrom_col <- chrom_col[[1]]
    x_rows <- toupper(trimws(as.character(plink[[chrom_col]]))) == "X"
    if (any(x_rows, na.rm = TRUE)) {
      set(plink, which(x_rows), chrom_col, "23")
      message("Converted PLINK ", chrom_col, " X values to 23: ", sum(x_rows, na.rm = TRUE))
    }
  }

  plink[, ID := sub("^X([:_])", "23\\1", as.character(ID), ignore.case = TRUE)]
  changed_ids <- plink[ID != plink_ID_original, .N]
  if (changed_ids > 0) {
    message("Converted PLINK IDs starting with X: or X_ to 23: or 23_: ", format_n(changed_ids))
  }

  plink
}

normalize_chr <- function(x) {
  value <- toupper(trimws(as.character(x)))
  value <- sub("^CHR", "", value)
  value[value == "X"] <- "23"
  suppressWarnings(as.integer(value))
}

format_n <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_cor <- function(x) {
  if (!is.finite(x)) {
    return("NA")
  }
  sprintf("%.6f", x)
}

parse_max_plot_points <- function() {
  value <- Sys.getenv("MAX_PLOT_POINTS", "1000000")
  if (tolower(value) %in% c("all", "inf", "infinite")) {
    return(Inf)
  }
  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed) || parsed <= 0) {
    fail("MAX_PLOT_POINTS must be a positive number, all, inf, or infinite")
  }
  parsed
}

read_table <- function(path, label) {
  tryCatch(
    fread(path),
    error = function(err) {
      if (grepl("\\.gz$", path) && nzchar(Sys.which("gzip"))) {
        message(
          "Direct fread failed for ", label,
          "; retrying through gzip -dc. Original error: ", conditionMessage(err)
        )
        cmd_dt <- tryCatch(
          fread(cmd = paste("gzip -dc", shQuote(path))),
          error = function(cmd_err) {
            message("gzip -dc fread retry failed; using base R gzfile. Error: ", conditionMessage(cmd_err))
            NULL
          }
        )
        if (!is.null(cmd_dt)) {
          return(cmd_dt)
        }
      }

      if (grepl("\\.gz$", path)) {
        message(
          "Direct fread failed for ", label,
          "; retrying with base R gzfile. Original error: ", conditionMessage(err)
        )
        header <- readLines(gzfile(path, "rt"), n = 1)
        sep <- if (grepl(",", header, fixed = TRUE)) "," else ""
        return(as.data.table(read.table(
          gzfile(path, "rt"),
          header = TRUE,
          sep = sep,
          quote = "",
          comment.char = "",
          check.names = FALSE
        )))
      }

      stop(err)
    }
  )
}

write_table_gz <- function(dt, path, sep = "\t") {
  direct_ok <- tryCatch(
    {
      fwrite(dt, path, sep = sep)
      TRUE
    },
    error = function(err) {
      message("Direct fwrite to gzip failed; using fallback. Original error: ", conditionMessage(err))
      if (file.exists(path)) {
        unlink(path)
      }
      FALSE
    }
  )
  if (direct_ok) {
    return(invisible(TRUE))
  }

  if (nzchar(Sys.which("gzip"))) {
    tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
    fwrite(dt, tmp, sep = sep)
    status <- system2("gzip", c("-f", tmp))
    gz_tmp <- paste0(tmp, ".gz")
    if (identical(status, 0L) && file.exists(gz_tmp)) {
      if (file.exists(path)) {
        unlink(path)
      }
      if (file.rename(gz_tmp, path)) {
        return(invisible(TRUE))
      }
    }
    if (file.exists(tmp)) {
      unlink(tmp)
    }
    if (file.exists(gz_tmp)) {
      unlink(gz_tmp)
    }
  }

  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(dt, file = con, sep = sep, quote = FALSE, row.names = FALSE, col.names = TRUE)
  invisible(TRUE)
}

subset_for_plot <- function(dt, max_points) {
  if (is.infinite(max_points) || nrow(dt) <= max_points) {
    return(dt)
  }
  set.seed(1)
  dt[sample.int(nrow(dt), max_points)]
}

build_scatter_plot <- function(dt, x_col, y_col, x_label, y_label, title, subtitle, max_points,
                               color_col = NULL) {
  finite_dt <- dt[is.finite(get(x_col)) & is.finite(get(y_col))]
  if (nrow(finite_dt) == 0) {
    fail("No finite points available for plot: ", title)
  }

  plot_dt <- subset_for_plot(finite_dt, max_points)
  axis_range <- range(c(plot_dt[[x_col]], plot_dt[[y_col]]), finite = TRUE)
  padding <- diff(axis_range) * 0.02
  if (!is.finite(padding) || padding == 0) {
    padding <- max(abs(axis_range), 1) * 0.02
  }
  axis_range <- axis_range + c(-padding, padding)

  point_alpha <- 0.4
  point_size <- 1.2

  abline_layer <- if (packageVersion("ggplot2") >= "3.4.0") {
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 0.65)
  } else {
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 0.65)
  }

  p <- ggplot(plot_dt, aes(x = .data[[x_col]], y = .data[[y_col]]))
  if (!is.null(color_col) && color_col %in% names(plot_dt)) {
    p <- p +
      geom_point(aes(color = .data[[color_col]]), alpha = point_alpha, size = point_size, stroke = 0) +
      scale_color_manual(
        values = c(
          "OBS_CT=max" = "steelblue",
          "OBS_CT<max" = "darkorange",
          "OBS_CT missing" = "#7f7f7f"
        ),
        drop = FALSE,
        name = "PLINK sample size"
      )
  } else {
    p <- p + geom_point(alpha = point_alpha, size = point_size, stroke = 0, color = "steelblue")
  }

  p <- p +
    abline_layer +
    coord_equal(xlim = axis_range, ylim = axis_range) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      legend.position = "bottom"
    )
}

save_comparison_panel <- function(merged_group, label, output_file, max_points) {
  matched_group <- merged_group[
    allele_status %in% c("matched", "flipped") &
      is.finite(plink_beta_aligned) &
      is.finite(hase_beta) &
      is.finite(SE) &
      is.finite(hase_se)
  ]

  if (nrow(matched_group) == 0) {
    message("Warning: no finite allele-alignable variants available for ", label, "; skipping plot")
    return(invisible(FALSE))
  }

  n_merged_group <- nrow(merged_group)
  n_matched_group <- merged_group[allele_status %in% c("matched", "flipped"), .N]
  n_flipped_group <- merged_group[allele_status == "flipped", .N]
  n_mismatch_group <- merged_group[allele_status == "mismatch", .N]
  plot_n <- min(nrow(matched_group), max_points)
  if (is.infinite(max_points)) {
    plot_n <- nrow(matched_group)
  }

  beta_cor_group <- suppressWarnings(cor(
    matched_group$plink_beta_aligned,
    matched_group$hase_beta,
    use = "complete.obs"
  ))
  se_cor_group <- suppressWarnings(cor(
    matched_group$SE,
    matched_group$hase_se,
    use = "complete.obs"
  ))

  subtitle <- paste0(
    label,
    "; merged=", format_n(n_merged_group),
    "; allele alignable=", format_n(n_matched_group),
    "; allele flipped=", format_n(n_flipped_group),
    "; mismatch=", format_n(n_mismatch_group),
    "; plotted=", format_n(plot_n)
  )

  beta_plot <- build_scatter_plot(
    matched_group,
    "plink_beta_aligned",
    "hase_beta",
    "PLINK beta aligned to HASE beta allele",
    "HASE beta",
    paste0("Beta; cor=", format_cor(beta_cor_group)),
    subtitle,
    max_points
  )

  se_plot <- build_scatter_plot(
    matched_group,
    "SE",
    "hase_se",
    "PLINK SE",
    "HASE standard_error",
    paste0("SE; cor=", format_cor(se_cor_group)),
    subtitle,
    max_points
  )

  png(output_file, width = 16, height = 7, units = "in", res = 300)
  device_open <- TRUE
  on.exit({
    if (device_open) {
      dev.off()
    }
  }, add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(beta_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(se_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  dev.off()
  device_open <- FALSE

  message(label, " beta correlation: ", format_cor(beta_cor_group))
  message(label, " SE correlation: ", format_cor(se_cor_group))
  invisible(TRUE)
}

message("Reading HASE results: ", hase_file)
hase <- read_table(hase_file, "HASE file")
message("HASE input rows: ", format_n(nrow(hase)))
message("Reading PLINK results: ", plink_file)
plink <- read_table(plink_file, "PLINK file")
message("PLINK input rows: ", format_n(nrow(plink)))

require_columns(
  hase,
  c("ID", "effect_allele", "non_effect_allele", "beta", "standard_error"),
  "HASE file"
)
require_columns(
  plink,
  c("ID", "A1", "TEST", "BETA", "SE"),
  "PLINK file"
)
message("Using HASE non_effect_allele as the beta allele for PLINK alignment")

plink <- plink[TEST == "ADD"]
message("PLINK rows after TEST == 'ADD': ", format_n(nrow(plink)))
if (nrow(plink) == 0) {
  fail("No PLINK rows remain after filtering TEST == 'ADD'")
}
plink <- normalize_plink_chromosome(plink)

hase_dups <- hase[, .N, by = ID][N > 1, .N]
plink_dups <- plink[, .N, by = ID][N > 1, .N]
if (hase_dups > 0) {
  message("Warning: HASE has duplicated IDs: ", format_n(hase_dups))
}
if (plink_dups > 0) {
  message("Warning: PLINK has duplicated IDs after TEST == 'ADD': ", format_n(plink_dups))
}

hase[, hase_beta := beta]
hase[, hase_se := standard_error]

merged <- merge(
  hase,
  plink,
  by = "ID",
  all = FALSE,
  suffixes = c("_hase", "_plink")
)

if (nrow(merged) == 0) {
  fail("No variants overlapped between HASE and PLINK by ID")
}

merged[, A1_clean := clean_allele(A1)]
merged[, hase_beta_allele := non_effect_allele]
merged[, hase_other_allele := effect_allele]
merged[, hase_beta_allele_clean := clean_allele(hase_beta_allele)]
merged[, hase_other_allele_clean := clean_allele(hase_other_allele)]

merged[, allele_status := fifelse(
  A1_clean == hase_beta_allele_clean,
  "matched",
  fifelse(A1_clean == hase_other_allele_clean, "flipped", "mismatch")
)]

merged[, plink_beta_aligned := fifelse(
  allele_status == "matched",
  BETA,
  fifelse(allele_status == "flipped", -BETA, NA_real_)
)]

chrom_col <- intersect(c("CHR", "CHR_hase", "chr", "chromosome", "#CHROM", "CHROM"), names(merged))
if (length(chrom_col) == 0) {
  fail("Merged data does not contain a chromosome column for chr1-22/chr23 plotting")
}
chrom_col <- chrom_col[[1]]
merged[, comparison_chr := normalize_chr(merged[[chrom_col]])]
missing_chr <- merged[is.na(comparison_chr), .N]
if (missing_chr > 0) {
  message("Warning: rows with unparseable chromosome values in ", chrom_col, ": ", format_n(missing_chr))
}

if ("OBS_CT" %in% names(merged)) {
  merged[, plink_obsct := suppressWarnings(as.numeric(OBS_CT))]
  finite_obsct <- merged[is.finite(plink_obsct), plink_obsct]
  if (length(finite_obsct) > 0) {
    max_obsct <- max(finite_obsct)
    min_obsct <- min(finite_obsct)
    merged[, plink_obsct_group := fifelse(
      !is.finite(plink_obsct),
      "OBS_CT missing",
      fifelse(plink_obsct == max_obsct, "OBS_CT=max", "OBS_CT<max")
    )]
    message("Merged PLINK OBS_CT min: ", format_n(min_obsct))
    message("Merged PLINK OBS_CT max: ", format_n(max_obsct))
    obsct_counts <- merged[, .N, by = plink_obsct_group][order(plink_obsct_group)]
    for (i in seq_len(nrow(obsct_counts))) {
      message("Merged ", obsct_counts$plink_obsct_group[i], " rows: ", format_n(obsct_counts$N[i]))
    }
  } else {
    merged[, plink_obsct_group := "OBS_CT missing"]
    fail("PLINK OBS_CT column exists but contains no finite values after merge")
  }
} else {
  fail("PLINK OBS_CT column not found; four-panel sample-size split requires OBS_CT")
}

matched_dt <- merged[
  allele_status %in% c("matched", "flipped") &
    is.finite(plink_beta_aligned) &
    is.finite(hase_beta) &
    is.finite(SE) &
    is.finite(hase_se)
]

if (nrow(matched_dt) == 0) {
  fail("No finite allele-alignable variants available after alignment")
}

if ("plink_obsct_group" %in% names(matched_dt)) {
  obsct_diagnostic <- matched_dt[, .(
    N = .N,
    beta_cor = suppressWarnings(cor(plink_beta_aligned, hase_beta, use = "complete.obs")),
    se_cor = suppressWarnings(cor(SE, hase_se, use = "complete.obs")),
    mean_abs_beta_diff = mean(abs(plink_beta_aligned - hase_beta), na.rm = TRUE)
  ), by = plink_obsct_group][order(plink_obsct_group)]
  for (i in seq_len(nrow(obsct_diagnostic))) {
    message(
      "Allele-alignable ", obsct_diagnostic$plink_obsct_group[i],
      ": N=", format_n(obsct_diagnostic$N[i]),
      "; beta cor=", format_cor(obsct_diagnostic$beta_cor[i]),
      "; SE cor=", format_cor(obsct_diagnostic$se_cor[i]),
      "; mean abs beta diff=", sprintf("%.8g", obsct_diagnostic$mean_abs_beta_diff[i])
    )
  }
}

n_merged <- nrow(merged)
n_matched <- merged[allele_status %in% c("matched", "flipped"), .N]
n_mismatch <- merged[allele_status == "mismatch", .N]
n_flipped <- merged[allele_status == "flipped", .N]

beta_cor <- suppressWarnings(cor(matched_dt$plink_beta_aligned, matched_dt$hase_beta, use = "complete.obs"))
se_cor <- suppressWarnings(cor(matched_dt$SE, matched_dt$hase_se, use = "complete.obs"))

merged_out <- file.path(out_dir, paste0(prefix, ".merged.tsv.gz"))
chr1_22_obsct_max_png <- file.path(out_dir, paste0(prefix, ".chr1_22.obsct_max.beta_se_scatter.png"))
chr1_22_obsct_lt_max_png <- file.path(out_dir, paste0(prefix, ".chr1_22.obsct_lt_max.beta_se_scatter.png"))
chr23_obsct_max_png <- file.path(out_dir, paste0(prefix, ".chr23.obsct_max.beta_se_scatter.png"))
chr23_obsct_lt_max_png <- file.path(out_dir, paste0(prefix, ".chr23.obsct_lt_max.beta_se_scatter.png"))

message("Writing merged comparison table: ", merged_out)
write_table_gz(merged, merged_out, sep = "\t")

max_plot_points <- parse_max_plot_points()

if (!"plink_obsct_group" %in% names(merged)) {
  fail("PLINK OBS_CT grouping is required for the requested four-panel split")
}

message("Writing chr1-22 OBS_CT=max beta/SE scatter panel: ", chr1_22_obsct_max_png)
save_comparison_panel(
  merged[comparison_chr %in% 1:22 & plink_obsct_group == "OBS_CT=max"],
  "chr1-22; OBS_CT=max",
  chr1_22_obsct_max_png,
  max_plot_points
)

message("Writing chr1-22 OBS_CT<max beta/SE scatter panel: ", chr1_22_obsct_lt_max_png)
save_comparison_panel(
  merged[comparison_chr %in% 1:22 & plink_obsct_group == "OBS_CT<max"],
  "chr1-22; OBS_CT<max",
  chr1_22_obsct_lt_max_png,
  max_plot_points
)

message("Writing chr23 OBS_CT=max beta/SE scatter panel: ", chr23_obsct_max_png)
save_comparison_panel(
  merged[comparison_chr == 23 & plink_obsct_group == "OBS_CT=max"],
  "chr23; OBS_CT=max",
  chr23_obsct_max_png,
  max_plot_points
)

message("Writing chr23 OBS_CT<max beta/SE scatter panel: ", chr23_obsct_lt_max_png)
save_comparison_panel(
  merged[comparison_chr == 23 & plink_obsct_group == "OBS_CT<max"],
  "chr23; OBS_CT<max",
  chr23_obsct_lt_max_png,
  max_plot_points
)

message("Done")
message("Merged rows: ", format_n(n_merged))
message("Allele alignable rows: ", format_n(n_matched))
message("PLINK beta flipped rows: ", format_n(n_flipped))
message("Allele mismatch rows: ", format_n(n_mismatch))
message("Beta correlation: ", format_cor(beta_cor))
message("SE correlation: ", format_cor(se_cor))
