#------------------------------------------------------------------------------
# Respiration Data Workup
# Purpose:
# Preprocess MO2 (metabolic rate) data from intermittent flow respirometry trials
# Author:
# HM-Adkins
# Updates:
# 2025-12-19, 2026-11-02 (removed tww conversion),
# 2026-11-10 migrted to github for version control

# AI Assistance:
#   Used to enhance code functionality, reproducibility, and testing.
#   All AI-generated code was reviewed and validated by HM-Adkins.

# Citations:
# Anthropic. (2025). Claude 3.5 Sonnet [Large language model].
#     https://www.anthropic.com/claude
# GitHub. (2021). GitHub Copilot [AI-powered code completion tool].
#     https://github.com/features/copilot
#------------------------------------------------------------------------------

library(tidyverse)
library(lattice)
library(writexl)

#------------------------------------------------------------------------------
# PLOT CONFIGURATION
#------------------------------------------------------------------------------

# Configuration object for plot themes and settings
# Users can modify this before running analysis or create custom configs
RESPIRATION_CONFIG <- list(
  colors = list(
    raw = "#D55E00",
    corrected = "#0072B2",
    control = "#009E73",
    gray = "gray",
    red = "red",
    blue = "blue"
  ),
  plot = list(
    point_size = 0.75,
    point_cex = 0.7,
    line_width = 2,
    control_line_width = 2
  ),
  thresholds = list(
    default_r2 = 0.9,
    default_control_r2 = 0.85,
    min_good_control_points = 3,
    control_regression_r2 = 0.6,
    high_background_pct = 50
  )
)

#------------------------------------------------------------------------------
# PLOT CONFIGURATION HELPER
#------------------------------------------------------------------------------

# Helper function to create custom configuration

create_custom_config <- function(colors = NULL, plot = NULL, thresholds = NULL) {
  cfg <- RESPIRATION_CONFIG # Start with defaults

  # Override colors if provided
  if (!is.null(colors)) {
    for (name in names(colors)) {
      cfg$colors[[name]] <- colors[[name]]
    }
  }

  # Override plot settings if provided
  if (!is.null(plot)) {
    for (name in names(plot)) {
      cfg$plot[[name]] <- plot[[name]]
    }
  }

  # Override thresholds if provided
  if (!is.null(thresholds)) {
    for (name in names(thresholds)) {
      cfg$thresholds[[name]] <- thresholds[[name]]
    }
  }

  cfg
}

#------------------------------------------------------------------------------
# CORE HELPER FUNCTIONS: get_column_info, validate_columns, plot_device, setup_analysis
#------------------------------------------------------------------------------

# 0.1 Get column indices and names
get_column_info <- function(data) {
  patterns <- list(
    # Descriptive columns
    chamber = "^chamber$|^parameter.*chamber",
    www_g = "wet.*weight(?!.*tissue)|www",
    resp_vol_ml = "resp.*vol",
    chamber_notes = "chamber.*note|note",
    tww_g = "wet.*tissue|tww",
    ch_vol_ml = "chamber.*vol(?!.*resp)",

    # Time columns
    datetime = "date.*time",
    time_sec = "time.*sec",

    # Chamber measurement columns
    mo2 = "CH.*MO2",
    r2 = "CH.*R.*2",
    slope = "CH.*slope",
    temp = "AMB.*temp.*avg.*temp"
  )

  # Define which columns use only the first match (a single column per field,
  # not one column per chamber). Includes descriptive columns (one value per
  # individual) and time columns (shared across all individuals).
  scalar_cols <- c(
    "chamber", "www_g", "resp_vol_ml", "chamber_notes", "tww_g", "ch_vol_ml",
    "datetime", "time_sec"
  )

  result <- list(idx = list(), names = list())

  for (nm in names(patterns)) {
    matches <- grep(patterns[[nm]], names(data), ignore.case = TRUE, perl = TRUE)

    if (length(matches) > 0) {
      # For scalar columns, take only the first match
      if (nm %in% scalar_cols) {
        result$idx[[nm]] <- matches[1]
        result$names[[nm]] <- names(data)[matches[1]]
      } else {
        # For vector columns (mo2, r2, slope, temp), keep all matches
        result$idx[[nm]] <- matches
        result$names[[nm]] <- names(data)[matches]
      }
    } else {
      result$idx[[nm]] <- NA
      result$names[[nm]] <- NA_character_
    }
  }

  result
}

# 0.2 Validate required columns exist
validate_columns <- function(col_info, required) {
  missing <- character()

  for (col in required) {
    idx <- col_info$idx[[col]]
    if (is.na(idx[1])) {
      missing <- c(missing, col)
    }
  }

  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
}

# 0.3 Setup plot device helper - plots to screen and PDF
plot_device <- function(output_file, width = 10, height = 8, expr) {
  # Plot to screen first
  eval.parent(substitute(expr))

  # Then save to PDF if output_file is provided
  if (!is.null(output_file)) {
    # Create directory if it doesn't exist
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    pdf(output_file, width = width, height = height)
    eval.parent(substitute(expr))
    dev.off()
  }

  invisible(NULL)
}

# 0.4 Setup function to validate inputs and create output directory
setup_analysis <- function(file_path, export_path = NULL, verbose = TRUE) {
  # Check file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  # Create output directory
  if (!is.null(export_path)) {
    output_dir <- dirname(export_path)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      if (verbose) message("Created output directory: ", output_dir)
    }
  }

  invisible(TRUE)
}

#------------------------------------------------------------------------------
# PART 1: DATA LOADING AND PREPROCESSING
# extract_descriptive_data, prepare_time_series, build_chamber_list,
# plot_raw_data, process_part1_data_loading
#------------------------------------------------------------------------------

source("Part1_Load_And_Prepare_Data.R")

#------------------------------------------------------------------------------
# PART 2A: CONTROL CHAMBER PROCESSING
# analyze_control_chamber, plot_control_regression, process_control_chambers
#------------------------------------------------------------------------------

source("Part2A_Control_Chamber_Analysis.R")

#------------------------------------------------------------------------------
# PART 2B: BACKGROUND CORRECTION
# correct_background, plot_background_correction, process_background_correction
#------------------------------------------------------------------------------

source("Part2B_Background_Correction.R")

#------------------------------------------------------------------------------
# PART 3: QUALITY CONTROL
# apply_quality_filters, plot_rmr_by_temperature, plot_quality_diagnostics,
# process_part3_quality_control
#------------------------------------------------------------------------------

source("Part3_Quality_Control_Filtering.R")

#------------------------------------------------------------------------------
# EXPORT HELPER
#------------------------------------------------------------------------------

export_results <- function(results, export_path, verbose = TRUE) {
  if (verbose) cat("\n=== EXPORTING ===\n")

  # Create directory if it doesn't exist
  output_dir <- dirname(export_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Metadata
  metadata_df <- tibble(
    parameter = names(results$metadata),
    value = as.character(results$metadata)
  )

  # Config (flatten nested list for Excel)
  config_df <- tibble(
    category = c(
      rep("colors", length(results$config$colors)),
      rep("plot", length(results$config$plot)),
      rep("thresholds", length(results$config$thresholds))
    ),
    parameter = c(
      names(results$config$colors),
      names(results$config$plot),
      names(results$config$thresholds)
    ),
    value = as.character(c(
      unlist(results$config$colors),
      unlist(results$config$plot),
      unlist(results$config$thresholds)
    ))
  )

  # Control data
  ctrl_combined <- bind_rows(lapply(results$control_analysis$individual_controls, function(x) {
    select(x$data, experiment, chamber, id, time, slope, slope_filtered, slope_corrected, corrected_mo2_ch)
  }))

  # Excel workbook
  sheets <- list(
    descriptive_data = results$descriptive_data,
    control_analysis = results$control_analysis$regression,
    control_chambers = ctrl_combined,
    control_mean = results$control_analysis$data,
    corrected_expt_data = results$corrected_expt_data,
    sample_sizes = results$sample_sizes,
    exclusion_summary = results$exclusion_summary,
    metadata = metadata_df,
    config = config_df
  )

  data_summary <- results$descriptive_data %>%
    left_join(
      results$sample_sizes %>%
        select(experiment, chamber, species, id, n_raw:pct_retained),
      by = c("experiment", "chamber", "species", "id")
    )

  write_xlsx(sheets[!sapply(sheets, is.null)], paste0(export_path, "_workup.xlsx"))
  write_csv(results$final_data, paste0(export_path, "_final.csv"))
  write_csv(data_summary, paste0(export_path, "_data_summary.csv"))

  if (verbose) cat("  ✓ Files exported\n")
}

#------------------------------------------------------------------------------
# MAIN WORKFLOW
#------------------------------------------------------------------------------

process_respiration_data <- function(
  file_path,
  experiment_start_time,
  experiment,
  tz,
  date_format = "%m/%d/%Y/%I:%M:%S %p",
  id_delimiter = "-",
  r2_threshold = RESPIRATION_CONFIG$thresholds$default_r2,
  control_r2_threshold = RESPIRATION_CONFIG$thresholds$default_control_r2,
  export_path = NULL,
  plot_title_prefix = NULL,
  generate_raw_data_plot = FALSE,
  generate_control_plot = FALSE,
  generate_correction_plot = FALSE,
  generate_rmr_plot = FALSE,
  generate_diagnostics_plot = FALSE,
  rmr_plot_title = NULL,
  verbose = TRUE
) {
  if (verbose) cat("\n", strrep("=", 70), "\n")
  if (verbose) cat("RESPIRATION DATA WORKUP\n")
  if (verbose) cat("File:", file_path, "\n")
  if (verbose) cat(strrep("=", 70), "\n")

  # Setup and validation
  setup_analysis(file_path, export_path, verbose)

  # Part 1
  part1 <- process_part1_data_loading(
    file_path = file_path,
    experiment_start_time = experiment_start_time,
    experiment = experiment,
    tz = tz,
    id_delimiter = id_delimiter,
    plot_title_prefix = plot_title_prefix,
    generate_raw_data_plot = generate_raw_data_plot,
    export_path = export_path,
    verbose = verbose,
    date_format = date_format
  )

  # Part 2a: Process control chambers
  control_analysis <- process_control_chambers(
    part1_results = part1,
    control_r2_threshold = control_r2_threshold,
    generate_control_plot = generate_control_plot,
    plot_title_prefix = plot_title_prefix,
    export_path = export_path,
    verbose = verbose
  )

  # Part 2b: Correct background using control data
  part2 <- process_background_correction(
    part1_results = part1,
    control_analysis = control_analysis,
    generate_correction_plot = generate_correction_plot,
    plot_title_prefix = plot_title_prefix,
    export_path = export_path,
    verbose = verbose
  )

  # Part 3
  part3 <- process_part3_quality_control(
    part1, part2, r2_threshold, generate_rmr_plot, generate_diagnostics_plot,
    rmr_plot_title, plot_title_prefix, export_path, verbose
  )

  # Compile results
  results <- list(
    descriptive_data = part1$descriptive_data,
    raw_data = part1$raw_data,
    all_respiration_data = part1$time_series$all,
    acclimation_data = part1$time_series$acclimation,
    experimental_data = part1$time_series$experimental,
    control_analysis = part2$control_analysis,
    corrected_expt_data = part3$raw_corrected,
    final_data = part3$data,
    sample_sizes = part3$sample_sizes,
    exclusion_summary = part3$exclusion_summary,
    config = RESPIRATION_CONFIG,
    metadata = c(part1$metadata, list(
      r2_threshold = r2_threshold,
      control_r2_threshold = control_r2_threshold,
      id_delimiter = id_delimiter,
      analysis_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      analyst = Sys.info()["user"]
    ))
  )

  # Export
  if (!is.null(export_path)) {
    export_results(results, export_path, verbose)
  }

  if (verbose) cat("\n", strrep("=", 70), "\n")
  if (verbose) cat("✓✓✓ COMPLETE ✓✓✓\n")
  if (verbose) cat(strrep("=", 70), "\n\n")

  results
}