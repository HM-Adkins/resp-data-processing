library(tidyverse)


local({
  source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
  script_path <- if (length(source_files) > 0L) {
    tail(source_files, 1L)[[1]]
  } else {
file_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_args)) sub("^--file=", "", file_args[[1]]) else
      file.path(getwd(), "Raw_Resp_Data_Processing.R")
  }
  script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
  target <- parent.env(environment())
  source(file.path(script_dir, "Raw_Part1_Load_And_Prepare_Data.R"), local = target)
  source(file.path(script_dir, "Raw_Part2_Regression.R"), local = target)
})


# CORE HELPERS
save_raw_checkpoint <- function(results, path, stage = "raw_part1") {
  if (!is.list(results)) stop("results must be a list.")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(schema_version = 1L, stage = stage, results = results),
    path
  )
  invisible(path)
}

load_raw_checkpoint <- function(path, expected_stage = "raw_part1") {
  if (!file.exists(path)) stop("Checkpoint file not found: ", path)
  checkpoint <- readRDS(path)
  if (!identical(checkpoint$schema_version, 1L) ||
      !identical(checkpoint$stage, expected_stage) ||
      !is.list(checkpoint$results)) {
    stop("Checkpoint has an unsupported schema or stage.")
  }
  checkpoint$results
}


process_raw_data <- function(
    file_path,
    experiment_start_time,
    experiment_end_time,
    tz,
    summary_file = NULL,
    experiment = NULL,
    export_path = NULL,
    plot_title_prefix = NULL,
    generate_raw_data_plot = FALSE,
    chamber_numbers = NULL,
    checkpoint_file = NULL,
    verbose = TRUE) {
  if (!is.null(checkpoint_file) && is.null(summary_file)) {
    stop("checkpoint_file requires summary_file to build Part 2 input.")
  }
  prepared <- process_raw_part1(
    file_path, experiment_start_time, experiment_end_time, tz,
    chamber_numbers = chamber_numbers, export_path = export_path,
    plot_title_prefix = plot_title_prefix,
    generate_raw_data_plot = generate_raw_data_plot, verbose = verbose
  )
  regression <- process_raw_part2(
    prepared, summary_file = summary_file, experiment = experiment,
    export_path = export_path, verbose = verbose  
  )
  if (!is.null(checkpoint_file)) {
    save_raw_checkpoint(regression$part1_results, checkpoint_file)
  }
list(
    raw_data = prepared$raw_data,
    time_series = prepared$time_series,
    measure_stats = regression$measure_stats,
    col_info = prepared$col_info,
    part1_results = regression$part1_results
  )
}