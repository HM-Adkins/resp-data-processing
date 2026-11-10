# Edit the SETTINGS block, save, then run the whole file.

# ---- SETTINGS: edit these ----------------------------------------------
stages <- c("raw", "control", "correction")  # any subset, runs in this order

settings <- list(
  file_path = "data/T.macrodon_FMR_18C_raw.csv",
  summary_file = "data/saw_FMR_TM_Data_summary.csv",
  experiment = 18,
  experiment_start_time = "2023-07-26 09:00:00",
  experiment_end_time = "2023-07-26 16:30:00",
  tz = "America/Chicago",
  chamber_numbers = c(1, 2, 3, 4, 5, 6, 7, 8),  # must include a control chamber
  plot_title_prefix = "T. macrodon 18C",
  export_path = "output/FMR_18C_TM_raw",
  checkpoint_file = "output/raw_part1.rds"
)
# ---- END SETTINGS -------------------------------------------------------

source("Resp_Data_Processing_Library.R")
source("Raw_Resp_Data_Processing.R")

dir.create(dirname(settings$export_path), recursive = TRUE, showWarnings = FALSE)
control_file <- paste0(settings$export_path, "_control_analysis.rds")
correction_file <- paste0(settings$export_path, "_background_correction.rds")

if ("raw" %in% stages) {
  process_raw_data(
    file_path = settings$file_path,
    experiment_start_time = settings$experiment_start_time,
    experiment_end_time = settings$experiment_end_time,
    tz = settings$tz,
    summary_file = settings$summary_file,
    experiment = settings$experiment,
    chamber_numbers = settings$chamber_numbers,
    export_path = settings$export_path,
    plot_title_prefix = settings$plot_title_prefix,
    generate_raw_data_plot = TRUE,
    checkpoint_file = settings$checkpoint_file
  )
}

# Later stages reload the checkpoint, so they work without "raw"
part1 <- load_raw_checkpoint(settings$checkpoint_file)

if ("control" %in% stages) {
  control_analysis <- process_control_chambers(
    part1,
    generate_control_plot = TRUE,
    plot_title_prefix = settings$plot_title_prefix,
    export_path = settings$export_path
  )
  saveRDS(control_analysis, control_file)
}

if ("correction" %in% stages) {
  if (!exists("control_analysis")) control_analysis <- readRDS(control_file)
  part2b <- process_background_correction(
    part1,
    control_analysis,
    generate_correction_plot = TRUE,
    plot_title_prefix = settings$plot_title_prefix,
    export_path = settings$export_path
  )
  saveRDS(part2b, correction_file)
}