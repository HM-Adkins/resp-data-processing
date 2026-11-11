#------------------------------------------------------------------------------
# Respiration Data Workup
# Purpose: Preprocess MO2 (metabolic rate) data from intermittent flow respirometry trials
# Author: HM-Adkins
# Last Updated: 2025-12-30
#------------------------------------------------------------------------------

# Requires packages: tidyverse, lattice and writexl

#------------------------------------------------------------------------------
# Process respiration data
#------------------------------------------------------------------------------

source("Resp_Data_Processing_Library.R")

#------------------------------------------------------------------------------

FMR_18C_TM_results <- process_respiration_data(
  file_path = "T.macrodon_FMR_18C.csv",
  experiment_start_time = "2023-07-26 9:00:00",
  experiment = 18,
  tz = "America/Chicago",
  export_path = "output/FMR_18C_TM",
  plot_title_prefix = "T. macrodon 18°C",
  generate_raw_data_plot = TRUE,
  generate_control_plot = TRUE,
  generate_correction_plot = TRUE,
  generate_rmr_plot = TRUE,
  generate_diagnostics_plot = TRUE,
  r2_threshold = 0.9,
  control_r2_threshold = 0.85,
  verbose = TRUE
)

#------------------------------------------------------------------------------

FMR_22C_R1_TM_results <- process_respiration_data(
  file_path = "T.macrodon_FMR_22C_RUN1.csv",
  experiment_start_time = "2023-06-21 9:00:00",
  experiment = 22,
  tz = "America/Chicago",
  export_path = "output/FMR_22C_R1_TM",
  plot_title_prefix = "T. macrodon 22°C",
  generate_raw_data_plot = TRUE,
  generate_control_plot = TRUE,
  generate_correction_plot = TRUE,
  generate_rmr_plot = TRUE,
  generate_diagnostics_plot = TRUE,
  r2_threshold = 0.9,
  control_r2_threshold = 0.85,
  verbose = TRUE
)

#------------------------------------------------------------------------------

FMR_26C_R1_TM_results <- process_respiration_data(
  file_path = "T.macrodon_FMR_26C_RUN1.csv",
  experiment_start_time = "2023-06-28 9:00:00",
  experiment = 26,
  tz = "America/Chicago",
  export_path = "output/FMR_26C_R1_TM",
  plot_title_prefix = "T. macrodon 26°C",
  generate_raw_data_plot = TRUE,
  generate_control_plot = TRUE,
  generate_correction_plot = TRUE,
  generate_rmr_plot = TRUE,
  generate_diagnostics_plot = TRUE,
  r2_threshold = 0.9,
  control_r2_threshold = 0.85,
  verbose = TRUE
)

#------------------------------------------------------------------------------

FMR_30C_R1_TM_results <- process_respiration_data(
  file_path = "T.macrodon_FMR_30C_RUN1.csv",
  experiment_start_time = "2023-07-12 9:00:00",
  experiment = 30,
  tz = "America/Chicago",
  export_path = "output/FMR_30C_R1_TM",
  plot_title_prefix = "T. macrodon 30°C",
  generate_raw_data_plot = TRUE,
  generate_control_plot = TRUE,
  generate_correction_plot = TRUE,
  generate_rmr_plot = TRUE,
  generate_diagnostics_plot = TRUE,
  r2_threshold = 0.9,
  control_r2_threshold = 0.85,
  verbose = TRUE
)

#------------------------------------------------------------------------------

FMR_34C_TM_results <- process_respiration_data(
  file_path = "T.macrodon_FMR_34C.csv",
  experiment_start_time = "2023-07-19 9:00:00",
  experiment = 34,
  tz = "America/Chicago",
  export_path = "output/FMR_34C_TM",
  plot_title_prefix = "T. macrodon 34°C",
  generate_raw_data_plot = TRUE,
  generate_control_plot = TRUE,
  generate_correction_plot = TRUE,
  generate_rmr_plot = TRUE,
  generate_diagnostics_plot = TRUE,
  r2_threshold = 0.9,
  control_r2_threshold = 0.85,
  verbose = TRUE
)

#------------------------------------------------------------------------------

setwd("C:/Users/ham0018/Box/2023 Truncilla Brazos Project/FMR Data/WORKUP 2026/output")

TM18 <- read.csv("FMR_18C_TM_final.csv")
TM22 <- read.csv("FMR_22C_R1_TM_final.csv")
TM26 <- read.csv("FMR_26C_R1_TM_final.csv")
TM30 <- read.csv("FMR_30C_R1_TM_final.csv")
TM34 <- read.csv("FMR_34C_TM_final.csv")
All_Data <- rbind(TM18, TM22, TM26, TM30, TM34)
str(All_Data)
unique(All_Data$species)

TM18ds <- read.csv("FMR_18C_TM_data_summary.csv")
TM22ds <- read.csv("FMR_22C_R1_TM_data_summary.csv")
TM26ds <- read.csv("FMR_26C_R1_TM_data_summary.csv")
TM30ds <- read.csv("FMR_30C_R1_TM_data_summary.csv")
TM34ds <- read.csv("FMR_34C_TM_data_summary.csv")
All_Data_summary <- rbind(TM18ds, TM22ds, TM26ds, TM30ds, TM34ds)
str(All_Data_summary)
unique(All_Data_summary$species)

# Export combined data to Excel
library(writexl)
write_xlsx(All_Data, "FMR_TM_Data_final.xlsx")
write_xlsx(All_Data_summary, "FMR_TM_Data_summary.xlsx")