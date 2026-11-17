#!/usr/bin/env Rscript

box::use(
  ./locate[...],
  ./extract[...],
  ./regress[...]
)

data_type <- "concentration"
data <- read.csv("data/T.macrodon_FMR_18C_raw.csv", stringsAsFactors = FALSE)

located <- match_chamber_columns(
  data,
  data_type = data_type,
  chamber_numbers = c(1:8)
)

extracted <- extract_descriptive_data(data, located)

descriptive <- format_descriptive_data(
  extracted$descriptive, located, experiment = 18
)

formatted <- format_time_series(
  extracted$time_series, located,
  tz = "America/Chicago",
  date_format = "%m/%d/%Y/%I:%M:%S %p"
)

experimental <- extract_experimental_data(
  formatted,
  located,
  data_type = data_type,
  start_time = "2023-07-26 09:00:00",
  end_time = "2023-07-26 16:30:00",
  tz = "America/Chicago"
)

chambers <- build_chamber_list(
  experimental,
  descriptive,
  located,
  data_type = data_type
)

regressions <- process_regressions(
  chambers,
  descriptive,
  bin_size = NULL
)

regression_chambers <- regressions$chamber_list

comparison <- compare_regressions(
  chambers,
  descriptive,
  bin_sizes = c(3, 6)
)

comparison$comparison