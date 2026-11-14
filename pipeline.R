#!/usr/bin/env Rscript

box::use(
  ./locate[match_chamber_columns],
  ./extract[
    extract_descriptive_data,
    format_time_series,
    extract_experimental_data,
    build_chamber_list
  ]
)

data_type <- "concentration"
data <- read.csv("data/T.macrodon_FMR_18C_raw.csv", stringsAsFactors = FALSE)

located <- match_chamber_columns(data, data_type, chamber_numbers = 1:8)
descriptive_data <- extract_descriptive_data(data, located, experiment = 1)
formatted <- format_time_series(
  data, located,
  tz = "America/Chicago",
  date_format = "%m/%d/%Y/%I:%M:%S %p"
)
experimental <- extract_experimental_data(
  formatted, located, data_type,
  start_time = "2023-07-26 09:00:00",
  end_time = "2023-07-26 16:30:00",
  tz = "America/Chicago"
)

chambers <- build_chamber_list(experimental, descriptive_data, located, data_type)
