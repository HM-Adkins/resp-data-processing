#!/usr/bin/env Rscript

box::use(./locate[match_chamber_columns])

data <- read.csv("data/T.macrodon_FMR_18C.csv", stringsAsFactors = FALSE)

located <- match_chamber_columns(
  data,
  data_type = "rate", # Use "rate" for rate data
  chamber_numbers = 1:8
)

col_info <- located$col_info
chamber_numbers <- located$chamber_numbers
chamber_columns <- located$columns