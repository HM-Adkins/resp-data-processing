raw_script_candidates <- c(
  "Raw_Resp_Data_Processing.R",
  file.path("..", "Raw_Resp_Data_Processing.R"),
  file.path("..", "..", "Raw_Resp_Data_Processing.R")
)
raw_script <- raw_script_candidates[file.exists(raw_script_candidates)][1]
if (is.na(raw_script)) stop("Could not locate Raw_Resp_Data_Processing.R.")
source(raw_script, local = TRUE)

testthat::test_that("raw chamber selection keeps only requested channels", {
  data <- data.frame(
    DateTime = "7/26/2023 9:00:00 AM",
    Phase = "M1",
    `CH1 O2 mgO2` = 1,
    `CH1 temp` = 20,
    `CH2 O2 mgO2` = 2,
    `CH2 temp` = 21,
    check.names = FALSE
  )

  selected <- select_raw_chambers(data, get_raw_column_info(data), 2)

  testthat::expect_identical(selected$chamber_numbers, 2L)
  testthat::expect_identical(
    names(selected$data),
    c("DateTime", "Phase", "CH2 O2 mgO2", "CH2 temp")
  )
})

testthat::test_that("raw chamber selection rejects missing chamber numbers", {
  data <- data.frame(
    DateTime = "7/26/2023 9:00:00 AM",
    Phase = "M1",
    `CH1 O2 mgO2` = 1,
    check.names = FALSE
  )

  testthat::expect_error(
    select_raw_chambers(data, get_raw_column_info(data), 2),
    "not found"
  )
})

testthat::test_that("raw Part 1 checkpoints round-trip", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  results <- list(chamber_list = list(), descriptive_data = data.frame())

  save_raw_checkpoint(results, path)

  testthat::expect_identical(load_raw_checkpoint(path), results)
  testthat::expect_error(
    load_raw_checkpoint(path, expected_stage = "raw_part2b"),
    "unsupported schema or stage"
  )
})

testthat::test_that("raw Part 1 output has the Part 2B chamber-list contract", {
  measure <- tibble::tibble(
    Phase = c("M1", "M1", "M2", "M2"),
    time_sec = c(0, 36, 3600, 3636),
    time_hour = c(0, 0.01, 1, 1.01),
    `CH1 O2 mgO2` = c(10, 9, 8, 7),
    `CH1 temp` = c(20, 20, 21, 21),
    `CH2 O2 mgO2` = c(5, 6, 7, 8),
    `CH2 temp` = c(20, 20, 21, 21)
  )
  stats <- tibble::tibble(
    chamber = c(1L, 1L, 2L, 2L),
    phase = c("M", "M", "M", "M"),
    loop = c(1L, 2L, 1L, 2L),
    r2 = rep(0.95, 4),
    slope_per_hr = c(-1, -1, 2, 2),
    mo2_mg_kg_hr = c(10, 10, NA_real_, NA_real_)
  )
  summary <- data.frame(
    experiment = 18L,
    chamber = c(1L, 2L),
    species = c("Fish", "Control"),
    id = c("A", "Control"),
    www_kg = c(0.1, NA_real_),
    resp_vol_ml = c(100, NA_real_),
    ch_vol_ml = c(500, 1000)
  )
  summary_path <- tempfile(fileext = ".csv")
  on.exit(unlink(summary_path), add = TRUE)
  utils::write.csv(summary, summary_path, row.names = FALSE)

  result <- build_raw_part1_results(
    raw_data = measure,
    time_series = list(measure = measure),
    measure_stats = stats,
    col_info = get_raw_column_info(measure),
    summary_file = summary_path,
    experiment = 18,
    chamber_numbers = c(1L, 2L),
    file_path = "raw.csv",
    experiment_start_time = "2023-07-26 09:00:00",
    experiment_end_time = "2023-07-26 16:30:00"
  )

  testthat::expect_equal(length(result$chamber_list), 2L)
  testthat::expect_true(all(c(
    "time", "r2", "slope", "raw_mo2_ch", "raw_mo2_kg_www"
  ) %in% names(result$chamber_list[[1]])))
  testthat::expect_equal(result$chamber_list[[1]]$raw_mo2_ch, c(1, 1))
  testthat::expect_equal(result$chamber_list[[2]]$raw_mo2_ch, c(2, 2))
  testthat::expect_true(all(is.na(result$chamber_list[[2]]$raw_mo2_kg_www)))
})
