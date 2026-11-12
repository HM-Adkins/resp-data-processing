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

testthat::test_that("bins use complete groups of valid observations within each chamber and phase", {
  measure <- tibble::tibble(
    Phase = rep(c("M1", "M2"), each = 7),
    time_sec = 1:14,
    time_hour = (1:14) / 3600,
    `CH1 O2 mgO2` = c(1, NA, 3, 5, 7, 9, 100, 11:17),
    `CH2 O2 mgO2` = c(21:27, NA, 29:34)
  )
  info <- get_raw_column_info(measure)
  points <- prepare_regression_points(measure, info, 2)
  m1_ch1 <- dplyr::filter(points, phase == "M1", chamber == 1)
  testthat::expect_equal(m1_ch1$o2, c(2, 6, 54.5))
  testthat::expect_equal(m1_ch1$time_sec, c(2, 4.5, 6.5))
  testthat::expect_equal(nrow(points), 12)
  testthat::expect_equal(
    nrow(dplyr::filter(points, phase == "M2", chamber == 2)), 3
  )
  for (size in list(0, -1, 1.5, NA_real_, Inf, c(2, 3), "2")) {
    testthat::expect_error(prepare_regression_points(measure, info, size), "bin_size")
  }
})

testthat::test_that("regression comparisons preserve missing and insufficient periods", {
  measure <- tibble::tibble(
    Phase = rep(c("M1", "M2"), each = 7),
    time_sec = 1:14,
    time_hour = (1:14) / 3600,
    `CH1 O2 mgO2` = c(1, 2, 4, 5, 6, 7, 100, rep(NA_real_, 7)),
    `CH2 O2 mgO2` = c(2, 4, 5, 6, 8, 9, 10, 10, 12, 13, 14, 16, 17, 18)
  )
  prepared <- list(
    raw_data = measure, time_series = list(measure = measure),
    col_info = get_raw_column_info(measure), metadata = list(chambers = c(1L, 2L))
  )
  compared <- compare_raw_regressions(prepared, c(2, 4), verbose = FALSE)
  testthat::expect_identical(names(compared$variants), c("unbinned", "bin_2", "bin_4"))
  testthat::expect_equal(nrow(compared$comparison), 12)
  subset <- dplyr::filter(compared$comparison, chamber == 1, phase == "M1")
  testthat::expect_equal(subset$n_points, c(7L, 3L, 1L))
  x <- c(1.5, 3.5, 5.5) / 3600
  expected <- summary(lm(c(1.5, 4.5, 6.5) ~ x))
  testthat::expect_equal(subset$slope_per_hr[2], unname(coef(expected)[2, 1]))
  testthat::expect_equal(subset$r2[2], expected$r.squared)
  testthat::expect_true(is.na(subset$slope_per_hr[3]))
  missing <- dplyr::filter(compared$comparison, chamber == 1, phase == "M2")
  testthat::expect_equal(missing$n_points, c(0L, 0L, 0L))
  testthat::expect_true(all(is.na(missing$r2)))
  testthat::expect_error(compare_raw_regressions(prepared, c(2, 2)), "unique")
  testthat::expect_error(compare_raw_regressions(prepared, numeric()), "bin_sizes")
  empty_bins <- calculate_measure_stats(measure, get_raw_column_info(measure), 100)
  testthat::expect_equal(empty_bins$n_points, rep(0L, 4))
  testthat::expect_true(all(is.na(empty_bins$slope_per_hr)))
  single_bins <- calculate_measure_stats(measure, get_raw_column_info(measure), 1)
  unbinned <- calculate_measure_stats(measure, get_raw_column_info(measure))
  testthat::expect_equal(single_bins, unbinned)
})

testthat::test_that("raw stages and legacy wrapper produce downstream-compatible results", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  seconds <- 0:601
  start <- as.POSIXct("2023-07-26 09:00:00", tz = "UTC")
  input <- data.frame(
    DateTime = format(start + seconds, "%m/%d/%Y/%I:%M:%S %p", tz = "UTC"),
    Phase = c("F1", rep("M1", 300), "F2", rep("M2", 300)),
    `CH1 O2 mgO2` = 10 - seconds / 1000 + sin(seconds) / 1000,
    `CH1 temp` = 20,
    `CH2 O2 mgO2` = 8 - seconds / 2000 + cos(seconds) / 1000,
    `CH2 temp` = 21,
    check.names = FALSE
  )
  summary <- data.frame(
    experiment = 18L, chamber = c(1L, 2L), species = c("Fish", "Control"),
    id = c("A", "Control"), www_kg = c(.1, NA),
    resp_vol_ml = c(100, NA), ch_vol_ml = c(500, 1000)
  )
  input_path <- file.path(dir, "raw.csv")
  summary_path <- file.path(dir, "summary.csv")
  write.csv(input, input_path, row.names = FALSE)
  write.csv(summary, summary_path, row.names = FALSE, na = "")
  prefix <- file.path(dir, "results", "sample")
  prepared <- process_raw_part1(
    input_path, "2023-07-26 09:00:00", "2023-07-26 09:10:01",
    tz = "UTC", chamber_numbers = c(1, 2), export_path = prefix, verbose = FALSE
  )
  checkpoint <- file.path(dir, "raw-part1.rds")
  saveRDS(prepared, checkpoint)
  prepared <- readRDS(checkpoint)
  compared <- compare_raw_regressions(
    prepared, c(10, 25), summary_path, 18, export_path = prefix, verbose = FALSE
  )
  legacy <- process_raw_data(
    input_path, "2023-07-26 09:00:00", "2023-07-26 09:10:01",
    tz = "UTC", summary_file = summary_path, experiment = 18,
    chamber_numbers = c(1, 2),
    checkpoint_file = file.path(dir, "mo2-input.rds"), verbose = FALSE
  )
  testthat::expect_equal(compared$variants$unbinned$measure_stats, legacy$measure_stats)
  testthat::expect_equal(
    load_raw_checkpoint(file.path(dir, "mo2-input.rds")), legacy$part1_results
  )
  testthat::expect_equal(
    dplyr::filter(compared$comparison, variant == "bin_10")$n_points,
    rep(30L, 4)
  )
  testthat::expect_true(file.exists(paste0(prefix, "_measure.csv")))
  testthat::expect_true(file.exists(paste0(prefix, "_regression_comparison.csv")))
  testthat::expect_true(file.exists(paste0(prefix, "_bin_25_measure_stats.csv")))
  testthat::expect_warning(
    process_raw_part2(prepared, summary_path, 18, bin_size = 301, verbose = FALSE),
    "No usable regressions for control chamber"
  )

  root <- dirname(normalizePath(raw_script))
  env <- new.env(parent = environment(process_raw_part2))
  env$resp_config <- list(thresholds = list(
    default_r2 = .9, default_control_r2 = .85,
    min_good_control_points = 3, control_regression_r2 = .6,
    high_background_pct = 50
  ))
  for (name in c("Part2A_Control_Chamber_Analysis.R",
                 "Part2B_Background_Correction.R", "Part3_Quality_Control_Filtering.R")) {
    sys.source(file.path(root, name), envir = env)
  }
  for (variant in compared$variants) {
    result <- variant$part1_results
    testthat::expect_equal(length(result$chamber_list), 2)
    testthat::expect_equal(nrow(result$chamber_list[[1]]), 2)
    qc <- env$process_part3_quality_control(result, verbose = FALSE)
    control <- suppressMessages(env$process_control_chambers(result, verbose = FALSE))
    corrected <- env$process_background_correction(result, control, verbose = FALSE)
    testthat::expect_equal(length(qc$part1_results$chamber_list), 2)
    testthat::expect_equal(length(corrected$corrected_chambers), 2)
    post_qc <- env$process_part3_quality_control(result, corrected, verbose = FALSE)
    testthat::expect_equal(nrow(post_qc$raw_corrected), 2)
  }
})

testthat::test_that("the raw entry point loads stage files outside the repository", {
  script_path <- normalizePath(raw_script)
  dir <- tempfile()
  dir.create(dir)
  old_dir <- getwd()
  on.exit({
    setwd(old_dir)
    unlink(dir, recursive = TRUE)
  }, add = TRUE)
  setwd(dir)
  env <- new.env(parent = globalenv())
  source(script_path, local = env)
  testthat::expect_true(is.function(env$process_raw_part1))
  testthat::expect_true(is.function(env$process_raw_part2))
  testthat::expect_true(is.function(env$compare_raw_regressions))
})