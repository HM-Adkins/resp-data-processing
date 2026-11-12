library(dplyr)
library(tibble)
library(forcats)

find_script <- function(filename) {
  candidates <- c(filename, file.path("..", filename),
                  file.path("..", "..", filename))
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) stop("Could not locate ", filename)
  path
}

source(find_script("Part3_Quality_Control_Filtering.R"), local = TRUE)
source(find_script("Part2A_Control_Chamber_Analysis.R"), local = TRUE)
source(find_script("Part2B_Background_Correction.R"), local = TRUE)

resp_config <- list(thresholds = list(
  default_r2 = 0.9,
  min_good_control_points = 3,
  control_regression_r2 = 0.6,
  high_background_pct = 50
))

make_part1 <- function() {
  desc <- data.frame(
    experiment = c(18L, 18L), chamber = c(1L, 2L),
    species = c("Fish", "Control"), id = c("A", "Control"),
    www_kg = c(0.1, NA_real_), ch_vol_ml = c(500, 500)
  )
  chamber <- function(i, r2) {
    tibble(
      experiment = desc$experiment[i], chamber = desc$chamber[i],
      species = desc$species[i], id = desc$id[i], time = 0:3,
      avg_temp = rep(20, 4), r2 = r2,
      raw_mo2_kg_www = if (i == 1) rep(10, 4) else rep(NA_real_, 4),
      raw_mo2_ch = rep(1, 4)
    )
  }
  control <- chamber(2, c(0.8, 0.9, NA, 0.6))
  control$slope <- c(-0.2, -0.3, -0.4, -0.5)
  list(
    descriptive_data = desc,
    chamber_list = list(chamber(1, c(0.89, 0.9, NA, 0.95)), control),
    metadata = list(file_path = "example.csv")
  )
}

testthat::test_that("Part 1 quality checks flag every row without altering inputs", {
  part1 <- make_part1()
  original <- part1
  qc <- process_part3_quality_control(part1, r2_threshold = 0.9, verbose = FALSE)

  testthat::expect_identical(part1, original)
  testthat::expect_equal(qc$part1_results$chamber_list[[1]]$r2_quality,
                         c("low_r2", "pass", "missing_r2", "pass"))
  testthat::expect_equal(qc$part1_results$chamber_list[[2]]$r2_quality,
                         c("low_r2", "pass", "missing_r2", "low_r2"))
  testthat::expect_equal(qc$summary$n_raw, c(4, 4))
  testthat::expect_equal(qc$summary$n_pass_r2, c(2, 1))
  testthat::expect_equal(qc$summary$n_low_r2, c(1, 2))
  testthat::expect_equal(qc$summary$n_missing_r2, c(1, 1))
  for (i in seq_along(part1$chamber_list)) {
    flagged <- qc$part1_results$chamber_list[[i]]
    flagged$r2_quality <- NULL
    testthat::expect_identical(flagged, part1$chamber_list[[i]])
  }
  testthat::expect_identical(qc$part1_results$descriptive_data, original$descriptive_data)
})

testthat::test_that("quality flags do not change control or correction calculations", {
  part1 <- make_part1()
  qc <- process_part3_quality_control(part1, r2_threshold = 0.9, verbose = FALSE)
  baseline_control <- suppressMessages(process_control_chambers(
    part1, control_r2_threshold = 0.85, verbose = FALSE
  ))
  flagged_control <- suppressMessages(process_control_chambers(
    qc$part1_results, control_r2_threshold = 0.85, verbose = FALSE
  ))
  testthat::expect_equal(flagged_control$data, baseline_control$data)
  baseline <- process_background_correction(part1, baseline_control, verbose = FALSE)
  flagged <- process_background_correction(qc$part1_results, flagged_control, verbose = FALSE)
  for (i in seq_along(part1$chamber_list)) {
    testthat::expect_equal(
      flagged$corrected_chambers[[i]]$corrected_mo2_ch,
      baseline$corrected_chambers[[i]]$corrected_mo2_ch
    )
  }
  result <- process_part3_quality_control(part1, baseline, r2_threshold = 0.9,
                                          verbose = FALSE)
  testthat::expect_equal(result$data$exclusion_reason,
                         c("low_r2", "included", "missing_r2", "included"))
  testthat::expect_equal(result$sample_sizes$n_final, 2)
})

testthat::test_that("invalid inputs and corrected-only plots fail explicitly", {
  part1 <- make_part1()
  testthat::expect_error(process_part3_quality_control(part1, r2_threshold = NA_real_),
                         "r2_threshold")
  testthat::expect_error(process_part3_quality_control(part1, r2_threshold = 1.1),
                         "r2_threshold")
  testthat::expect_error(process_part3_quality_control(part1, generate_rmr_plot = TRUE),
                         "require Part 2")
  testthat::expect_error(process_part3_quality_control(part1, generate_diagnostics_plot = TRUE),
                         "require Part 2")
  part1$chamber_list[[2]]$id <- "wrong"
  testthat::expect_error(process_part3_quality_control(part1), "does not match")
})