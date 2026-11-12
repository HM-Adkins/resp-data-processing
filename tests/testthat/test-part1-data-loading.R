library(dplyr)
library(stringr)
library(tibble)

part1_candidates <- c(
  "Part1_Load_And_Prepare_Data.R",
  file.path("..", "Part1_Load_And_Prepare_Data.R"),
  file.path("..", "..", "Part1_Load_And_Prepare_Data.R")
)
part1_path <- part1_candidates[file.exists(part1_candidates)][1]
if (is.na(part1_path)) stop("Could not locate Part1_Load_And_Prepare_Data.R.")
source(part1_path, local = TRUE)

part1_col_info <- list(
  idx = list(datetime = 7L, time_sec = 8L),
  names = list(
    chamber = "Chamber",
    www_g = "Wet weight",
    resp_vol_ml = "Resp volume",
    chamber_notes = "Chamber notes",
    tww_g = "Tissue weight",
    ch_vol_ml = "Chamber volume",
    datetime = "DateTime",
    time_sec = "Time sec",
    mo2 = c("CH2 MO2", "CH1 MO2"),
    r2 = c("CH2 R2", "CH1 R2"),
    slope = "CH2 slope",
    temp = character()
  )
)

testthat::test_that("metadata uses the explicit experiment and checks chamber notes", {
  raw <- data.frame(
    Chamber = c("Chamber 1", "Chamber 2"),
    `Wet weight` = c(100, NA),
    `Resp volume` = c(250, NA),
    `Chamber notes` = c("Fish-A1", "Control-Control"),
    `Tissue weight` = c(80, NA),
    `Chamber volume` = c(500, 500),
    check.names = FALSE
  )

  testthat::expect_no_message(
    desc <- extract_descriptive_data(raw, part1_col_info, 123, verbose = FALSE)
  )
  testthat::expect_equal(desc$experiment, c(123L, 123L))
  testthat::expect_equal(desc$species, c("Fish", "Control"))
  testthat::expect_equal(desc$id, c("A1", "Control"))
  testthat::expect_equal(desc$www_kg[1], 0.1)
  testthat::expect_equal(desc$tww_g[1], 80)
  testthat::expect_equal(desc$tww_kg[1], 0.08)

  raw$`Chamber notes`[1] <- "missing delimiter"
  testthat::expect_error(
    extract_descriptive_data(raw, part1_col_info, 123, verbose = FALSE),
    "Could not parse species and ID"
  )
})

testthat::test_that("tissue wet weight is identified using the tww abbreviation", {
  input <- data.frame(
    Chamber = 1,
    www_g = 100,
    `Resp volume` = 250,
    `Chamber notes` = "Fish-A1",
    tww_g = 80,
    `Chamber volume` = 500,
    check.names = FALSE
  )

  info <- get_column_info(input)
  testthat::expect_identical(info$names$tww_g, "tww_g")
  testthat::expect_silent(validate_columns(info, c("chamber", "tww_g")))
  testthat::expect_error(
    validate_columns(info, c("datetime", "slope")),
    "Missing required columns: datetime, slope"
  )
  desc <- extract_descriptive_data(input, info, experiment = 1, verbose = FALSE)
  testthat::expect_equal(desc$tww_g, 80)
  testthat::expect_equal(desc$tww_kg, 0.08)
})

testthat::test_that("time series parsing uses the supplied time zone", {
  data <- data.frame(
    DateTime = c(
      "07/26/2023 08:00:00 AM",
      "07/26/2023 09:00:00 AM",
      "07/26/2023 10:00:00 AM"
    ),
    `Time sec` = c(0, 3600, 7200),
    check.names = FALSE
  )
  info <- list(
    names = list(datetime = "DateTime", time_sec = "Time sec")
  )

  series <- prepare_time_series(
    data, info, "2023-07-26 09:00:00", tz = "America/Chicago",
    date_format = "%m/%d/%Y %I:%M:%S %p"
  )

  testthat::expect_equal(series$all$time_hour, c(-1, 0, 1))
  testthat::expect_equal(
    as.numeric(series$all$DateTime[2]),
    as.numeric(as.POSIXct("2023-07-26 09:00:00", tz = "America/Chicago"))
  )
  testthat::expect_error(
    prepare_time_series(data, info, "2023-07-26 09:00:00", tz = "invalid"),
    "valid time-zone"
  )

  data$DateTime <- sub("2023 ", "2023/", data$DateTime)
  default_series <- prepare_time_series(
    data, info, "2023-07-26 09:00:00", tz = "America/Chicago"
  )
  testthat::expect_equal(default_series$all$time_hour, c(-1, 0, 1))
})

testthat::test_that("chamber measurements are matched by chamber number", {
  desc <- data.frame(
    experiment = c(1L, 1L),
    chamber = c(1L, 2L),
    species = c("Fish", "Control"),
    id = c("A1", "Control"),
    www_kg = c(0.1, NA_real_),
    ch_vol_ml = c(500, 500)
  )
  data <- data.frame(
    time_hour = c(0, 1),
    `CH2 MO2` = c(20, 21),
    `CH1 MO2` = c(10, 11),
    `CH2 R2` = c(0.9, 0.9),
    `CH1 R2` = c(0.8, 0.8),
    `CH.2.slope..mgO2.L.hr.` = c(-0.2, -0.3),
    `CH.1.slope..mgO2.L.hr.` = c(-1, -2),
    check.names = FALSE
  )
  info <- list(
    names = list(
      mo2 = c("CH2 MO2", "CH1 MO2"),
      r2 = c("CH2 R2", "CH1 R2"),
      slope = c("CH.2.slope..mgO2.L.hr.", "CH.1.slope..mgO2.L.hr."),
      temp = character()
    )
  )

  chambers <- build_chamber_list(data, desc, info)

  testthat::expect_equal(chambers[[1]]$raw_mo2_kg_www, c(10, 11))
  testthat::expect_equal(chambers[[1]]$r2, c(0.8, 0.8))
  testthat::expect_equal(chambers[[1]]$slope, c(-1, -2))
  testthat::expect_equal(chambers[[2]]$slope, c(-0.2, -0.3))
  testthat::expect_identical(
    match_chamber_columns(info$names$mo2, "MO2", desc$chamber),
    c("CH1 MO2", "CH2 MO2")
  )
  testthat::expect_identical(
    match_chamber_columns(c("CH.2.MO2", "CH.1.MO2"), "MO2", desc$chamber),
    c("CH.1.MO2", "CH.2.MO2")
  )
  testthat::expect_error(
    build_chamber_list(data, desc, utils::modifyList(info, list(names = list(
      mo2 = info$names$mo2,
      r2 = info$names$r2,
      slope = character(),
      temp = character()
    )))),
    "No slope columns"
  )
  missing_slope <- info
  missing_slope$names$slope <- "CH.2.slope..mgO2.L.hr."
  testthat::expect_error(
    build_chamber_list(data, desc, missing_slope),
    "Missing slope column\\(s\\) for chamber\\(s\\): 1"
  )
})

testthat::test_that("raw plots use the same chamber mapping as the data", {
  data <- data.frame(
    time_hour = c(0, 1),
    `CH2 MO2` = c(20, 21),
    `CH1 MO2` = c(10, 11),
    check.names = FALSE
  )
  desc <- data.frame(
    chamber = c(1L, 2L),
    species = c("Fish", "Fish"),
    id = c("A", "B")
  )
  info <- list(names = list(mo2 = c("CH2 MO2", "CH1 MO2")))
  recorded <- list()
  plot_env <- new.env(parent = environment(plot_raw_data))
  plot_env$resp_config <- list(
    plot = list(point_size = 1),
    colors = list(blue = "blue", red = "red")
  )
  plot_env$plot_device <- function(output_file, width, height, expr) {
    eval.parent(substitute(expr))
  }
  plot_env$plot <- function(x, y, ...) {
    recorded[[length(recorded) + 1L]] <<- y
  }
  plot_env$par <- function(...) c(0, 1, 0, 1)
  plot_env$polygon <- function(...) NULL
  plot_env$abline <- function(...) NULL
  plot_function <- plot_raw_data
  environment(plot_function) <- plot_env

  plot_function(data, desc, info)

  testthat::expect_equal(recorded, list(c(10, 11), c(20, 21)))
})

testthat::test_that("raw control plots match chambers and label per-chamber units", {
  data <- data.frame(
    time_hour = c(-1, 0, 1),
    `CH.3.MO2` = c(3, 4, 5),
    `CH.1.MO2` = c(10, 11, 12),
    `CH.2.MO2` = c(-1, 0, 1),
    check.names = FALSE
  )
  desc <- data.frame(chamber = c(1L, 2L, 3L),
                     id = c("Fish", "Control", "Control"))
  info <- list(names = list(mo2 = c("CH.3.MO2", "CH.1.MO2", "CH.2.MO2")))
  recorded <- list()
  plot_env <- new.env(parent = environment(plot_raw_control_data))
  plot_env$resp_config <- list(
    plot = list(point_size = 1), colors = list(blue = "blue", red = "red")
  )
  plot_env$plot_device <- function(output_file, width, height, expr) {
    eval.parent(substitute(expr))
  }
  plot_env$plot <- function(x, y, ...) {
    recorded[[length(recorded) + 1L]] <<- list(y = y, args = list(...))
  }
  plot_env$par <- function(...) c(0, 1, 0, 1)
  plot_env$polygon <- function(...) NULL
  plot_env$abline <- function(...) NULL
  plot_function <- plot_raw_control_data
  environment(plot_function) <- plot_env

  plot_function(data, desc, info)

  testthat::expect_equal(length(recorded), 2)
  testthat::expect_equal(recorded[[1]]$y, c(-1, 0, 1))
  testthat::expect_equal(recorded[[2]]$y, c(3, 4, 5))
  testthat::expect_equal(recorded[[1]]$args$main, "Control - Chamber 2")
  testthat::expect_equal(recorded[[2]]$args$main, "Control - Chamber 3")
  testthat::expect_equal(recorded[[1]]$args$ylab, "MO2 (mg O2/chamber/hr)")
  testthat::expect_equal(recorded[[1]]$args$ylim, c(-1, 5))
  testthat::expect_equal(recorded[[2]]$args$ylim, c(-1, 5))

  desc$id <- "Fish"
  testthat::expect_error(plot_function(data, desc, info), "No control chamber")
  desc$id[2] <- "Control"
  data$`CH.2.MO2` <- NA_real_
  testthat::expect_error(plot_function(data, desc, info), "no finite time and MO2 pairs")
})

testthat::test_that("raw control plots handle constant and single-control data", {
  data <- data.frame(time_hour = c(-Inf, 0, 1), `CH.2.MO2` = c(100, 0, 0),
                     check.names = FALSE)
  desc <- data.frame(chamber = 2L, id = "Control")
  info <- list(names = list(mo2 = "CH.2.MO2"))
  plot_env <- new.env(parent = environment(plot_raw_control_data))
  plot_env$resp_config <- list(
    plot = list(point_size = 1), colors = list(blue = "blue", red = "red")
  )
  plot_env$plot_device <- function(output_file, width, height, expr) {
    eval.parent(substitute(expr))
  }
  plot_env$plot <- function(x, y, ...) {
    testthat::expect_equal(list(...)$ylim, c(-1, 1))
  }
  plot_env$par <- function(...) c(0, 1, 0, 1)
  plot_env$polygon <- function(x, ...) {
    testthat::expect_true(all(is.finite(x)))
  }
  plot_env$abline <- function(...) NULL
  plot_function <- plot_raw_control_data
  environment(plot_function) <- plot_env
  plot_function(data, desc, info)
})

testthat::test_that("Part 1 writes a separate control PDF when plotting is enabled", {
  library_candidates <- c(
    "Resp_Data_Processing_Library.R",
    file.path("..", "Resp_Data_Processing_Library.R"),
    file.path("..", "..", "Resp_Data_Processing_Library.R")
  )
  library_path <- library_candidates[file.exists(library_candidates)][1]
  if (is.na(library_path)) stop("Could not locate Resp_Data_Processing_Library.R.")
  expressions <- parse(library_path)
  device_expression <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("plot_device"))
  }, expressions)[[1]]
  env <- new.env(parent = environment(process_part1_data_loading))
  sys.source(part1_path, envir = env)
  eval(device_expression, envir = env)
  env$resp_config <- list(
    plot = list(point_size = 1), colors = list(blue = "blue", red = "red")
  )

  dir <- tempfile()
  dir.create(dir)
  old_dir <- getwd()
  on.exit({
    setwd(old_dir)
    unlink(dir, recursive = TRUE)
  }, add = TRUE)
  setwd(dir)
  input <- data.frame(
    Chamber = c("CH1", "CH2", "CH3"),
    www_g = c(100, NA, NA),
    `Resp vol` = c(250, NA, NA),
    `Chamber notes` = c("Fish-A", "Control-Control", "Control-Control"),
    tww_g = c(80, NA, NA),
    `Chamber vol` = c(500, 500, 500),
    DateTime = c("07/26/2023/08:00:00 AM", "07/26/2023/09:00:00 AM",
                 "07/26/2023/10:00:00 AM"),
    `Time sec` = c(0, 3600, 7200),
    `CH1 MO2` = c(8, 9, 10), `CH2 MO2` = c(1, 2, 3),
    `CH3 MO2` = c(2, 3, 4), `CH1 R2` = c(.9, .9, .9),
    `CH2 R2` = c(.9, .9, .9), `CH3 R2` = c(.9, .9, .9),
    `CH1 slope` = c(-1, -1, -1), `CH2 slope` = c(-2, -2, -2),
    `CH3 slope` = c(-3, -3, -3), check.names = FALSE
  )
  utils::write.csv(input, "input.csv", row.names = FALSE, na = "")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressMessages(env$process_part1_data_loading(
    "input.csv", "2023-07-26 09:00:00", experiment = 18,
    tz = "America/Chicago", generate_raw_data_plot = TRUE,
    export_path = "plots/example", verbose = FALSE
  ))
  testthat::expect_true(file.exists("plots/example_raw_data_plots.pdf"))
  testthat::expect_true(file.exists("plots/example_raw_control_plots.pdf"))
  testthat::expect_gt(file.info("plots/example_raw_control_plots.pdf")$size, 0)
})