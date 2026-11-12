#!/usr/bin/env Rscript

library(tidyverse)


prepare_regression_points <- function(measure_data, col_info, bin_size = NULL) {
  if (!is.null(bin_size) &&
      (length(bin_size) != 1L || !is.numeric(bin_size) ||
       is.na(bin_size) || !is.finite(bin_size) || bin_size < 1 ||
       bin_size != floor(bin_size) || bin_size > .Machine$integer.max)) {
    stop("bin_size must be NULL or a single positive whole number.")
  }
  phase_col <- col_info$names$phase
  o2_cols <- col_info$names$o2

  points <- measure_data %>%
    select(all_of(c(phase_col, "time_sec", "time_hour", o2_cols))) %>%
    pivot_longer(
      cols = all_of(o2_cols),
      names_to = "chamber",
      values_to = "o2"
    ) %>%
    mutate(
      chamber = as.integer(str_extract(chamber, "(?i)(?<=CH)\\d+")),
      phase = toupper(trimws(as.character(.data[[phase_col]])))
    ) %>%
    filter(is.finite(o2), is.finite(time_sec), is.finite(time_hour)) %>%
    arrange(chamber, phase, time_sec)

  if (!is.null(bin_size)) {
    points <- points %>%
      group_by(phase, chamber) %>%
      mutate(bin = (row_number() - 1L) %/% bin_size + 1L) %>%
      group_by(phase, chamber, bin) %>%
      filter(n() == bin_size) %>%
      summarise(
        time_sec = mean(time_sec),
        time_hour = mean(time_hour),
        o2 = mean(o2),
        .groups = "drop"
      )
  }
  points
}


calculate_measure_stats <- function(measure_data, col_info, bin_size = NULL) {
  points <- prepare_regression_points(measure_data, col_info, bin_size)
  summarise_regression_points(points, measure_data, col_info)
}


summarise_regression_points <- function(points, measure_data, col_info) {
  phases <- unique(toupper(trimws(as.character(
    measure_data[[col_info$names$phase]]
  ))))
  chambers <- as.integer(str_extract(col_info$names$o2, "(?i)(?<=CH)\\d+"))
  if (anyNA(chambers) || anyDuplicated(chambers)) {
    stop("O2 column names must identify unique chamber numbers.")
  }
  combinations <- tidyr::expand_grid(phase = phases, chamber = chambers)
  stats <- points %>%
    group_by(phase, chamber) %>%
    summarise(
      n_points = n(),
      start_time_sec = first(time_sec),
      end_time_sec = last(time_sec),
      start_time_hr = first(time_hour),
      end_time_hr = last(time_hour),
      r2 = if (n() >= 3 && n_distinct(time_hour) >= 2)
        summary(lm(o2 ~ time_hour))$r.squared else NA_real_,
      slope_per_hr = if (n() >= 3 && n_distinct(time_hour) >= 2)
        coef(lm(o2 ~ time_hour))[[2]] else NA_real_,
      mean_o2 = mean(o2),
      first_o2 = first(o2),
      last_o2 = last(o2),
      .groups = "drop"
    )
  combinations %>%
    left_join(stats, by = c("phase", "chamber")) %>%
    mutate(n_points = coalesce(n_points, 0L)) %>%
    mutate(
      measure_number = as.integer(str_extract(phase, "\\d+"))
    ) %>%
    arrange(measure_number, chamber) %>%
    select(-measure_number)
}


add_mo2 <- function(measure_stats, summary_file, experiment = NULL) {
  summary_data <- read.csv(summary_file, stringsAsFactors = FALSE)

  if (!is.null(experiment)) {
    if (!"experiment" %in% names(summary_data)) {
      stop("The summary file must contain an experiment column when experiment is set.")
    }
    summary_data <- summary_data %>%
      filter(.data$experiment == .env$experiment)
  }

  if (nrow(summary_data) == 0) {
    stop("No rows found in the summary file for the selected experiment.")
  }

if (!all(c("chamber", "resp_vol_ml", "www_kg") %in% names(summary_data))) {
    stop("Summary file must contain chamber, resp_vol_ml, and www_kg columns.")
  }
  summary_data$chamber <- as.integer(
    stringr::str_extract(as.character(summary_data$chamber), "\\d+")
  )
  if (anyNA(summary_data$chamber)) {
    stop("Summary file contains invalid chamber identifiers.")
  }
  summary_data <- summary_data %>%
    select(chamber, resp_vol_ml, www_kg)

  if (anyDuplicated(summary_data$chamber)) {
    stop("Summary file has more than one row per chamber. Filter to a single experiment.")
  }

  missing_chambers <- setdiff(unique(measure_stats$chamber), summary_data$chamber)
  if (length(missing_chambers) > 0) {
    stop("No summary data for chamber(s): ", paste(missing_chambers, collapse = ", "))
  }

  measure_stats %>%
    left_join(summary_data, by = "chamber", relationship = "many-to-one") %>%
    mutate(
      loop = as.integer(stringr::str_extract(as.character(phase), "\\d+")),
      phase = stringr::str_extract(as.character(phase), "^[A-Za-z]+"),
      mo2_mg_kg_hr = -1 * slope_per_hr * (resp_vol_ml / 1000) / www_kg
    ) %>%
    relocate(phase, loop)
}

# Convert raw measure regressions into the chamber-list shape expected by
# Part 2A/2B/3. The summary file supplies the chamber metadata used downstream.
build_raw_part1_results <- function(
    raw_data,
    time_series,
    measure_stats,
    col_info,
    summary_file,
    experiment = NULL,
    chamber_numbers,
    file_path,
    experiment_start_time,
    experiment_end_time) {
  if (is.null(summary_file) || !file.exists(summary_file)) {
    stop("A valid summary_file is required to prepare Part 2A/2B input.")
  }

  summary_data <- read.csv(summary_file, stringsAsFactors = FALSE)
  required <- c(
    "chamber", "species", "id", "www_kg", "resp_vol_ml", "ch_vol_ml"
  )
  missing <- setdiff(required, names(summary_data))
  if (length(missing) > 0L) {
    stop("Summary file is missing required columns: ",
         paste(missing, collapse = ", "))
  }
  if (!is.null(experiment) && "experiment" %in% names(summary_data)) {
    summary_data <- dplyr::filter(
      summary_data,
      .data$experiment == .env$experiment
    )
  }


  summary_data$chamber <- as.integer(
    stringr::str_extract(as.character(summary_data$chamber), "\\d+")
  )
  if (anyNA(summary_data$chamber) || anyDuplicated(summary_data$chamber)) {
    stop("Summary data must contain exactly one valid row per chamber.")
  }

  desc_data <- summary_data[
    match(chamber_numbers, summary_data$chamber),
    ,
    drop = FALSE
  ]
  if (anyNA(desc_data$chamber)) {
    absent <- chamber_numbers[is.na(desc_data$chamber)]
    stop("Summary data is missing chamber(s): ", paste(absent, collapse = ", "))
  }
  if (!"experiment" %in% names(desc_data)) {
    desc_data$experiment <- if (is.null(experiment)) NA_integer_ else experiment
  }

  phase_col <- col_info$names$phase
  phase <- toupper(trimws(as.character(time_series$measure[[phase_col]])))
  temp_cols <- col_info$names$temp
  if (length(temp_cols) == 1L && is.na(temp_cols)) {
    avg_temp <- rep(NA_real_, nrow(time_series$measure))
  } else {
    avg_temp <- rowMeans(
      time_series$measure[, temp_cols, drop = FALSE],
      na.rm = TRUE
    )
    avg_temp[!is.finite(avg_temp)] <- NA_real_
  }
  phase_times <- tibble::tibble(
    phase_key = phase,
    time = time_series$measure$time_hour,
    avg_temp = avg_temp
  ) %>%
    dplyr::filter(grepl("^M[0-9]+$", .data$phase_key)) %>%
    dplyr::group_by(.data$phase_key) %>%
    dplyr::summarise(
      time = mean(.data$time, na.rm = TRUE),
      avg_temp = mean(.data$avg_temp, na.rm = TRUE),
      .groups = "drop"
    )
  phase_times$avg_temp[!is.finite(phase_times$avg_temp)] <- NA_real_

  measure_stats <- measure_stats %>%
    dplyr::mutate(
      chamber = as.integer(.data$chamber),
      phase_key = paste0(.data$phase, .data$loop)
    )
  phase_keys <- phase_times$phase_key
  if (length(phase_keys) == 0L) {
    stop("No measure phases are available for the selected experiment.")
  }

  chamber_list <- lapply(seq_len(nrow(desc_data)), function(i) {
    chamber_info <- desc_data[i, , drop = FALSE]
    chamber_number <- chamber_info$chamber[[1]]
    experiment_value <- if ("experiment" %in% names(chamber_info)) {
      chamber_info$experiment[[1]]
    } else if (!is.null(experiment)) {
      experiment
    } else {
      NA_integer_
    }
    chamber_stats <- measure_stats %>%
      dplyr::filter(.data$chamber == chamber_number) %>%
      dplyr::select(
        phase_key,
        r2,
        slope = slope_per_hr,
        raw_mo2_kg_www = mo2_mg_kg_hr
      )

    chamber_data <- tibble::tibble(phase_key = phase_keys) %>%
      dplyr::left_join(phase_times, by = "phase_key") %>%
      dplyr::left_join(chamber_stats, by = "phase_key") %>%
      dplyr::mutate(
        experiment = experiment_value,
        chamber = chamber_number,
        species = chamber_info$species[[1]],
        id = as.character(chamber_info$id[[1]]),
        raw_mo2_ch = .data$raw_mo2_kg_www * chamber_info$www_kg[[1]]
      )

    if (identical(chamber_data$id[[1]], "Control")) {
      chamber_data <- chamber_data %>%
        dplyr::mutate(
          raw_mo2_ch = abs(.data$slope * chamber_info$ch_vol_ml[[1]] / 1000),
          raw_mo2_kg_www = NA_real_
        )
    }

    dplyr::select(
      chamber_data,
      experiment, chamber, species, id, time, avg_temp, r2, slope,
      raw_mo2_kg_www, raw_mo2_ch
    )
  })

  list(
    raw_data = raw_data,
    descriptive_data = desc_data,
    time_series = time_series,
    chamber_list = chamber_list,
    col_info = col_info,
    metadata = list(
      file_path = file_path,
      experiment_start = experiment_start_time,
      experiment_end = experiment_end_time,
      chambers = chamber_numbers,
      n_chambers = nrow(desc_data),
      n_controls = sum(desc_data$id == "Control", na.rm = TRUE),
      n_experimental = sum(desc_data$id != "Control", na.rm = TRUE)
    )
  )
}


process_raw_part2 <- function(
    raw_part1,
    summary_file = NULL,
    experiment = NULL,
    bin_size = NULL,
    export_path = NULL,
    verbose = TRUE) {
  measure <- raw_part1$time_series$measure
  col_info <- raw_part1$col_info
  if (is.null(measure) || is.null(col_info) || is.null(raw_part1$metadata)) {
    stop("raw_part1 must be the result of process_raw_part1().")
  }
  points <- prepare_regression_points(measure, col_info, bin_size)
  stats <- summarise_regression_points(points, measure, col_info)
  if (!is.null(summary_file)) {
    stats <- add_mo2(stats, summary_file, experiment)
  }
  part1_results <- NULL
  if (!is.null(summary_file)) {
    part1_results <- build_raw_part1_results(
      raw_data = raw_part1$raw_data,
      time_series = raw_part1$time_series,
      measure_stats = stats,
      col_info = col_info,
      summary_file = summary_file,
      experiment = experiment,
      chamber_numbers = raw_part1$metadata$chambers,
      file_path = raw_part1$metadata$file_path,
      experiment_start_time = raw_part1$metadata$experiment_start,
      experiment_end_time = raw_part1$metadata$experiment_end
    )
    part1_results$metadata["bin_size"] <- list(bin_size)
    control_idx <- which(part1_results$descriptive_data$id == "Control")
    missing_controls <- vapply(control_idx, function(i) {
      chamber <- part1_results$chamber_list[[i]]
      !any(is.finite(chamber$slope) & is.finite(chamber$r2))
    }, logical(1))
    if (any(missing_controls)) {
      warning(
        "No usable regressions for control chamber(s): ",
        paste(part1_results$descriptive_data$chamber[control_idx[missing_controls]],
              collapse = ", "),
        ". Main Part 2A treats insufficient control data as negligible background; ",
        "inspect these results before background correction.",
        call. = FALSE
      )
    }
  }
  if (!is.null(export_path)) {
    dir.create(dirname(export_path), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(stats, paste0(export_path, "_measure_stats.csv"))
  }
  if (verbose) {
    message("Calculated regressions for ", nrow(stats), " chamber-periods.")
  }
  list(
    regression_points = points,
    measure_stats = stats,
    part1_results = part1_results,
    bin_size = bin_size
  )
}


compare_raw_regressions <- function(
    raw_part1,
    bin_sizes,
    summary_file = NULL,
    experiment = NULL,
    export_path = NULL,
    verbose = TRUE) {
  if (!is.numeric(bin_sizes) || length(bin_sizes) == 0L || anyNA(bin_sizes) ||
      any(!is.finite(bin_sizes)) || any(bin_sizes < 1) ||
      any(bin_sizes != floor(bin_sizes)) ||
      any(bin_sizes > .Machine$integer.max) || anyDuplicated(bin_sizes)) {
    stop("bin_sizes must contain unique positive whole numbers.")
  }
  sizes <- c(list(NULL), as.list(bin_sizes))
  labels <- c("unbinned", paste0("bin_", bin_sizes))
  variants <- setNames(lapply(seq_along(sizes), function(i) {
    process_raw_part2(
      raw_part1, summary_file, experiment, bin_size = sizes[[i]],
      export_path = if (is.null(export_path)) NULL else
        paste0(export_path, "_", labels[[i]]),
      verbose = verbose
    )
  }), labels)
  comparison <- bind_rows(lapply(seq_along(variants), function(i) {
    stats <- variants[[i]]$measure_stats
    phase_key <- if ("loop" %in% names(stats)) {
      paste0(stats$phase, stats$loop)
    } else {
      stats$phase
    }
    tibble(
      variant = labels[[i]],
      bin_size = if (is.null(sizes[[i]])) NA_integer_ else as.integer(sizes[[i]]),
      chamber = stats$chamber,
      phase = phase_key,
      n_points = stats$n_points,
      slope_per_hr = stats$slope_per_hr,
      r2 = stats$r2
    )
  }))
  if (!is.null(export_path)) {
    readr::write_csv(comparison, paste0(export_path, "_regression_comparison.csv"))
  }
  list(variants = variants, comparison = comparison)
}