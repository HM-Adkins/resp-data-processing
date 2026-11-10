library(tidyverse)

#------------------------------------------------------------------------------
# CORE HELPERS
#------------------------------------------------------------------------------

get_raw_column_info <- function(data) {
  patterns <- list(
    datetime = "date.*time",
    phase = "^phase$",
    o2 = "CH\\d+.*O2.*mgO2",
    temp = "CH\\d+.*temp"
  )
  scalar_cols <- c("datetime", "phase")
  result <- list(idx = list(), names = list())

  for (nm in names(patterns)) {
    matches <- grep(patterns[[nm]], names(data), ignore.case = TRUE, perl = TRUE)
    if (length(matches) > 0) {
      if (nm %in% scalar_cols) matches <- matches[1]
      result$idx[[nm]] <- matches
      result$names[[nm]] <- names(data)[matches]
    } else {
      result$idx[[nm]] <- NA
      result$names[[nm]] <- NA_character_
    }
  }
  result
}

plot_raw_device <- function(output_file, width = 10, height = 8, expr) {
  eval.parent(substitute(expr))

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    pdf(output_file, width = width, height = height)
    eval.parent(substitute(expr))
    dev.off()
  }

  invisible(NULL)
}

#------------------------------------------------------------------------------
# PLOTTING
#------------------------------------------------------------------------------

plot_raw_instrument_data <- function(data, col_info, title_prefix = NULL, output_file = NULL) {
  plot_raw_device(output_file, 12, 8, {
    o2_cols <- col_info$names$o2
    n <- length(o2_cols)
    if (n == 0L || all(is.na(o2_cols))) {
      stop("No chamber O2 columns are available to plot.")
    }
    nr <- ceiling(sqrt(n))
    nc <- ceiling(n / nr)
    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1), oma = c(1, 1, 3, 0))

    # Shared y-axis across chambers, fitted to the data
    all_o2 <- unlist(data[o2_cols], use.names = FALSE)
    y_lim <- range(all_o2, na.rm = TRUE)
    if (any(!is.finite(y_lim))) {
      stop("Cannot plot chamber O2 data because all values are missing.")
    }
    if (diff(y_lim) == 0) {
      y_lim <- y_lim + c(-1, 1) * max(abs(y_lim[1]) * 0.05, 1)
    }

    for (i in seq_len(n)) {
      plot(data$time_hour, data[[o2_cols[i]]],
        main = paste("Chamber",
                     stringr::str_extract(o2_cols[i], "(?i)(?<=CH)\\d+")),
        xlab = ifelse(i > n - nc, "Time (hrs)", ""),
        ylab = "O2 (mg O2/L)",
        ylim = y_lim,
        pch = 16, cex = 0.5
      )
    }

    if (!is.null(title_prefix)) {
      mtext(title_prefix, outer = TRUE, font = 2, cex = 1.2)
    }
  })

  if (!is.null(output_file)) message("Raw plot: ", output_file)
}

#------------------------------------------------------------------------------
# DATA PREPARATION
#------------------------------------------------------------------------------

prepare_raw_data <- function(
    data,
    col_info,
    experiment_start_time,
    experiment_end_time,
    tz,
    date_format = "%m/%d/%Y/%I:%M:%S %p"
) {
  dt_col <- col_info$names$datetime
  data_end_indices <- c(col_info$idx$o2, col_info$idx$temp)
  data_end_indices <- data_end_indices[!is.na(data_end_indices)]
  if (length(data_end_indices) == 0L) {
    stop("No chamber O2 or temperature columns were identified.")
  }
  end_idx <- max(data_end_indices)
  start_time <- as.POSIXct(experiment_start_time, tz = tz)
  end_time <- as.POSIXct(experiment_end_time, tz = tz)
  
  if (is.na(start_time)) stop("Invalid experiment start time.")
  if (is.na(end_time)) stop("Invalid experiment end time.")
  if (start_time >= end_time) {
    stop("Experiment start time must be before end time.")
  }
  
  start_idx <- col_info$idx$datetime
  data <- data[, seq.int(start_idx, end_idx), drop = FALSE] %>%
    mutate(
      !!sym(dt_col) := as.POSIXct(
        .data[[dt_col]],
        format = date_format,
        tz = tz
      )
    )
  
  if (anyNA(data[[dt_col]])) {
    stop("Some Date & Time values could not be parsed.")
  }
  
  time_filtered <- filter(
    data,
    .data[[dt_col]] >= start_time,
    .data[[dt_col]] <= end_time
  )
  
  phase_col <- grep("^phase$", names(time_filtered), ignore.case = TRUE)
  if (length(phase_col) != 1) {
    stop("Could not identify exactly one Phase column.")
  }
  
  phase <- toupper(trimws(as.character(time_filtered[[phase_col]])))
  fwm_rows <- grepl("^[FWM][0-9]+$", phase)
  
  if (!any(grepl("^M[0-9]+$", phase))) {
    stop("No M-phase values found within the experiment time window.")
  }
  
  # Remove entire F/W/M loops whose M phase has fewer than 300 rows
  phase_number <- suppressWarnings(as.integer(sub("^[FWM]", "", phase)))
  loops <- unique(phase_number[fwm_rows])
  
  m_counts <- vapply(
    loops,
    function(loop_number) sum(phase == paste0("M", loop_number)),
    integer(1)
  )
  
  incomplete_loops <- loops[m_counts < 300]
  keep_loop <- !(fwm_rows & phase_number %in% incomplete_loops)
  
  time_filtered <- time_filtered[keep_loop, , drop = FALSE]
  phase <- toupper(trimws(as.character(time_filtered[[phase_col]])))
  m_rows <- grepl("^M[0-9]+$", phase)
  
  if (!any(m_rows)) {
    stop("No M phases with at least 300 rows remain after filtering.")
  }
  
  # Keep F/W phases only up to the highest remaining M phase
  highest_m <- max(as.integer(sub("^M", "", phase[m_rows])))
  fw_rows <- grepl("^[FW][0-9]+$", phase)
  fw_numbers <- as.integer(sub("^[FW]", "", phase[fw_rows]))
  
  keep <- rep(TRUE, nrow(time_filtered))
  keep[fw_rows] <- fw_numbers <= highest_m
  
  time_filtered <- time_filtered[keep, , drop = FALSE]
  phase <- toupper(trimws(as.character(time_filtered[[phase_col]])))
  m_rows <- grepl("^M[0-9]+$", phase)
  
  f1_rows <- phase == "F1"
  if (!any(f1_rows)) {
    stop("No F1 phase found within the retained experiment data.")
  }
  
  time_zero <- min(time_filtered[[dt_col]][f1_rows])
  
  add_time <- function(df) {
    df %>%
      mutate(
        time_sec = as.numeric(
          difftime(.data[[dt_col]], time_zero, units = "secs")
        ),
        time_hour = time_sec / 3600
      )
  }
  
  experimental <- add_time(time_filtered)
  measure <- add_time(time_filtered[m_rows, , drop = FALSE])
  
  list(
    all = data,
    acclimation = filter(data, .data[[dt_col]] < start_time),
    experimental = experimental,
    measure = measure,
    post_experiment = filter(data, .data[[dt_col]] > end_time)
  )
}

#-------------------------------------------------------------------------------

select_raw_chambers <- function(data, col_info, chamber_numbers = NULL) {
  o2_names <- col_info$names$o2
  if (length(o2_names) == 1L && is.na(o2_names)) {
    stop("No chamber O2 columns were identified.")
  }

  available <- as.integer(
    stringr::str_extract(o2_names, "(?i)(?<=CH)\\d+")
  )
  if (anyNA(available) || anyDuplicated(available)) {
    stop("O2 column names must identify unique chamber numbers.")
  }

  selected <- if (is.null(chamber_numbers)) {
    available
  } else {
    if (!is.numeric(chamber_numbers) ||
        length(chamber_numbers) == 0L ||
        anyNA(chamber_numbers) ||
        any(!is.finite(chamber_numbers)) ||
        any(chamber_numbers != as.integer(chamber_numbers)) ||
        anyDuplicated(chamber_numbers)) {
      stop("chamber_numbers must be unique, non-missing whole numbers.")
    }
    as.integer(chamber_numbers)
  }

  missing <- setdiff(selected, available)
  if (length(missing) > 0L) {
    stop(
      "Requested chamber(s) not found in the raw input: ",
      paste(missing, collapse = ", ")
    )
  }

  keep_o2 <- o2_names[available %in% selected]
  temp_names <- col_info$names$temp
  if (length(temp_names) == 1L && is.na(temp_names)) {
    temp_names <- character()
  } else {
    temp_chambers <- as.integer(
      stringr::str_extract(temp_names, "(?i)(?<=CH)\\d+")
    )
    temp_names <- temp_names[!is.na(temp_chambers) & temp_chambers %in% selected]
  }
  shared_names <- c(col_info$names$datetime, col_info$names$phase)
  keep_names <- names(data)[names(data) %in% c(shared_names, keep_o2, temp_names)]
  selected_data <- dplyr::select(data, dplyr::all_of(keep_names))

  list(
    data = selected_data,
    chamber_numbers = selected,
    col_info = get_raw_column_info(selected_data)
  )
}

calculate_measure_stats <- function(measure_data, col_info) {
  phase_col <- col_info$names$phase
  o2_cols <- col_info$names$o2

  measure_data %>%
    select(all_of(c(phase_col, "time_sec", "time_hour", o2_cols))) %>%
    pivot_longer(
      cols = all_of(o2_cols),
      names_to = "chamber",
      values_to = "o2"
    ) %>%
    mutate(
      chamber = as.integer(str_extract(chamber, "(?<=CH)\\d+")),
      phase = toupper(trimws(as.character(.data[[phase_col]])))
    ) %>%
    filter(!is.na(o2)) %>%
    arrange(chamber, time_sec) %>%
    group_by(phase, chamber) %>%
    summarise(
      n_points = n(),
      start_time_sec = first(time_sec),
      end_time_sec = last(time_sec),
      start_time_hr = first(time_hour),
      end_time_hr = last(time_hour),
      r2 = if (n() >= 3) summary(lm(o2 ~ time_hour))$r.squared else NA_real_,
      slope_per_hr = if (n() >= 3) coef(lm(o2 ~ time_hour))[[2]] else NA_real_,
      mean_o2 = mean(o2),
      first_o2 = first(o2),
      last_o2 = last(o2),
      .groups = "drop"
    ) %>%
    mutate(
      measure_number = as.integer(str_extract(phase, "\\d+"))
    ) %>%
    arrange(measure_number, chamber) %>%
    select(-measure_number)
}

#------------------------------------------------------------------------------

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
# Part 2A/2B. The summary file supplies the chamber metadata used downstream.
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

save_raw_checkpoint <- function(results, path, stage = "raw_part1") {
  if (!is.list(results)) stop("results must be a list.")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(schema_version = 1L, stage = stage, results = results),
    path
  )
  invisible(path)
}

load_raw_checkpoint <- function(path, expected_stage = "raw_part1") {
  if (!file.exists(path)) stop("Checkpoint file not found: ", path)
  checkpoint <- readRDS(path)
  if (!identical(checkpoint$schema_version, 1L) ||
      !identical(checkpoint$stage, expected_stage) ||
      !is.list(checkpoint$results)) {
    stop("Checkpoint has an unsupported schema or stage.")
  }
  checkpoint$results
}

#------------------------------------------------------------------------------
# MAIN WORKFLOW
#------------------------------------------------------------------------------

process_raw_data <- function(
  file_path,
  experiment_start_time,
  experiment_end_time,
  tz,
  summary_file = NULL,
  experiment = NULL,
  export_path = NULL,
  plot_title_prefix = NULL,
  generate_raw_data_plot = FALSE,
  chamber_numbers = NULL,
  checkpoint_file = NULL,
  verbose = TRUE
) {
  if (verbose) cat("\n=== RAW DATA PREPROCESSING ===\nFile:", file_path, "\n\n")
  if (!file.exists(file_path)) stop("File not found: ", file_path)

  raw_data <- read.csv(file_path, stringsAsFactors = FALSE)
  if (verbose) cat("[1] Loaded", nrow(raw_data), "rows x", ncol(raw_data), "columns\n")

  col_info <- get_raw_column_info(raw_data)
  input_data <- raw_data
  selected <- select_raw_chambers(raw_data, col_info, chamber_numbers)
  selected_data <- selected$data
  col_info <- selected$col_info

  if (verbose) cat("[2] Filtering to experiment window and measure phases...\n")
  time_series <- prepare_raw_data(
    selected_data, col_info,
    experiment_start_time, experiment_end_time,
    tz = tz
  )
  if (verbose) {
    cat(
      "    Before start:", nrow(time_series$acclimation),
      "| All phases:", nrow(time_series$experimental),
      "| Measure only:", nrow(time_series$measure),
      "| After end:", nrow(time_series$post_experiment), "\n"
    )
  }

  if (generate_raw_data_plot && !is.null(export_path)) {
    if (verbose) cat("[3] Plotting filtered data...\n")

    title_all <- paste(c(plot_title_prefix, "All Phases"), collapse = " - ")
    title_measure <- paste(c(plot_title_prefix, "Measure Phases"), collapse = " - ")

    plot_raw_instrument_data(
      time_series$experimental,
      get_raw_column_info(time_series$experimental),
      title_all,
      paste0(export_path, "_all_phases_plots.pdf")
    )

    plot_raw_instrument_data(
      time_series$measure,
      get_raw_column_info(time_series$measure),
      title_measure,
      paste0(export_path, "_measure_plots.pdf")
    )
  }

  # Regression stats for each measure phase and chamber
  if (verbose) cat("[4] Calculating measure-phase regressions...\n")
  measure_stats <- calculate_measure_stats(
    time_series$measure,
    get_raw_column_info(time_series$measure)
  )
  if (!is.null(summary_file)) { # <<< NEW
    measure_stats <- add_mo2(measure_stats, summary_file, experiment)
  }

  if (!is.null(export_path)) {
    write_csv(measure_stats, paste0(export_path, "_measure_stats.csv"))
    if (verbose) cat("    Exported:", paste0(export_path, "_measure_stats.csv"), "\n")
  }

  # Export only the measure-phase data
  if (!is.null(export_path)) {
    dir.create(dirname(export_path), recursive = TRUE, showWarnings = FALSE)
    write_csv(time_series$measure, paste0(export_path, "_measure.csv"))
    if (verbose) cat("[5] Exported:", paste0(export_path, "_measure.csv"), "\n")
  }

  if (verbose) cat("✓ Complete\n")

  results <- list(
    raw_data = input_data,
    time_series = time_series,
    measure_stats = measure_stats,
    col_info = col_info,
    part1_results = NULL
  )
  if (!is.null(summary_file)) {
    results$part1_results <- build_raw_part1_results(
      raw_data = raw_data,
      time_series = time_series,
      measure_stats = measure_stats,
      col_info = get_raw_column_info(time_series$measure),
      summary_file = summary_file,
      experiment = experiment,
      chamber_numbers = selected$chamber_numbers,
      file_path = file_path,
      experiment_start_time = experiment_start_time,
      experiment_end_time = experiment_end_time
    )
  }
  if (!is.null(checkpoint_file)) {
    if (is.null(results$part1_results)) {
      stop("checkpoint_file requires summary_file to build Part 2 input.")
    }
    save_raw_checkpoint(results$part1_results, checkpoint_file)
  }
  results
}