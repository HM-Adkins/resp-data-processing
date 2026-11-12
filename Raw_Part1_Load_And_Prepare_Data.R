#!/usr/bin/env Rscript


# PART 1: DATA LOADING AND PROCESSING
# get_raw_column_info


library(tidyverse)


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


# Get named column indices from the data
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


process_raw_part1 <- function(
    file_path,
    experiment_start_time,
    experiment_end_time,
    tz,
    chamber_numbers = NULL,
    export_path = NULL,
    plot_title_prefix = NULL,
    generate_raw_data_plot = FALSE,
    verbose = TRUE) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  raw_data <- read.csv(file_path, stringsAsFactors = FALSE)
  selected <- select_raw_chambers(
    raw_data, get_raw_column_info(raw_data), chamber_numbers
  )
  time_series <- prepare_raw_data(
    selected$data, selected$col_info,
    experiment_start_time, experiment_end_time, tz = tz
  )

  if (generate_raw_data_plot && !is.null(export_path)) {
    plot_raw_instrument_data(
      time_series$experimental, selected$col_info,
      paste(c(plot_title_prefix, "All Phases"), collapse = " - "),
      paste0(export_path, "_all_phases_plots.pdf")
    )
    plot_raw_instrument_data(
      time_series$measure, selected$col_info,
      paste(c(plot_title_prefix, "Measure Phases"), collapse = " - "),
      paste0(export_path, "_measure_plots.pdf")
    )
  }
  if (!is.null(export_path)) {
    dir.create(dirname(export_path), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(time_series$measure, paste0(export_path, "_measure.csv"))
  }
  if (verbose) {
    message("Prepared ", nrow(time_series$measure), " measurement rows for ",
            length(selected$chamber_numbers), " chambers.")
  }
  list(
    raw_data = raw_data,
    time_series = time_series,
    col_info = selected$col_info,
    metadata = list(
      file_path = file_path,
      experiment_start = experiment_start_time,
      experiment_end = experiment_end_time,
      tz = tz,
      chambers = selected$chamber_numbers
    )
  )
}