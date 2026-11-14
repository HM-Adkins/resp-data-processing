#!/usr/bin/env Rscript


# 1.1 Extract descriptive metadata
extract_descriptive_data <- function(
  raw_data,
  col_info,
  experiment,
  id_delimiter = "-",
  verbose = TRUE
) {
  if (length(experiment) != 1L || is.na(experiment) ||
    !grepl("^[0-9]+$", as.character(experiment))) {
    stop("experiment must be a single non-negative integer.")
  }
  experiment <- as.integer(experiment)
  if (is.na(experiment)) {
    stop("experiment is outside the supported integer range.")
  }
  if (length(id_delimiter) != 1L || is.na(id_delimiter) ||
    !nzchar(id_delimiter)) {
    stop("id_delimiter must be a non-empty string.")
  }

  chamber_notes_col <- col_info$names$chamber_notes
  descriptive_cols <- unlist(
    col_info$names[c(
      "chamber", "www_g", "resp_vol_ml", "chamber_notes",
      "tww_g", "ch_vol_ml"
    )],
    use.names = FALSE
  )
  desc <- raw_data %>%
    select(all_of(descriptive_cols)) %>%
    filter(!is.na(.data[[chamber_notes_col]]), .data[[chamber_notes_col]] != "")
  if (nrow(desc) == 0L) stop("No chamber metadata rows were found.")
  rename_map <- c(
    "chamber" = col_info$names$chamber,
    "www_g" = col_info$names$www_g,
    "resp_vol_ml" = col_info$names$resp_vol_ml,
    "chamber_notes" = col_info$names$chamber_notes,
    "tww_g" = col_info$names$tww_g,
    "ch_vol_ml" = col_info$names$ch_vol_ml
  )

  desc <- desc %>%
    rename(!!!rename_map)

  notes <- str_split_fixed(as.character(desc$chamber_notes),
    fixed(id_delimiter),
    n = 2
  )
  species <- str_trim(notes[, 1])
  id <- str_trim(notes[, 2])
  invalid_notes <- which(species == "" | id == "")
  if (length(invalid_notes) > 0L) {
    stop(
      "Could not parse species and ID from chamber notes in row(s): ",
      paste(invalid_notes, collapse = ", ")
    )
  }

  chamber <- as.integer(str_extract(as.character(desc$chamber), "\\d+"))
  if (anyNA(chamber)) {
    stop(
      "Invalid chamber number(s) in descriptive data row(s): ",
      paste(which(is.na(chamber)), collapse = ", ")
    )
  }
  if (anyDuplicated(chamber)) {
    stop("Descriptive data must contain exactly one row per chamber.")
  }

  numeric_fields <- c("www_g", "tww_g", "resp_vol_ml", "ch_vol_ml")
  for (field in numeric_fields) {
    values <- as.character(desc[[field]])
    values[trimws(values) == ""] <- NA_character_
    converted <- suppressWarnings(as.numeric(values))
    invalid <- which(!is.na(values) & is.na(converted))
    if (length(invalid) > 0L) {
      stop(
        "Non-numeric ", field, " value(s) in descriptive row(s): ",
        paste(invalid, collapse = ", ")
      )
    }
    desc[[field]] <- converted
  }

  desc <- desc %>%
    select(-chamber_notes) %>%
    transmute(
      experiment = .env$experiment,
      chamber = .env$chamber,
      species = .env$species,
      id = .env$id,
      www_g = www_g,
      tww_g = tww_g,
      resp_vol_ml = resp_vol_ml,
      ch_vol_ml = ch_vol_ml,
      www_kg = www_g / 1000,
      tww_kg = tww_g / 1000
    )

  if (verbose) {
    message(sprintf(
      "Extracted %d chambers (%d experimental, %d control)",
      nrow(desc),
      sum(desc$id != "Control", na.rm = TRUE),
      sum(desc$id == "Control", na.rm = TRUE)
    ))
  }

  desc
}


# 1.2 Prepare time series
prepare_time_series <- function(
  data,
  col_info,
  experiment_start_time,
  tz,
  date_format = "%m/%d/%Y/%I:%M:%S %p"
) {
  if (length(tz) != 1L || is.na(tz) || !tz %in% OlsonNames()) {
    stop("tz must be a valid time-zone name.")
  }

  dt_col <- col_info$names$datetime
  ts_col <- col_info$names$time_sec
  start_time <- as.POSIXct(
    experiment_start_time,
    format = "%Y-%m-%d %H:%M:%S",
    tz = tz
  )
  if (length(start_time) != 1L || is.na(start_time)) {
    stop("Invalid experiment_start_time: ", experiment_start_time)
  }

  parsed_datetime <- suppressWarnings(
    as.POSIXct(data[[dt_col]], format = date_format, tz = tz)
  )
  invalid_datetime <- which(is.na(parsed_datetime))
  if (length(invalid_datetime) > 0L) {
    stop(
      "Could not parse datetime value(s) in row(s): ",
      paste(invalid_datetime, collapse = ", ")
    )
  }

  time_seconds <- suppressWarnings(as.numeric(as.character(data[[ts_col]])))
  invalid_time <- which(!is.finite(time_seconds))
  if (length(invalid_time) > 0L) {
    stop(
      "Invalid time-in-seconds value(s) in row(s): ",
      paste(invalid_time, collapse = ", ")
    )
  }

  data[[dt_col]] <- parsed_datetime
  data$time_hour <- time_seconds / 3600

  exp_start_idx <- which(data[[dt_col]] >= start_time)[1]
  if (length(exp_start_idx) == 0L || is.na(exp_start_idx)) {
    stop(
      "No data found at or after experiment start time: ",
      experiment_start_time
    )
  }

  offset <- data$time_hour[exp_start_idx]
  data <- mutate(data, time_hour = time_hour - offset)

  list(
    all = data,
    acclimation = filter(data, .data[[dt_col]] < start_time),
    experimental = filter(data, .data[[dt_col]] >= start_time)
  )
}


# 1.3 Build list of chamber-specific data sets
build_chamber_list <- function(resp_data, desc_data, col_info) {
  n_chambers <- nrow(desc_data)
  time_vec <- resp_data$time_hour

  chamber_numbers <- desc_data$chamber
  if (anyNA(chamber_numbers) || anyDuplicated(chamber_numbers)) {
    stop("Descriptive data must contain one valid row per chamber.")
  }

  mo2_cols <- match_chamber_columns(col_info$names$mo2, "MO2", chamber_numbers)
  r2_cols <- match_chamber_columns(col_info$names$r2, "R2", chamber_numbers)
  slope_cols <- match_chamber_columns(
    col_info$names$slope, "slope", chamber_numbers
  )

  # Calculate average temperature
  temp_cols <- col_info$names$temp[!is.na(col_info$names$temp)]
  avg_temp <- if (length(temp_cols) > 0) {
    rowMeans(resp_data[, temp_cols, drop = FALSE], na.rm = TRUE)
  } else {
    rep(NA_real_, nrow(resp_data))
  }
  avg_temp[!is.finite(avg_temp)] <- NA_real_

  chamber_list <- lapply(seq_len(n_chambers), function(i) {
    chamber_data <- tibble(
      experiment = desc_data$experiment[i],
      chamber = desc_data$chamber[i],
      species = desc_data$species[i],
      id = desc_data$id[i],
      time = time_vec,
      avg_temp = avg_temp,
      r2 = resp_data[[r2_cols[i]]],
      slope = resp_data[[slope_cols[i]]],
      raw_mo2_kg_www = resp_data[[mo2_cols[i]]]
    ) %>%
      mutate(raw_mo2_ch = raw_mo2_kg_www * desc_data$www_kg[i])

    if (desc_data$id[i] == "Control") {
      chamber_data <- chamber_data %>%
        mutate(
          raw_mo2_ch = raw_mo2_kg_www,
          raw_mo2_kg_www = NA_real_
        )
    }

    chamber_data
  })

  chamber_list
}


# 1.4 Raw MO2 data plot
plot_raw_data <- function(
  data,
  desc_data,
  col_info,
  title_prefix = NULL,
  output_file = NULL
) {
  cfg <- resp_config
  mo2_cols <- match_chamber_columns(
    col_info$names$mo2, "MO2", desc_data$chamber
  )

  plot_device(output_file, 10, 8, {
    exp_idx <- which(desc_data$id != "Control")
    n <- length(exp_idx)
    nr <- ceiling(sqrt(n))
    nc <- ceiling(n / nr)
    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1), oma = c(1, 1, 3, 0))
    time_zero <- min(data$time_hour, na.rm = TRUE)

    # Calculate global y-axis range across all experimental chambers
    all_mo2 <- unlist(lapply(exp_idx, function(idx) {
      data[[mo2_cols[idx]]]
    }))
    y_max <- max(all_mo2, na.rm = TRUE)
    y_min <- min(all_mo2, na.rm = TRUE)

    # Set y_lim: start at 0 unless minimum is below 0
    y_lim <- c(ifelse(y_min < 0, y_min, 0), y_max)
    x_label <- c(rep("", 4), rep("Time (hrs)", 2))

    for (i in seq_along(exp_idx)) {
      idx <- exp_idx[i]
      mo2 <- data[[mo2_cols[idx]]]

      plot(data$time_hour, mo2,
        main = paste(desc_data$species[idx], "-", desc_data$id[idx]),
        xlab = x_label[i], ylab = "MO2 (mg O2/kg WWW/hr)",
        ylim = y_lim,
        pch = 16, cex = cfg$plot$point_size
      )

      # Shade acclimation
      usr <- par("usr") # Get actual plot region coordinates
      polygon(c(time_zero, 0, 0, time_zero),
        c(usr[3], usr[3], usr[4], usr[4]),
        col = rgb(0.4, 0.4, 0.4, 0.2), border = NA
      )
      abline(v = 0, col = cfg$colors$blue, lty = 2)
      abline(h = 0, col = cfg$colors$red, lty = 2)
    }

    if (!is.null(title_prefix)) {
      mtext(title_prefix, outer = TRUE, font = 2, cex = 1.2)
    }
  })

  if (!is.null(output_file)) message("Raw plot: ", output_file)
}


# Plot raw control MO2 separately because it is not mass-specific.
plot_raw_control_data <- function(
  data,
  desc_data,
  col_info,
  title_prefix = NULL,
  output_file = NULL
) {
  control_idx <- which(desc_data$id == "Control")
  if (length(control_idx) == 0L) {
    stop("No control chamber found for the raw control plot.")
  }
  control_numbers <- desc_data$chamber[control_idx]
  mo2_cols <- match_chamber_columns(
    col_info$names$mo2, "MO2", control_numbers
  )
  if (!is.numeric(data$time_hour) || !any(is.finite(data$time_hour))) {
    stop("Cannot plot controls without finite time values.")
  }
  values <- lapply(mo2_cols, function(name) {
    mo2 <- data[[name]]
    if (!is.numeric(mo2) ||
      !any(is.finite(mo2) & is.finite(data$time_hour))) {
      stop("Cannot plot control column ", name, ": no finite time and MO2 pairs.")
    }
    mo2
  })
  finite_values <- unlist(lapply(values, function(mo2) {
    mo2[is.finite(mo2) & is.finite(data$time_hour)]
  }), use.names = FALSE)
  limits <- range(c(0, finite_values[is.finite(finite_values)]))
  if (limits[1] == limits[2]) {
    limits <- limits + c(-1, 1)
  }
  cfg <- resp_config

  plot_device(output_file, 10, 8, {
    n <- length(control_idx)
    nr <- ceiling(sqrt(n))
    nc <- ceiling(n / nr)
    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1), oma = c(1, 1, 3, 0))
    time_zero <- min(data$time_hour[is.finite(data$time_hour)])

    for (i in seq_along(control_idx)) {
      plot(
        data$time_hour, values[[i]],
        main = paste("Control - Chamber", control_numbers[i]),
        xlab = if (i > (nr - 1) * nc) "Time (hrs)" else "",
        ylab = "MO2 (mg O2/chamber/hr)",
        ylim = limits,
        pch = 16, cex = cfg$plot$point_size
      )
      usr <- par("usr")
      polygon(c(time_zero, 0, 0, time_zero),
        c(usr[3], usr[3], usr[4], usr[4]),
        col = rgb(0.4, 0.4, 0.4, 0.2), border = NA
      )
      abline(v = 0, col = cfg$colors$blue, lty = 2)
      abline(h = 0, col = cfg$colors$red, lty = 2)
    }
    if (!is.null(title_prefix)) {
      mtext(paste(title_prefix, "- Raw Controls"),
        outer = TRUE,
        font = 2, cex = 1.2
      )
    }
  })

  if (!is.null(output_file)) message("Raw control plot: ", output_file)
}


# 1.5 RUN PART 1 FUNCTIONS AND PLOT RAW DATA
process_part1_data_loading <- function(
  file_path,
  experiment_start_time,
  experiment,
  tz,
  date_format = "%m/%d/%Y/%I:%M:%S %p",
  id_delimiter = "-",
  plot_title_prefix = NULL,
  generate_raw_data_plot = FALSE,
  export_path = NULL,
  verbose = TRUE
) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)

  if (verbose) cat("\n=== PART 1: DATA LOADING & PREPROCESSING ===\n\n")

  if (verbose) cat("Loading data...\n")
  raw_data <- read.csv(file_path, stringsAsFactors = FALSE)
  if (verbose) {
    cat(
      "      ", nrow(raw_data), "rows x", ncol(raw_data), "columns\n\n"
    )
  }

  col_info <- get_column_info(raw_data)
  validate_columns(col_info, c(
    "chamber", "www_g", "resp_vol_ml", "chamber_notes", "tww_g", "ch_vol_ml",
    "datetime", "time_sec", "mo2", "r2"
  ))

  if (verbose) cat("[1.1] Extracting metadata...\n")
  desc_data <- extract_descriptive_data(
    raw_data, col_info, experiment, id_delimiter, verbose
  )

  if (verbose) cat("[1.2] Processing time series...\n")
  resp_data <- raw_data %>%
    select(all_of(names(raw_data)[col_info$idx$datetime:ncol(raw_data)]))

  col_info_resp <- get_column_info(resp_data)
  time_series <- prepare_time_series(
    resp_data, col_info_resp, experiment_start_time, tz, date_format
  )
  if (verbose) {
    cat(
      "      Acclimation:", nrow(time_series$acclimation),
      "| Experimental:", nrow(time_series$experimental), "\n\n"
    )
  }

  if (verbose) cat("[1.3] Building list of chamber datasets...\n")
  # Update col_info_resp to include the newly created time_hour column
  col_info_resp <- get_column_info(time_series$experimental)
  chamber_list <- build_chamber_list(
    time_series$experimental, desc_data,
    col_info_resp
  )
  if (verbose) cat("      ", length(chamber_list), "chambers ready\n\n")

  if (generate_raw_data_plot && !is.null(export_path)) {
    if (verbose) cat("[1.4] Plotting raw data...\n")
    col_info_all <- get_column_info(time_series$all)
    plot_raw_data(
      time_series$all, desc_data, col_info_all, plot_title_prefix,
      paste0(export_path, "_raw_data_plots.pdf")
    )
    if (any(desc_data$id == "Control")) {
      plot_raw_control_data(
        time_series$all, desc_data, col_info_all, plot_title_prefix,
        paste0(export_path, "_raw_control_plots.pdf")
      )
    } else if (verbose) {
      message("No control chambers found; skipping raw control plot.")
    }
  }

  if (verbose) cat("✓ Part 1 Complete\n")

  list(
    raw_data = raw_data,
    descriptive_data = desc_data,
    time_series = time_series,
    chamber_list = chamber_list,
    col_info = lapply(col_info_resp, function(fields) {
      fields[setdiff(
        names(fields),
        c(
          "chamber", "www_g", "resp_vol_ml",
          "chamber_notes", "tww_g", "ch_vol_ml"
        )
      )]
    }),
    metadata = list(
      species = unique(desc_data$species[desc_data$id != "Control"]),
      file_path = file_path,
      experiment = as.integer(experiment),
      experiment_start = experiment_start_time,
      tz = tz,
      n_chambers = nrow(desc_data),
      n_controls = sum(desc_data$id == "Control"),
      n_experimental = sum(desc_data$id != "Control")
    )
  )
}
