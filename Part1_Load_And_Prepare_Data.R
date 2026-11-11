#------------------------------------------------------------------------------
# PART 1: DATA LOADING AND PREPROCESSING
# extract_descriptive_data, prepare_time_series, build_chamber_list,
# plot_raw_data, process_part1_data_loading
#------------------------------------------------------------------------------

# 1. 1 Extract descriptive metadata
extract_descriptive_data <- function(
  raw_data, col_info, id_delimiter = "-", plot_title_prefix = NULL
) {
  chamber_notes_col <- col_info$names$chamber_notes
  desc <- raw_data %>%
    select(1:(col_info$idx$datetime - 1)) %>%
    filter(!is.na(.data[[chamber_notes_col]]), .data[[chamber_notes_col]] != "")

  # Rename using col_info mappings
  rename_map <- c(
    "chamber" = col_info$names$chamber,
    "www_g" = col_info$names$www_g,
    "resp_vol_ml" = col_info$names$resp_vol_ml,
    "chamber_notes" = col_info$names$chamber_notes,
    "wtw_g" = col_info$names$wtw_g,
    "ch_vol_ml" = col_info$names$ch_vol_ml
  )

  desc <- desc %>%
    rename(!!!rename_map)

  # Parse species-id from notes
  parsed <- str_match(
    desc$chamber_notes,
    paste0("^\\s*(.+?)\\s*", id_delimiter, "\\s*(.+?)\\s*$")
  )

  # Restructure descriptive data and calculate additional fields
  desc <- desc %>%
    select(-chamber_notes) %>%
    transmute(
      experiment = if (!is.null(plot_title_prefix)) {
        as.integer(str_extract(plot_title_prefix, "(?i)\\d+"))
      } else {
        NA_integer_
      },
      chamber = as.integer(str_extract(chamber, "\\d+")),
      species = parsed[, 2],
      id = parsed[, 3],
      www_g = as.numeric(www_g),
      wtw_g = as.numeric(wtw_g),
      resp_vol_ml = as.numeric(resp_vol_ml),
      ch_vol_ml = as.numeric(ch_vol_ml),
      www_kg = www_g / 1000,
      wtw_kg = wtw_g / 1000
    )

  message(sprintf(
    "Extracted %d chambers (%d experimental, %d control)",
    nrow(desc),
    sum(desc$id != "Control", na.rm = TRUE),
    sum(desc$id == "Control", na.rm = TRUE)
  ))

  desc
}

# 1.2 Prepare time series data
prepare_time_series <- function(data, col_info, experiment_start_time, date_format = "%m/%d/%Y/%I:%M:%S %p") {
  dt_col <- col_info$names$datetime
  ts_col <- col_info$names$time_sec

  start_time <- as.POSIXct(experiment_start_time)

  # Select only datetime column onwards (no descriptive columns)
  data <- data %>%
    select(col_info$idx$datetime:ncol(data))

  # Format datetime and time_hour columns
  data <- data %>%
    mutate(
      !!sym(dt_col) := as.POSIXct(.data[[dt_col]], format = date_format, tz = Sys.timezone()),
      time_hour = .data[[ts_col]] / 3600
    )

  # Scale time_hour to experiment start
  exp_start_idx <- which(data[[dt_col]] >= start_time)[1]

  if (is.na(exp_start_idx)) {
    stop("No data found at or after experiment start time: ", experiment_start_time)
  }

  offset <- data$time_hour[exp_start_idx]

  data <- mutate(data, time_hour = time_hour - offset)

  # Split into acclimation and experimental periods
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

  # Calculate average temperature
  temp_cols <- col_info$names$temp[!is.na(col_info$names$temp)]
  avg_temp <- if (length(temp_cols) > 0) {
    rowMeans(resp_data[, temp_cols, drop = FALSE], na.rm = TRUE)
  } else {
    rep(NA_real_, nrow(resp_data))
  }

  chamber_list <- lapply(seq_len(n_chambers), function(i) {
    chamber_data <- tibble(
      experiment = desc_data$experiment[i],
      chamber = desc_data$chamber[i],
      species = desc_data$species[i],
      id = desc_data$id[i],
      time = time_vec,
      avg_temp = avg_temp,
      r2 = resp_data[[col_info$names$r2[i]]],
      raw_mo2_kg_www = resp_data[[col_info$names$mo2[i]]]
    ) %>%
      mutate(raw_mo2_ch = raw_mo2_kg_www * desc_data$www_kg[i])

    # Handle control chambers differently
    if (desc_data$id[i] == "Control") {
      # Add slope column for this control chamber (use same index i)
      chamber_data$slope <- resp_data[[col_info$names$slope[i]]]

      # For control: set raw_mo2_ch = raw_mo2_kg_www, then set raw_mo2_kg_www = NA
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

# 1.4 Raw data plot
plot_raw_data <- function(data, desc_data, col_info, title_prefix = NULL, output_file = NULL) {
  cfg <- RESPIRATION_CONFIG

  plot_device(output_file, 10, 8, {
    exp_idx <- which(desc_data$id != "Control")
    n <- length(exp_idx)

    nr <- ceiling(sqrt(n))
    nc <- ceiling(n / nr)
    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1), oma = c(1, 1, 3, 0))

    time_zero <- min(data$time_hour, na.rm = TRUE)

    # Calculate global y-axis range across all experimental chambers
    all_mo2 <- unlist(lapply(exp_idx, function(idx) {
      data[[col_info$names$mo2[idx]]]
    }))
    y_max <- max(all_mo2, na.rm = TRUE)
    y_min <- min(all_mo2, na.rm = TRUE)

    # Set y_lim: start at 0 unless minimum is below 0
    y_lim <- c(ifelse(y_min < 0, y_min, 0), y_max)
    x_label <- c(rep("", 4), rep("Time (hrs)", 2))

    for (i in seq_along(exp_idx)) {
      idx <- exp_idx[i]
      mo2 <- data[[col_info$names$mo2[idx]]]

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

# RUN PART 1 FUNCTIONS AND PLOT RAW DATA: 1.5 Load and preprocess
process_part1_data_loading <- function(
  file_path,
  experiment_start_time,
  id_delimiter = "-",
  plot_title_prefix = NULL,
  generate_raw_data_plot = FALSE,
  export_path = NULL,
  verbose = TRUE
) {
  if (verbose) cat("\n=== PART 1: DATA LOADING & PREPROCESSING ===\n\n")

  # Load
  if (verbose) cat("[1. 1] Loading data...\n")
  raw_data <- read.csv(file_path, stringsAsFactors = FALSE)
  if (verbose) cat("      ", nrow(raw_data), "rows ×", ncol(raw_data), "columns\n\n")

  # Get column info
  col_info <- get_column_info(raw_data)

  # Validate required columns are present (descriptive + measurement)
  validate_columns(col_info, c(
    "chamber", "www_g", "resp_vol_ml", "chamber_notes", "wtw_g", "ch_vol_ml",
    "datetime", "time_sec", "mo2", "r2"
  ))

  # Extract metadata
  if (verbose) cat("[1.2] Extracting metadata...\n")
  desc_data <- extract_descriptive_data(raw_data, col_info, id_delimiter, plot_title_prefix)

  # Prepare respiration data
  if (verbose) cat("[1.3] Processing time series...\n")
  resp_data <- raw_data %>%
    select(-c(1:(col_info$idx$datetime - 1)))

  col_info_resp <- get_column_info(resp_data)
  time_series <- prepare_time_series(resp_data, col_info_resp, experiment_start_time)
  if (verbose) cat("      Acclimation:", nrow(time_series$acclimation), "| Experimental:", nrow(time_series$experimental), "\n\n")

  # Build chamber structures
  if (verbose) cat("[1.4] Building list of chamber datasets...\n")
  # Update col_info_resp to include the newly created time_hour column
  col_info_resp <- get_column_info(time_series$experimental)
  chamber_list <- build_chamber_list(time_series$experimental, desc_data, col_info_resp)
  if (verbose) cat("      ", length(chamber_list), "chambers ready\n\n")

  # Plot raw data
  if (generate_raw_data_plot && !is.null(export_path)) {
    if (verbose) cat("[1.5] Plotting raw data...\n")
    col_info_all <- get_column_info(time_series$all)
    plot_raw_data(
      time_series$all, desc_data, col_info_all, plot_title_prefix,
      paste0(export_path, "_raw_data_plots.pdf")
    )
  }

  if (verbose) cat("✓ Part 1 Complete\n")

  list(
    raw_data = raw_data,
    descriptive_data = desc_data,
    time_series = time_series,
    chamber_list = chamber_list,
    col_info = col_info_resp,
    metadata = list(
      species = unique(desc_data$species[desc_data$id != "Control"]),
      file_path = file_path,
      experiment_start = experiment_start_time,
      n_chambers = nrow(desc_data),
      n_controls = sum(desc_data$id == "Control"),
      n_experimental = sum(desc_data$id != "Control")
    )
  )
}
