#' Extract descriptive and time-series data
#'
#' Separates descriptive columns from the remaining input columns.
#'
#' @param data Original input data frame.
#' @param located Result of `locate$match_chamber_columns()`.
#' @returns A list with `descriptive` and `time_series` data frames. The
#'   descriptive columns are removed from `time_series`.
#' @export
extract_descriptive_data <- function(data, located) {
  fields <- c(
    "chamber", "chamber_notes", "www_g", "tww_g",
    "resp_vol_ml", "ch_vol_ml"
  )
  columns <- unlist(located$col_info$names[fields], use.names = FALSE)

  list(
    descriptive = data[, columns, drop = FALSE],
    time_series = data[, !names(data) %in% columns, drop = FALSE]
  )
}


#' Format descriptive chamber data
#'
#' Parses chamber numbers, species and IDs; converts measurements to numeric;
#' calculates weights in kilograms; and orders rows by chamber.
#'
#' @param data Descriptive data frame from [extract_descriptive_data()].
#' @param located Result of `locate$match_chamber_columns()`.
#' @param experiment Experiment identifier.
#' @param id_delimiter Separator between species and ID in chamber notes.
#' @returns A formatted data frame with one row per selected chamber.
#' @export
format_descriptive_data <- function(data, located, experiment,
                                    id_delimiter = "-") {
  fields <- c(
    "chamber", "chamber_notes", "www_g", "tww_g",
    "resp_vol_ml", "ch_vol_ml"
  )
  names(data) <- fields

  notes <- as.character(data$chamber_notes)
  data <- data[!is.na(notes) & nzchar(trimws(notes)), , drop = FALSE]
  if (nrow(data) == 0L) {
    stop("No chamber metadata rows found.", call. = FALSE)
  }

  chamber_text <- as.character(data$chamber)
  chamber <- suppressWarnings(as.integer(
    regmatches(chamber_text, regexpr("[0-9]+", chamber_text))
  ))
  if (anyNA(chamber) || anyDuplicated(chamber)) {
    stop("Metadata must identify each chamber exactly once.", call. = FALSE)
  }

  notes <- strsplit(as.character(data$chamber_notes), id_delimiter, fixed = TRUE)
  valid_notes <- vapply(notes, function(x) {
    length(x) >= 2L &&
      nzchar(trimws(x[1])) &&
      nzchar(trimws(paste(x[-1], collapse = id_delimiter)))
  }, logical(1))
  if (any(!valid_notes)) {
    stop("Chamber notes must contain species and ID.", call. = FALSE)
  }

  species <- vapply(notes, function(x) trimws(x[1]), character(1))
  id <- vapply(notes, function(x) {
    trimws(paste(x[-1], collapse = id_delimiter))
  }, character(1))

  numeric_fields <- c("www_g", "tww_g", "resp_vol_ml", "ch_vol_ml")
  for (field in numeric_fields) {
    values <- as.character(data[[field]])
    values[!is.na(values) & !nzchar(trimws(values))] <- NA_character_
    converted <- suppressWarnings(as.numeric(values))

    if (any(!is.na(values) & !is.finite(converted))) {
      stop("Invalid values in ", field, ".", call. = FALSE)
    }
    data[[field]] <- converted
  }

  result <- data.frame(
    experiment = experiment,
    chamber = chamber,
    species = species,
    id = id,
    www_g = data$www_g,
    tww_g = data$tww_g,
    resp_vol_ml = data$resp_vol_ml,
    ch_vol_ml = data$ch_vol_ml,
    www_kg = data$www_g / 1000,
    tww_kg = data$tww_g / 1000
  )

  positions <- match(located$chamber_numbers, result$chamber)
  if (anyNA(positions)) {
    stop(
      "Missing metadata for chamber(s): ",
      paste(located$chamber_numbers[is.na(positions)], collapse = ", "),
      call. = FALSE
    )
  }

  result[positions, , drop = FALSE]
}


#' Format the input time series
#'
#' Parses the located date-time column in the specified time zone.
#'
#' @param data Time-series data frame, such as `time_series` from
#'   [extract_descriptive_data()].
#' @param located Result of `locate$match_chamber_columns()`.
#' @param tz An Olson time-zone name, such as `"UTC"`.
#' @param date_format Format of the input date-time strings, passed to
#'   `as.POSIXct()`.
#' @returns The input data frame with its date-time column parsed as POSIXct.
#' @export
format_time_series <- function(data, located, tz, date_format) {
  column <- located$col_info$names$datetime
  parsed <- suppressWarnings(as.POSIXct(
    as.character(data[[column]]),
    format = date_format,
    tz = tz
  ))

  if (anyNA(parsed)) {
    stop(
      "Could not parse date-time in row(s): ",
      paste(which(is.na(parsed)), collapse = ", "),
      call. = FALSE
    )
  }

  data[[column]] <- parsed
  data
}


#' Extract the experimental time window
#'
#' Selects observations from start through end, inclusive. For concentration
#' data, excludes rows before the first F1 and F/W/M loops with fewer than
#' 300 M-phase observations in the window. F/W phases after the last complete
#' M phase are also excluded. Rate data retain all rows in the window.
#'
#' @param data Data returned by [format_time_series()].
#' @param located Result of `locate$match_chamber_columns()`.
#' @param data_type Either `"concentration"` or `"rate"`.
#' @param start_time,end_time Boundaries as `"YYYY-MM-DD HH:MM:SS"` strings.
#' @param tz Olson time-zone name used to interpret the boundaries.
#' @returns Selected rows with `time_sec` and `time_hour` elapsed-time columns.
#' @export
extract_experimental_data <- function(data, located, data_type, start_time,
                                      end_time, tz) {
  data_type <- match.arg(data_type, c("concentration", "rate"))

  start <- as.POSIXct(start_time, format = "%Y-%m-%d %H:%M:%S", tz = tz)
  end <- as.POSIXct(end_time, format = "%Y-%m-%d %H:%M:%S", tz = tz)
  if (length(start) != 1L || length(end) != 1L ||
      is.na(start) || is.na(end) || start >= end) {
    stop("Experiment boundaries must be valid and start before end.",
         call. = FALSE)
  }

  datetime <- data[[located$col_info$names$datetime]]
  if (!inherits(datetime, "POSIXct") || anyNA(datetime)) {
    stop("Format date-time values before selecting the experiment.",
         call. = FALSE)
  }

  experimental <- data[datetime >= start & datetime <= end, , drop = FALSE]
  if (nrow(experimental) == 0L) {
    stop("No observations found in the experiment window.", call. = FALSE)
  }

  if (data_type == "concentration") {
    phase_col <- located$col_info$names$phase
    phases <- toupper(trimws(as.character(experimental[[phase_col]])))
    f1_rows <- which(!is.na(phases) & phases == "F1")
    if (length(f1_rows) == 0L) {
      stop("No F1 phase found in the experiment window.", call. = FALSE)
    }

    origin <- experimental[[located$col_info$names$datetime]][f1_rows[1]]
    experimental <- experimental[
      seq.int(f1_rows[1], nrow(experimental)), , drop = FALSE
    ]
    phases <- toupper(trimws(as.character(experimental[[phase_col]])))

    loop_rows <- grepl("^[FWM][0-9]+$", phases)
    phase_numbers <- suppressWarnings(
      as.integer(sub("^[FWM]", "", phases))
    )
    loops <- unique(phase_numbers[loop_rows])
    m_counts <- vapply(loops, function(loop_number) {
      sum(phases == paste0("M", loop_number), na.rm = TRUE)
    }, integer(1))

    incomplete_loops <- loops[m_counts < 300L]
    keep <- !(loop_rows & phase_numbers %in% incomplete_loops)
    experimental <- experimental[keep, , drop = FALSE]

    phases <- toupper(trimws(as.character(experimental[[phase_col]])))
    m_rows <- grepl("^M[0-9]+$", phases)
    if (!any(m_rows)) {
      stop("No M phases with at least 300 observations remain.",
           call. = FALSE)
    }

    last_complete_m <- max(as.integer(sub("^M", "", phases[m_rows])))
    fw_rows <- grepl("^[FW][0-9]+$", phases)
    fw_numbers <- suppressWarnings(
      as.integer(sub("^[FW]", "", phases[fw_rows]))
    )
    keep <- rep(TRUE, nrow(experimental))
    keep[fw_rows] <- fw_numbers <= last_complete_m
    experimental <- experimental[keep, , drop = FALSE]

    dates <- experimental[[located$col_info$names$datetime]]
    seconds <- as.numeric(difftime(dates, origin, units = "secs"))
  } else {
    seconds <- suppressWarnings(as.numeric(as.character(
      experimental[[located$col_info$names$time_sec]]
    )))
    if (any(!is.finite(seconds))) {
      stop("Rate input contains invalid time-in-seconds values.",
           call. = FALSE)
    }
    seconds <- seconds - seconds[1]
  }

  experimental$time_sec <- seconds
  experimental$time_hour <- seconds / 3600
  experimental
}


#' Build a chamber list
#'
#' Combines each chamber's metadata and matched measurements with experimental
#' timing and mean ambient temperature. Control MO2 is retained per chamber.
#'
#' @param data Data returned by [extract_experimental_data()].
#' @param descriptive_data Data returned by [format_descriptive_data()].
#' @param located Result of `locate$match_chamber_columns()`.
#' @param data_type Either `"concentration"` or `"rate"`.
#' @returns A named list of chamber data frames.
#' @export
build_chamber_list <- function(data, descriptive_data, located, data_type) {
  data_type <- match.arg(data_type, c("concentration", "rate"))
  fields <- if (data_type == "concentration") {
    "o2"
  } else {
    c("mo2", "slope", "r2")
  }

  if (!all(fields %in% names(located$columns))) {
    stop("Located columns do not match data_type.", call. = FALSE)
  }

  positions <- match(located$chamber_numbers, descriptive_data$chamber)
  if (anyNA(positions) || anyDuplicated(descriptive_data$chamber)) {
    stop("Descriptive data must contain exactly one row for each chamber.",
         call. = FALSE)
  }

  temps <- located$col_info$names$temp
  if (is.null(temps) || length(temps) == 0L) {
    avg_temp <- rep(NA_real_, nrow(data))
  } else {
    if (!all(vapply(data[temps], is.numeric, logical(1)))) {
      stop("Ambient temperature columns must be numeric.", call. = FALSE)
    }
    avg_temp <- rowMeans(data[, temps, drop = FALSE], na.rm = TRUE)
    avg_temp[is.nan(avg_temp)] <- NA_real_
  }

  chambers <- lapply(seq_along(located$chamber_numbers), function(i) {
    meta <- descriptive_data[positions[i], , drop = FALSE]

    result <- data.frame(
      experiment = meta$experiment,
      chamber = meta$chamber,
      species = meta$species,
      id = meta$id,
      datetime = data[[located$col_info$names$datetime]],
      time_sec = data$time_sec,
      time_hour = data$time_hour,
      time = data$time_hour,
      avg_temp = avg_temp
    )

    if (data_type == "concentration") {
      result$phase <- data[[located$col_info$names$phase]]
      result$o2 <- data[[located$columns$o2[i]]]
    } else {
      result$r2 <- data[[located$columns$r2[i]]]
      result$slope <- data[[located$columns$slope[i]]]
      mo2 <- data[[located$columns$mo2[i]]]

      if (!is.numeric(mo2)) {
        stop("MO2 columns must be numeric.", call. = FALSE)
      }

      if (meta$id == "Control") {
        result$raw_mo2_kg_www <- NA_real_
        result$raw_mo2_ch <- mo2
      } else {
        result$raw_mo2_kg_www <- mo2
        result$raw_mo2_ch <- mo2 * meta$www_kg
      }
    }

    result
  })

  names(chambers) <- paste0("CH", located$chamber_numbers)
  chambers
}