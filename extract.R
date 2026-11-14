#' Extract chamber metadata
#'
#' Extracts and standardizes metadata for the selected chambers. Chamber notes
#' must contain species and ID separated by `id_delimiter`. Missing weights
#' and volumes are retained as NA, including for control chambers.
#'
#' @param data Original input data frame.
#' @param located Result of `locate$match_chamber_columns()` on `data`.
#' @param experiment A single non-negative integer experiment identifier.
#' @param id_delimiter Separator between species and ID in chamber notes.
#' @returns A data frame ordered by selected chamber with experiment, chamber,
#'   species, ID, weights in grams and kilograms, and volumes in millilitres.
#' @export
extract_descriptive_data <- function(data, located, experiment,
                                     id_delimiter = "-") {
  if (length(experiment) != 1L || !is.numeric(experiment) ||
      is.na(experiment) || !is.finite(experiment) || experiment < 0 ||
      experiment != floor(experiment) || experiment > .Machine$integer.max) {
    stop("experiment must be a non-negative whole number.", call. = FALSE)
  }

  if (!is.character(id_delimiter) || length(id_delimiter) != 1L ||
      is.na(id_delimiter) || !nzchar(id_delimiter)) {
    stop("id_delimiter must be a non-empty string.", call. = FALSE)
  }

  fields <- c(
    "chamber", "chamber_notes", "www_g", "tww_g",
    "resp_vol_ml", "ch_vol_ml"
  )
  column_names <- unlist(located$col_info$names[fields], use.names = FALSE)
  desc <- data[, column_names, drop = FALSE]
  names(desc) <- fields

  notes <- as.character(desc$chamber_notes)
  desc <- desc[!is.na(notes) & nzchar(trimws(notes)), , drop = FALSE]

  if (nrow(desc) == 0L) {
    stop("No chamber metadata rows found.", call. = FALSE)
  }

  chamber <- as.character(desc$chamber)
  digits <- regmatches(chamber, regexpr("[0-9]+", chamber))

  if (length(digits) != nrow(desc)) {
    stop("Some metadata rows have invalid chamber numbers.", call. = FALSE)
  }

  chamber <- suppressWarnings(as.integer(digits))
  if (anyNA(chamber) || any(chamber < 1L) || anyDuplicated(chamber)) {
    stop("Metadata must identify each chamber exactly once.", call. = FALSE)
  }

  notes <- strsplit(as.character(desc$chamber_notes), id_delimiter, fixed = TRUE)
  valid_notes <- vapply(notes, function(x) {
    length(x) >= 2L && nzchar(trimws(x[1])) &&
      nzchar(trimws(paste(x[-1], collapse = id_delimiter)))
  }, logical(1))

  if (any(!valid_notes)) {
    stop(
      "Could not parse species and ID from chamber notes in metadata row(s): ",
      paste(which(!valid_notes), collapse = ", "),
      call. = FALSE
    )
  }

  species <- vapply(notes, function(x) trimws(x[1]), character(1))
  id <- vapply(notes, function(x) {
    trimws(paste(x[-1], collapse = id_delimiter))
  }, character(1))

  for (field in c("www_g", "tww_g", "resp_vol_ml", "ch_vol_ml")) {
    values <- as.character(desc[[field]])
    values[!is.na(values) & !nzchar(trimws(values))] <- NA_character_
    converted <- suppressWarnings(as.numeric(values))
    invalid <- !is.na(values) & !is.finite(converted)

    if (any(invalid)) {
      stop(
        "Invalid ", field, " in metadata row(s): ",
        paste(which(invalid), collapse = ", "),
        call. = FALSE
      )
    }

    desc[[field]] <- converted
  }

  result <- data.frame(
    experiment = as.integer(experiment),
    chamber = chamber,
    species = species,
    id = id,
    www_g = desc$www_g,
    tww_g = desc$tww_g,
    resp_vol_ml = desc$resp_vol_ml,
    ch_vol_ml = desc$ch_vol_ml,
    www_kg = desc$www_g / 1000,
    tww_kg = desc$tww_g / 1000
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
#' Parses the located date-time column in the specified time zone. The
#' original column is replaced with POSIXct values; all other columns remain.
#'
#' @param data Original input data frame.
#' @param located Result of `locate$match_chamber_columns()` on `data`.
#' @param tz An Olson time-zone name, such as `"UTC"`.
#' @param date_format Format of the input date-time strings, passed to
#'   `as.POSIXct()`. For example, `"%d-%m-%Y/%H:%M:%S"`.
#' @returns The input data frame with its date-time column parsed.
#' @export
format_time_series <- function(data, located, tz, date_format) {
  if (!is.character(tz) || length(tz) != 1L || is.na(tz) ||
      !tz %in% OlsonNames()) {
    stop("tz must be a valid time-zone name.", call. = FALSE)
  }

  if (!is.character(date_format) || length(date_format) != 1L ||
      is.na(date_format) || !nzchar(date_format)) {
    stop("date_format must be a non-empty format string.", call. = FALSE)
  }

  column <- located$col_info$names$datetime
  parsed <- suppressWarnings(as.POSIXct(
    as.character(data[[column]]), format = date_format, tz = tz
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
#' data, rows before the first F1 phase are excluded. F/W/M loops whose M phase
#' contains fewer than 300 observations within the window are excluded, as are
#' F/W phases after the last complete M phase. Rate data retain every row in
#' the time window. Concentration elapsed time starts at the first retained F1
#' timestamp; rate elapsed time starts at the first in-window observation.
#'
#' @param data Data returned by [format_time_series()].
#' @param located Result of `locate$match_chamber_columns()` on the input.
#' @param data_type Either `"concentration"` or `"rate"`.
#' @param start_time,end_time Experiment boundaries as
#'   `"YYYY-MM-DD HH:MM:SS"` strings interpreted in `tz`.
#' @param tz The Olson time zone used by [format_time_series()].
#' @returns The selected rows with `time_sec` and `time_hour` elapsed-time
#'   columns.
#' @export
extract_experimental_data <- function(data, located, data_type, start_time,
                                      end_time, tz) {
  data_type <- match.arg(data_type, c("concentration", "rate"))
  if (!is.character(tz) || length(tz) != 1L || is.na(tz) ||
      !tz %in% OlsonNames()) {
    stop("tz must be a valid time-zone name.", call. = FALSE)
  }

  parse_boundary <- function(value) {
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
      stop("Experiment boundaries must be date-time strings.", call. = FALSE)
    }
    suppressWarnings(as.POSIXct(
      value, format = "%Y-%m-%d %H:%M:%S", tz = tz
    ))
  }
  start <- parse_boundary(start_time)
  end <- parse_boundary(end_time)
  if (is.na(start) || is.na(end) || start >= end) {
    stop(
      "Experiment boundaries must be valid and start before end.",
      call. = FALSE
    )
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

    # Discard all rows before the first F1.
    experimental <- experimental[seq.int(f1_rows[1], nrow(experimental)), ,
                                  drop = FALSE]
    phases <- toupper(trimws(as.character(experimental[[phase_col]])))

    # Remove loops whose M phase has fewer than 300 observations.
    fwm_rows <- grepl("^[FWM][0-9]+$", phases)
    phase_numbers <- suppressWarnings(
      as.integer(sub("^[FWM]", "", phases))
    )
    loops <- unique(phase_numbers[fwm_rows])
    m_counts <- vapply(loops, function(loop_number) {
      sum(phases == paste0("M", loop_number), na.rm = TRUE)
    }, integer(1))

    incomplete_loops <- loops[m_counts < 300L]
    keep <- !(fwm_rows & phase_numbers %in% incomplete_loops)
    experimental <- experimental[keep, , drop = FALSE]

    phases <- toupper(trimws(as.character(experimental[[phase_col]])))
    m_rows <- grepl("^M[0-9]+$", phases)
    if (!any(m_rows)) {
      stop("No M phases with at least 300 observations remain.",
           call. = FALSE)
    }

    # Exclude F/W phases after the last complete M phase.
    last_complete_m <- max(as.integer(sub("^M", "", phases[m_rows])))
    fw_rows <- grepl("^[FW][0-9]+$", phases)
    fw_numbers <- suppressWarnings(
      as.integer(sub("^[FW]", "", phases[fw_rows]))
    )
    keep <- rep(TRUE, nrow(experimental))
    keep[fw_rows] <- fw_numbers <= last_complete_m
    experimental <- experimental[keep, , drop = FALSE]

    phases <- toupper(trimws(as.character(experimental[[phase_col]])))
    f1_rows <- which(!is.na(phases) & phases == "F1")
    dates <- experimental[[located$col_info$names$datetime]]
    seconds <- as.numeric(difftime(
      dates, min(dates[f1_rows]), units = "secs"
    ))
  } else {
    values <- experimental[[located$col_info$names$time_sec]]
    seconds <- suppressWarnings(as.numeric(as.character(values)))
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
#' Combines each selected chamber's metadata and matched measurements with the
#' shared experimental timing and mean ambient temperature. Concentration
#' chambers contain `o2` and `phase`; rate chambers contain `r2`, `slope`,
#' `raw_mo2_kg_www`, and `raw_mo2_ch`. Control chamber MO2 is retained as a
#' per-chamber value rather than multiplied by body weight.
#'
#' @param data Data returned by [extract_experimental_data()].
#' @param descriptive_data Data returned by [extract_descriptive_data()].
#' @param located Result of `locate$match_chamber_columns()`.
#' @param data_type Either `"concentration"` or `"rate"`.
#' @returns A named list in selected chamber order. Each entry is a data frame
#'   of metadata, elapsed time, average temperature, and measurements.
#' @export
build_chamber_list <- function(data, descriptive_data, located, data_type) {
  data_type <- match.arg(data_type, c("concentration", "rate"))
  fields <- if (data_type == "concentration") "o2" else c("mo2", "slope", "r2")

  if (!all(fields %in% names(located$columns))) {
    stop("Located columns do not match data_type.", call. = FALSE)
  }

  positions <- match(located$chamber_numbers, descriptive_data$chamber)
  if (anyNA(positions) || anyDuplicated(descriptive_data$chamber)) {
    stop("Descriptive data must contain exactly one row for each chamber.",
         call. = FALSE)
  }

  temps <- located$col_info$names$temp

  if (length(temps) == 0L) {
    avg_temp <- rep(NA_real_, nrow(data))
  } else {
    if (!all(vapply(data[temps], is.numeric, logical(1)))) {
      stop("Ambient temperature columns must be numeric.", call. = FALSE)
    }

    if (any(vapply(data[temps], function(x) any(is.infinite(x)), logical(1)))) {
      stop("Ambient temperature contains infinite values.", call. = FALSE)
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