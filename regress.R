#' Prepare regression points
#'
#' Collects M-phase oxygen measurements from the chamber list. Optionally
#' averages consecutive observations into bins within each chamber and phase.
#' Incomplete bins and observations with missing or non-finite values are
#' excluded with a warning.
#'
#' @param chambers Concentration chamber list from `extract$build_chamber_list()`.
#' @param bin_size Number of observations per bin, or `NULL` for unbinned data.
#' @returns A data frame of regression points by chamber and phase.
#' @export
prepare_regression_points <- function(chambers, bin_size = NULL) {
  if (!is.null(bin_size) &&
      (length(bin_size) != 1 || !is.numeric(bin_size) ||
       is.na(bin_size) || !is.finite(bin_size) ||
       bin_size < 1 || bin_size != floor(bin_size))) {
    stop("bin_size must be NULL or a positive whole number.", call. = FALSE)
  }

  points <- dplyr::bind_rows(lapply(chambers, function(chamber_data) {
    if (!all(vapply(
      chamber_data[c("o2", "time_sec", "time_hour")],
      is.numeric, logical(1)
    ))) {
      stop("O2 and time columns must be numeric.", call. = FALSE)
    }

    dplyr::transmute(
      chamber_data,
      chamber = .data$chamber,
      phase = toupper(trimws(as.character(.data$phase))),
      time_sec = .data$time_sec,
      time_hour = .data$time_hour,
      o2 = .data$o2
    )
  }))

  points <- dplyr::filter(points, grepl("^M[0-9]+$", .data$phase))
  if (nrow(points) == 0) {
    stop("No measurement phases found.", call. = FALSE)
  }

  losses <- points |>
    dplyr::group_by(.data$chamber, .data$phase) |>
    dplyr::summarise(
      excluded = sum(
        !is.finite(.data$o2) |
          !is.finite(.data$time_sec) |
          !is.finite(.data$time_hour)
      ),
      valid = dplyr::n() - .data$excluded,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      discarded = if (is.null(bin_size)) 0 else .data$valid %% bin_size
    ) |>
    dplyr::filter(.data$excluded > 0 | .data$discarded > 0)

  if (nrow(losses) > 0) {
    details <- paste0(
      "chamber ", losses$chamber, " phase ", losses$phase, ": ",
      losses$excluded, " invalid observation(s), ",
      losses$discarded, " observation(s) in incomplete bins"
    )
    warning(paste(details, collapse = "; "), call. = FALSE)
  }

  points <- points |>
    dplyr::filter(
      is.finite(.data$o2),
      is.finite(.data$time_sec),
      is.finite(.data$time_hour)
    ) |>
    dplyr::arrange(.data$chamber, .data$phase, .data$time_sec)

  if (!is.null(bin_size)) {
    points <- points |>
      dplyr::group_by(.data$chamber, .data$phase) |>
      dplyr::mutate(bin = (dplyr::row_number() - 1) %/% bin_size + 1) |>
      dplyr::group_by(.data$chamber, .data$phase, .data$bin) |>
      dplyr::filter(dplyr::n() == bin_size) |>
      dplyr::summarise(
        time_sec = mean(.data$time_sec),
        time_hour = mean(.data$time_hour),
        o2 = mean(.data$o2),
        .groups = "drop"
      )
  }

  points
}


#' Summarise regression points
#'
#' Calculates regression statistics for each chamber and M phase. Phases
#' without usable points are retained with `n_points = 0` and missing
#' statistics. At least three points and two distinct times are required
#' for a regression.
#'
#' @param points Data returned by [prepare_regression_points()].
#' @param chambers Original concentration chamber list.
#' @returns A data frame of chamber-phase regression statistics.
#' @export
summarise_regression_points <- function(points, chambers) {
  phase_times <- dplyr::bind_rows(lapply(chambers, function(chamber_data) {
    dplyr::transmute(
      chamber_data,
      chamber = .data$chamber,
      phase = toupper(trimws(as.character(.data$phase))),
      time_hour = .data$time_hour,
      avg_temp = .data$avg_temp
    )
  })) |>
    dplyr::filter(grepl("^M[0-9]+$", .data$phase)) |>
    dplyr::group_by(.data$chamber, .data$phase) |>
    dplyr::summarise(
      time = if (any(is.finite(.data$time_hour))) {
        mean(.data$time_hour[is.finite(.data$time_hour)])
      } else {
        NA_real_
      },
      avg_temp = if (any(is.finite(.data$avg_temp))) {
        mean(.data$avg_temp[is.finite(.data$avg_temp)])
      } else {
        NA_real_
      },
      .groups = "drop"
    )

  stats <- points |>
    dplyr::group_by(.data$chamber, .data$phase) |>
    dplyr::summarise(
      n_points = dplyr::n(),
      start_time_sec = dplyr::first(.data$time_sec),
      end_time_sec = dplyr::last(.data$time_sec),
      start_time_hr = dplyr::first(.data$time_hour),
      end_time_hr = dplyr::last(.data$time_hour),
      r2 = if (dplyr::n() >= 3 && dplyr::n_distinct(.data$time_hour) >= 2) {
        summary(stats::lm(o2 ~ time_hour))$r.squared
      } else {
        NA_real_
      },
      slope_per_hr = if (
        dplyr::n() >= 3 && dplyr::n_distinct(.data$time_hour) >= 2
      ) {
        stats::coef(stats::lm(o2 ~ time_hour))[[2]]
      } else {
        NA_real_
      },
      mean_o2 = mean(.data$o2),
      first_o2 = dplyr::first(.data$o2),
      last_o2 = dplyr::last(.data$o2),
      .groups = "drop"
    )

  phase_times |>
    dplyr::left_join(stats, by = c("chamber", "phase")) |>
    dplyr::mutate(
      n_points = dplyr::coalesce(.data$n_points, 0L),
      loop = as.integer(sub("^M", "", .data$phase))
    ) |>
    dplyr::arrange(.data$loop, .data$chamber) |>
    dplyr::select(-.data$loop)
}


#' Add oxygen consumption rates
#'
#' Uses regression slopes and chamber metadata to calculate oxygen consumption.
#' Experimental chambers use respiration volume and whole wet weight; controls
#' use chamber volume and retain a per-chamber rate.
#'
#' @param measure_stats Data returned by [summarise_regression_points()].
#' @param descriptive_data Data returned by `extract$format_descriptive_data()`.
#' @returns Regression statistics with `phase`, `loop`, and MO2 columns.
#' @export
add_mo2 <- function(measure_stats, descriptive_data) {
  if (anyDuplicated(descriptive_data$chamber)) {
    stop("Descriptive data contains duplicate chambers.", call. = FALSE)
  }

  missing <- setdiff(unique(measure_stats$chamber), descriptive_data$chamber)
  if (length(missing) > 0) {
    stop(
      "Missing metadata for chamber(s): ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  metadata <- dplyr::select(
    descriptive_data,
    .data$chamber, .data$id, .data$www_kg,
    .data$resp_vol_ml, .data$ch_vol_ml
  )

  result <- dplyr::left_join(measure_stats, metadata, by = "chamber")
  control <- result$id == "Control"

  if (anyNA(control) ||
      any(!control & (!is.finite(result$www_kg) | result$www_kg <= 0)) ||
      any(!control & (!is.finite(result$resp_vol_ml) | result$resp_vol_ml <= 0)) ||
      any(control & (!is.finite(result$ch_vol_ml) | result$ch_vol_ml <= 0))) {
    stop("Missing or invalid weight or volume in chamber metadata.",
         call. = FALSE)
  }

  result |>
    dplyr::mutate(
      loop = as.integer(sub("^M", "", .data$phase)),
      phase = "M",
      raw_mo2_kg_www = ifelse(
        control, NA_real_,
        -.data$slope_per_hr * (.data$resp_vol_ml / 1000) / .data$www_kg
      ),
      raw_mo2_ch = ifelse(
        control,
        abs(.data$slope_per_hr * .data$ch_vol_ml / 1000),
        .data$raw_mo2_kg_www * .data$www_kg
      )
    ) |>
    dplyr::relocate(.data$phase, .data$loop)
}


#' Build a regression chamber list
#'
#' Organizes regression results into one data frame per chamber for subsequent
#' correction and filtering steps.
#'
#' @param measure_stats Data returned by [add_mo2()].
#' @param descriptive_data Data returned by `extract$format_descriptive_data()`.
#' @returns A named list of chamber regression data frames.
#' @export
build_regression_chamber_list <- function(measure_stats, descriptive_data) {
  chambers <- lapply(seq_len(nrow(descriptive_data)), function(i) {
    meta <- descriptive_data[i, , drop = FALSE]

    measure_stats |>
      dplyr::filter(.data$chamber == meta$chamber) |>
      dplyr::transmute(
        experiment = meta$experiment,
        chamber = .data$chamber,
        species = meta$species,
        id = meta$id,
        time = .data$time,
        avg_temp = .data$avg_temp,
        r2 = .data$r2,
        slope = .data$slope_per_hr,
        raw_mo2_kg_www = .data$raw_mo2_kg_www,
        raw_mo2_ch = .data$raw_mo2_ch
      )
  })

  names(chambers) <- paste0("CH", descriptive_data$chamber)
  chambers
}


#' Process concentration regressions
#'
#' Prepares points, calculates regression statistics and MO2, and builds the
#' chamber list used by subsequent pipeline steps.
#'
#' @param chambers Concentration chamber list.
#' @param descriptive_data Formatted descriptive chamber data.
#' @param bin_size Number of observations per bin, or `NULL` for unbinned data.
#' @returns A list containing `regression_points`, `measure_stats`,
#'   `chamber_list`, and `bin_size`.
#' @export
process_regressions <- function(chambers, descriptive_data, bin_size = NULL) {
  points <- prepare_regression_points(chambers, bin_size)
  stats <- summarise_regression_points(points, chambers)
  stats <- add_mo2(stats, descriptive_data)

  list(
    regression_points = points,
    measure_stats = stats,
    chamber_list = build_regression_chamber_list(stats, descriptive_data),
    bin_size = bin_size
  )
}


#' Compare regression bin sizes
#'
#' Runs unbinned regressions and regressions for each requested bin size.
#'
#' @param chambers Concentration chamber list.
#' @param descriptive_data Formatted descriptive chamber data.
#' @param bin_sizes Positive whole-number bin sizes to compare.
#' @returns A list containing full `variants` and a `comparison` data frame.
#' @export
compare_regressions <- function(chambers, descriptive_data, bin_sizes) {
  if (!is.numeric(bin_sizes) || length(bin_sizes) == 0 ||
      anyNA(bin_sizes) || any(!is.finite(bin_sizes)) ||
      any(bin_sizes < 1 | bin_sizes != floor(bin_sizes)) ||
      anyDuplicated(bin_sizes)) {
    stop("bin_sizes must contain unique positive whole numbers.",
         call. = FALSE)
  }

  sizes <- c(list(NULL), as.list(bin_sizes))
  labels <- c("unbinned", paste0("bin_", bin_sizes))

  variants <- stats::setNames(lapply(sizes, function(size) {
    process_regressions(chambers, descriptive_data, bin_size = size)
  }), labels)

  comparison <- dplyr::bind_rows(lapply(seq_along(variants), function(i) {
    stats <- variants[[i]]$measure_stats

    dplyr::transmute(
      stats,
      variant = labels[i],
      bin_size = if (is.null(sizes[[i]])) NA_integer_ else sizes[[i]],
      chamber = .data$chamber,
      phase = paste0(.data$phase, .data$loop),
      n_points = .data$n_points,
      slope_per_hr = .data$slope_per_hr,
      r2 = .data$r2
    )
  }))

  list(variants = variants, comparison = comparison)
}