#'
#' Locate recognized data columns
#'
#' Searches the input data's column names for the metadata, timing, oxygen
#' concentration, respiration-rate, regression, and temperature fields
#' recognized by the pipeline. Matching is case-insensitive. Fields with no
#' matching columns are omitted; fields with multiple matches return every
#' match in input-column order.
#'
#' @param data A data frame whose column names will be searched.
#' @returns A list with two named components. `idx` maps each located field to
#'   its column position or positions in `data`; `names` maps each field to
#'   the corresponding column name or names.
#' @examples
#' data <- data.frame(`CH1 O2 [mgO2/L]` = 1, check.names = FALSE)
#' get_column_info(data)$names$o2
#' @export
get_column_info <- function(data) {
  patterns <- list(
    chamber = "^chamber$|^parameter.*chamber|^parameter.*channel",
    chamber_notes = "note",
    www_g = "^(?!.*tissue)(?:(?=.*wet)(?=.*weight)|.*www)",
    tww_g = "wet.*tissue|tww|wtw",
    resp_vol_ml = "resp.*vol",
    ch_vol_ml = "chamber.*vol",
    datetime = "date.*time",
    time_sec = "time.*sec",
    phase = "^phase$",
    temp = "^amb[.\\s]+temp[.\\s]+n[^a-z0-9]*[12][.\\s]+(?:avg[.\\s]+temp[.\\s]+)?[^a-z0-9]*c[^a-z0-9]*$",
    o2 = "ch[.\\s]*\\d+[.\\s]+o2.*mgo2",
    mo2 = "ch.*mo2",
    slope = "ch.*slope",
    r2 = "ch.*r.*2"
  )

  result <- list(idx = list(), names = list())

  for (nm in names(patterns)) {
    matches <- grep(
      patterns[[nm]], names(data),
      ignore.case = TRUE, perl = TRUE
    )
    if (length(matches) == 0L) next
    result$idx[[nm]] <- matches
    result$names[[nm]] <- names(data)[matches]
  }
  result
}

#'
#' Validate located columns
#'
#' Checks that a located-column result contains all fields required for the
#' selected input type and that fields expected to identify one input column
#' do not have multiple matches. Multiple columns are permitted for
#' chamber-specific or trough-specific fields (`temp`, `o2`, `mo2`, `slope`,
#' and `r2`).
#'
#' Concentration input requires `chamber`, `chamber_notes`, `www_g`, `tww_g`,
#' `resp_vol_ml`, `ch_vol_ml`, `datetime`, `phase`, and `o2`. Rate input
#' requires `chamber`, `chamber_notes`, `www_g`, `tww_g`, `resp_vol_ml`,
#' `ch_vol_ml`, `datetime`, `time_sec`, `mo2`, `slope`, and `r2`.
#'
#' @param col_info A list returned by [get_column_info()].
#' @param data_type Input type: `"concentration"` for oxygen-concentration
#'   observations or `"rate"` for precomputed respiration-rate data.
#' @returns Invisibly returns `NULL` when validation succeeds. Otherwise,
#'   throws an error describing missing fields or ambiguous multiple matches.
validate_columns <- function(col_info, data_type) {
  required <- switch(
    data_type,
    concentration = c(
      "chamber", "chamber_notes", "www_g", "tww_g",
      "resp_vol_ml", "ch_vol_ml", "datetime", "phase", "o2"
      ),
    rate = c(
      "chamber", "chamber_notes", "www_g", "tww_g",
      "resp_vol_ml", "ch_vol_ml", "datetime", "time_sec",
      "mo2", "slope", "r2"
    )
  )

  multi_column_fields <- c("temp", "o2", "mo2", "slope", "r2")
  flagged_fields <- names(col_info$idx)[
    lengths(col_info$idx) > 1L &
      !names(col_info$idx) %in% multi_column_fields
  ]

  if (length(flagged_fields) > 0L) {
    details <- vapply(flagged_fields, function(field) {
      paste0(field, ": ", paste(col_info$names[[field]], collapse = ", "))
    }, character(1))

    stop(
      "Expected one column per field. Multiple matches found: ",
      paste(details, collapse = "; "),
      call. = FALSE
    )
  }

  missing <- required[
    vapply(required, function(field) {
      idx <- col_info$idx[[field]]
      is.null(idx) || length(idx) == 0L || all(is.na(idx))
    }, logical(1))
  ]

  if (length(missing) > 0L) {
    stop("Missing required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  invisible(NULL)
}

#'
#' Locate and match input columns to chambers
#'
#' Identifies columns, validates the fields required for the selected input
#' type, and matches chamber-specific measurements to the requested chamber
#' numbers. Matched column names are returned in `chamber_numbers` order,
#' regardless of their order in the input data.
#'
#' Concentration input returns matched `o2` columns. Rate input returns matched
#' `mo2`, `slope`, and `r2` columns. Every requested chamber must have a
#' corresponding column for each matched field. Ambient temperature columns
#' are included in `col_info` but are not chamber-matched here.
#'
#' @param data A data frame containing concentration or rate input.
#' @param data_type Either `"concentration"` (the default) or `"rate"`.
#' @param chamber_numbers Positive, unique whole-number chamber identifiers in
#'   the desired output order. Defaults to `c(1, 2, 3, 4, 5, 6, 7, 8)`.
#' @returns A list containing `col_info` (from [get_column_info()]),
#'   `chamber_numbers` (the requested chamber order as integers), and
#'   `columns` (a named list of matched column-name vectors).
#' @details The function errors if required fields are missing, a field
#'   expected to match one column matches multiple columns, chamber numbers
#'   are invalid or duplicated, or a requested chamber has no matching
#'   measurement column.
#' @examples
#' data <- data.frame(
#'   `Parameter/Chamber` = c("Chamber 1", "Chamber 2"),
#'   `Chamber notes` = c("Fish-A", "Fish-B"),
#'   `Wet weight [g]` = c(100, 120),
#'   `Wet tissue weight [g]` = c(80, 90),
#'   `Resp. volume [mL]` = c(250, 250),
#'   `Chamber vol [ml]` = c(500, 500),
#'   `Date & Time` = c("2023-01-01", "2023-01-01"),
#'   `Phase` = c("M1", "M1"),
#'   `CH2 O2 [mgO2/L]` = c(1, 2),
#'   `CH1 O2 [mgO2/L]` = c(2, 1),
#'   check.names = FALSE
#' )
#' match_chamber_columns(data, "concentration", c(1, 2))$columns$o2
#' @export
match_chamber_columns <- function(
    data,
    data_type = c("concentration", "rate"),
    chamber_numbers = 1:8) {
  data_type <- match.arg(data_type)

  if (!is.numeric(chamber_numbers) || !length(chamber_numbers) ||
      anyNA(chamber_numbers) || any(!is.finite(chamber_numbers)) ||
      any(chamber_numbers < 1 | chamber_numbers != floor(chamber_numbers) |
          chamber_numbers > .Machine$integer.max) ||
      anyDuplicated(chamber_numbers)) {
    stop("chamber_numbers must be unique positive whole numbers.",
         call. = FALSE)
  }

  col_info <- get_column_info(data)
  validate_columns(col_info, data_type)

  fields <- switch(
    data_type,
    concentration = "o2",
    rate = c("mo2", "slope", "r2")
  )

  columns <- lapply(fields, function(field) {
    column_names <- col_info$names[[field]]
    numbers <- stringr::str_match(
      column_names,
      stringr::regex("ch[.\\s]*(\\d+)", ignore_case = TRUE)
    )[, 2]

    numbers <- suppressWarnings(as.integer(numbers))
    if (anyNA(numbers) || anyDuplicated(numbers)) {
      stop("Invalid or duplicate chamber numbers in ", field, " columns.",
           call. = FALSE)
    }

    positions <- match(chamber_numbers, numbers)
    if (anyNA(positions)) {
      stop(
        "Missing ", field, " columns for chambers: ",
        paste(chamber_numbers[is.na(positions)], collapse = ", "),
        call. = FALSE
      )
    }

    column_names[positions]
  })
  names(columns) <- fields

  list(
    col_info = col_info,
    chamber_numbers = as.integer(chamber_numbers),
    columns = columns
  )
}