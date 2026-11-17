#' Locate specific columns
#'
#' Searches the input data's column names for the fields specified in the
#' patterns list. Matching is not case-sensitive. Fields without matches
#' are omitted. Fields with multiple matches return every match in input-column
#' order.
#'
#' @param data A data frame whose column names will be searched.
#' @returns A list with two named components. `idx` maps each located field to
#'   its column position or positions in `data`; `names` maps each field to
#'   the corresponding column name or names.
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

#' Validate located columns
#'
#' Checks for required fields and rejects multiple matches for fields that
#' should identify a single column. Multiple columns are permitted for
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
#'   observations or `"rate"` for respiration-rate data.
#' @returns Invisibly returns `NULL` when validation succeeds. Otherwise,
#'   gives an error.
validate_columns <- function(col_info, data_type) {
  common <- c(
    "chamber", "chamber_notes", "www_g", "tww_g",
    "resp_vol_ml", "ch_vol_ml", "datetime"
  )

  required <- c(
    common,
    switch(
      data_type,
      concentration = c("phase", "o2"),
      rate = c("time_sec", "mo2", "slope", "r2")
    )
  )

  missing <- setdiff(required, names(col_info$idx))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  multiple_allowed <- c("temp", "o2", "mo2", "slope", "r2")
  multiple <- names(col_info$idx)[lengths(col_info$idx) > 1]
  ambiguous <- setdiff(multiple, multiple_allowed)

  if (length(ambiguous) > 0) {
    details <- vapply(ambiguous, function(field) {
      paste0(field, ": ", paste(col_info$names[[field]], collapse = ", "))
    }, character(1))

    stop("Multiple matches found: ", paste(details, collapse = "; "),
         call. = FALSE)
  }

  invisible(NULL)
}

#' Order columns by chamber
#'
#' Extracts chamber numbers from column names and returns the names in the
#' requested chamber order.
#'
#' @param column_names Column names containing chamber identifiers.
#' @param field Measurement field name, used in error messages.
#' @param chamber_numbers Chamber numbers in the desired order.
#' @returns A character vector of column names ordered by chamber.
order_columns_by_chamber <- function(column_names, field, chamber_numbers) {
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
}

#' Locate and match input columns to chambers
#'
#' Locates and validates input columns, then matches measurement columns to
#' chambers. Ambient temperature columns remain in `col_info` and are not
#' matched to chambers.
#'
#' @param data A concentration or rate data frame.
#' @param data_type Either `"concentration"` or `"rate"`.
#' @param chamber_numbers Unique positive whole-number chamber identifiers
#'   in the desired output order. Defaults to `1:8`.
#' @returns A list containing `col_info`, `chamber_numbers`, and `columns`.
#'   The `columns` list contains `o2` for concentration data or `mo2`,
#'   `slope`, and `r2` for rate data.
#' @export
match_chamber_columns <- function(data,
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
    order_columns_by_chamber(col_info$names[[field]], field, chamber_numbers)
  })
  names(columns) <- fields

  list(
    col_info = col_info,
    chamber_numbers = as.integer(chamber_numbers),
    columns = columns
  )
}