#' Validate required columns are present
#'
#' Checks for required fields. Concentration input requires `chamber`,
#' `chamber_notes`, `www_g`, `tww_g`, `resp_vol_ml`, `ch_vol_ml`, `datetime`,
#' `phase`, and `o2`. Rate input requires `chamber`, `chamber_notes`, `www_g`,
#' `tww_g`, `resp_vol_ml`, `ch_vol_ml`, `datetime`, `time_sec`, `mo2`,
#' `slope`, and `r2`.
#'
#' @param col_info A list returned by [locate_columns()].
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

  invisible(NULL)
}

#' Locate data columns
#'
#' Identifies input columns for the metadata, timing, phase, temperature, and
#' chamber-specific measurement fields recognized by the pipeline.
#'
#' @param data A data frame whose column names will be searched.
#' @param data_type Either `"concentration"` or `"rate"`. Determines which
#'   columns are required.
#' @returns A list with two named components. `idx` maps each located field to
#'   its column position or positions in `data`; `names` maps each field to
#'   the corresponding column name or names.
#' @export
locate_columns <- function(data, data_type = c("concentration", "rate")) {
  data_type <- match.arg(data_type)
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

  if (nm %in% names(patterns)[1:9] && length(matches) > 1) {
      stop(
        "Multiple columns matched ", nm, ": ",
        paste(names(data)[matches], collapse = ", "),
        call. = FALSE
      )
    }

    if (length(matches) == 0) next
    result$idx[[nm]] <- matches
    result$names[[nm]] <- names(data)[matches]
  }

  validate_columns(result, data_type)
  result
}

#' Validate required columns are present
#'
#' Checks for required fields. Concentration input requires `chamber`,
#' `chamber_notes`, `www_g`, `tww_g`, `resp_vol_ml`, `ch_vol_ml`, `datetime`,
#' `phase`, and `o2`. Rate input requires `chamber`, `chamber_notes`, `www_g`,
#' `tww_g`, `resp_vol_ml`, `ch_vol_ml`, `datetime`, `time_sec`, `mo2`,
#' `slope`, and `r2`.
#'
#' @param col_info A list returned by [locate_columns()].
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

  invisible(NULL)
}

#' Select a subset of chambers
#'
#' Use only when analyzing fewer than all available chambers. Selected
#' measurement columns retain their original order.
#'
#' @param col_info Result of [locate_columns()].
#' @param data_type Either `"concentration"` or `"rate"`.
#' @param chamber_numbers Chamber numbers to include.
#' @returns A list containing `col_info`, selected `chamber_numbers`, and
#'   selected measurement `columns`.
#' @export
select_chamber_subset <- function(col_info, data_type, chamber_numbers) {
  data_type <- match.arg(data_type, c("concentration", "rate"))
  fields <- switch(
    data_type,
    concentration = "o2",
    rate = c("mo2", "slope", "r2")
  )

  selected <- lapply(fields, function(field) {
    column_names <- col_info$names[[field]]
    numbers <- as.integer(stringr::str_match(
      column_names,
      stringr::regex("ch[.\\s]*(\\d+)", ignore_case = TRUE)
    )[, 2])

    missing <- setdiff(chamber_numbers, numbers)
    if (length(missing) > 0) {
      stop(
        "Missing ", field, " columns for chambers: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }

    keep <- numbers %in% chamber_numbers
    list(numbers = numbers[keep], columns = column_names[keep])
  })

  selected_numbers <- selected[[1]]$numbers
  if (any(!vapply(selected, function(x) {
    identical(x$numbers, selected_numbers)
  }, logical(1)))) {
    stop("Selected measurement columns have inconsistent chamber order.",
         call. = FALSE)
  }

  list(
    col_info = col_info,
    chamber_numbers = selected_numbers,
    columns = setNames(lapply(selected, `[[`, "columns"), fields)
  )
}