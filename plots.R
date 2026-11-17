resp_config <- list(
  colors = list(
    raw = "#D55E00",
    corrected = "#0072B2",
    control = "#009E73",
    gray = "gray",
    red = "red",
    blue = "blue"
  ),
  plot = list(
    point_size = 0.75,
    point_cex = 0.7,
    line_width = 2,
    control_line_width = 2
  ),
  thresholds = list(
    default_r2 = 0.9,
    default_control_r2 = 0.85,
    min_good_control_points = 3,
    control_regression_r2 = 0.6,
    high_background_pct = 50
  )
)


# Helper function to create custom configuration
create_custom_config <- function(
  colors = NULL, plot = NULL, thresholds = NULL
) {
  cfg <- resp_config

  # Override colors if provided
  if (!is.null(colors)) {
    for (name in names(colors)) {
      cfg$colors[[name]] <- colors[[name]]
    }
  }

  # Override plot settings if provided
  if (!is.null(plot)) {
    for (name in names(plot)) {
      cfg$plot[[name]] <- plot[[name]]
    }
  }

  # Override thresholds if provided
  if (!is.null(thresholds)) {
    for (name in names(thresholds)) {
      cfg$thresholds[[name]] <- thresholds[[name]]
    }
  }

  cfg
}


# Plot device helper (plots to screen and PDF)
plot_device <- function(output_file, width = 10, height = 8, expr) {
  # Plot to screen
  eval.parent(substitute(expr))

  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    pdf(output_file, width = width, height = height)
    eval.parent(substitute(expr))
    dev.off()
  }

  invisible(NULL)
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
