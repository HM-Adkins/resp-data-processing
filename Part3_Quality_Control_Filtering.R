#------------------------------------------------------------------------------
# PART 3: QUALITY CONTROL
# flag_part1_quality, apply_quality_filters, plot_rmr_by_temperature,
# plot_quality_diagnostics, process_part3_quality_control
#------------------------------------------------------------------------------

# 3.1 Flag Part 1 quality
flag_part1_quality <- function(part1_results, r2_threshold, verbose = TRUE) {
  if (!is.numeric(r2_threshold) || length(r2_threshold) != 1L ||
        !is.finite(r2_threshold) || r2_threshold < 0 || r2_threshold > 1) {
    stop("r2_threshold must be a single finite number between 0 and 1.")
  }
  if (!is.list(part1_results)) {
    stop("part1_results must be a list of results from Part 1.")
  }
  chambers <- part1_results$chamber_list
  desc <- part1_results$descriptive_data
  if (!is.list(chambers) || !is.data.frame(desc) ||
        length(chambers) != nrow(desc) || length(chambers) == 0L) {
    stop("Part 1 results must contain a chamber list aligned with descriptive data.")
  }
  required_desc <- c("experiment", "chamber", "species", "id")
  if (!all(required_desc %in% names(desc))) {
    stop("Part 1 descriptive data is missing chamber identification columns.")
  }
  if (anyNA(desc$chamber) || anyNA(desc$id) || anyDuplicated(desc$chamber)) {
    stop("Part 1 descriptive data must identify each chamber exactly once.")
  }

  summaries <- vector("list", length(chambers))
  for (i in seq_along(chambers)) {
    chamber <- chambers[[i]]
    if (!is.data.frame(chamber) ||
          !all(c("chamber", "id", "r2") %in% names(chamber)) ||
          !is.numeric(chamber$r2)) {
      stop("Chamber ", i, " must contain chamber, id, and numeric r2 columns.")
    }
    if (anyNA(chamber$chamber) || anyNA(chamber$id) ||
          any(chamber$chamber != desc$chamber[i]) ||
          any(as.character(chamber$id) != as.character(desc$id[i]))) {
      stop("Chamber ", i, " does not match its descriptive data row.")
    }
    if (any(!is.finite(chamber$r2) & !is.na(chamber$r2))) {
      stop("Chamber ", i, " contains non-finite R2 values.")
    }

    chamber$r2_quality <- ifelse(
      is.na(chamber$r2), "missing_r2",
      ifelse(chamber$r2 < r2_threshold, "low_r2", "pass")
    )
    chambers[[i]] <- chamber
    counts <- table(factor(
      chamber$r2_quality,
      levels = c("pass", "low_r2", "missing_r2")
    ))
    n_raw <- nrow(chamber)
    summaries[[i]] <- tibble(
      experiment = desc$experiment[i],
      chamber = desc$chamber[i],
      species = desc$species[i],
      id = desc$id[i],
      n_raw = .env$n_raw,
      n_pass_r2 = unname(counts[1]),
      n_low_r2 = unname(counts[2]),
      n_missing_r2 = unname(counts[3]),
      perc_low_r2 = round(100 * n_low_r2 / n_raw, 2)
    )
  }
  part1_results$chamber_list <- chambers
  summary <- bind_rows(summaries)
  if (verbose) {
    cat("\n=== Part 1 R2 quality summary (rows retained) ===\n")
    print(summary, n = Inf)
  }
  list(part1_results = part1_results,
       summary = summary,
       r2_threshold = r2_threshold)
}


# 3.1 Apply filters
apply_quality_filters <- function(chamber_list, desc_data, r2_threshold, verbose = TRUE) {
  # Combine experimental chambers
  exp_idx <- which(desc_data$id != "Control")
  combined <- bind_rows(chamber_list[exp_idx]) %>%
    mutate(
      species = fct_inorder(as.factor(species)),
      id = fct_inorder(as.factor(id)),
      filtered_mo2_kg_www = case_when(
        is.na(r2) | r2 < r2_threshold ~ NA_real_,
        is.na(corrected_mo2_ch) | corrected_mo2_ch < 0 ~ NA_real_,
        TRUE ~ corrected_mo2_kg_www
      ),
      exclusion_reason = case_when(
        !is.na(filtered_mo2_kg_www) ~ "included",
        is.na(r2) ~ "missing_r2",
        r2 < r2_threshold ~ "low_r2",
        corrected_mo2_ch < 0 ~ "negative_mo2",
        TRUE ~ "other"
      )
    )

  # Calculate stats
  sample_sizes <- combined %>%
    group_by(experiment, chamber, species, id) %>%
    summarise(
      n_raw = n(),
      n_final = sum(!is.na(filtered_mo2_kg_www)),
      n_r2_low = sum(r2 < r2_threshold, na.rm = TRUE),
      n_negative = sum(corrected_mo2_ch < 0, na.rm = TRUE),
      mean_r2 = mean(r2, na.rm = TRUE),
      mean_bg_pct = mean(background_pct, na.rm = TRUE),
      pct_retained = round(100 * n_final / n_raw, 1),
      .groups = "drop"
    )

  exclusion_summary <- count(combined, exclusion_reason) %>% arrange(desc(n))

  if (verbose) {
    cat(sprintf(
      "\n=== Quality Summary ===\nTotal: %d | Valid: %d | Retention: %.1f%%\n\n",
      nrow(combined),
      sum(!is.na(combined$filtered_mo2_kg_www)),
      100 * mean(!is.na(combined$filtered_mo2_kg_www))
    ))
    print(sample_sizes, n = Inf)
  }

  list(
    data = select(combined, experiment, chamber, species, id, time, avg_temp, r2,
      rmr = filtered_mo2_kg_www, correction_flag, exclusion_reason
    ),
    sample_sizes = sample_sizes,
    exclusion_summary = exclusion_summary,
    raw_corrected = combined
  )
}

# 3.2 RMR plot
plot_rmr_by_temperature <- function(data, title, output_file = NULL) {
  cfg <- resp_config

  plot_data <- data %>%
    filter(!is.na(rmr)) %>%
    mutate(panel_id = paste(species, id, sep = " - "))

  # Calculate y-axis range: start at 0 unless minimum is below 0
  y_max <- max(plot_data$rmr, na.rm = TRUE)
  y_min <- min(plot_data$rmr, na.rm = TRUE)
  y_lim <- c(ifelse(y_min < 0, y_min, 0), y_max)

  p <- xyplot(rmr ~ avg_temp | panel_id,
    data = plot_data,
    main = title,
    xlab = "Temperature (°C)",
    ylab = "RMR (mg O2/kg WWW/hr)",
    ylim = y_lim,
    panel = function(x, y, ...) {
      panel.xyplot(x, y, pch = 16, cex = cfg$plot$point_cex, ...)
      if (length(x) > 5) panel.loess(x, y, lwd = cfg$plot$line_width)
      lims <- current.panel.limits()
      panel.text(lims$xlim[2], lims$ylim[2],
        paste0("n=", length(y)),
        adj = c(1, 1), offset = 0.5
      )
    }
  )

  # Print to screen first
  print(p)

  # Then save to PDF if output_file is provided
  if (!is.null(output_file)) {
    pdf(output_file, width = 12, height = 8)
    print(p)
    dev.off()
    message("RMR plot: ", output_file)
  }

  invisible(NULL)
}

# 3. 3 Diagnostics
plot_quality_diagnostics <- function(results, output_file = NULL) {
  plot_device(output_file, 12, 8, {
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(1, 1, 0, 0))

    hist(results$raw_corrected$r2, main = "R² Distribution", xlab = "R²", col = "lightblue", breaks = 20)
    abline(v = 0.9, col = "red", lwd = 2, lty = 2)

    hist(results$raw_corrected$corrected_mo2_kg_www,
      main = "Corrected MO2",
      xlab = "MO2", col = "lightgreen", breaks = 30
    )
    abline(v = 0, col = "red", lwd = 2, lty = 2)

    labels <- paste(results$sample_sizes$species, results$sample_sizes$id, sep = "\n")
    barplot(results$sample_sizes$n_final,
      names.arg = labels, main = "Sample Sizes",
      ylab = "N", col = "coral", las = 2, cex.names = 0.7
    )

    barplot(results$sample_sizes$pct_retained,
      names.arg = labels, main = "Retention Rate",
      ylab = "% Retained", col = "gold", las = 2, ylim = c(0, 100), cex.names = 0.7
    )
    abline(h = 50, col = "red", lty = 2)
  })

  if (!is.null(output_file)) message("Diagnostics: ", output_file)
}

# RUN PART 3 FUNCTION AND PLOT DATA: 3.4 Quality control
process_part3_quality_control <- function(
  part1_results,
  part2_results = NULL,
  r2_threshold = resp_config$thresholds$default_r2,
  generate_rmr_plot = FALSE,
  generate_diagnostics_plot = FALSE,
  rmr_plot_title = NULL,
  plot_title_prefix = NULL,
  export_path = NULL,
  verbose = TRUE) {
  if (is.null(part2_results)) {
    if (generate_rmr_plot || generate_diagnostics_plot) {
      stop("RMR and diagnostics plots require Part 2 background-corrected results.")
    }
    return(flag_part1_quality(part1_results, r2_threshold, verbose))
  }

  if (verbose) cat("\n=== PART 3: QUALITY CONTROL ===\n\n")

  if (verbose) cat("[3.1] Filtering data...\n")
  results <- apply_quality_filters(
    part2_results$corrected_chambers,
    part1_results$descriptive_data,
    r2_threshold,
    verbose
  )

  if (generate_rmr_plot && !is.null(export_path)) {
    if (verbose) cat("\n[3.2] Plotting RMR...\n")
    title <- if (!is.null(rmr_plot_title)) {
      rmr_plot_title
    } else if (!is.null(plot_title_prefix)) {
      paste(plot_title_prefix, "- Corrected RMR (MO2)")
    } else {
      "Corrected RMR (MO2)"
    }
    plot_rmr_by_temperature(results$data, title, paste0(export_path, "_experiment_plots.pdf"))
  }

  if (generate_diagnostics_plot && !is.null(export_path)) {
    if (verbose) cat("[3.3] Plotting diagnostics...\n")
    plot_quality_diagnostics(results, paste0(export_path, "_diagnostic_plots.pdf"))
  }

  if (verbose) cat("✓ Part 3 Complete\n")

  results
}
