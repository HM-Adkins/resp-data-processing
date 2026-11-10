#------------------------------------------------------------------------------
# PART 2A: CONTROL CHAMBER PROCESSING
# analyze_control_chamber, plot_control_regression, process_control_chambers
#------------------------------------------------------------------------------

# 2A.1 Analyze control chamber(s)
analyze_control_chamber <- function(chamber_list, desc_data, r2_threshold, 
                                    min_good = RESPIRATION_CONFIG$thresholds$min_good_control_points) {
  
  ctrl_idx <- which(desc_data$id == "Control")
  if (length(ctrl_idx) == 0) stop("No control chamber found")
  
  cfg <- RESPIRATION_CONFIG
  
  # Process each control
  ctrl_results <- lapply(ctrl_idx, function(i) {
    ctrl <- chamber_list[[i]]
    chamber_num <- desc_data$chamber[i]
    
    # Filter good slopes
    good <- !is.na(ctrl$r2) & ctrl$r2 > r2_threshold & !is.na(ctrl$slope)
    n_good <- sum(good, na.rm = TRUE)
    ctrl$slope_filtered <- ifelse(good, ctrl$slope, 0)
    
    if (n_good < min_good) {
      # Negligible
      message(sprintf("Control %d: Negligible background (%d points)", chamber_num, n_good))
      return(list(
        data = mutate(ctrl, slope_corrected = 0, corrected_mo2_ch = 0),
        method = "negligible",
        regression = NULL
      ))
    }
    
    # Fit model
    model_data <- ctrl[good, ]
    fit <- lm(slope ~ time, data = model_data)
    summ <- summary(fit)
    
    if (summ$r.squared < cfg$thresholds$control_regression_r2) {
      # Use mean
      message(sprintf("Control %d: Mean method (R²=%.3f)", chamber_num, summ$r.squared))
      ctrl$slope_corrected <- mean(ctrl$slope_filtered)
    } else {
      # Use regression
      message(sprintf("Control %d: Regression (R²=%.3f, p=%.4f)", 
                      chamber_num, summ$r.squared, summ$coefficients[2, 4]))
      ctrl$slope_corrected <- predict(fit, newdata = ctrl)
      # Zero out predictions outside time range
      time_range <- range(model_data$time)
      ctrl$slope_corrected[ctrl$time < time_range[1] | ctrl$time > time_range[2]] <- 0
    }
    
    # Convert to MO2
    volume_L <- desc_data$ch_vol_ml[i] / 1000
    ctrl$corrected_mo2_ch <- abs(ctrl$slope_corrected * volume_L)
    
    list(
      data = ctrl,
      method = if (summ$r.squared < cfg$thresholds$control_regression_r2) "mean" else "regression",
      regression = tibble(
        chamber = chamber_num,
        method = if (summ$r.squared < cfg$thresholds$control_regression_r2) "mean" else "regression",
        intercept = coef(fit)[1],
        slope = coef(fit)[2],
        r_squared = summ$r.squared,
        p_value = summ$coefficients[2, 4],
        n_points = nrow(model_data)
      )
    )
  })
  
  # Average across controls
  all_mo2 <- do.call(cbind, lapply(ctrl_results, function(x) x$data$corrected_mo2_ch))
  avg_mo2 <- if (length(ctrl_idx) > 1) rowMeans(all_mo2) else as.vector(all_mo2)
  
  list(
    data = tibble(time = chamber_list[[ctrl_idx[1]]]$time, mean_control_corrected_mo2_ch = avg_mo2),
    individual_controls = ctrl_results,
    regression = bind_rows(lapply(ctrl_results, function(x) x$regression)),
    method = unique(sapply(ctrl_results, function(x) x$method)),
    n_controls = length(ctrl_idx)
  )
}

# 2A.2 Control plot
plot_control_regression <- function(ctrl_analysis, title_prefix = NULL, output_file = NULL) {
  
  cfg <- RESPIRATION_CONFIG
  
  plot_device(output_file, 10, 6, {
    n <- ctrl_analysis$n_controls
    par(mfrow = c(2, n), mar = c(4, 4, 3, 1), oma = c(1, 1, 0, 0))
    
    # First row: Slope plots
    for (i in seq_len(n)) {
      ctrl <- ctrl_analysis$individual_controls[[i]]
      
      good <- !is.na(ctrl$data$slope) & ctrl$data$r2 > 0.85
      
      plot(ctrl$data$time, ctrl$data$slope,
           ylim = c(min(ctrl$data$slope, na.rm = TRUE) - 0.2, max(ctrl$data$slope, na.rm = TRUE) + 0.4),
           main = paste0(title_prefix, " - Control ", ctrl$data$chamber[1], " (Slope)"),
           xlab = "Time (hrs)", ylab = "Slope",
           pch = 19, col = ifelse(good, cfg$colors$control, cfg$colors$gray))
      
      if (!is.null(ctrl$regression)) {
        abline(a = ctrl$regression$intercept, b = ctrl$regression$slope, 
               col = cfg$colors$control, lty = 2, lwd = cfg$plot$control_line_width)
        legend("topleft", 
               legend = sprintf("R²=%.3f\np=%.4f", ctrl$regression$r_squared, ctrl$regression$p_value),
               bty = "n")
      }
    }
    
    # Second row: raw_mo2_ch plots
    for (i in seq_len(n)) {
      ctrl <- ctrl_analysis$individual_controls[[i]]
      
      good <- !is.na(ctrl$data$slope) & ctrl$data$r2 > 0.85
      
      # Calculate y-axis range: start at 0 unless minimum is below 0
      y_max <- max(ctrl$data$raw_mo2_ch, na.rm = TRUE)
      y_min <- min(ctrl$data$raw_mo2_ch, na.rm = TRUE)
      y_lim <- c(ifelse(y_min < 0, y_min, 0), y_max)
      
      plot(ctrl$data$time, ctrl$data$raw_mo2_ch,
           main = paste0(title_prefix, " - Control ", ctrl$data$chamber[1], " (Raw MO2)"),
           xlab = "Time (hrs)", ylab = "MO2 (mg O2/chamber/hr)",
           ylim = y_lim,
           pch = 19, col = ifelse(good, cfg$colors$control, cfg$colors$gray))
      
      abline(h = 0, col = cfg$colors$red, lty = 2)
    }
  })
  
  if (!is.null(output_file)) message("Control plot: ", output_file)
}


# RUN PART 2A FUNCTIONS: 2A.3
process_control_chambers <- function(
    part1_results,
    control_r2_threshold = RESPIRATION_CONFIG$thresholds$default_control_r2,
    generate_control_plot = FALSE,
    plot_title_prefix = NULL,
    export_path = NULL,
    verbose = TRUE) {
  
  if (verbose) cat("\n=== PART 2A: CONTROL CHAMBER PROCESSING ===\n\n")
  
  if (verbose) cat("[1] Analyzing control chamber(s)...\n")
  control_analysis <- analyze_control_chamber(
    part1_results$chamber_list,
    part1_results$descriptive_data,
    control_r2_threshold
  )
  
  if (generate_control_plot && !is.null(export_path)) {
    if (verbose) cat("[2] Plotting control chambers...\n")
    plot_control_regression(
      control_analysis,
      plot_title_prefix,
      paste0(export_path, "_control_plots.pdf")
    )
  }

  if (verbose) cat("✓ Control chamber processing complete\n")
  
  control_analysis
}
