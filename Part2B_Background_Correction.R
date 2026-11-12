#------------------------------------------------------------------------------
# PART 2B: BACKGROUND CORRECTION
# correct_background, plot_background_correction, process_background_correction
#------------------------------------------------------------------------------

# 2B.1 Correct experimental chamber data using the control results
correct_background <- function(chamber_list, ctrl_analysis, desc_data) {
  
  bg_mo2 <- ctrl_analysis$data$mean_control_corrected_mo2_ch
  bg_threshold <- resp_config$thresholds$high_background_pct
  
  lapply(seq_along(chamber_list), function(i) {
    chamber_list[[i]] %>%
      mutate(
        corrected_mo2_ch = raw_mo2_ch - bg_mo2,
        corrected_mo2_kg_www = corrected_mo2_ch / desc_data$www_kg[i],
        background_pct = abs(bg_mo2 / raw_mo2_ch) * 100,
        correction_flag = case_when(
          background_pct > bg_threshold ~ "high_background",
          corrected_mo2_ch < 0 ~ "negative",
          TRUE ~ "ok"
        )
      )
  })
}

# 2B.2 Background correction comparison plot
plot_background_correction <- function(corrected_chambers, ctrl_analysis, desc_data, title_prefix = NULL, output_file = NULL) {
  
  cfg <- resp_config
  
  plot_device(output_file, 12, 10, {
    exp_idx <- which(desc_data$id != "Control")
    n <- length(exp_idx)
    
    nr <- ceiling(sqrt(n))
    nc <- ceiling(n / nr)
    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1), oma = c(4, 1, 3, 0))
    
    control_mo2 <- ctrl_analysis$data$mean_control_corrected_mo2_ch
    time_vec <- ctrl_analysis$data$time
    
    # Calculate global y-axis range across all chambers
    all_values <- unlist(lapply(exp_idx, function(i) {
      c(corrected_chambers[[i]]$raw_mo2_ch, 
        corrected_chambers[[i]]$corrected_mo2_ch,
        control_mo2)
    }))
    y_max <- max(all_values, na.rm = TRUE)
    y_min <- min(all_values, na.rm = TRUE)
    
    # Set y_lim: start at 0 unless minimum is below 0
    y_lim <- c(ifelse(y_min < 0, y_min, 0), y_max)
    x_label <- c(rep("",4), rep("Time (hrs)", 2))
    
    for (i in seq_along(exp_idx)) {
      idx <- exp_idx[i]
      chamber_data <- corrected_chambers[[idx]]
      
      # Plot raw data
      plot(chamber_data$time, chamber_data$raw_mo2_ch,
           main = paste(desc_data$species[idx], "-", desc_data$id[idx]),
           xlab = x_label[i], 
           ylab = "MO2 (mg O2/chamber/hr)",
           ylim = y_lim,
           pch = 16, cex = cfg$plot$point_size, col = cfg$colors$raw)
      
      # Add corrected data
      points(chamber_data$time, chamber_data$corrected_mo2_ch,
             pch = 16, cex = cfg$plot$point_size, col = cfg$colors$corrected)
      
      # Add control
      points(time_vec, control_mo2, pch = 19, cex = cfg$plot$point_size, col = cfg$colors$control)
      
      # Add zero line
      abline(h = 0, col = cfg$colors$red, lty = 2)
    }
    
    if (! is.null(title_prefix)) {
      mtext(paste(title_prefix, "- Background Correction"), 
            outer = TRUE, font = 2, cex = 1.2)
    }
    
    # Add horizontal legend at the bottom
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend("bottom", 
           legend = c("Raw", "Corrected", "Control"),
           col = c(cfg$colors$raw, cfg$colors$corrected, cfg$colors$control),
           pch = c(16, 16, 19),
           bty = "n",
           cex = 1.2,
           pt.cex = 1.5,
           horiz = TRUE,
           xpd = TRUE)
  })
  
  if (!is.null(output_file)) message("Background correction plot: ", output_file)
}

# RUN PART 2B FUNCTIONS: 2B.3
process_background_correction <- function(
    part1_results,
    control_analysis,
    generate_correction_plot = FALSE,
    plot_title_prefix = NULL,
    export_path = NULL,
    verbose = TRUE) {
  
  if (verbose) cat("\n=== PART 2B: BACKGROUND CORRECTION ===\n\n")
  
  if (verbose) cat("[1] Applying control-based correction...\n")
  corrected_chambers <- correct_background(
    part1_results$chamber_list,
    control_analysis,
    part1_results$descriptive_data
  )
  
  if (generate_correction_plot && !is.null(export_path)) {
    if (verbose) cat("[2] Plotting background correction...\n")
    plot_background_correction(
      corrected_chambers,
      control_analysis,
      part1_results$descriptive_data,
      plot_title_prefix,
      paste0(export_path, "_correction_plots.pdf")
    )
  }
  
  if (verbose) cat("✓ Background correction complete\n")
  
  list(
    control_analysis = control_analysis,
    corrected_chambers = corrected_chambers
  )
}
