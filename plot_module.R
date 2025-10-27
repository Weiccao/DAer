library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

#' @title load_and_summarize_results
#' @description Load raw simulation results from .Rdata files, unnest replication results,
#'              and compute mean MSE and average computation time.
#' @param base_dir Directory containing the raw .Rdata files
#' @param scenarios Vector of scenario names, e.g., c("sim1", "sim2")
#' @param n_obs_list Vector of sample sizes
#' @param censor_rate_list Vector of censoring rates
#' @param target_Tau Vector of expectile levels
#' @return A list containing:
#'         beta_MSE - Tidy data frame of mean squared errors by coefficient, tau, method
#'         time     - Tidy data frame of average execution time by scenario and method
load_and_summarize_results <- function(base_dir, scenarios, n_obs_list, censor_rate_list, target_Tau) {
  
  all_results_list <- list() # Initialize list to store all results
  
  # Loop over all scenarios, sample sizes, and censor rates
  for (scen in scenarios) {
    for (n in n_obs_list) {
      for (cr in censor_rate_list) {
        
        # Construct file path for the current experimental setting
        setting_name <- sprintf("%s_n%d_cr%.1f", scen, n, cr)
        file_path <- file.path(base_dir, paste0(setting_name, ".Rdata"))
        
        # Skip if the file does not exist
        if (!file.exists(file_path)) {
          warning(paste("File not found, skipped:", file_path))
          next
        }
        
        # Load the .Rdata file into a local environment to avoid overwriting variables
        env <- local({
          load(file_path)
          environment()
        })
        
        replication_results_list <- env$replication_results_list
        K <- length(replication_results_list) # Number of replications
        
        # Process each replication
        tidy_data_for_file <- lapply(1:K, function(k) {
          
          rep_res <- replication_results_list[[k]]
          methods_to_process <- c("DAer", "IPW", "FULL")
          
          # Process each estimation method
          results_per_method <- lapply(methods_to_process, function(method_name) {
            
            method_res <- rep_res[[method_name]]
            if (is.null(method_res)) return(NULL)
            
            # Compute squared error of estimated coefficients
            mse_matrix <- (method_res$finalbeta - rep_res$betatrue)^2
            colnames(mse_matrix) <- paste0("beta_", 0:(ncol(mse_matrix) - 1))
            
            # Convert matrix to data frame and add replication info
            mse_df <- as.data.frame(mse_matrix)
            mse_df$replication <- k
            mse_df$Method <- method_name 
            mse_df$tau <- target_Tau
            mse_df$time <- method_res$time
            
            mse_df
          })
          
          # Combine results from all methods for this replication
          bind_rows(results_per_method)
          
        }) %>% bind_rows() # Combine all replications for this file
        
        # Add scenario, sample size, and censor rate columns
        tidy_data_for_file$scenario <- scen
        tidy_data_for_file$n_obs <- n
        tidy_data_for_file$censor_rate <- cr
        
        all_results_list[[setting_name]] <- tidy_data_for_file
      }
    }
  }
  
  # Combine all results across settings into one data frame
  full_tidy_data <- bind_rows(all_results_list)
  
  # Pivot beta columns into long format and calculate mean MSE
  final_summary_df <- full_tidy_data %>%
    pivot_longer(
      cols = starts_with("beta_"),
      names_to = "coefficient",
      values_to = "squared_error"
    ) %>%
    group_by(scenario, n_obs, censor_rate, Method, coefficient, tau) %>%
    summarise(MSE = mean(squared_error, na.rm = TRUE), .groups = 'drop')
  
  # Calculate mean computation time for each scenario and method
  mean_time_summary <- full_tidy_data %>%
    distinct(scenario, n_obs, censor_rate, Method, replication, .keep_all = TRUE) %>%
    group_by(scenario, n_obs, censor_rate, Method) %>%
    summarise(Mean_Time_sec = mean(time, na.rm = TRUE), .groups = 'drop')
  
  return(list(beta_MSE = final_summary_df, time = mean_time_summary))
}


#' @title create_single_panel
#' @description Create a single panel plot of MSE versus tau for a specific coefficient
#'              and censoring rates.
#' @param data Tidy data frame containing MSE by tau and method
#' @param coeff_label String label of the coefficient (for title)
#' @param caption_text Optional caption text for the plot
#' @return A ggplot object of MSE vs tau for the specified coefficient
create_single_panel <- function(data, coeff_label, caption_text = NULL) {
  
  # Define custom colors and shapes for different censor rates
  custom_colors <- c("10" = "#d62728", "20" = "#2ca02c", "30" = "#1f77b4")
  custom_shapes <- c("10" = 16, "20" = 17, "30" = 15)
  
  # Create the base ggplot
  p <- ggplot(data, aes(x = tau, y = MSE, color = rate_label, shape = rate_label, group = rate_label)) +
    geom_line(linewidth = 0.8) +  # Line connecting points for each censor rate
    geom_point(size = 2.5) +      # Scatter points
    scale_y_continuous(labels = function(x) sprintf("%.3f", x)) + # Format Y-axis
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +             # X-axis ticks
    scale_color_manual(values = custom_colors, name = "Censor Rate (%)") +
    scale_shape_manual(values = custom_shapes, name = "Censor Rate (%)") +
    labs(
      title = parse(text = coeff_label),  # Use parsed expression for title
      x = NULL, 
      y = "MSE"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"  # No legend for individual panels
    )
  
  # Add optional caption if provided
  if (!is.null(caption_text)) {
    p <- p + 
      labs(caption = caption_text) +
      theme(plot.caption = element_text(hjust = 1, size = 12, face = "bold", vjust = 1))
  }
  
  return(p)
}


#' @title create_stacked_plot
#' @description Combine individual coefficient MSE panels into a stacked layout
#'              for a given scenario, adding captions for each row of sample sizes.
#' @param summary_df Tidy data frame of MSE results
#' @param scenario_to_plot Scenario name to plot
#' @param output_dir Directory to save the final stacked plot
#' @return Saves the combined plot as an EPS file; returns nothing
create_stacked_plot <- function(summary_df, scenario_to_plot, output_dir = "Figures/MSE") {
  
  # Create output directory if it does not exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Filter data for the selected scenario and prepare labels
  base_data <- summary_df %>%
    filter(scenario == scenario_to_plot) %>%
    mutate(
      coefficient_label = str_replace(coefficient, "beta_", "beta[") %>% paste0("](tau)"), # Convert to plot expression
      rate_label = factor(censor_rate * 100)
    )
  
  n_obs_values <- sort(unique(base_data$n_obs))
  coefficient_labels <- sort(unique(base_data$coefficient_label))
  
  # Loop over each sample size to create rows
  plots_by_n_row <- lapply(seq_along(n_obs_values), function(i) {
    n_val <- n_obs_values[i]
    
    # Create panels for each coefficient
    beta_plots_row <- lapply(coefficient_labels, function(coeff_lab) {
      plot_data_subset <- base_data %>%
        filter(n_obs == n_val, coefficient_label == coeff_lab)
      create_single_panel(plot_data_subset, coeff_lab)
    })
    
    # Arrange coefficients horizontally
    row_plot <- wrap_plots(beta_plots_row, nrow = 1)
    
    # Add caption for the row
    caption_expression <- sprintf("(%s)~~italic(n)==%d", letters[i], n_val)
    caption_plot <- ggplot() + 
      annotate("text", x = 0, y = 0, label = caption_expression, parse = TRUE, size = 5) + 
      theme_void()
    
    # Combine row with caption
    row_plot_with_caption <- row_plot / caption_plot + 
      plot_layout(heights = c(1, 0.1))
    final_row_group <- row_plot_with_caption +
      theme(plot.margin = margin(t = 5, r = 5, b = 20, l = 5, unit = "pt"))
    
    return(final_row_group)
  })
  
  # Stack all rows vertically
  final_plot_layout <- wrap_plots(plots_by_n_row, ncol = 1)
  
  # Add shared legend at bottom
  final_plot <- final_plot_layout + 
    plot_layout(guides = 'collect') &
    theme(
      legend.position = 'bottom',
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    )
  
  # Save plot as EPS file
  output_filename <- file.path(output_dir, sprintf("mse_%s.eps", scenario_to_plot))
  ggsave(output_filename, plot = final_plot, device = "eps", width = 10, height = 10)
  
  cat(sprintf("Stacked plot for scenario '%s' saved to %s\n", scenario_to_plot, output_filename))
}


#' @title plot_surv
#' @description Generate and save a Kaplan-Meier survival plot.
#' @param fit_object A survival fit object from the survival package
#' @param legend_labs Character vector for legend labels corresponding to strata
#' @param x_axis_label Label for the x-axis
#' @param filename File name to save the plot (under "Figures" folder)
#' @return Saves the plot as a file; prints a message to console
plot_surv <- function(fit_object, legend_labs, x_axis_label, filename) {
  
  # Create survival plot using ggsurvplot
  p <- ggsurvplot(
    fit = fit_object,           
    pval = TRUE,                   # Display p-value for log-rank test
    conf.int = TRUE,               # Show confidence intervals
    font.legend = c(20),           # Font size for legend text
    font.tickslab = c(15),         # Font size for axis tick labels
    font.x = c(20),                # Font size for x-axis label
    font.y = c(20),                # Font size for y-axis label
    pval.size = 8,                 # Font size for p-value text
    legend.title = " ",            # Blank legend title
    legend.labs = legend_labs,     # Custom labels for strata
    xlab = x_axis_label,           # x-axis label
    linetype = "strata",           # Different line type per strata
    surv.median.line = "hv",       # Add median survival lines (horizontal/vertical)
    ggtheme = theme_bw(),           # Use black-and-white theme
    palette = "lancet"             # Use lancet color palette
  )
  
  # Save the plot to the "Figures" folder
  full_path <- file.path("Figures", filename)
  ggsave(
    filename = full_path, 
    plot = p$plot,
    width = 8, 
    height = 6, 
    units = "in"
  )
  
  # Print message to console
  cat(paste("Survival plot saved to:", full_path, "\n"))
}


#' @title plot_CI
#' @description Generate coefficient plots with confidence intervals for multiple variables
#'              and save each as a separate file.
#' @param summary_df Tidy data frame containing Mean, Lower_CI, Upper_CI, Tau, Method, Variable
#' @param vars_to_plot Character vector of variable names to plot
#' @param file_format File format to save plots (default "eps")
#' @return Saves one plot per variable in "Figures/fig-turnover/beta"
plot_CI <- function(summary_df, vars_to_plot, file_format = "eps") {
  
  # Create output directory if it does not exist
  if (!dir.exists("Figures/fig-turnover/beta")) dir.create("Figures/fig-turnover/beta")
  
  # Loop over each variable to plot
  for (variable in vars_to_plot) {
    
    # Filter data for current variable
    plot_data <- summary_df %>%
      filter(Variable == variable)
    
    if (nrow(plot_data) == 0) {
      warning(paste("Variable not found in summary data:", variable, "- skipped"))
      next 
    }
    
    # Determine y-axis limits dynamically based on CI range
    y_min_val <- min(plot_data$Lower_CI, na.rm = TRUE) # Minimum lower bound
    y_max_val <- max(plot_data$Upper_CI, na.rm = TRUE) # Maximum upper bound
    range_buffer <- (y_max_val - y_min_val) * 0.05     # 5% buffer for margin
    y_lower_limit <- y_min_val - range_buffer
    y_upper_limit <- y_max_val + range_buffer
    final_y_limits <- c(min(y_lower_limit, 0), max(y_upper_limit, 0)) # Ensure zero included
    
    # Create ggplot for current variable
    p <- ggplot(plot_data, aes(x = Tau, y = Mean, color = Method, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") + # Reference line at 0
      geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI, linetype = Method), alpha = 0.25) + # CI ribbon
      geom_line(linewidth = 1, aes(linetype = Method)) +  # Mean line
      geom_point(size = 3.5, aes(shape = Method)) +        # Points for mean estimates
      scale_y_continuous(labels = scales::label_number(accuracy = 0.01), limits = final_y_limits) + # Y-axis formatting
      scale_x_continuous(breaks = seq(0.1, 0.9, 0.1), limits = c(0.1, 0.9)) + # X-axis
      scale_color_lancet() +  # Custom color palette
      scale_fill_lancet() +
      scale_shape_manual(values = c(16, 17)) +  # Custom shapes for methods
      labs(
        x = expression(tau),
        y = expression(hat(beta)(tau))
      ) +
      theme_bw() +
      theme(
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = 'top',
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.key.size = unit(35, "pt")
      )
    
    # Save plot with variable name as filename
    filename <- paste0(variable, ".", file_format)
    full_path <- file.path("Figures/fig-turnover/beta", filename)
    ggsave(full_path, p, width = 8, height = 6, device = file_format)
    
    # Print confirmation
    cat(paste("Plot saved successfully to:", normalizePath(full_path), "\n"))
  }
}
