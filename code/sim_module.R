#' @title run_single_simulation
#' @description Run a single simulation replication for a specific setting
#' @param n_obs Number of observations
#' @param scenario Scenario name
#' @param target_Tau Target quantiles
#' @param censor_rate Target censoring rate
#' @param H_values Number of imputation cycles for DAer
#' @param m_Tau Number of imputation quantiles
#' @param ifAuto Logical flag for automatic censoring
#' @param seed Integer, random seed for reproducibility.
#' @return A list containing estimation results, true beta values, Y quantiles, and actual censoring rate
run_single_simulation <- function(n_obs, scenario, target_Tau, censor_rate, 
                                  H_values, m_Tau, ifAuto = FALSE, seed) {
  
  # Generate data
  true_data <- generate_data(n = n_obs, scenario = scenario, target_Tau = target_Tau, seed = seed)
  censor_type <- true_data$censor_params$censor_type
  
  if(ifAuto){
    dataset <- apply_censoring_auto(true_data, censor_rate = censor_rate, censor_type = censor_type)
  } else {
    dataset <- apply_censoring_fixed(true_data, scenario = scenario, censor_rate = censor_rate)
  }

  act_cr <- dataset$actual_censor_rate
  
  # Set formula
  predictor_names <- colnames(dataset$X)
  base_formula <- reformulate(termlabels = predictor_names, response = "Y")
  
  # Get m_Tau
  m_Tau <- max(floor(sqrt(n_obs)), m_Tau)
  imp_Tau <- seq(1, m_Tau, 1) / (1 + m_Tau)  # Imputation quantiles
  
  # Run DAer estimation 
  daer_res <- DAer(dataset, H = H_values, imptau = imp_Tau, target_Tau = target_Tau, 
                   base_formula = base_formula, censor_type = censor_type)
  
  # Run IPW and Full estimation if right censoring
  if(censor_type == 'right'){
    full_res <- Full_est(dataset, target_Tau = target_Tau, base_formula = base_formula)
    ipw_res <- IPW_est(dataset, target_Tau = target_Tau, base_formula = base_formula)
  }
  
  if(censor_type == 'right'){
    results <- list(DAer = daer_res, IPW = ipw_res, Full = full_res,
                    betatrue = dataset$betatrue, Ytau = dataset$true_Y_tau, Realrate = act_cr)
  } else {
    results <- list(DAer = daer_res, 
                    betatrue = dataset$betatrue, Ytau = dataset$true_Y_tau, Realrate = act_cr)
  }
  
  return(results)
}



#' @title summarize_scenario_results
#' @description Aggregate results across replications, calculate mean bias, MSE, time, and actual censoring rates, and save to Excel
#' @param scenario Scenario name
#' @param n_obs_list Vector of sample sizes
#' @param censor_rate_list Vector of target censoring rates
#' @param target_Tau Vector of target quantiles
summarize_scenario_results <- function(scenario, n_obs_list, censor_rate, target_Tau,
                                       raw_results_dir, summary_dir) {
  
  scenario_summary_list <- list()
  scenario_times_list <- list()
  cr <- censor_rate
  
  for (n_obs in n_obs_list) {
    setting_name <- sprintf("%s_n%d_cr%.1f", scenario, n_obs, cr)
    file_path <- file.path(raw_results_dir, paste0(setting_name, ".Rdata"))
      
    if (!file.exists(file_path)) {
      warning(paste("File not found, skipping:", file_path))
      next
    }
      
    # Load K replication results
    load(file_path) # Loads 'replication_results_list'
    K <- length(replication_results_list)
      
    methods <- names(replication_results_list[[1]])
    methods <- methods[!methods %in% c("betatrue", "Ytau", "Realrate")] # Filter methods
      
    # Aggregate bias and MSE
    summed_metrics <- lapply(methods, function(m) {
      errors_per_replication <- lapply(replication_results_list, function(res) {
        bias_beta_k <- res[[m]]$finalbeta - res$betatrue
        bias_y_k <- res[[m]]$final_y_tau - res$Ytau
        list(
          mse_b = bias_beta_k^2,
          bias_y = bias_y_k,
          mse_y = bias_y_k^2
        )
      })
        
      # Sum across K replications
      sum_mse_beta <- Reduce("+", lapply(errors_per_replication, function(err) err$mse_b))
      sum_bias_y <- Reduce("+", lapply(errors_per_replication, function(err) err$bias_y))
      sum_mse_y <- Reduce("+", lapply(errors_per_replication, function(err) err$mse_y))
        
      if(scenario %in% c('sim1', 'sim2', 'sim3', 'sim4')){
        list(
          mean_mse_beta = round((sum_mse_beta / K) * 100, 3),
          mean_bias_y = round(colMeans(sum_bias_y / K), 3),
          mean_mse_y = round(colMeans(sum_mse_y / K),3))
        }else{
        list(
          mean_mse_beta = round((sum_mse_beta / K), 3),
          mean_bias_y = round(colMeans(sum_bias_y / K), 3),
          mean_mse_y = round(colMeans(sum_mse_y / K),3))
        }
      
    })
    names(summed_metrics) <- methods
      
    # Store metrics for each method 
    for(m in methods) {
      temp_df <- data.frame(
        n_obs = n_obs,
        method = m,
        tau = target_Tau,
        MSE_beta = summed_metrics[[m]]$mean_mse_beta,
        BIAS_y = summed_metrics[[m]]$mean_bias_y,
        MSE_y = summed_metrics[[m]]$mean_mse_y
      )
      scenario_summary_list[[length(scenario_summary_list) + 1]] <- temp_df
    }
      
    # Compute average computation time for each method
    avg_times <- sapply(methods, function(m) {
      times_k <- sapply(replication_results_list, function(res) res[[m]]$time)
      mean(unlist(times_k), na.rm = TRUE)
    })
      
    times_df <- data.frame(
      n_obs = n_obs,
      method = names(avg_times),
      average_time_sec = avg_times
    )
    scenario_times_list[[length(scenario_times_list) + 1]] <- times_df
  }
  
  # Combine all results
  final_df <- dplyr::bind_rows(scenario_summary_list)
  final_times_df <- dplyr::bind_rows(scenario_times_list)
  if(scenario %in% c('sim1', 'sim2', 'sim3', 'sim4')){
    average_rep_df <- final_df %>%
      dplyr::filter(tau %in% c(0.2, 0.5, 0.8)) %>%
      dplyr::mutate(method = factor(method, levels = c('DAer', 'IPW', 'Full'))) %>%
      dplyr::select(tau, n_obs, method, starts_with("MSE_beta."), BIAS_y, MSE_y) %>%
      dplyr::arrange(tau, n_obs, method)
  }else{
    average_rep_df <- final_df %>%
      dplyr::filter(tau %in% c(0.2, 0.5, 0.8)) %>%
      dplyr::filter(method == 'DAer')%>%
      dplyr::select(tau, n_obs, method, starts_with("MSE_beta."), BIAS_y, MSE_y) %>%
      dplyr::arrange(tau, n_obs, method)
  }
  
  # Save 
  wb <- createWorkbook()
  addWorksheet(wb, "average_rep")
  writeData(wb, "average_rep", average_rep_df)
  
  addWorksheet(wb, "Average_Time")
  writeData(wb, "Average_Time", final_times_df)
  
  output_excel_file <- file.path(summary_dir, paste0("summary_", scenario, ".xlsx"))
  openxlsx::saveWorkbook(wb, output_excel_file, overwrite = TRUE)
  
  cat(sprintf("Saved Excel summary for %s to: %s\n", scenario, output_excel_file))
}
