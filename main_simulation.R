library(openxlsx)
library(expectreg)
library(dirttee)
library(dplyr)
source("estimation_module.R")
source("datagen_module.R")
source("plot_module.R")

# --- Simulation parameters ---
K <- 1000  # Number of replications per setting
scenarios_to_run <- paste0("sim", 1:8)  # List of scenarios to simulate
n_obs_list <- c(100, 200, 400)  # Sample sizes to test
censor_rate_list <- c(0.1, 0.2, 0.3)  # Target censoring rates
target_Tau <- seq(0.1, 0.9, 0.1)  # Quantiles for estimation
H_val <- 20  # Number of imputation cycles
m_Tau_val <- 99  # Number of imputation quantiles


#' @title run_single_simulation
#' @description Run a single simulation replication for a specific setting
#' @param n_obs Number of observations
#' @param scenario Scenario name
#' @param target_Tau Target quantiles
#' @param censor_rate Target censoring rate
#' @param H_values Number of imputation cycles for DAer
#' @param m_Tau Number of imputation quantiles
#' @param ifAuto Logical flag for automatic censoring
#' @return A list containing estimation results, true beta values, Y quantiles, and actual censoring rate
run_single_simulation <- function(n_obs, scenario, target_Tau, censor_rate, 
                                  H_values, m_Tau, ifAuto = FALSE) {
  
  # --- Step 1: Generate true dataset ---
  true_data <- generate_data(n = n_obs, scenario = scenario, target_Tau = target_Tau)
  censor_type <- true_data$censor_params$censor_type
  
  # --- Step 2: Apply censoring ---
  if(ifAuto){
    dataset <- apply_censoring_auto(true_data, censor_rate = censor_rate, censor_type = censor_type)
  } else {
    dataset <- apply_censoring_fixed(true_data, scenario = scenario, censor_rate = censor_rate)
  }
  
  # Store actual censoring rate
  act_cr <- dataset$actual_censor_rate
  
  # --- Step 3: Define regression formula and imputation quantiles ---
  predictor_names <- colnames(dataset$X)
  base_formula <- reformulate(termlabels = predictor_names, response = "Y")
  
  # Ensure m_Tau does not exceed max(sqrt(n_obs), m_Tau)
  m_Tau <- max(floor(sqrt(n_obs)), m_Tau)
  imp_Tau <- seq(1, m_Tau, 1) / (1 + m_Tau)  # Imputation quantiles
  
  # --- Step 4: Run DAer estimation ---
  daer_res <- DAer(dataset, H = H_values, imptau = imp_Tau, target_Tau = target_Tau, 
                   base_formula = base_formula, censor_type = censor_type)
  full_res <- Full_est(dataset, target_Tau = target_Tau, base_formula = base_formula)
  
  # --- Step 5: Run IPW and Full estimation if right censoring ---
  if(censor_type == 'right'){
    ipw_res <- IPW_est(dataset, target_Tau = target_Tau, base_formula = base_formula)
  }
  
  # --- Step 6: Package results into a list ---
  if(censor_type == 'right'){
    results <- list(DAer = daer_res, IPW = ipw_res, Full = full_res,
                    betatrue = dataset$betatrue, Ytau = dataset$true_Y_tau, Realrate = act_cr)
  } else {
    results <- list(DAer = daer_res, Full = full_res,
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
summarize_scenario_results <- function(scenario, n_obs_list, censor_rate_list, target_Tau) {
  
  scenario_summary_list <- list()
  scenario_times_list <- list()
  scenario_rate_list <- list()
  
  # --- Loop over each combination of n_obs and censor_rate ---
  for (n_obs in n_obs_list) {
    for (cr in censor_rate_list) {
      
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
      
      # --- Aggregate bias and MSE across replications ---
      summed_metrics <- lapply(methods, function(m) {
        errors_per_replication <- lapply(replication_results_list, function(res) {
          bias_beta_k <- res[[m]]$finalbeta - res$betatrue
          bias_y_k <- res[[m]]$final_y_tau - res$Ytau
          list(
            bias_b = bias_beta_k,
            mse_b = bias_beta_k^2,
            bias_y = bias_y_k,
            mse_y = bias_y_k^2
          )
        })
        
        # Sum across K replications
        sum_bias_beta <- Reduce("+", lapply(errors_per_replication, function(err) err$bias_b))
        sum_mse_beta <- Reduce("+", lapply(errors_per_replication, function(err) err$mse_b))
        sum_bias_y <- Reduce("+", lapply(errors_per_replication, function(err) err$bias_y))
        sum_mse_y <- Reduce("+", lapply(errors_per_replication, function(err) err$mse_y))
        
        list(
          mean_bias_beta = sum_bias_beta / K,
          mean_mse_beta = sum_mse_beta / K,
          mean_bias_y = colMeans(sum_bias_y / K),
          mean_mse_y = colMeans(sum_mse_y / K)
        )
      })
      names(summed_metrics) <- methods
      
      # --- Store metrics for each method ---
      for(m in methods) {
        temp_df <- data.frame(
          n_obs = n_obs,
          censor_rate = cr,
          method = m,
          tau = target_Tau,
          BIAS_beta = summed_metrics[[m]]$mean_bias_beta,
          MSE_beta = summed_metrics[[m]]$mean_mse_beta,
          BIAS_y = summed_metrics[[m]]$mean_bias_y,
          MSE_y = summed_metrics[[m]]$mean_mse_y
        )
        scenario_summary_list[[length(scenario_summary_list) + 1]] <- temp_df
      }
      
      # --- Compute average computation time for each method ---
      avg_times <- sapply(methods, function(m) {
        times_k <- sapply(replication_results_list, function(res) res[[m]]$time)
        mean(unlist(times_k), na.rm = TRUE)
      })
      
      times_df <- data.frame(
        n_obs = n_obs,
        censor_rate = cr,
        method = names(avg_times),
        average_time_sec = avg_times
      )
      scenario_times_list[[length(scenario_times_list) + 1]] <- times_df
      
      # --- Compute mean actual censoring rate ---
      all_real_rates <- sapply(replication_results_list, function(res) res$Realrate)
      avg_real_rate <- mean(all_real_rates, na.rm = TRUE)
      rate_df <- data.frame(
        n_obs = n_obs,
        censor_rate_target = cr,
        censor_rate_actual_mean = avg_real_rate
      )
      scenario_rate_list[[length(scenario_rate_list) + 1]] <- rate_df
    }
  }
  
  # --- Combine all results into data frames ---
  final_df <- dplyr::bind_rows(scenario_summary_list)
  final_times_df <- dplyr::bind_rows(scenario_times_list)
  final_rates_df <- dplyr::bind_rows(scenario_rate_list)
  
  # --- Save to Excel ---
  wb <- createWorkbook()
  addWorksheet(wb, "BIAS_beta"); writeData(wb, "BIAS_beta", final_df %>% select(n_obs, censor_rate, method, tau, starts_with("BIAS_beta.")))
  addWorksheet(wb, "MSE_beta"); writeData(wb, "MSE_beta", final_df %>% select(n_obs, censor_rate, method, tau, starts_with("MSE_beta.")))
  addWorksheet(wb, "BIAS_y"); writeData(wb, "BIAS_y", final_df %>% select(n_obs, censor_rate, method, tau, BIAS_y))
  addWorksheet(wb, "MSE_y"); writeData(wb, "MSE_y", final_df %>% select(n_obs, censor_rate, method, tau, MSE_y))
  addWorksheet(wb, "Average_Time"); writeData(wb, "Average_Time", final_times_df)
  addWorksheet(wb, "Censor_Rates"); writeData(wb, "Censor_Rates", final_rates_df)
  
  output_excel_file <- file.path(summary_dir, paste0("summary_", scenario, ".xlsx"))
  openxlsx::saveWorkbook(wb, output_excel_file, overwrite = TRUE)
  
  cat(sprintf("Saved Excel summary for %s to: %s\n", scenario, output_excel_file))
}


# --- EXECUTE SIMULATIONS AND SAVE RAW DATA ---
raw_results_dir <- "results/raw_data"
dir.create(raw_results_dir, recursive = TRUE, showWarnings = FALSE)

for (scenario in scenarios_to_run) {
  for (n_obs in n_obs_list) {
    for (cr in censor_rate_list) {
      
      setting_name <- sprintf("%s_n%d_cr%.1f", scenario, n_obs, cr)
      cat(sprintf("\n--- Starting: %s (K=%d) ---\n", setting_name, K))
      
      # Run K replications
      replication_results_list <- replicate(K, {
        run_single_simulation(
          n_obs = n_obs,
          scenario = scenario,
          target_Tau = target_Tau,
          censor_rate = cr,
          H_values = H_val,
          m_Tau = m_Tau_val,
          ifAuto = FALSE
        )
      }, simplify = FALSE)
      
      # Save results for this setting
      output_file <- file.path(raw_results_dir, paste0(setting_name, ".Rdata"))
      save(replication_results_list, file = output_file)
      
      cat(sprintf("--- Finished and saved results to: %s ---\n", output_file))
    }
  }
}

cat("\n\nAll simulations are complete. Raw data saved in 'results/raw_data'.\n")


# --- RUN THE SUMMARY PROCESS FOR ALL SCENARIOS ---
summary_dir <- "results/summary"
dir.create(summary_dir, showWarnings = FALSE)

for (scenario in scenarios_to_run) {
  summarize_scenario_results(scenario, n_obs_list, censor_rate_list, target_Tau)
}


# --- PLOT FIGURES ---
# Load and summarize all raw data into one tidy data frame
full_summary_data <- load_and_summarize_results(
  base_dir = raw_results_dir,
  scenarios = scenarios_to_run,
  n_obs_list = n_obs_list,
  censor_rate_list = censor_rate_list,
  target_Tau = target_Tau
)

# Plot MSE for scenarios 1-4
scenarios_mse <- paste0("sim", 1:4)
MSE <- filter(full_summary_data$beta_MSE, Method == 'DAer')
for (scen in scenarios_mse) {
    create_stacked_plot(
      summary_df = MSE,
      scenario_to_plot = scen
    )
}

