# This script should be executed from the project root directory.
# All paths are specified relative to this root to ensure portability.
library(stats)
library(utils)
library(openxlsx)
library(expectreg)
library(survival)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggsci) 
library(patchwork)

source("code/estimation_module.R")
source("code/datagen_module.R")
source("code/plot_module.R")
source("code/sim_module.R")

# ---- 0. Readme text ----
# This script runs all simulation studies reported in Section 3.
# Note that the package 'dirttee' has been removed from the CRAN repository.
# It can be downloaded manually from:
# https://cran.r-project.org/src/contrib/Archive/dirttee/
if (!requireNamespace("dirttee", quietly = TRUE)) {
  stop(
    "The package 'dirttee' is required for reproducing",
    "Please install it manually from the CRAN Archive:\n",
    "https://cran.r-project.org/src/contrib/Archive/dirttee/\n"
  )
}
library(dirttee)

# All outputs are stored under the directory 'results/simulation'.
# Specifically, the following items are generated:
# 1. Raw simulation results (RData format) saved in the 'raw_data' directory;
# 2. Summary results aggregated over 1,000 replications, saved in the
#    'summary_tabs' directory;
# 3. Figures, including Figure 2 in Section 3 and Figures S1–S3 in the
#    Supplementary Materials, saved in the 'figs' directory.


# ---- 1. Simulation parameters and path file setup ----
K <- 1000  # Number of replications per setting
scenarios_to_run <- paste0("sim", 1:8)  # List of scenarios to simulate
n_obs_list <- c(100, 200, 400)  # Sample sizes to test
censor_rate_list <- c(0.1, 0.2, 0.3)  # Target censoring rates
target_Tau <- seq(0.1, 0.9, 0.1)  # Target expectiles for estimation
H_val <- 20  # Number of imputation cycles
m_Tau_val <- 99  # Number of imputation expectiles
set_seed <- 1234 # Set seed = 1234

raw_results_dir <- "results/simulation/raw_data"
summary_dir <- "results/simulation/summary_tabs"

# ---- 2. EXECUTE SIMULATIONS AND SAVE RAW DATA ----
dir.create(raw_results_dir, recursive = TRUE, showWarnings = FALSE)

for (scenario in scenarios_to_run) {
  for (n_obs in n_obs_list) {
    for (cr in censor_rate_list) {
      
      setting_name <- sprintf("%s_n%d_cr%.1f", scenario, n_obs, cr)
      cat(sprintf("\n--- Starting: %s (K=%d) ---\n", setting_name, K))
      
      
      # Run K replications
      replication_results_list <- vector("list", K)
      
      for (k in 1:K) {
        replication_results_list[[k]] <- run_single_simulation(
          n_obs = n_obs,
          scenario = scenario,
          target_Tau = target_Tau,
          censor_rate = cr,
          H_values = H_val,
          m_Tau = m_Tau_val,
          ifAuto = FALSE,
          seed = set_seed + k # for each Monte Carlo simulation, seed = 1234 + k
        )
      }
      
      # Save results for this setting
      output_file <- file.path(raw_results_dir, paste0(setting_name, ".Rdata"))
      save(replication_results_list, file = output_file)
      
      cat(sprintf("--- Finished and saved results to: %s ---\n", output_file))
    }
  }
}

cat("\n\nAll simulations are complete. Raw data saved in 'results/raw_data'.\n")


# --- 3. RUN THE SUMMARY PROCESS FOR ALL SCENARIOS ---
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

for (scenario in scenarios_to_run) {
  summarize_scenario_results(scenario, n_obs_list, 
                             censor_rate = 0.2, target_Tau, 
                             raw_results_dir, summary_dir)
}


# --- 4. PLOT FIGURES ---
# Load and summarize all raw data into one tidy data frame
scenarios_mse <- paste0("sim", 1:4)
full_summary_data <- load_and_summarize_results(
  base_dir = raw_results_dir,
  scenarios = scenarios_mse,
  n_obs_list = n_obs_list,
  censor_rate_list = censor_rate_list,
  target_Tau = target_Tau
)

# Plot MSE for scenarios 1-4
MSE <- filter(full_summary_data$beta_MSE, Method == 'DAer')
for (scen in scenarios_mse) {
    create_stacked_plot(
      summary_df = MSE,
      scenario_to_plot = scen,
      output_dir = "results/simulation/figs",
      set_device = 'pdf'
    )
}

