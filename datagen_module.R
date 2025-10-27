#' @title generate_data
#' @description Simulate covariates, event times, and true expectile regression
#'              coefficients for a given scenario.
#' @param n Number of observations to generate
#' @param scenario Simulation scenario, e.g., "sim1", "sim2", ..., "sim8"
#' @param target_Tau Vector of tau values for which expectiles are calculated (default seq(0.1, 0.9, 0.1))
#' @return A list containing:
#'         X - Covariate matrix (n x 2)
#'         true_Y - True event times
#'         true_Y_tau - True expectile values at each tau for each observation
#'         betatrue - True expectile regression coefficients (intercept + slopes) for each tau
#'         censor_params - Parameters used in the scenario (for reference)
generate_data <- function(n, scenario = "sim1", target_Tau = seq(0.1, 0.9, 0.1)) {
  
  # 1. Set scenario-specific parameters
  params <- switch(
    scenario,
    "sim1" = list(beta0 = 1, beta1 = c(2, -1), gamma0 = 1, gamma1 = c(2, -1), x_dist = "unif", err_dist = "norm", sd = 0, censor_type = "right"),
    "sim2" = list(beta0 = 1, beta1 = c(2, -1), gamma0 = 1, gamma1 = c(2, -1), x_dist = "unif", err_dist = "t", sd = 0, censor_type = "right"),
    "sim3" = list(beta0 = 1, beta1 = c(1, 1), gamma0 = 1, gamma1 = c(0, 1), x_dist = "norm", err_dist = "norm", sd = 0.5, censor_type = "right"),
    "sim4" = list(beta0 = 1, beta1 = c(1, 1), gamma0 = 1, gamma1 = c(0, 1), x_dist = "norm", err_dist = "t", sd = 0.5, censor_type = "right"),
    "sim5" = list(beta0 = 0, beta1 = c(2, 2), gamma0 = 1, gamma1 = c(0, 0), x_dist = "norm", err_dist = "t", sd = 0.5, censor_type = "left"),
    "sim6" = list(beta0 = 0.2, beta1 = c(3, -2), gamma0 = 0.5, gamma1 = c(0.5, -1), x_dist = "chisq_unif", err_dist = "norm", sd = 0, censor_type = "left"),
    "sim7" = list(beta0 = 1, beta1 = c(2, -1), gamma0 = 1, gamma1 = c(2, -1), x_dist = "unif", err_dist = "norm", sd = 0, censor_type = "interval", L_fixed = -0.3),
    "sim8" = list(beta0 = 1, beta1 = c(2, -1), gamma0 = 1, gamma1 = c(2, -1), x_dist = "unif", err_dist = "t", sd = 0, censor_type = "interval", L_fixed = -0.3),
    stop("Unknown scenario specified")
  )
  
  SNR <- 4 # Signal-to-noise ratio

  # 2. Generate covariates X
  if (params$x_dist == "unif") {
    X <- cbind(x1 = runif(n, 0, 1), x2 = runif(n, 0, 1))
  } else if (params$x_dist == "norm") {
    sd.Z <- matrix(c(1, params$sd, params$sd, 1), 2, 2)
    Z <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = sd.Z)
    X <- cbind(x1 = Z[, 1], x2 = pnorm(Z[, 2]))
  } else if (params$x_dist == "chisq_unif") {
    X <- cbind(x1 = rchisq(n, df = 1), x2 = runif(n))
  }
  X <- scale(X, center = TRUE, scale = FALSE) # Center covariates
  colnames(X) <- c("x1", "x2")

  # 3. Generate error term and true expectiles
  eps <- if (params$err_dist == "norm") rnorm(n) else rt(n, df = 3)
  mu.tau <- if (params$err_dist == "norm") enorm(target_Tau) else et(target_Tau, df = 3)
  
  t_loc <- params$beta0 + X %*% params$beta1
  sigma.X <- params$gamma0 * eps + X %*% params$gamma1 * eps
  p <- sqrt(var(t_loc) / (SNR * var(sigma.X))) %>% as.numeric()
  
  true_Y <- t_loc + sigma.X * p # True event times

  beta0_true <- params$beta0 + p * mu.tau
  beta1_true <- matrix(rep(params$beta1, each = length(target_Tau)), nrow = length(target_Tau)) + p * mu.tau %*% t(params$gamma1)
  betatrue <- cbind(beta0_true, beta1_true)

  true_Y_tau <- cbind(1, X) %*% t(betatrue)
  
  return(list(
    X = X,
    true_Y = true_Y,
    true_Y_tau = true_Y_tau,
    betatrue = betatrue,
    censor_params = params
  ))
}


#' @title apply_censoring_auto
#' @description Automatically apply censoring (right, left, interval) based on
#'              a target censor rate.
#' @param data A list containing true_Y and X
#' @param censor_rate Desired censoring proportion (0 < censor_rate < 1)
#' @param censor_type Type of censoring: "right", "left", or "interval"
#' @param L_fixed Fixed left boundary for interval censoring (optional, default -0.3)
#' @return The data list updated with:
#'         Y - Observed/censored outcome
#'         delta - Indicator of observation (1=observed, 0=censored)
#'         L, R - Censoring bounds (for interval censoring)
apply_censoring_auto <- function(data, censor_rate, censor_type, L_fixed = -0.3) {
  
  if (!(censor_rate > 0 && censor_rate < 1)) stop("censor_rate must be between 0 and 1.")
  
  true_Y <- data$true_Y
  n <- length(true_Y)
  
  if (censor_type == "right") {
    censor_point <- quantile(true_Y, probs = 1 - censor_rate, na.rm = TRUE)
    C <- runif(n, min = censor_point * 0.5, max = censor_point * 1.5)
    
    data$Y <- pmin(true_Y, C)
    data$delta <- as.integer(true_Y <= C)
    data$L <- data$Y
    data$R <- C
    
  } else if (censor_type == "left") {
    censor_point <- quantile(true_Y, probs = censor_rate, na.rm = TRUE)
    C <- runif(n, min = censor_point * 0.5, max = censor_point * 1.5)
    
    data$Y <- pmax(true_Y, C)
    data$delta <- as.integer(true_Y >= C)
    data$L <- C
    data$R <- data$Y
    
  } else if (censor_type == "interval") {
    if (is.null(L_fixed)) stop("For interval censoring, L_fixed must be provided.")
    if (length(L_fixed) == 1) L_fixed <- rep(L_fixed, n)
    
    # Determine upper bounds to achieve target censor rate
    find_rmax_error <- function(r_max, t_vector, l_vector, target) {
      if (r_max <= max(l_vector)) return(1)
      r_random <- runif(n, min = l_vector, max = r_max)
      current_rate <- mean(t_vector >= l_vector & t_vector <= r_random)
      return(current_rate - target)
    }
    
    search_interval <- c(max(L_fixed) + 1e-4, max(true_Y) * 2)
    optimal_rmax_result <- uniroot(
      f = find_rmax_error, interval = search_interval,
      t_vector = true_Y, l_vector = L_fixed, target = censor_rate, tol = 0.001
    )
    optimal_rmax <- optimal_rmax_result$root
    
    L <- L_fixed
    R <- runif(n, min = L, max = optimal_rmax)
    is_censored <- (true_Y >= L & true_Y <= R)
    
    data$delta <- as.integer(is_censored)
    data$Y <- ifelse(data$delta == 0, (L + R) / 2, true_Y)
    data$L <- ifelse(data$delta == 1, data$Y, L)
    data$R <- ifelse(data$delta == 1, data$Y, R)
    
  } else {
    stop("Unknown censor_type. Must be 'right', 'left', or 'interval'.")
  }
  
  return(data)
}


#' @title apply_censoring_fixed
#' @description Apply censoring using pre-defined parameters for each
#'              scenario and target censor rate.
#' @param data A list containing true_Y
#' @param scenario Name of the simulation scenario (e.g., "sim1")
#' @param censor_rate Desired censoring proportion
#' @return The data list updated with:
#'         Y - Observed/censored outcome
#'         delta - Indicator of observation (1=observed, 0=censored)
#'         L, R - Censoring bounds for interval or left/right censoring
#'         actual_censor_rate - Actual proportion of censored observations
apply_censoring_fixed <- function(data, scenario, censor_rate) {
  param_grid <- tribble(
    ~scenario, ~censor_rate, ~censor_type, ~param1, ~param2,
    # --- Scenario 1-4 (Right Censoring) ---
    # For right censoring, param2 is R_max
    "sim1",      0.1,          "right",      NA,     10.0,
    "sim1",      0.2,          "right",      NA,     5.0,
    "sim1",      0.3,          "right",      NA,     4.0,
    "sim2",      0.1,          "right",      NA,     8.0,
    "sim2",      0.2,          "right",      NA,     5.0,
    "sim2",      0.3,          "right",      NA,     3.0,
    "sim3",      0.1,          "right",      NA,     10.0,
    "sim3",      0.2,          "right",      NA,     5.0,
    "sim3",      0.3,          "right",      NA,     3.5,
    "sim4",      0.1,          "right",      NA,     9.0,
    "sim4",      0.2,          "right",      NA,     5.5,
    "sim4",      0.3,          "right",      NA,     3.5,
    # --- Scenario 5-6 (Left Censoring) ---
    # For left censoring, param1 is a fixed C value
    "sim5",      0.1,          "left",       -3.0,   NA,
    "sim5",      0.2,          "left",       -2.0,   NA,
    "sim5",      0.3,          "left",       -1.5,   NA,
    "sim6",      0.1,          "left",       -3.0,   NA, 
    "sim6",      0.2,          "left",       -2.5,   NA,
    "sim6",      0.3,          "left",       -2.0,   NA,
    # --- Scenario 7-8 (Interval Censoring) ---
    # For interval, param1 is L_fixed, param2 is R_max
    "sim7",      0.1,          "interval",   -0.1,   0.5,
    "sim7",      0.2,          "interval",   -0.3,   1.0,
    "sim7",      0.3,          "interval",   -0.3,   1.5,
    "sim8",      0.1,          "interval",   -0.3,   0.5,
    "sim8",      0.2,          "interval",   -0.3,   1.0, 
    "sim8",      0.3,          "interval",   -0.3,   1.5
  )

  true_Y <- data$true_Y
  n <- length(true_Y)
  
  # Lookup parameter set for the scenario
  params <- param_grid %>%
    filter(scenario == !!scenario, censor_rate == !!censor_rate)
  if (nrow(params) == 0) stop(sprintf("No parameters found for scenario '%s' and censor_rate '%s'", scenario, censor_rate))
  if (nrow(params) > 1) stop(sprintf("Multiple parameter sets found for scenario '%s' and censor_rate '%s'", scenario, censor_rate))
  
  censor_type <- params$censor_type
  
  # Apply censoring based on type and fixed parameters
  if (censor_type == "right") {
    R_max <- params$param2
    C <- runif(n, min = 0, max = R_max)
    data$Y <- pmin(true_Y, C)
    data$delta <- as.integer(true_Y <= C)
    data$L <- NA; data$R <- C
    
  } else if (censor_type == "left") {
    L_fixed <- params$param1
    C <- rep(L_fixed, n)
    data$Y <- pmax(true_Y, C)
    data$delta <- as.integer(true_Y >= C)
    data$L <- C; data$R <- NA
    
  } else if (censor_type == "interval") {
    L_fixed <- params$param1
    R_max <- params$param2
    L <- rep(L_fixed, n)
    R <- runif(n, min = L_fixed, max = R_max)
    
    is_censored <- (true_Y >= L & true_Y <= R)
    data$delta <- as.integer(!is_censored)
    data$Y <- ifelse(is_censored, (L + R) / 2, true_Y)
    data$L <- L; data$R <- R
    
  } else {
    stop("Unknown censor_type in parameter grid.")
  }
  
  data$actual_censor_rate <- mean(data$delta == 0)
  return(data)
}
