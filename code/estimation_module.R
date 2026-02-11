#' @title intBeta
#' @description Estimate initial beta coefficients using only fully
#'              observed (uncensored) data.
#' @param dataset A list containing:
#'                Y     : observed or censored outcomes
#'                X     : covariate matrix
#'                delta : censoring indicator (1=observed)
#' @param imptau Vector of taus for which initial betas are computed
#' @param base_formula Formula for the linear/expectile regression
#' @return A matrix of initial beta estimates:
#'         Column 1 = intercepts for each tau
#'         Columns 2+ = slopes (replicated for each tau)
intBeta <- function(dataset, imptau, base_formula) {
  complete_cases <- dataset[['delta']] == 1
  dat0 <- data.frame(Y = dataset[['Y']][complete_cases], dataset[['X']][complete_cases, ])
  
  m <- lm(base_formula, dat0)
  coeffs <- m$coefficients
  
  beta1 <- matrix(rep(coeffs[-1], each = length(imptau)), nrow = length(imptau))
  beta0 <- coeffs[1] + expectile(m$residuals, imptau)
  
  cbind(beta0, beta1)
}


#' @title impy
#' @description Impute censored Y values based on current beta estimates.
#' @param dataset A list containing: Y, X, delta, L, R
#' @param beta Matrix of beta coefficients (intercepts + slopes)
#' @param censor_type Type of censoring: "right", "left", or "interval"
#' @return Vector of Y values with imputed values for censored observations
impy <- function(dataset, beta, censor_type) {
  Y <- dataset[['Y']]
  X <- dataset[['X']]
  delta <- dataset[['delta']]
  L <- dataset[['L']]
  R <- dataset[['R']]
  
  cen_ind <- which(delta == 0)
  if (length(cen_ind) == 0) return(Y)
  
  X_cen <- X[cen_ind, , drop = FALSE]
  yj_matrix <- beta[, 1] + X_cen %*% t(beta[, -1, drop = FALSE])
  
  for (j in 1:length(cen_ind)) {
    ind <- cen_ind[j]
    yj_preds <- yj_matrix[j, ]
    
    imp_ind <- switch(
      censor_type,
      'right' = which(yj_preds > R[ind]),
      'left' = which(yj_preds < L[ind]),
      'interval' = which(yj_preds > L[ind] & yj_preds < R[ind])
    )
    
    if (length(imp_ind) > 0) {
      Y[ind] <- yj_preds[sample(imp_ind, 1)]
    }
  }
  return(Y)
}


#' @title DAer_est
#' @description Estimate beta coefficients using imputed Y values for
#'              both target and imputation taus.
#' @param dataset A list containing covariates X
#' @param impY Vector of imputed Y values
#' @param target_Tau Vector of taus for final estimation
#' @param imptau Vector of taus used for imputation
#' @param base_formula Regression formula
#' @param update Logical, whether to update beta for imputation taus
#' @return A list containing:
#'         bhat    - Beta estimates for target_Tau
#'         impbeta - Beta estimates for imptau (used for next imputation step)
DAer_est <- function(dataset, impY, target_Tau, imptau, base_formula, update = TRUE) {
  dat0 <- data.frame(Y = impY, dataset[['X']])
  
  m_eval <- expectreg.ls(base_formula, dat0, expectiles = target_Tau, quietly = TRUE)
  bhat <- cbind(m_eval$intercepts, matrix(unlist(m_eval$coefficients), nrow = length(target_Tau)))
  predictor_names <- names(m_eval$coefficients)
  colnames(bhat) <- c("Intercept", predictor_names)
  
  impbeta <- NULL
  if(update){
    m_imp <- expectreg.ls(base_formula, dat0, expectiles = imptau, quietly = TRUE)
    impbeta <- cbind(m_imp$intercepts, matrix(unlist(m_imp$coefficients), nrow = length(imptau)))
    colnames(impbeta) <- c("Intercept", predictor_names)
  }
  
  list(bhat = bhat, impbeta = impbeta)
}


#' @title DAer
#' @description Perform iterative DAer procedure: initial beta estimation,
#'              iterative imputation, and final beta estimation.
#' @param dataset A list containing Y, X, delta, L, R
#' @param H Number of iterations
#' @param imptau Taus used for imputation
#' @param target_Tau Taus for final estimation
#' @param base_formula Regression formula
#' @param censor_type Censoring type ("right", "left", "interval")
#' @return A list containing:
#'         finalbeta   - Final beta estimates for target_Tau
#'         final_y_tau - Predicted Y for each tau
#'         time        - Execution time in seconds
DAer <- function(dataset, H, imptau, target_Tau, base_formula, censor_type) {
  start_time <- Sys.time()
  
  impbeta <- intBeta(dataset, imptau, base_formula)
  imp.y <- NULL
  
  for (h in 1:H) {
    imp.y <- impy(dataset, beta = impbeta, censor_type = censor_type)
    result <- DAer_est(dataset, imp.y, target_Tau, imptau, base_formula, update = TRUE)
    impbeta <- result$impbeta
  }
  
  final_est <- DAer_est(dataset, imp.y, target_Tau, imptau, base_formula, update = FALSE)
  end_time <- Sys.time()
  
  bhat <- final_est$bhat
  final_y_tau <- cbind(1, dataset$X) %*% t(bhat)
  
  list(finalbeta = bhat, final_y_tau = final_y_tau, time = as.numeric(end_time - start_time))
}


#' @title IPW_est
#' @description Estimate beta coefficients using inverse probability
#'              weighting (IPW) to adjust for censoring.
#' @param dataset A list containing Y, X, delta
#' @param target_Tau Vector of taus for estimation
#' @param base_formula Regression formula
#' @return A list containing:
#'         finalbeta   - Beta estimates for target_Tau
#'         final_y_tau - Predicted Y for each tau
#'         time        - Execution time in seconds
IPW_est <- function(dataset, target_Tau, base_formula) {
  dat0 <- data.frame(Y = dataset[['Y']], dataset[['X']], delta = dataset[['delta']])
  start_time <- Sys.time()
  
  ipc_formula <- update(base_formula, Surv(Y, delta) ~ .)
  m <- expectreg.ipc(ipc_formula, data = dat0, expectiles = target_Tau)
  
  end_time <- Sys.time()
  
  bhat <- cbind(m$intercepts, matrix(unlist(m$coefficients), nrow = length(target_Tau)))
  predictor_names <- names(m$coefficients)
  colnames(bhat) <- c("Intercept", predictor_names)
  
  final_y_tau <- cbind(1, dataset$X) %*% t(bhat)
  list(finalbeta = bhat, final_y_tau = final_y_tau, time = as.numeric(end_time - start_time))
}


#' @title Full_est
#' @description Estimate beta coefficients ignoring censoring (full data),
#'              i.e., simple expectile regression.
#' @param dataset A list containing Y, X
#' @param target_Tau Vector of taus for estimation
#' @param base_formula Regression formula
#' @return A list containing:
#'         finalbeta   - Beta estimates for target_Tau
#'         final_y_tau - Predicted Y for each tau
#'         time        - Execution time in seconds
Full_est <- function(dataset, target_Tau, base_formula) {
  dat_full <- data.frame(Y = dataset[['Y']], dataset[['X']])
  start_time <- Sys.time()
  
  m <- expectreg.ls(base_formula, data = dat_full, expectiles = target_Tau)
  end_time <- Sys.time()
  
  bhat <- cbind(m$intercepts, matrix(unlist(m$coefficients), nrow = length(target_Tau)))
  predictor_names <- names(m$coefficients)
  colnames(bhat) <- c("Intercept", predictor_names)
  
  final_y_tau <- cbind(1, dataset$X) %*% t(bhat)
  list(finalbeta = bhat, final_y_tau = final_y_tau, time = as.numeric(end_time - start_time))
}
