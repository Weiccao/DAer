#' Initializes the imputation parameters and returns them.
#'
#' Return initialized imputation parameters, the expectile parameter is estimated under common slope assumption.
#'
#' @importFrom expectreg expectile
#'
#'
#' @param dataset A censoring dataset
#' @param imptau The imputation expectile. The default value isn\{0.01,...0.09\}
#' @param formula A formula expression for regression model
#'
#' @return Initialized imputation parameters
#' @export

intBeta <- function(dataset, imptau = seq(0.01,0.99,0.01), formula){
  ind <- dataset[['delta']]
  Y <- dataset[['Y']][which(ind==1)]
  X <- dataset[['X']][which(ind==1),]
  #X <- as.matrix(X)
  formula <- as.formula(formula)

  dat0 <- data.frame(y = Y, X)

  m <- lm(formula,dat0)
  beta1 <- m$coefficients[-1]
  beta1 <- matrix(rep(beta1,each = length(imptau)),
                  length(imptau),length(beta1))
  beta0 <- m$coefficients[1]
  epsilon <- m$residuals

  beta0 <- beta0+expectile(epsilon,imptau)
  intbeta <- cbind(beta0,beta1)
  return(intbeta)
}


#' Computes and returns the imputed values for censored samples.
#'
#' Return the imputation value for censoring samples
#'
#' @import expectreg
#' @importFrom("dplyr", "%>%")
#'
#' @param dataset A censoring dataset
#' @param beta The imputation parameters
#' @param censor.type The censoring mechanism. Options include right-censoring, left-censoring, and interval-censoring.
#'
#' @return Imputed value for censoring samples
#' @export
#'
#'
impy <- function (dataset, beta, censor.type = c('right', ' left', 'interval')) {
  cen.ind <- dataset[['delta']]
  cen.ind <- which(cen.ind == 0)
  X <-dataset[['X']] %>% t()
  Y <- dataset[['Y']]

  if(censor.type == 'right'){
    if (length(cen.ind)) {
      for (j in 1:length(cen.ind)) {
        ind <- cen.ind[j]; yj <- beta[,1]+beta[,-1]%*%X[,ind]
        imp.ind <- which(yj > Y[ind])
        if(length(imp.ind) !=0){
          if (length(imp.ind) != 1)
            Y[ind] <- yj[sample(imp.ind, 1)]
          else Y[ind] <- yj[imp.ind]
        }
        else Y[ind] <- Y[ind]
      }
    }
  }

  else if(censor.type == 'left'){
    l0 <- dataset[['l0']]
    if (length(cen.ind)) {
      for (j in 1:length(cen.ind)) {
        ind <- cen.ind[j]
        yj <- beta[,1]+beta[,-1]%*%X[,ind]
        imp.ind <- which(yj < l0)
        if(length(imp.ind) !=0){
          if (length(imp.ind) != 1)
            Y[ind] <- yj[sample(imp.ind, 1)]
          else Y[ind] <- yj[imp.ind]
        }
        else Y[ind] <- Y[ind]
      }
    }
  }

  else{
    l0 <- dataset[['l0']]
    if (length(cen.ind)) {
      for (j in 1:length(cen.ind)) {
        ind <- cen.ind[j]
        yj <- beta[,1]+beta[,-1]%*%X[,ind]
        imp.ind <- which(yj > l0 &  yj < Y[ind])
        if(length(imp.ind) !=0){
          if (length(imp.ind) != 1)
            Y[ind] <- yj[sample(imp.ind, 1)]

          else Y[ind] <- yj[imp.ind]
        }
        else Y[ind] <- Y[ind]
      }
    }
  }

  impY <- Y
  return(impY)
}


#' Function for Data Augmentation
#'
#' @import expectreg
#' @import dirttee
#' @import dplyr
#'
#' @param dataset A censoring dataset
#' @param impY Dependent value after imputation
#' @param tau The expectile level, representing the target point on the response distribution.
#' @param imptau TA sequence of expectile levels used for imputing censored values. The default is a dense sequence from 0.01 to 0.99.
#' @param formula A formula expression for regression model
#'
#' @return A list object, which is basically a list consisting of:
#' \item{bhat}{Parameter estimation of interest}
#' \item{yhat}{Fitted value}
#' \item{intbeta}{The updata imputation parameters}

DAer_est <- function(dataset, impY, tau = seq(0.1,0.9,0.1),
                     imptau = seq(0.01,0.99,0.01), formula){

  formula <- as.formula(formula)
  dat0 <- data.frame(y = impY,
                     dataset[['X']])

  m <- expectreg.ls(formula, dat0,
                    expectiles = tau, quietly=TRUE)
  beta1 <- matrix(unlist(m$coefficients),length(tau),
                  length(m$coefficients))
  beta0 <- m$intercepts
  bhat <- cbind(beta0,beta1)
  yhat <- m$fitted

  m <- expectreg.ls(formula, dat0,
                    expectiles = imptau, quietly=TRUE)

  beta1 <- matrix(unlist(m$coefficients),length(imptau),
                  length(m$coefficients))
  beta0 <- m$intercepts
  intbeta <- cbind(beta0,beta1)

  return(result = list(bhat = bhat, yhat = yhat, intbeta = intbeta))
}

#' Combines the estimators from multiple iterations to produce a robust and aggregated estimate.
#'
#' @param betahat estimators
#' @param yhat fitted value
#' @param H iteration time
#'

DAer_ensemble <- function(betahat, yhat, H){
  finalbeta <- 0
  finaly <- 0
  for (h in 1:H) {
    finalbeta <- finalbeta+betahat[[h]]
    finaly <- finaly+yhat[[h]]
  }
  finalbeta <- finalbeta/H
  finaly <- finaly/H

  return(result = list(finalbeta = finalbeta, finalyhat = finaly))
}


#' Implements the data augmentation steps specific to expectile regression.
#'
#' @import expectreg
#' @import dirttee
#' @import dplyr
#'
#' @param dataset A censoring dataset
#' @param tau The expectile level, representing the target point on the response distribution.
#' @param imptau A sequence of expectile levels used for imputing censored values. The default is a dense sequence from 0.01 to 0.99.
#' @param H The total number of iterations for the algorithm. The default value is set to 20.
#' @param censortype Censoring  mechanism, 'right' for right censoring, 'left' for left censoring, 'interval' for interval censoring.
#' @param formula The regression formula, consisting of the response variable, '~' and the sum of all effects that should be taken into consideration.
#'
#' @return A list object, which is basically a list consisting of:
#' \item{finalbeta}{Parameter estimation at interest $tau$ level}
#' \item{finalyhat}{Fitted value}
#' @export
#'
#' @example
#' library(dirttee)
#' library(dplyr)
#' data("colcancer")
#' dat <- list()
#' dat[['delta']] <- colcancer$death
#' dat[['Y']] <- colcancer$logfollowup
#' dat[['X']] <- as.matrix(dplyr::select(colcancer, LNE, age))
#' f <- 'y ~ LNE+age'
#' res <- DAer(dat, censortype = 'right', formula = f)

DAer <- function(dataset, imptau = seq(0.01,0.99,0.01), tau = seq(0.1,0.9,0.1),
                 H = 20, censortype = c('right', ' left', 'interval'), formula){

  intbeta  <- intBeta(dataset, imptau,formula = formula)

  betahat <- NULL
  yhat <- NULL

  for(h in 1:H){
    imp.y <- impy(dataset, beta = intbeta, censor.type = censortype) # 插补后数据
    result <- DAer_est(dataset, imp.y, tau = tau, formula = formula)

    betahat[[h]] <- result$bhat
    yhat[[h]] <- result$yhat

    intbeta <- result$intbeta
  }
  result <- DAer_ensemble(betahat, yhat, H=H)
  return(result)
}
