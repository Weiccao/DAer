#' Function for parameters initialization
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


#' Function for imputation censoring samples
#'
#' Return the imputaion value for censoring samples
#'
#' @import expectreg
#' @importFrom dplyr %>%
#'
#' @param dataset A censoring dataset
#' @param beta The imputation parameters
#' @param censor.type Censoring  mechanism, 'right' for right censoring, 'left' for left censoring, 'interval' for interval censoring
#'
#' @return Imputation value for censoring samples
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
#' @param tau A set of expectile levels of interest
#' @param imptau The imputation expectile. The default value isn\{0.01,...0.09\}
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

#' Resemble iteration
#'
#' @param betahat estimate
#' @param yhat fitted value
#' @param H itertation time
#'

DAer_resemble <- function(betahat, yhat, H){
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


#' Data Augmentation for censored expectile regression
#'
#' @import expectreg
#' @import dirttee
#' @import dplyr
#'
#' @param dataset A censoring dataset
#' @param tau A set of expectile levels of interest
#' @param imptau The imputation expectile. The default value isn\{0.01,...0.09\}
#' @param H Iteration time, default H = 20
#' @param censortype Censoring  mechanism, 'right' for right censoring, 'left' for left censoring, 'interval' for interval censoring.
#' @param formula A formula expression for regression model
#'
#' @return A list object, which is basically a list consisting of:
#' \item{finalbeta}{Parameter estimation of interest}
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
  result <- DAer_resemble(betahat, yhat, H=H)
  return(result)
}
