#' Title Function for generating censored data
#'
#' @param n Sample size
#' @param xdis Distribution of X
#' @param edis Distribution of error term
#' @param sd Stand error, when X simulate from Gaussian distribution
#' @param tau The interested expectile level
#' @param beta0 True value of intercept
#' @param beta1 True value of slop
#' @param gamma0 True value of intercept for heterogeneous scenario
#' @param gamma1 True value of slope  for heterogeneous scenario
#' @param censor.type Censoring mechanism, 'right' for right censoring, 'left' for left censoring, 'interval' for interval censoring.
#' @param SNR Signal-to-noise ratio
#' @param c0 Parameter to generate right censoring time C
#' @param l0 Parameter to generate left fixed censoring time L
#' @param rate Parameter to control censored rate
#'
#' @return A censored dataset, a list object
#' @export

data.generate <- function(n, xdis, edis, sd = 0.5, tau = seq(0.1,0.9,0.1),
                          beta0 , beta1, gamma0, gamma1,
                          censor.type, SNR, c0 = 0, l0 = -0.3, rate){
  library(dplyr)
  if(xdis == 'norm'){
    sd.Z <- matrix(c(1,sd, sd,1),2,2)
    Z <- MASS::mvrnorm(n=n, mu=c(0,0), Sigma=sd.Z)
    x1 <- Z[,1]; x2 <- pnorm(Z[,2])
    x1  <- x1-mean(x1); x2 <- x2-mean(x2)
  }
  else{
    x1 <- runif(n,0,1)
    x2 <- runif(n,0,1)
    x1  <- x1-mean(x1); x2 <- x2-mean(x2)
  }
  X <- cbind(x1,x2)

  if(edis=='norm'){
    eps <- rnorm(n)
    mu.tau <- expectreg::enorm(tau)
  }
  else{
    eps <- rt(n,3)
    mu.tau <- expectreg::et(tau,3)
  }

  t <- beta0+X%*%beta1
  sigma.X <- gamma0*eps+X%*%gamma1*eps
  p  <- sqrt(var(t)/(SNR*var(sigma.X))) %>% as.numeric()
  t <- t + sigma.X%*%p
  beta0 <- beta0+p*mu.tau
  beta1 <- matrix(rep(beta1,each=9),9)+p*mu.tau%*%t(gamma1)
  betatrue <- cbind(beta0, beta1)

  if(censor.type == 'right'){
    c <- runif(n,0,c0)
    delta <- ifelse(c <= t, 0, 1)
    censor.p <- length(which(delta == 0))/n %>% round(2)

    if(censor.p >= rate){
      c <- runif(n,0,(c0+0.5))
      delta <- ifelse(c <= t, 0, 1)
      censor.p <- length(which(delta == 0))/n %>% round(2)
      Y <- pmin(t,c)

    }
    else{
      Y <- pmin(t,c)
    }
  }

  else if(censor.type == 'left'){
    c <- l0
    delta <- ifelse(c >= t, 0, 1)
    censor.p <- length(which(delta == 0))/n %>% round(2)

    if(censor.p >= rate){
      c <- l0-1
      delta <- ifelse(c >= t, 0, 1)
      censor.p <- length(which(delta == 0))/n %>% round(2)
      Y <- pmax(t,c)

    }
    else{
      Y <- pmax(t,c)
    }
  }

  else{
    Y <- NULL
    c <- runif(n,0,c0)
    delta <- ifelse(t >= l0 & t <= c, 0, 1)
    censor.p <- length(which(delta == 0))/n %>% round(2)


    if(censor.p >= rate){
      c <- runif(n,0,(c0-0.2))
      delta <- ifelse(c <= t, 0, 1)
      censor.p <- length(which(delta == 0))/n %>% round(2)
      for (i in 1:n) {
        if(t[i] <= l0){
          t[i]=t[i]
        }
        else{
          if(t[i] > c[i]){
            t[i]=t[i]
          }
          else{
            t[i] = c[i]
          }
        }
        Y <- c(Y,t[i])
      }
    }
    else{
      for (i in 1:n) {
        if(t[i] <= l0){
          t[i]=t[i]
        }
        else{
          if(t[i] > c[i]){
            t[i]=t[i]
          }
          else{
            t[i] = c[i]
          }
        }
        Y <- c(Y,t[i])
      }
    }
  }

  dat <- list(Y=Y, X=X, delta = delta,
              Time=t, C=c, betatrue = betatrue, l0=l0)

  return(dat)
}

