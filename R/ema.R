#' EMA filter for mixed moment matrix of E(R_i R_j)
#'
#' @param R the data-set of arithmetic returns, multivariate time-series
#' @param lambda the smoothing parameter for the EMA filter
#'
#' @description {Estimate the mixed moment matrix using an EMA filter
#' }
#' @return matrix
#' @export ema_mmm
ema_mmm <- function(R, lambda)
{
  N <- nrow(R)
  # Compute weights
  ws <- (1-lambda)^(0:(N-1))
  ws <- ws*(lambda/(1-(1-lambda)^N))
  # Compute weighted sample covariance in matrix form
  Sigma <- t(rev(ws)*R)%*%R
  # Tidy and return
  colnames(Sigma) <- colnames(R)
  rownames(Sigma) <- colnames(R)
  return(Sigma)
}

#' EMA filter for alphas and betas of a set of stocks (arithmetic returns)
#'
#' @param R the data-set of arithmetic returns, multivariate time-series, assumes the
#' first column is a market index like SPY etc.
#' @param lambda the smoothing parameter for the EMA filter
#'
#' @description {Estimate the betas of a set of stocks using an EMA filter
#' }
#' @return vector
#' @export ema_regression
ema_regression <- function(R, lambda=0.95)
{
  w <- findistr::ewmc(R, lambda, h=1)
  betas <- (w/w[1,1])[1,-1]
  y_hat <- findistr::ema(R, lambda, 1)[-1]
  x_hat <- as.numeric(findistr::ema(R[,1], lambda, 1))
  alphas <- y_hat-betas*x_hat
  return(rbind(alphas, betas))
}

#' EMA filter for DAT matrix of E((1+R_j)/(1+R_i))
#'
#' @param R the data-set of arithmetic returns, multivariate time-series
#' @param lambda the smoothing parameter for the EMA filter
#'
#' @description {Estimate the mixed DAT matrix using an EMA filter
#' }
#' @return matrix
#' @export ema_dat_mat
ema_dat_mat <- function(R, lambda = 0.94)
{
  N <- nrow(R)
  # Center data
  U <- 1+R
  V <- 1/(1+R)
  # Compute weights
  ws <- (1-lambda)^(0:(N-1))
  ws <- ws*(lambda/(1-(1-lambda)^N))
  # Compute weighted sample covariance in matrix form
  Sigma <- t(t(rev(ws)*U)%*%V)
  # Tidy and return
  colnames(Sigma) <- colnames(R)
  rownames(Sigma) <- colnames(R)
  return(Sigma)
}

#' EMA filter for DAT algorithm
#'
#' @param R the data-set of arithmetic returns, multivariate time-series
#' @param lambda the smoothing parameter for the EMA filter
#'
#' @description {Estimate the dominant asset using an EMA filter
#' }
#' @return matrix
#' @export ema_dat
ema_dat <- function(R, lambda = 0.94)
{
  n <- ncol(R)
  # EMA DAT
  datMat <- ema_dat_mat(R, lambda)
  dom_index <- which(apply(datMat<=1, 1, all))
  a <- rep(0, n)
  if(length(dom_index) == 0)
  {
    warning("No single dominant asset, underbetting")
    approxDom <- datMat-1
    approxDom <- which.min(apply(approxDom, 1, sum))
    dom_index <- approxDom
    a[dom_index] <- 1/n
  } else
  {
    a[dom_index] <- 1
  }


  return(list(dom_index = dom_index, allocations = a))
}
