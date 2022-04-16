#' The equal-weighted strategy
#'
#' @param X data-set of arithmetic returns or log-returns
#' @description {The equally weighted diversified strategy.}
#'
#' @return vector
#' @export equal_weight
equal_weight <- function(X)
{
  d <- ncol(X)

  return(t(matrix(rep(1/d, d))))
}


#' The Kelly-criterion strategy assuming a GBM model on log-returns
#'
#' @param X data-set of periodic log-returns
#' @param lambda smoothing parameter for optional EMA filter
#' @param h timescale to use
#' @description {Maximize log growth by solving a QP problem.}
#'
#' @return vector
#' @export kelly_gbm
kelly_gbm <- function(X, lambda = 0, h=1/252)
{
  gbm <- NULL
  if(lambda == 0)
  {
    gbm <- findistr::fitGBMs(X)
  } else
  {
    gbm <- findistr::fit_ema_gbm(X, lambda, h)
  }
  gbmk <- mv_solver(gbm$drift, gbm$Sigma)
  return(gbmk)
}


#' Minimum variance portfolio
#'
#' @param R data-set of periodic arithmetic returns
#' @param lambda optional smoothing parameter for EMA, use 0 for naive estimates (default)
#' @param restraint the percentage of wealth to restrain to for investing
#' @description {Minimum variance portfolio.}
#' @return list of allocations, growth-rate and volatility
#' @export min_var
min_var <- function(R, lambda=0, restraint = 1)
{
  Sigma <- stats::cov(R)
  if(lambda > 0)
  {
    Sigma <- findistr::ewmc(R, lambda, h=1)
  }

  drifts <- rep(0, ncol(R))
  return(mv_solver(drifts, Sigma, restraint = restraint))

}

#' Approximate log-optimal portfolio in discrete time with beta-hedging.
#'
#' @param R data-set of periodic arithmetic returns
#' @param lambda optional smoothing parameter for EMA, use 0 for naive estimates (default)
#' @param betas regression coefficients from a linear regression against market-index, or NULL
#' @param restraint the percentage of wealth to restrain to for investing
#' @description {The same Taylor-approximation to the log-growth problem but with
#' beta neutral hedging as a constraint, as an option.}
#' @return list of allocations, growth-rate and volatility
#' @export kelly_taylor
kelly_taylor <- function(R, lambda=0, betas=NULL, restraint=1)
{
  N <- nrow(R)
  mu <- NULL
  M <- NULL
  # Naive point-estimates of mean and mixed product moments
  # of arithmetic returns
  if(lambda == 0)
  {
    mu <- rep(1, N)%*%R/N
    M <- t(R)%*%R/N
  } else
  {
    mu <- findistr::ema(R, lambda, h=1)
    M <- ema_mmm(R, lambda)
  }

  # Solve the QP problem
  q <- mv_solver(mu, M, betas, restraint)
  return(q)
}

#' Optimal log growth under linear regression dynamics to market index
#'
#' @param R first column is assumed to be a market index, like SPY, QQQ, etc
#' @param lambda the smoothing parameter for the EMA filter
#' @param restraint portfolio restraint to pick
#'
#' @description {This uses the Taylor approximation and substitutes the linear
#' regression equations against a market index for the vector of returns.
#' }
#' @return list
#' @export kelly_index_reg
kelly_index_reg <- function(R, lambda=0, restraint=1)
{
  M1 <- mean(R[,1])
  M2 <- mean(R[,1]^2)
  if(lambda==0)
  {
    betas <- trader::computeBetas(R)
    alphas <- trader::computeAlphas(R)
  } else if (lambda >0)
  {
    reg <- ema_regression(R, lambda)
    betas <- reg[2,]
    alphas <- reg[1,]
  }

  RR <- alphas%*%t(alphas)+M1*(alphas%*%t(betas)+betas%*%t(alphas))+M2*betas%*%t(betas)+stats::cov(R[,-1]-alphas-R[,1]%*%(betas))
  m <- alphas+M1*betas
  rownames(RR) <- colnames(R[,-1])
  betas1 <- betas
  if(all(betas > 0) || all(betas < 0))
  {
    warning("Beta-neutral constraints have been ignored since all betas have
            the same sign so the constraint is inconsistent/infeasible.")
    betas1 <- rep(0, ncol(R)-1)
  }
  w1 <- mv_solver(m, RR, betas1, restraint = restraint)
  return(w1)
}
