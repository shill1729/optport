#' QP solver for min 0.5 x^T A x-mu^t x subject to Gx >= b
#'
#' @param drift Vector of means
#' @param Sigma Covariance matrix
#' @param betas optional betas
#' @param restraint restraint for portfolio investment
#'
#' @description {The Kelly-criterion in continuous time. A quadratic programming
#' problem with optional constraints for maximizing log-growth assuming a GBM model.
#' }
#' @return list
#' @export mv_solver
mv_solver <- function(drift, Sigma, betas=NULL, restraint = 1.0)
{
  # Function body starts here
  d <- length(drift)
  # For maximizing growth rate
  objective_function_sign <- -1
  # If a zero-vector drift is passed, then we are minimizing variance.
  if(all(drift == rep(0, d)))
  {
    objective_function_sign <- 1
  }
  num_eq_constr <- 1

  # The representation of constraints for long only is:
  # - sum x_i = -restraint (with meq = 1), and x_i >0
  if(!is.null(betas))
  {
    num_eq_constr <- 2
    # betas^T x = 0
    bvec <- c(restraint, 0, rep(0, d))
    A <- cbind(rep(1, d), betas, diag(1, d))
  } else
  {
    bvec <- c(restraint, rep(0, d))
    A <- cbind(rep(1, d), diag(1, d))
  }

  # Now pass to solve.QP to solve it!
  w <- quadprog::solve.QP(Dmat=Sigma, dvec=drift, Amat=A, bvec=bvec, meq = num_eq_constr)
  # Tidy up output
  allocations <- matrix(round(w$solution, 6))
  rownames(allocations) <- rownames(Sigma)
  colnames(allocations) <- "weight"
  expected_growth <- w$value*objective_function_sign
  volat <- as.numeric(sqrt(t(allocations)%*%Sigma%*%allocations))

  output <- list(allocations = allocations,
                 max_growth = expected_growth,
                 volatility = volat
                 )
  if(!is.null(betas))
  {
    output$beta = t(betas)%*%allocations
  }
  return(output)
}

