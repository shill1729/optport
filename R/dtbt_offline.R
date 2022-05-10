#' Discrete-time offline backtest
#'
#' @param stock the stock prices
#' @param w0 initial wealth to invest
#' @param p fraction determining training set
#' @param lambda EMA filter smoothing parameter
#'
#' @description {An offline backtest that fits models once and then
#' backtests it on out-of-sample data.
#' }
#' @return NULL
#' @export dtbt_offline
dtbt_offline <- function(stock, w0=1000, p=0.5, lambda=0.5)
{
  R <- trader::stockReturns(stock, "arithmetic")
  # Function body begins here
  n <- nrow(R)
  d <- ncol(R)
  ntrain <- ceiling(p*n)
  ntest <- n-ntrain
  print("Train/Test/Sample:")
  print(c(ntrain, ntest, n))
  train_region <- 1:ntrain
  betas <- ema_regression(R[train_region, ], lambda)[2, ]
  # Wrap up all strategies into a list
  allocators <- list(taylor2 = function(x) kelly_taylor(x[train_region,-1]),
                     marketReg = kelly_index_reg,
                     min_var = function(x) min_var(x[train_region, -1]),
                     ema_taylor = function(x) kelly_taylor(x[train_region,-1], lambda),
                     ema_beta = function(x) kelly_index_reg(x[train_region,], lambda),
                     ema_min_var = function(x) min_var(x[train_region, -1], lambda),
                     ema_t2_beta = function(x) kelly_taylor(x[train_region, -1], lambda, betas),
                     emat_dat = function(x) ema_dat(x[train_region, -1], lambda)
  )
  # Compute the optimal allocations
  w <- lapply(allocators, function(Y) Y(R)$allocations)
  strategies <- append(list(equal_weight = t(equal_weight(R[train_region,-1]))), w)
  nstrat <- length(strategies)
  backtest_region <- (ntrain+1):n

  # Forward backtest with no updating
  portfolios <- lapply(strategies, function(x) w0*c(1, cumprod(1+R[backtest_region,-1]%*%x)))
  portfolios <- do.call(cbind, portfolios)

  # Add market index returns
  portfolios <- cbind(portfolios, w0*c(1, cumprod(1+R[backtest_region, 1])))
  colnames(portfolios)[nstrat+1] <- "spy"
  print(utils::head(portfolios))

  log_port <- apply(portfolios, 2, function(x) log(x/w0))
  max_v <- max(log_port)
  min_v <- min(log_port)
  tt <- stats::time(stock[c(ntrain+1, backtest_region+1),])

  measures <- list(mean = mean,
                   sd = stats::sd,
                   sum = sum,
                   max = max,
                   min = min,
                   probLoss = function(x) mean(x<0)
  )
  y <- apply(portfolios, 2, function(x) diff(log(x)))
  backtest_stats <- lapply(measures, function(x) apply(y, 2, x))
  backtest_stats$return <- exp(backtest_stats$sum)-1
  backtest_stats$sharpe <- backtest_stats$mean/backtest_stats$sd
  backtest_stats <- t(do.call(rbind, backtest_stats))
  print(backtest_stats)


  # Plot the log growths log(W_t/W0)
  plot(tt,log_port[, nstrat+1], type = "l", ylim = c(min_v, max_v), lty = "dashed")
  for(i in 1:nstrat)
  {
    graphics::lines(tt, log_port[, i], col = i)
  }
  graphics::legend(x="topleft",
         legend = c("spy", names(strategies)),
         # fill = c(1, 1:nstrat),
         lty = c(2, rep(1, nstrat)),
         col = c(1, 1:nstrat),
         cex = 0.6
  )
  # Find out the best performing strategies out of the measures
  strat_names <- rownames(backtest_stats)
  max_winners <- apply(backtest_stats, 2, which.max)
  max_winners <- strat_names[max_winners]
  min_winners <- apply(backtest_stats, 2, which.min)
  min_winners <- strat_names[min_winners]
  winners <- rbind(max_winners, min_winners)
  rownames(winners) <- c("max", "min")
  colnames(winners) <- colnames(backtest_stats)
  winners <- t(winners)
  print(winners)
}

