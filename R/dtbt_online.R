#' Discrete-time online backtest
#'
#' @param stock data-set of stock prices
#' @param w0 initial wealth to invest
#' @param p fraction determining training set
#' @param lambda EMA filter smoothing parameter
#'
#' @description {An online backtest that fits models initially and then
#' updates it as it trades on out-of-sample data.
#' }
#' @return NULL
#' @export dtbt_online
dtbt_online <- function(stock, w0=1000, p=0.5, lambda = 0.5)
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
  noise <- sqrt(diag(stats::cov(R)))


  # Compute betas
  # betas <- ema_regression(R[train_region, ], lambda)[2, ]
  betas <- trader::computeBetas(R[train_region, ])
  # Wrap up all strategies into a list
  allocators <- list(taylor2 = function(x) kelly_taylor(x),
                     min_var = function(x) min_var(x),
                     ema_taylor = function(x) kelly_taylor(x, lambda),
                     ema_min_var = function(x) min_var(x, lambda),
                     ema_t2_beta = function(x) kelly_taylor(x, lambda, betas),
                     emat_dat = function(x) ema_dat(x, lambda)
  )
  # Compute the optimal allocations
  w <- lapply(allocators, function(Y) Y(R[train_region,-1])$allocations)
  w$index_reg <- function(x) kelly_index_reg(R)$allocations
  w$ema_index_reg <- function(x) kelly_index_reg(R, lambda)$allocations
  strategies <- append(list(equal_weight = t(equal_weight(R[train_region,-1]))), w)
  nstrat <- length(strategies)

  backtest_region <- (ntrain+1):n
  portfolios <- matrix(0, nrow=ntest+1, ncol=nstrat)
  portfolios[1, ]<- w0
  colnames(portfolios) <- c("ew", names(w))
  time1 <- Sys.time()
  print("Starting backtest....")
  # Forward backtest with no updating
  for(i in 2:(ntest+1))
  {
    # Compute the optimal allocations
    # Adding some noise/slippage to last return for training set
    # We train up to the current close return+noise
    # and trade on the current close return.
    training_dat <- R[1:(ntrain+i-1), ]
    training_dat[ntrain+i-1,] <- training_dat[ntrain+i-1,]+matrix(stats::rnorm(d, 0, noise), ncol=d)
    # Compute EMA betas or rolling betas
    betas <- ema_regression(training_dat, lambda)[2, ]
    betas <- trader::computeBetas(training_dat)

    w <- lapply(allocators, function(Y) Y(training_dat[,-1])$allocations)
    w$index_reg <- kelly_index_reg(training_dat)$allocations
    w$ema_index_reg <- kelly_index_reg(training_dat, lambda)$allocations

    w$ema_t2_beta <- kelly_taylor(training_dat[,-1], lambda, betas)$allocations
    # print(Sys.time()-time1)
    # Combine all strategies into a list
    strategies <- append(list(equal_weight = t(equal_weight(R[train_region,-1]))), w)
    for(j in 1:nstrat)
    {
      portfolios[i, j] <- portfolios[i-1, j]*(1+R[ntrain+i-1,-1]%*%strategies[[j]])
    }
    if(i == ceiling(ntest/2))
    {
      print("Halfway done...")
    }
  }

  # Add market index returns
  dim(portfolios)


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
  backtest_stats <- data.frame(backtest_stats)
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
  print(Sys.time()-time1)
}
