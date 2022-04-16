#' Continuous-time online backtest
#'
#' @param stock the stock prices
#' @param w0 initial wealth to invest
#' @param p fraction determining training set
#' @param lambda EMA filter smoothing parameter
#' @param timescale timescale to scale the estimates by
#'
#' @description {An online backtest that fits models initially and then
#' updates it as it trades on out-of-sample data.
#' }
#' @return NULL
#' @export ctbt_online
ctbt_online <- function(stock, w0=1000, p=0.5, lambda = 0.5, timescale=1/252)
{
  X <- trader::stockReturns(stock)
  # Use slippage std's from daily log-returns std's
  noise <- (sqrt(diag(findistr::fitGBMs(X, timeScale=1)$Sigma)))[-1]
  dS <- apply(stock, 2, diff)
  # Function body begins here
  n <- nrow(X)
  d <- ncol(X)
  ntrain <- ceiling(p*n)
  ntest <- n-ntrain
  print("Train/Test/Sample:")
  print(c(ntrain, ntest, n))
  train_region <- 1:ntrain

  # Wrap up all strategies into a list
  allocators <- list(gbm = function(x) kelly_gbm(x,h=timescale),
                     min_var = function(x) min_var(x),
                     ema_gbm = function(x) kelly_gbm(x, lambda, h=timescale),
                     ema_min_var = function(x) min_var(x, lambda)
  )
  # Compute the optimal allocations
  w <- lapply(allocators, function(Y) Y(X)$allocations)
  strategies <- append(list(equal_weight = t(equal_weight(X[train_region,]))), w)
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
    training_dat <- X[1:(ntrain+i-1), ]
    training_dat[ntrain+i-1,] <- training_dat[ntrain+i-1,]+matrix(stats::rnorm(d, 0, noise), ncol=d)



    w <- lapply(allocators, function(Y) Y(training_dat)$allocations)
    # Combine all strategies into a list
    strategies <- append(list(equal_weight = t(equal_weight(X[train_region,]))), w)
    s <- as.numeric(stock[ntrain+i,])
    for(j in 1:nstrat)
    {


      shares <- strategies[[j]]*portfolios[i-1, j]/s
      portfolios[i, j] <- portfolios[i-1, j]+t(shares)%*%dS[ntrain+i-1,]
    }

    if(i == ceiling(ntest/2))
    {
      print("Halfway done...")
    }
  }

  # Add market index returns
  dim(portfolios)


  portfolios <- cbind(portfolios, w0*c(1, cumprod(1+exp(X[backtest_region, 1])-1)))
  colnames(portfolios)[nstrat+1] <- colnames(X)[1]
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
         legend = c(colnames(X)[1], names(strategies)),
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
