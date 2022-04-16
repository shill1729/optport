# library(ravapi)
# tickers <- c("SPY", "QQQ", "VXX", "CGC", "NTDOY", "SNE", "GE",
#              "KO")
# period <- "daily"
# interval <- NULL
# stock <- getAssets(tickers, period, interval)
# h1 <- timescale(period, interval, FALSE)
# w0 <- 1000
# p <- 0.1
# lambda <- 0.5
# par(mfrow=c(2,2))
# dtbt_offline(stock, w0, p, lambda)
# dtbt_online(stock, w0, p, lambda)
# ctbt_online(stock[,-1], w0, p, lambda, h1)
# ctbt_online(stock[,-1], w0, p, 0.01, h1)
#
