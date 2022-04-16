# library(ravapi)
# library(optport)
# period <- "intraday"
# interval <- "60min"
# h <- timescale(period, interval, TRUE)
# crypto <- getAssets(c("BTC", "ETH", "ETC", "LTC", "DOGE"), period,
#                     interval)
# print(tail(crypto))
# w0 <- 1000
# p <- 0.1
# lambda <- 0.1
# # par(mfrow=c(1,1))
# # ctbt_online(crypto, w0, p, lambda, h)
# w <- kelly_gbm(trader::stockReturns(crypto), lambda, h=h)
# print(w)
#
#
