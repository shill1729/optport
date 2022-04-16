
# optport

<!-- badges: start -->
<!-- badges: end -->

The goal of optport is to provide portfolio optimization functions for discrete and continuous time models of price dynamics. Additionally functions for basic backtests are available.

## Installation

You can install the latest dev version of optport from github with:

``` r
devtools::install_github("shill1729/optport")
```

## Example

A crypto backtest and allocation:

``` r
library(ravapi)
library(optport)
period <- "intraday"
interval <- "60min"
h <- timescale(period, interval, TRUE)
crypto <- getAssets(c("BTC", "ETH", "ETC", "LTC", "DOGE"), period,
                    interval)
print(tail(crypto))
w0 <- 1000
p <- 0.1
lambda <- 0.1
# par(mfrow=c(1,1))
ctbt_online(crypto, w0, p, lambda, h)
w <- kelly_gbm(trader::stockReturns(crypto), lambda, h=h)
print(w)
```

A stock backtest and allocation:

``` r
library(ravapi)
library(optport)
tickers <- c("SPY", "QQQ", "VXX", "CGC", "NTDOY", "SNE", "GE",
             "KO")
period <- "daily"
interval <- NULL
stock <- getAssets(tickers, period, interval)
h1 <- timescale(period, interval, FALSE)
w0 <- 1000
p <- 0.1
lambda <- 0.5
par(mfrow=c(2,2))
dtbt_offline(stock, w0, p, lambda)
dtbt_online(stock, w0, p, lambda)
ctbt_online(stock[,-1], w0, p, lambda, h1)
ctbt_online(stock[,-1], w0, p, 0.01, h1)

```



