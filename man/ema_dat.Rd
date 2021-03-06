% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ema.R
\name{ema_dat}
\alias{ema_dat}
\title{EMA filter for DAT algorithm}
\usage{
ema_dat(R, lambda = 0.94)
}
\arguments{
\item{R}{the data-set of arithmetic returns, multivariate time-series}

\item{lambda}{the smoothing parameter for the EMA filter}
}
\value{
matrix
}
\description{
{Estimate the dominant asset using an EMA filter
}
}
