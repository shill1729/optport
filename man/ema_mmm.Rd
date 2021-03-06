% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ema.R
\name{ema_mmm}
\alias{ema_mmm}
\title{EMA filter for mixed moment matrix of E(R_i R_j)}
\usage{
ema_mmm(R, lambda)
}
\arguments{
\item{R}{the data-set of arithmetic returns, multivariate time-series}

\item{lambda}{the smoothing parameter for the EMA filter}
}
\value{
matrix
}
\description{
{Estimate the mixed moment matrix using an EMA filter
}
}
