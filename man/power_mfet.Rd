% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_power.R
\name{power_mfet}
\alias{power_mfet}
\title{Find the power of the MFET -------------------------------------------------}
\usage{
power_mfet(
  p,
  .gamma0,
  .odds_ratio,
  .m,
  .n,
  .df,
  .alpha,
  .precision,
  .superiority
)
}
\arguments{
\item{p}{Vector containing (pi1, pi2).}

\item{.gamma0}{Some randomisation probability gamma0.}

\item{.odds_ratio}{The null hypothesis odds ratio being tested. No default.}

\item{.m}{Integer input responses and sample sizes. Tests u/m versus v/n. No default.}

\item{.n}{Integer input responses and sample sizes. Tests u/m versus v/n. No default.}

\item{.df}{Testing frame (data frame) generated as part as construct_test_frame().}

\item{.alpha}{The nominal significance level α. Defaults to 0.05.}

\item{.precision}{Defines the precision by which confidence limits, p-values, and size is determined. Defaults to 1E-03.}

\item{.superiority}{A logical. Defaults to FALSE. Setting this to TRUE will calculate the power for testing superiority.}
}
\description{
Define power of the MFET as a function of z = (u, u+v) for the null OR=1.
}
\keyword{find}
\keyword{power}
\keyword{test}
