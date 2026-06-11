#' Local size of the SAS Proc FREQ exact test
#'
#' Computes the unconditional rejection probability (local size) of the SAS
#' Proc FREQ exact test for \eqn{H_0}: OR = \eqn{\theta_0} at nuisance
#' parameter \eqn{p_0}, where
#' \eqn{p_1 = p_0 / (p_0 + \theta_0 (1 - p_0))}. For each possible table
#' \eqn{(i, j)}, rejection is determined by \code{sas_procfreq_pvalue()}: the
#' sum of Fisher non-central hypergeometric probabilities (given
#' \eqn{T = i + j} and OR = \eqn{\theta_0}) that are no greater than the
#' observed probability, compared against \eqn{\alpha}.
#'
#' @param nuisance The nuisance parameter \eqn{p_0}: success probability in
#'   group 1 under \eqn{H_0}. Must be in \eqn{[0, 1]}.
#' @param .m Number of trials in group 1.
#' @param .n Number of trials in group 2.
#' @param .odds_ratio The null hypothesis odds ratio \eqn{\theta_0}. No default.
#' @param .alpha Nominal significance level \eqn{\alpha}. No default.
#'
#' @keywords find local size hypothesis test sas procfreq

# Size:
local_size_sas_freq <- function(nuisance, .odds_ratio, .m, .n, .alpha) {

  p0 <- min(max(0, nuisance), 1)
  p1 <- p0 / (p0 + .odds_ratio * (1 - p0))

  size <- 0

  for (i in 0:.m) {
    for (j in 0:.n) {

      t      <- i + j
      p_val  <- sas_procfreq_pvalue(i, t, .m, .n, .odds_ratio)
      reject <- as.numeric(p_val <= .alpha)

      size <- size +
        reject *
        stats::dbinom(i, size = .m, prob = p0) *
        stats::dbinom(j, size = .n, prob = p1)

    }
  }

  return(size)

}
