#' Local size of Woolf's asymptotic test
#'
#' Computes the unconditional rejection probability (local size) of Woolf's
#' asymptotic Wald test for \eqn{H_0}: OR = \eqn{\theta_0} at nuisance
#' parameter \eqn{p_0}, where
#' \eqn{p_1 = p_0 / (p_0 + \theta_0 (1 - p_0))}. Rejects \eqn{H_0} when
#' \eqn{|\log \hat{\theta} - \log \theta_0| / SE(\log \hat{\theta}) > z_{\alpha/2}},
#' where \eqn{\hat{\theta} = u(n-v) / ((m-u)v)}. Uses the Haldane correction
#' (replacing zero cells with 0.5) to handle tables where one or more cells
#' are zero, matching the SAS macro.
#'
#' @param nuisance The nuisance parameter \eqn{p_0}: success probability in
#'   group 1 under \eqn{H_0}. Must be in \eqn{[0, 1]}.
#' @param .m Number of trials in group 1.
#' @param .n Number of trials in group 2.
#' @param .odds_ratio The null hypothesis odds ratio \eqn{\theta_0}. No default.
#' @param .alpha Nominal significance level \eqn{\alpha}. No default.
#'
#' @keywords find local size hypothesis test woolf asymptotic

local_size_woolf <- function(nuisance, .odds_ratio, .m, .n, .alpha) {

  p0 <- min(max(0, nuisance), 1)
  p1 <- p0 / (p0 + .odds_ratio * (1 - p0))

  z_crit <- stats::qnorm(1 - .alpha / 2)
  size   <- 0

  for (i in 0:.m) {
    for (j in 0:.n) {

      # Haldane correction: replace zeros with 0.5 (matches SAS macro)
      a <-   max(i,      0.5)
      b <-   max(.m - i, 0.5)
      c <-   max(j,      0.5)
      d <-   max(.n - j, 0.5)

      log_or <- log(a * d / (b * c))
      se     <- sqrt(1/a + 1/b + 1/c + 1/d)

      reject <- as.numeric(abs(log_or - log(.odds_ratio)) / se > z_crit)

      size <- size +
        reject *
        stats::dbinom(i, size = .m, prob = p0) *
        stats::dbinom(j, size = .n, prob = p1)
    }
  }

  return(size)

}
