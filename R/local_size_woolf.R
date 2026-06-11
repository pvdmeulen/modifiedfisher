# Find local size of the Woolf asymptotic test --------------------------------

#' Local size of Woolf's asymptotic (log-OR Wald) test
#
#' Rejects H0 when |log(OR_hat) / SE(log OR_hat)| > z_(alpha/2). Uses the
#' Haldane 0.5 correction for zero cells, exactly as the SAS macro does
#' (a=max(u, 0.5); b=max(m-u, 0.5); etc.). Because this is a purely asymptotic
#' test it does not need a test frame.

#' @param nuisance Nuisance parameter.
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .odds_ratio The null hypothesis odds ratio being tested. No default.
#' @param .alpha The nominal significance level α. Defaults to 0.05.
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
