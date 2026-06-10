# Find local size of SAS Proc Freq test ---------------------------------------

#' Local size of the SAS Proc FREQ "exact" test
#
#' The SAS Proc FREQ Fisher exact p-value sums all hypergeometric probabilities
#' (under H0: OR = 1, i.e. the central hypergeometric) that are no greater than
#' the observed table probability. Reject H0 if that sum <= alpha. This uses
#' only base-R dhyper(),  no extra packages needed. Because the null
#' distribution is always the central hypergeometric (OR = 1), no test frame is
#' required.

#' @param nuisance Nuisance parameter.
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .odds_ratio The null hypothesis odds ratio being tested. No default.
#' @param .alpha The nominal significance level α. Defaults to 0.05.
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
      p_val  <- sas_procfreq_pvalue(i, t, .m, .n)
      reject <- as.numeric(p_val <= .alpha)

      size <- size +
        reject *
        stats::dbinom(i, size = .m, prob = p0) *
        stats::dbinom(j, size = .n, prob = p1)

    }
  }

  return(size)

}
