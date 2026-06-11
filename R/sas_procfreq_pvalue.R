# Find p-value from SAS Proc Freq test ----------------------------------------

#' P-value from the SAS Proc FREQ exact test used in the local size calculation
#
#' The SAS Proc FREQ Fisher exact p-value sums all hypergeometric probabilities
#' (under H0: OR = 1, i.e. the central hypergeometric) that are no greater than
#' the observed table probability. Reject H0 if that sum <= alpha. This uses
#' only base-R dhyper(),  no extra packages needed. Because the null
#' distribution is always the central hypergeometric (OR = 1), no test frame is
#' required.

#' @param s Integer input. Number of successes in m. No default.
#' @param t Integer input. Number of successes in n. No default.
#' @param m Integer input. Total number of trials m. No default.
#' @param n Integer input. Total number of trials n. No default.
#'
#' @keywords find local size hypothesis test sas procfreq


# P-value:
sas_procfreq_pvalue <- function(s, t, m, n, odds_ratio = 1) {

  lower   <- max(0, t - n)
  upper   <- min(m, t)
  support <- lower:upper

  probs         <- BiasedUrn::dFNCHypergeo(support, m, n, t, odds_ratio)
  observed_prob <- BiasedUrn::dFNCHypergeo(s,       m, n, t, odds_ratio)

  # Sum all probabilities <= observed (small tolerance for floating point)
  pvalue <- sum(probs[probs <= observed_prob + 1e-10])

  return(pvalue)

}
