#' P-value from the SAS Proc FREQ exact test
#'
#' Computes the two-sided Fisher exact p-value for \eqn{H_0}: OR =
#' \code{odds_ratio} by summing all probabilities in the conditional
#' distribution of \eqn{U} given \eqn{T = t} that are no greater than the
#' observed probability \eqn{P_{\theta_0}(U = s \mid T = t)}. The conditional
#' distribution is the Fisher non-central hypergeometric with the specified odds
#' ratio, computed via \code{BiasedUrn::dFNCHypergeo()}. When
#' \code{odds_ratio = 1} this reduces to the standard central hypergeometric,
#' matching the p-value reported by SAS Proc FREQ for \eqn{H_0}: OR = 1.
#'
#' @param s Number of successes observed in group 1. No default.
#' @param t Total number of successes across both groups
#'   (\eqn{t = s + } successes in group 2); the conditioning variable in the
#'   hypergeometric distribution. No default.
#' @param m Number of trials in group 1. No default.
#' @param n Number of trials in group 2. No default.
#' @param odds_ratio Null hypothesis odds ratio \eqn{\theta_0}. Defaults to 1.
#'
#' @export
#' @keywords internal
#' @seealso [local_size_procfreq()] for the local size of the SAS Proc FREQ exact test; [power_procfreq()] for the power of the SAS Proc FREQ exact test; [modified_fisher_exact_test()] for the main user-facing function.
pvalue_procfreq <- function(s, t, m, n, odds_ratio = 1) {

  lower   <- max(0, t - n)
  upper   <- min(m, t)
  support <- lower:upper

  probs         <- BiasedUrn::dFNCHypergeo(support, m, n, t, odds_ratio)
  observed_prob <- BiasedUrn::dFNCHypergeo(s,       m, n, t, odds_ratio)

  # Sum all probabilities <= observed (small tolerance for floating point)
  pvalue <- sum(probs[probs <= observed_prob + 1e-10])

  return(pvalue)

}
