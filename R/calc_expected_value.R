#' Conditional expected value of S given T
#'
#' Computes \eqn{E_{\theta_0}(S \mid T = t)}, the conditional expected value
#' of \eqn{S} (number of successes in group 1) given total successes
#' \eqn{T = t}, under the Fisher non-central hypergeometric distribution with
#' odds ratio \eqn{\theta_0}. That is,
#' \eqn{\sum_s s \cdot P_{\theta_0}(S = s \mid T = t)}
#' over the support of \eqn{S}. Used by \code{.find_gamma12()} when solving for
#' the randomisation probabilities \eqn{\gamma_1} and \eqn{\gamma_2}.
#'
#' @param .odds_ratio The null hypothesis odds ratio \eqn{\theta_0}. No default.
#' @param .m Number of trials in group 1.
#' @param .n Number of trials in group 2.
#' @param .t Total number of successes \eqn{t = u + v}; the conditioning value
#'   for the non-central hypergeometric distribution. The expectation
#'   \eqn{E_{\theta_0}(S \mid T = t)} is returned.
#' @param .precision Numerical precision passed to
#'   \code{BiasedUrn::dFNCHypergeo()}. No default.
#'
#' @noRd
#' @importFrom BiasedUrn dFNCHypergeo
.calc_expected_value <- function(.odds_ratio, .m, .n, .t, .precision){

  lower <- max(.t-.n, 0)
  upper <- min(.m, .t)
  support <- lower:upper

  # If odds ratios are at extreme ends, return min/max a:

  if(.odds_ratio == 0)
    return(lower)

  if(.odds_ratio == Inf)
    return(upper)

  # Else, return the expected value given our null odds ratio and the
  # non-centric hypergeometric distribution:

  exp_val <- support * BiasedUrn::dFNCHypergeo(
    x = support, m1 = .m, m2 = .n, n = .t,
    odds = .odds_ratio, precision = .precision)

  exp_val <- sum(exp_val)

  return(exp_val)

}
