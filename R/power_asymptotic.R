#' Power of Woolf's asymptotic test
#'
#' Computes the unconditional power of Woolf's asymptotic Wald test for
#' \eqn{H_0}: OR = 1, at response rates \eqn{(\pi_1, \pi_2)}. Rejects
#' \eqn{H_0} when
#' \eqn{|\log \hat{\theta}| / SE(\log \hat{\theta}) > z_{\alpha/2}},
#' where \eqn{\hat{\theta} = u(n-v) / ((m-u)v)}. Uses the Haldane correction
#' (replacing zero cells with 0.5) for tables where one or more cells are zero.
#'
#' @param p Length-2 numeric vector \eqn{(\pi_1, \pi_2)}: success probability
#'   in group 1 (\eqn{\pi_1}) and group 2 (\eqn{\pi_2}) at which power is
#'   evaluated.
#' @param .m Number of trials in group 1.
#' @param .n Number of trials in group 2.
#' @param .alpha Nominal significance level \eqn{\alpha}. No default.
#' @param .superiority Logical. If \code{TRUE}, power is computed only over
#'   tables where the observed rate in group 2 exceeds that in group 1.
#'   Defaults to \code{FALSE}.
#'
#' @return A single numeric value: the power of Woolf's asymptotic test at the
#'   response rates \eqn{(\pi_1, \pi_2)}, in \eqn{[0, 1]}.
#' @examples
#' power_asymptotic(p = c(0.2, 0.6), .m = 6, .n = 4, .alpha = 0.05)
#' @export
#' @keywords internal
#' @family power
power_asymptotic <- function(p, .m, .n, .alpha, .superiority = FALSE) {

  p0    <- p[[1]]
  p1    <- p[[2]]
  power <- 0
  z_crit <- stats::qnorm(1 - .alpha / 2)

  for (u in 0:.m) {
    for (v in 0:.n) {

      a <- max(u,      0.5)
      b <- max(.m - u, 0.5)
      c <- max(v,      0.5)
      d <- max(.n - v, 0.5)

      reject <- as.numeric(abs(log(a*d/(b*c)) / sqrt(1/a+1/b+1/c+1/d)) > z_crit)

      contrib <- reject *
        stats::dbinom(u, size = .m, prob = p0) *
        stats::dbinom(v, size = .n, prob = p1)

      power <- power + if (.superiority) (v/.n > u/.m) * contrib else contrib

    }
  }

  return(power)

}
