#' Power of the SAS Proc FREQ exact test
#'
#' Computes the unconditional power of SAS Proc FREQ's exact test for
#' \eqn{H_0}: OR = 1, at response rates \eqn{(\pi_1, \pi_2)}. For each
#' possible table \eqn{(u, v)}, the rejection decision uses
#' \code{pvalue_procfreq()}: the sum of central hypergeometric
#' probabilities given \eqn{T = u + v} that are no greater than the observed
#' probability, compared against \eqn{\alpha}. This mirrors the rejection rule
#' used in \code{local_size_sas_freq()}.
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
#' @export
#' @keywords internal
#' @family power
power_procfreq <- function(p, .m, .n, .alpha, .superiority = FALSE) {

  p0    <- p[[1]]
  p1    <- p[[2]]
  power <- 0

  for (u in 0:.m) {
    for (v in 0:.n) {

      # Rejection rule: Fisher exact p-value under H0: OR = 1 (central hypergeometric)
      t      <- u + v
      p_val  <- pvalue_procfreq(u, t, .m, .n)
      reject <- as.numeric(p_val <= .alpha)

      contrib <- reject *
        stats::dbinom(u, size = .m, prob = p0) *
        stats::dbinom(v, size = .n, prob = p1)

      power <- power + if (.superiority) (v/.n > u/.m) * contrib else contrib

    }
  }

  return(power)

}
