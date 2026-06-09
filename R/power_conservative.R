#' Find the power of the conservative FE test ---------------------------------
#'
#' Identical to power_mfet(), but gamma0 is fixed at 1 so the
#' boundary is never included in the rejection region. The test frame must be
#' built for OR = 1 (the null).
#'
#' @param p Vector containing (pi1, pi2).
#' @param .gamma0 Some randomisation probability gamma0.
#' @param .odds_ratio The null hypothesis odds ratio being tested. No default.
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .df Testing frame (data frame) generated as part as construct_test_frame().
#' @param .alpha The nominal significance level α. Defaults to 0.05.
#' @param .precision Defines the precision by which confidence limits, p-values, and size is determined. Defaults to 1E-03.
#' @param .superiority A logical. Defaults to FALSE. Setting this to TRUE will calculate the power for testing superiority.
#'
#' @keywords find power test

power_conservative <- function(p, .m, .n, .df, .alpha, .precision,
                               .superiority = FALSE) {

  p0 <- p[[1]]
  p1 <- p[[2]]
  power <- 0

  for (u in 0:.m) {
    for (v in 0:.n) {
      z <- c(u, u + v)

      # gamma0 = 1 means reject only when strictly outside [c1, c2]
      reject <- mod_fe_test(
        z, .df, .gamma0 = 1, .odds_ratio = 1,
        .m = .m, .n = .n, .alpha = .alpha, .precision = .precision
      )

      contribution <- reject *
        stats::dbinom(u, size = .m, prob = p0) *
        stats::dbinom(v, size = .n, prob = p1)

      if (.superiority) {
        power <- power + (v / .n > u / .m) * contribution
      } else {
        power <- power + contribution
      }
    }
  }

  return(power)

}
