#' Find the power of the randomised FE test -----------------------------------
#'
#' Uses gamma1/gamma2 from the test frame directly as rejection probabilities at
#' the boundary, rather than applying a gamma0 threshold. The test frame must be
#' built for OR = 1 (the null).
#'
#' @param p Vector containing (pi1, pi2).
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .df Testing frame (data frame) generated as part as construct_test_frame().
#' @param .alpha The nominal significance level α. Defaults to 0.05.
#' @param .precision Defines the precision by which confidence limits, p-values, and size is determined. Defaults to 1E-03.
#' @param .superiority A logical. Defaults to FALSE. Setting this to TRUE will calculate the power for testing superiority.
#'
#' @keywords find power test

power_randomised <- function(p, .m, .n, .df, .alpha, .precision,
                             .superiority = FALSE) {

  p0 <- p[[1]]
  p1 <- p[[2]]
  power <- 0

  for (u in 0:.m) {
    for (v in 0:.n) {

      t  <- u + v
      c1 <- .df$c1[[t + 1]]
      c2 <- .df$c2[[t + 1]]
      g1 <- .df$gamma1[[t + 1]]
      g2 <- .df$gamma2[[t + 1]]

      if (u < c1 | u > c2) {
        phi <- 1
      } else if (c1 == c2 & u == c1) {
        phi <- g1 + g2
      } else if (u == c1) {
        phi <- g1
      } else if (u == c2) {
        phi <- g2
      } else {
        phi <- 0
      }

      contribution <- phi *
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
