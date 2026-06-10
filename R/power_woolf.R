#' Find the power of the Woolf asymptotic test --------------------------------
#'
#' Identical to power_mfet(), but gamma0 is fixed at 1 so the
#' boundary is never included in the rejection region. The test frame must be
#' built for OR = 1 (the null).
#'
#' @param p Vector containing (pi1, pi2).
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .alpha The nominal significance level α. No default.
#' @param .superiority A logical. Defaults to FALSE. Setting this to TRUE will calculate the power for testing superiority.
#'
#' @keywords find power test woolf asymptotic

power_woolf <- function(p, .m, .n, .alpha, .superiority = FALSE) {

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
