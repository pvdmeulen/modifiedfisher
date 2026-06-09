# Find local size of randomised test ------------------------------------------

#' Local size of the randomised (UMPU) Fisher Exact Test The existing
#' local_size() uses mod_fe_test(), which applies a gamma0 threshold and returns
#' 0 or 1. That works for the modified and conservative tests, but the fully
#' randomised test genuinely randomises at the boundary: it rejects with
#' probability gamma1 (lower) or gamma2 (upper), not deterministically. This
#' function computes that expected size directly. Find the local size of the
#' MFET as a function of the null OR theta_0.

#' @param nuisance Nuisance parameter.
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .df Testing frame (data frame) generated as part as construct_test_frame().
#' @param .odds_ratio The null hypothesis odds ratio being tested. No default.
#' @param .alpha The nominal significance level α. Defaults to 0.05.
#' @param .precision Defines the precision by which confidence limits, p-values, and size is determined. Defaults to 1E-03.
#'
#' @keywords find local size hypothesis test randomised randomized

local_size_randomised <- function(nuisance, .m, .n, .df, .odds_ratio,
                                  .alpha, .precision) {

  p0 <- min(max(0, nuisance), 1)
  p1 <- p0 / (p0 + .odds_ratio * (1 - p0))

  size <- 0

  for (i in 0:.m) {
    for (j in 0:.n) {

      t  <- i + j
      c1 <- .df$c1[[t + 1]]
      c2 <- .df$c2[[t + 1]]
      g1 <- .df$gamma1[[t + 1]]
      g2 <- .df$gamma2[[t + 1]]

      # Expected rejection probability for the randomised test
      if (i < c1 | i > c2) {

        phi <- 1                        # strictly outside: always reject

      } else if (c1 == c2 & i == c1) {

        phi <- g1 + g2                  # single boundary point

      } else if (i == c1) {

        phi <- g1                       # lower boundary: randomise with gamma1

      } else if (i == c2) {

        phi <- g2                       # upper boundary: randomise with gamma2

      } else {

        phi <- 0                        # strictly inside: never reject

      }

      size <- size +
        phi *
        stats::dbinom(i, size = .m, prob = p0) *
        stats::dbinom(j, size = .n, prob = p1)
    }
  }

  return(size)

}
