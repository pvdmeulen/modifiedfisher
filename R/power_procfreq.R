#' Find the power of the SAS Proc Freq exact test -----------------------------
#'
#' For each (u, v) pair the rejection decision is made using the SAS Proc FREQ
#' Fisher exact p-value (sum of central-hypergeometric probabilities, given the
#' marginal total t = u + v, that are no greater than the observed probability).
#' This mirrors the way local_size_procfreq() determines rejection, and is the
#' unconditional power of SAS Proc FREQ's exact test for H0: OR = 1.
#'
#' @param p Vector containing (pi1, pi2).
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .alpha The nominal significance level α. No default.
#' @param .superiority A logical. Defaults to FALSE. Setting this to TRUE will calculate the power for testing superiority.
#'
#' @keywords find power test sas procfreq

power_procfreq <- function(p, .m, .n, .alpha, .superiority = FALSE) {

  p0    <- p[[1]]
  p1    <- p[[2]]
  power <- 0

  for (u in 0:.m) {
    for (v in 0:.n) {

      # Rejection rule: Fisher exact p-value under H0: OR = 1 (central hypergeometric)
      t      <- u + v
      p_val  <- sas_procfreq_pvalue(u, t, .m, .n)
      reject <- as.numeric(p_val <= .alpha)

      contrib <- reject *
        stats::dbinom(u, size = .m, prob = p0) *
        stats::dbinom(v, size = .n, prob = p1)

      power <- power + if (.superiority) (v/.n > u/.m) * contrib else contrib

    }
  }

  return(power)

}
