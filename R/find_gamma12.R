#' Find randomisation probabilities \eqn{\gamma_1} and \eqn{\gamma_2}
#'
#' Solves the \eqn{2 \times 2} linear system for \eqn{(\gamma_1, \gamma_2)}
#' given critical values \eqn{(c_1, c_2)}, so that the test is conditionally
#' of size \eqn{\alpha} and unbiased (Neyman-structure conditions). Returns
#' \eqn{(-1, -1)} when no admissible solution exists, signalling
#' \code{construct_test_frame()} to expand the search square.
#'
#' @param c1 Lower critical value: \eqn{S \leq c_1} falls in the lower
#'   rejection region; \eqn{S = c_1} exactly is rejected with probability
#'   \eqn{\gamma_1}.
#' @param c2 Upper critical value: \eqn{S \geq c_2} falls in the upper
#'   rejection region; \eqn{S = c_2} exactly is rejected with probability
#'   \eqn{\gamma_2}.
#' @param .m Number of trials in group 1.
#' @param .n Number of trials in group 2.
#' @param .t Total number of successes \eqn{t = u + v}; the conditioning value
#'   in the non-central hypergeometric distribution.
#' @param .odds_ratio The null hypothesis odds ratio \eqn{\theta_0}. No default.
#' @param .alpha Nominal significance level \eqn{\alpha}. No default.
#' @param .precision Numerical precision. No default.
#'
#' @keywords construct randomisation values testing gamma
#' @importFrom BiasedUrn dFNCHypergeo

find_gamma <- function(c1, c2, .odds_ratio, .m, .n, .t, .alpha,
                       .precision){

  # requires .m = m, .n = n, .t (t), .odds_ratio = odds_ratio,
  # and critical values c1 and c2

  A <- matrix(0, nrow = 2, ncol = 2)
  b <- c(0, 0)
  #gamma <- c(-1, -1) # doesn't seem to be used below

  lower <- max(.t-.n, 0)
  upper <- min(.m, .t)

  if(.odds_ratio == 1){
    exp_value <- .m/(.m+.n)*.t
  } else {
    exp_value <- calc_expected_value(.odds_ratio, .m, .n,
                                     .t, .precision)
  }

  for(s in lower:upper){

    # If support (a) is leq c1, or geq c2, simply calculate
    # probability (except when OR is zero)

    if(s <= c1 | s >= c2){

      if(.odds_ratio < .precision){

        if(s < c1){ p <- 1 }
        if(s >= c1){ p <- 0 }

      } else {

        p <- BiasedUrn::dFNCHypergeo(n = .t, x = s, m1 = .m, m2 = .n,
                                     odds = .odds_ratio, precision = .precision)
      }

      if(s == c1){
        A[1, 1] <- A[1, 1] + p
        A[2, 1] <- A[2, 1] + s*p
      }

      if(s<c1 | s>c2){
        b[1] <- b[1] + p
        b[2] <- b[2] + s*p
      }

      if(s == c2){
        A[1, 2] <- A[1, 2] + p
        A[2, 2] <- A[2, 2] + s*p
      }

      c = c(.alpha, .alpha*exp_value)

      # When solution is found, det(A) is not equal to (c.a.) 0:

      if(det(A) <= -1e-07 | det(A) >= 1e-07){

        #solution <- solve(A) %*% (c-b) # same as below
        solution <- solve(A, c - b)

      } else {

        # Fail state when solution is not found:
        solution <- c(-1, -1)

      }

      # If this is the final solution, go to the next value up in
      # random_fe_test function.

      # Unless we're at the extreme values of .t (and no randomisation
      # between values can occur, just 'downrating' the last probabilities
      # to meet our particular probability mass .alpha/2):

      if(.t == 0 | .t == (.m + .n)){

        solution <- c(.alpha/2, .alpha/2)

      }

    } # End of s <= c1, s >= c2
  } # End of s for loop

  return(solution)

}
