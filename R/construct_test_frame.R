#' Construct testing frame and randomisation values
#'
#' Constructs a data frame of critical values (\eqn{c_1}, \eqn{c_2}) and
#' randomisation probabilities (\eqn{\gamma_1}, \eqn{\gamma_2}) for every
#' possible value of the total \eqn{T = 0, \ldots, m + n}, given the null
#' odds ratio, significance level \eqn{\alpha}, and precision. Starting from
#' the \eqn{\alpha/2} quantiles of the Fisher non-central hypergeometric
#' distribution, a spiral search over \eqn{(c_1, c_2)} is used whenever the
#' initial solution for \eqn{(\gamma_1, \gamma_2)} falls outside \eqn{[0, 1]}.
#'
#' @param .odds_ratio The null hypothesis odds ratio \eqn{\theta_0}. No default.
#' @param .m Number of trials in group 1.
#' @param .n Number of trials in group 2.
#' @param .alpha Nominal significance level \eqn{\alpha}. No default.
#' @param .precision Numerical precision for quantile calculations and
#'   \code{BiasedUrn::dFNCHypergeo()}. No default.
#' @param .message A logical. Defaults to \code{FALSE}. Setting this to
#'   \code{TRUE} will print progress messages; useful for debugging.
#'
#' @return A data frame with \code{m + n + 1} rows, one per possible total
#'   \eqn{T = 0, \ldots, m + n}, and columns \code{t} (the total), \code{c1}
#'   and \code{c2} (lower and upper critical values), \code{d1} and \code{d2}
#'   (the \eqn{\alpha/2} quantiles used as starting points), and \code{gamma1}
#'   and \code{gamma2} (the randomisation probabilities at \code{c1} and
#'   \code{c2}).
#' @examples
#' # Critical values and randomisation probabilities for m = 6, n = 4
#' # (reproduces Table 1 of van der Meulen et al., 2021):
#' construct_test_frame(.odds_ratio = 1, .m = 6, .n = 4,
#'                      .alpha = 0.05, .precision = 1e-3)
#' @export
#' @family modified
#' @seealso [modified_fisher_exact_test()] for the main user-facing function; [optimise_gamma0()] which uses this frame to find the optimal gamma0; [size_modified()] for the resulting test size.
#' @importFrom BiasedUrn qFNCHypergeo
construct_test_frame <- function(.odds_ratio, .m, .n, .alpha, .precision,
                                 .message = FALSE){

  # Set OR to arbitrarily small number if null is 0 (avoid errors)
  if(.odds_ratio < .precision){
    message(paste0("Odds ratio set to ", .precision, " for testing purposes"))
    or <- .precision
  } else { or <- .odds_ratio }

  # Test:
  #or <- 2

  df <- data.frame(
    "t" = rep(0, .m+.n+1),
    "c1" = rep(0, .m+.n+1),
    "d1" = rep(0, .m+.n+1),
    "gamma1" = rep(0, .m+.n+1),
    "c2" = rep(0, .m+.n+1),
    "d2" = rep(0, .m+.n+1),
    "gamma2" = rep(0, .m+.n+1)
  )

  # Loop version:

  # Test:
  #s <- 8

  for(s in 0:(.m+.n)){

    #message(paste0("Starting on s = ", s))
    df$t[[s+1]] <- s

    # If no successes, set gamma to chosen (two sided) .alpha:

    if(s == 0){

      #df$c1[[s+1]] <- s
      #df$c2[[s+1]] <- s

      df$gamma1[[s+1]] <- .alpha/2
      df$gamma2[[s+1]] <- .alpha/2

      if(.message){
        message(
          paste0("Solution for t = ", s, " found: ", .alpha/2, ", ", .alpha/2)
        )
      }

    }

    # If s > 0, set critical values to the relevant quantiles and find gamma
    # which mixes so that total probability mass is still 1-.alpha:

    if(s > 0){

      # Find critical values that are (conservatively) <= .alpha/2:

      d1 <- BiasedUrn::qFNCHypergeo(n = s, p = .alpha/2, m1 = .m, m2 = .n,
                                    odds = or, precision = .precision,
                                    lower.tail = TRUE)

      d2 <- BiasedUrn::qFNCHypergeo(n = s, p = .alpha/2, m1 = .m, m2 = .n,
                                    odds = or, precision = .precision,
                                    lower.tail = FALSE)

      df$d1[[s]] <- d1
      df$d2[[s]] <- d2

      # If total successes are 1, critical values are 0 and 1 and gamma is .alpha:

      # if(s == 1){
      #
      #   #df$c1[[s+1]] <- s-1
      #   #df$c2[[s+1]] <- s
      #
      #   #df$d1[[s]] <- 0
      #   #df$d2[[s]] <- 1
      #
      #   df$gamma1[[s+1]] <- .alpha
      #   df$gamma2[[s+1]] <- .alpha
      #
      #   if(message){
      #     message(
      #       paste0("Solution for t = ", s, " found: ", .alpha, ", ", .alpha)
      #     )
      #   }
      #
      # }

      # Find randomisation probabilities such that probability mass stays the
      # same:

      correct <- FALSE

      c1 <- d1
      c2 <- d2

      gamma <- .find_gamma12(.t = s, .odds_ratio = or, c1 = c1, c2 = c2,
                            .m = .m, .n = .n, .alpha = .alpha,
                            .precision = .precision)

      # Test if gamma is 0 < gamma < 1:

      correct <- (0 <= gamma[[1]] & gamma[[1]] <= 1 &
                    0 <= gamma[[2]] & gamma[[2]] <= 1)

      # If solution is admissible:

      if(correct){

        if(.message){
          message(
            paste0("Solution for t = ", s," found: ",
                   paste(round(gamma, 3), collapse = ", "))
          )
        }

        df$c1[[s+1]] <- c1
        df$c2[[s+1]] <- c2

        df$gamma1[[s+1]] <- round(gamma[[1]], digits = 7)
        df$gamma2[[s+1]] <- round(gamma[[2]], digits = 7)

        # If .alpha/2 quantiles are still both zero, set
        # gamma to .alpha:

        # } else if(correct == FALSE & d1 == d2){
        #
        #   df$c1[[s+1]] <- d1
        #   df$c2[[s+1]] <- d2+1
        #
        #   df$gamma1[[s+1]] <- .alpha
        #   df$gamma2[[s+1]] <- .alpha
        #
        #   if(message){
        #     message(
        #       paste0("Solution for t = ", s, " found: ", .alpha, ", ", .alpha)
        #     )
        #   }
        #
        #   correct <- TRUE

      } else {

        if(.message){message(paste0("Solution for t = ", s,
                                   " not found, starting iteration"))}

        # If not, search for one that is by looking around quantile square-wise
        # and expanding the square (kxk) until a solution is found:

        k <- 1

        while(correct == FALSE){

          for(i in max(0, d1-k):min(d1+k, min(s, .m))){

            for(j in max(0, d2-k):min(d2+k, min(s, .m))){

              if(max(d1-i, d2-j) == k | min(d1-i, d2-j) == -k){

                c1 <- i
                c2 <- j

                gamma <- .find_gamma12(.t = s, .odds_ratio = or,
                                      c1 = c1, c2 = c2, .m = .m, .n = .n,
                                      .alpha = .alpha, .precision = .precision)

                check <- (0 <= gamma[[1]] & gamma[[1]] <= 1 &
                            0 <= gamma[[2]] & gamma[[2]] <= 1)

                if(check == TRUE){

                  correct <- TRUE

                  df$c1[[s+1]] <- c1
                  df$c2[[s+1]] <- c2

                  df$gamma1[[s+1]] <- round(gamma[[1]], digits = 7)
                  df$gamma2[[s+1]] <- round(gamma[[2]], digits = 7)

                  if(.message){
                    message(paste0("Solution for t = ", s,
                                   " found at iteration ",
                                   k, ": ", paste(round(gamma, 3), collapse = ", "))
                    )
                  }

                } # End of check == TRUE

              } # End of min/max == k

            } # End of for j loop

          } # End of for i loop

          k <- k+1

          if(k == 1e4){
            message(paste0("Solution for t = ", s,
                           " not found, stopped iterating at k = ", k))
            correct <- TRUE
          }

        } # End of while correct == FALSE

        k <- 0

      } # End of iteration loop

    } # End of s > 1

  } # End of s loop

  df$d1[[.m+.n+1]] <- .m
  df$d2[[.m+.n+1]] <- .m

  return(df)

} # End of function
