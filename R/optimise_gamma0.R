#' Find the optimal \eqn{\gamma_0}
#'
#' Finds the optimal \eqn{\gamma_0} threshold: the largest value of
#' \eqn{\gamma_0} such that the actual size of the MFET (the maximum of the
#' local size over the nuisance parameter \eqn{p_0}) does not exceed
#' \eqn{\alpha}. Uses a bisection search over the \eqn{2(m + n + 1)} sorted
#' randomisation probability values from the test frame.
#'
#' @param .odds_ratio The null hypothesis odds ratio \eqn{\theta_0}. No default.
#' @param .m Number of trials in group 1.
#' @param .n Number of trials in group 2.
#' @param .alpha Nominal significance level \eqn{\alpha}. No default.
#' @param .precision Numerical precision. No default.
#' @param .method Numerical method for maximising the local size: \code{"zoom"}
#'   or \code{"trust"}.
#' @param .maze Number of grid points evaluated at each zoom iteration.
#' @param .zoom_iter Number of zoom iterations (\code{"zoom"} method only).
#'
#' @return A single numeric value: the optimal threshold \eqn{\gamma_0}.
#' @examples
#' optimise_gamma0(.odds_ratio = 1, .m = 6, .n = 4, .alpha = 0.05,
#'                 .precision = 1e-3, .method = "zoom", .maze = 10,
#'                 .zoom_iter = 6)
#' @export
#' @family modified
#' @seealso [modified_fisher_exact_test()] for the main user-facing function; [construct_test_frame()] for the test frame passed to this function; [size_modified()] for the size function being maximised.
optimise_gamma0 <- function(.odds_ratio, .m, .n, .alpha, .precision,
                            .method, .maze, .zoom_iter){

  df <- construct_test_frame(.odds_ratio, .m, .n, .alpha, .precision)

  gamind <- c(df$gamma1, df$gamma2)

  # Sort gammas:
  gamind <- sort(gamind)

  # Start bisecting on the index of gamind:
  g0 <- 0
  g1 <- 2*(.m+.n+1)
  g <- .m+.n+1

  sz0 <- size_modified(.c = gamind[[g]], .odds_ratio, .m, .n, .df = df,
                   .alpha, .precision, .method, .maze, .zoom_iter)

  # Initialise size_old before the while loop:
  size_old <- sz0

  while(as.integer(g1-g0) > 1){

    if(sz0 > .alpha){ g0 <- g } else { g1 <- g }
    if(sz0 < .alpha){

      size_old <- sz0
      # Only keep this old size if it is correct, old and new
      # can both be wrong

    }

    g <- as.integer((g0+g1)/2)
    sz0 <- size_modified(.c = gamind[[g]], .odds_ratio, .m, .n, .df = df,
                     .alpha, .precision, .method, .maze, .zoom_iter)

    if( g1-g0 == 1 ) {

      if( sz0 > .alpha ){
        g_opt <- g1
        size_opt <- size_old
      } else {
        g_opt <- g
        size_opt <- sz0
      }

      if( g == g0 & sz0 > .alpha ){

        g_opt <- g+1
        size_opt <- size_old

      }

    }

  }

  gamma0 <- gamind[[g_opt]]

  return(gamma0)

}
