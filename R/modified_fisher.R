#' The non-conservative, size-\eqn{\alpha} modified Fisher exact test
#'
#' Computes the non-conservative, size-\eqn{\alpha} modified Fisher exact test
#' for comparing two proportions \eqn{u/m} and \eqn{v/n}, testing
#' \eqn{H_0}: OR = \code{odds_ratio}. Returns an \code{htest} object
#' containing test-based two-sided p-values and \eqn{(1 - \alpha)} confidence
#' intervals for the odds ratio, which agree by construction. The test
#' maximises the actual size over the nuisance parameter \eqn{p_0} subject to
#' it remaining no greater than \eqn{\alpha}, making it less conservative than
#' the standard Fisher exact test while strictly controlling the Type I error
#' rate.
#'
#' @param u Number of successes observed in group 1 (out of \code{m} trials).
#'   No default.
#' @param m Number of trials in group 1. No default.
#' @param v Number of successes observed in group 2 (out of \code{n} trials).
#'   No default.
#' @param n Number of trials in group 2. No default.
#' @param odds_ratio The null hypothesis odds ratio \eqn{\theta_0}. No default.
#' @param alpha Nominal significance level \eqn{\alpha}. Defaults to 0.05.
#' @param precision Numerical precision for p-values, confidence limits, and
#'   size calculations. Defaults to 1e-03.
#' @param message Logical. If \code{TRUE}, prints progress messages during
#'   execution. Defaults to \code{FALSE}.
#' @param method Numerical method for maximising the local size over the
#'   nuisance parameter. Either \code{"zoom"} (default) or \code{"trust"}
#'   (trust-region method from the trust package).
#' @param maze Number of grid points evaluated at each zoom iteration.
#'   Defaults to 10.
#' @param zoom_iter Number of zoom iterations. Defaults to 6.
#' @param power Logical. If \code{FALSE}, skips the power calculation.
#'   Defaults to \code{TRUE}.
#' @param superiority Logical. If \code{TRUE}, power is computed only over
#'   tables where the observed rate in group 2 exceeds that in group 1.
#'   Defaults to \code{FALSE}.
#' @param power_at_pi1 Success probability in group 1 at which to evaluate
#'   power. Defaults to 0.5.
#' @param power_at_pi2 Success probability in group 2 at which to evaluate
#'   power. Defaults to 0.75.
#' @param conf_int Logical. If \code{FALSE}, skips the \eqn{(1 - \alpha)}
#'   two-sided confidence interval. Defaults to \code{TRUE}.
#' @param pvalue Logical. If \code{FALSE}, skips the two-sided test-based
#'   p-value. Defaults to \code{TRUE}.
#' @param local_size_data Logical. If \code{TRUE}, attaches a
#'   \code{local.size.data} data frame to the results for plotting the size of
#'   the test as a function of the nuisance parameter. Defaults to
#'   \code{FALSE}.
#'
#' @keywords non randomised randomized conservative fisher exact test modified
#' @rdname modified_fisher_exact_test
#' @export
#'
#' @return Returns an \code{htest} object.
#' @examples
#' \dontrun{
#' # Example here
#' }

modified_fisher_exact_test <- function(u, m, v, n, odds_ratio,
                                       alpha = 0.05,
                                       precision = 1e-03,
                                       message = FALSE,
                                       maze = 10, # Number of points at each iteration to select nuisance parameter with max size from
                                       method = "zoom", #c("zoom", "trust")
                                       zoom_iter = 6, # Number of iterations used for the zoom optimisation method.
                                       conf_int = TRUE,
                                       pvalue = TRUE,
                                       local_size_data = FALSE,
                                       power = TRUE,
                                       superiority = FALSE,
                                       power_at_pi1 = 0.5,
                                       power_at_pi2 = 0.75
){

  # RECREATE SAS PROGRAMME ======================================================

  # require:
  # biasurn
  # nlopt

  #require(BiasedUrn)
  #require(trust) # investigate this package for nloptr option

  t <- u+v

  # Data inputs:

  # u -- observed a (# of successes for column 1)
  # m -- observed column sum 1
  # v -- observed b (# of successes for column 2)
  # n -- observed column sum 2

  ## Check size of input table (original test) --------------------------------

  # Checks here

  # If size of table is > 2x2, there's an option to simulate p-value
  # (see original function)

  # If table is 2x2, proceed:

  ## Set some optional results to zero ------------------------------------------

  RESULT_pval <- NULL
  RESULT_estimate <- NULL
  RESULT_conf_int <- NULL

  plot_data <- NULL
  power_mfet <- NULL

  #DATANAME <- deparse(substitute(data))
  DATANAME <- paste0("u = ", u, ", v = ", v, ", m = ", m, ", n = ", n)
  METHOD <- paste0("Non-Conservative Size-\u03b1 Modified Fisher's Exact Test")

  ## Input values -------------------------------------------------------------

  ## Create support - which tables do we need to calculate probability for ----

  # Calculate all possible values for a, given the row and column totals:

  #lower <- max(.t-.n, 0)
  #upper <- min(.m, .t)
  #support <- lower:upper

  # HELPER FUNCTIONS ------------------------------------------------------------

  # Generic:
  # source("R/calc_exp_value.R")
  # source("R/find_gamma12.R")
  # source("R/construct_test_frame.R")
  # #source("R/random_fe_test.R")
  #
  # # Modified FE test:
  # source("R/mod_fe_test.R")
  # source("R/local_size.R")  # expand to also calc local size for other tests
  # source("R/mod_fe_size.R")
  # source("R/optimise_gamma0.R")
  # source("R/accept.R")      # expand to also output accept/reject for other tests
  # source("R/local_power.R") # expand to also calc local power for other tests

  # START MAIN FUNCTION =========================================================

  # Require:
  # m
  # n
  # alpha
  # u
  # t (u + v)

  # and:
  z <- c(u, t)

  df <- construct_test_frame(.odds_ratio = odds_ratio, .m = m, .n = n,
                             .alpha = alpha, .precision = precision,
                             .message = message)

  if(power == TRUE){

    or <- 1

    df_power <- construct_test_frame(.odds_ratio = or, .m = m, .n = n,
                                     .alpha = alpha, .precision = precision)

    gamma0 <- optimise_gamma0(.odds_ratio = or, .m = m, .n = n, .alpha = alpha,
                              .precision = precision, .method = method,
                              .maze = maze, .zoom_iter = zoom_iter)

    pi1 <- power_at_pi1
    pi2 <- power_at_pi2

    p <- c(pi1, pi2)

    power_mfet <- power_mfet(p, gamma0, or, .m = m, .n = n, .df = df_power,
                             .alpha = alpha, .precision = precision,
                             .superiority = superiority)*100

  }

  # Input values, with +0.5 adjustment if 0:

  a <- max(u, 0.5)
  b <- max(m-u, 0.5)
  c <- max(v, 0.5)
  d <- max(n-v, 0.5)

  or0 <- (a*d)/(b*c)
  RESULT_estimate <- or0
  #calc_expected_value(.odds_ratio = 1, m, n, t, precision)

  z_alpha <- -stats::qnorm(alpha/2)

  upper0 <- exp(log(or0) + z_alpha*sqrt(1/a+1/b+1/c+1/d))
  lower0 <- exp(log(or0) - z_alpha*sqrt(1/a+1/b+1/c+1/d))

  if(conf_int){

    # Starting bisecting upper limit

    f1 <- or0 + 0.75*(upper0-or0)
    f2 <- 1.25*upper0
    crit <- abs(f2-f1)

    if(b>0 & c>0){

      while(crit > precision){

        or <- (f1+f2)/2

        df2 <- construct_test_frame(.odds_ratio = or, .m = m, .n = n,
                                    .alpha = alpha, .precision = precision)

        answer <- accept(z, .odds_ratio = or, .m = m, .n = n, .df = df2,
                         .alpha = alpha, .precision = precision,
                         .method = method, .maze = maze, .zoom_iter = zoom_iter)

        if(answer == 1){ f1 <- or } else { f2 <- or }

        crit <- abs(f2-f1)

      }

      conf_int_upper <- f1

    } else {

      conf_int_upper <- Inf

    }

    # Starting bisecting lower limit:

    f1 <- 0.5*lower0
    f2 <- lower0+0.25*(or0-lower0)
    crit <- abs(f2-f1)

    precision_conf_int <- 0.1

    if(a>0 & d>0){

      while(crit > precision){

        or <- (f1+f2)/2

        df2 <- construct_test_frame(.odds_ratio = or, .m = m, .n = n,
                                    .alpha = alpha, .precision = precision)

        answer <- accept(z, .odds_ratio = or, .m = m, .n = n, .df = df2,
                         .alpha = alpha, .precision = precision,
                         .method = method, .maze = maze, .zoom_iter = zoom_iter)

        if(answer == 1){ f2 <- or } else { f1 <- or }

        crit <- abs(f2-f1)

      }

      conf_int_lower <- f2

    } else {

      conf_int_lower <- 0

    }

  } # End of conf_int == TRUE

  # Find optimal gamma0 for size data and for output later:

  opt_gamma0 <- optimise_gamma0(.odds_ratio = odds_ratio, .m = m, .n = n,
                                .alpha = alpha, .precision = precision,
                                .method = method, .maze = maze,
                                .zoom_iter = zoom_iter)

  if(local_size_data){

    plot_data <- data.frame(
      "pi1" = seq(0, 1, by = 1/100)
    )

    for(row in 1:101){

      point <- plot_data$pi1[[row]]

      plot_data$size[[row]] <- local_size(point, .gamma0 = opt_gamma0,
                                          .odds_ratio = odds_ratio,
                                          .m = m, .n = n, .df = df,
                                          .alpha = alpha,
                                          .precision = precision)*100

    }

    plot_data$size <- as.numeric(plot_data$size)
    plot_data$method <- method

  }

  if(pvalue){

    pval_lower <- 0
    pval_upper <- 1
    a0 <- 0.5
    crit <- abs(pval_upper-pval_lower)

    # p-value for H0 :

    while(crit > precision){

      df3 <- construct_test_frame(.odds_ratio = odds_ratio, .m = m, .n = n,
                                  .alpha = a0, .precision = precision)

      reject <- 1-accept(z, .odds_ratio = odds_ratio, .m = m, .n = n,
                         .df = df3, .alpha = a0, .precision = precision,
                         .method = method, .maze = maze, .zoom_iter = zoom_iter)

      if(reject == 1){
        pval_upper <- a0
        a0 <- (pval_lower+pval_upper)/2
      } else {
        pval_lower <- a0
        a0 <- (pval_lower+pval_upper)/2
      }

      crit <- pval_upper-pval_lower

    }

    RESULT_pval <- a0

  }

  ## Create output --------------------------------------------------------------

  RESULT_conf_int <- if(conf_int) c(conf_int_lower, conf_int_upper)
  attr(odds_ratio, "names") <- "odds ratio"
  attr(RESULT_estimate, "names") <- "odds ratio"
  attr(RESULT_conf_int, "conf.level") <- 1-alpha

  # Put results into list:

  RESULTS <- list(
    p.value = if(pvalue) RESULT_pval,
    estimate = RESULT_estimate,
    conf.int = if(conf_int) RESULT_conf_int,
    null.value = odds_ratio,
    alternative = "two.sided",
    method = METHOD,
    data.name = DATANAME,
    support.data = df,
    local.size.data = if(local_size_data) plot_data,
    gamma0 = opt_gamma0,
    power = if(power) power_mfet,
    fn_args = list(
      "alpha" = alpha,
      "precision" = precision,
      "maze" = maze,
      "method" = method,
      "zoom_iter" = zoom_iter,
      "superiority" = superiority,
      "power" = power,
      "power_at_pi1" = power_at_pi1,
      "power_at_pi2" = power_at_pi2
    )
  )

  attr(RESULTS, "class") <- "htest"

  return(RESULTS)

}
