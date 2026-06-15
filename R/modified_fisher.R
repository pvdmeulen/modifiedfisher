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
#' @rdname modified_fisher_exact_test
#' @return An object of class \code{htest}: a list with components
#'   \code{p.value} (the two-sided test-based p-value, if \code{pvalue = TRUE}),
#'   \code{estimate} (the sample odds ratio), \code{conf.int} (the two-sided
#'   \eqn{(1 - \alpha)} confidence interval for the odds ratio, if
#'   \code{conf_int = TRUE}), \code{null.value} (the null odds ratio),
#'   \code{alternative}, \code{method}, and \code{data.name}. It additionally
#'   carries \code{support.data} (the test frame from
#'   \code{construct_test_frame()}), \code{gamma0} (the optimal threshold),
#'   \code{power} (the power at \code{power_at_pi1} and \code{power_at_pi2}, if
#'   \code{power = TRUE}), \code{local.size.data} (a data frame for plotting the
#'   size as a function of the nuisance parameter, if \code{local_size_data =
#'   TRUE}), and \code{fn_args} (the arguments used).
#' @family modified
#' @seealso [construct_test_frame()] for the underlying test frame of critical values and randomisation probabilities; [optimise_gamma0()] for the gamma0 optimisation; [size_modified()] for the size of the test maximised over the nuisance parameter; [local_size_modified()] for the local size at a fixed nuisance parameter value; [power_modified()] for the power of the test.
#' @export
#' @examples
#' # Example 1 of van der Meulen, Raymond & van der Meulen (2021):
#' # 5/12 successes versus 7/11, testing H0: OR = 1 at alpha = 0.05.
#' modified_fisher_exact_test(u = 5, m = 12, v = 7, n = 11, odds_ratio = 1)
#'
#' \donttest{
#' # The size of the test as a function of the nuisance parameter can be
#' # returned for diagnostic plotting (slower, as the size is evaluated on a
#' # grid of 101 nuisance-parameter values):
#' res <- modified_fisher_exact_test(u = 5, m = 12, v = 7, n = 11,
#'                                   odds_ratio = 1, local_size_data = TRUE)
#' head(res$local.size.data)
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

  # INPUT VALIDATION ==========================================================

  ## u, v, m , n --------------------------------------------------------------
  for (arg_name in c("u", "m", "v", "n")) {
    val <- get(arg_name)
    if (!is.numeric(val) || length(val) != 1L || is.na(val)) {
      stop(sprintf("`%s` must be a single non-missing numeric value.", arg_name))
    }
    if (val != round(val)) {
      stop(sprintf("`%s` must be an integer (got %s).", arg_name, val))
    }
    if (val < 0) {
      stop(sprintf("`%s` must be non-negative (got %s).", arg_name, val))
    }
  }

  if (m == 0 || n == 0) {
    stop("`m` and `n` must both be at least 1 (zero-size groups are not supported).")
  }

  if (u > m) stop(sprintf("`u` (%s) cannot exceed `m` (%s).", u, m))
  if (v > n) stop(sprintf("`v` (%s) cannot exceed `n` (%s).", v, n))

  ## Odds ratio ---------------------------------------------------------------

  if (!is.numeric(odds_ratio) || length(odds_ratio) != 1L || is.na(odds_ratio)) {
    stop("`odds_ratio` must be a single non-missing numeric value.")
  }
  if (odds_ratio <= 0 || !is.finite(odds_ratio)) {
    stop(sprintf("`odds_ratio` must be a finite positive number (got %s).", odds_ratio))
  }

  ## Alpha --------------------------------------------------------------------

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha)) {
    stop("`alpha` must be a single non-missing numeric value.")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop(sprintf("`alpha` must be strictly between 0 and 1 (got %s).", alpha))
  }

  ## Precision ----------------------------------------------------------------

  if (!is.numeric(precision) || length(precision) != 1L || is.na(precision)) {
    stop("`precision` must be a single non-missing numeric value.")
  }
  if (precision <= 0 || precision >= 0.1) {
    stop(sprintf("`precision` must be a small positive number, typically <= 1e-2 (got %s).", precision))
  }

  ## Method -------------------------------------------------------------------

  method <- match.arg(method, choices = c("zoom", "trust"))

  ## Zoom method arguments ----------------------------------------------------

  for (arg_name in c("maze", "zoom_iter")) {
    val <- get(arg_name)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val != round(val) || val < 1) {
      stop(sprintf("`%s` must be a single positive integer (got %s).", arg_name, val))
    }
  }

  if (maze < 3) {
    warning("`maze` < 3 may produce an unreliable zoom-in optimisation; the SAS macro default is 10.")
  }

  ## Power at pi1/pi2 ---------------------------------------------------------

  for (arg_name in c("power_at_pi1", "power_at_pi2")) {
    val <- get(arg_name)
    if (!is.numeric(val) || length(val) != 1L || is.na(val)) {
      stop(sprintf("`%s` must be a single non-missing numeric value.", arg_name))
    }
    if (val < 0 || val > 1) {
      stop(sprintf("`%s` must be between 0 and 1 (got %s).", arg_name, val))
    }
  }

  ## Logical arguments --------------------------------------------------------

  logical_args <- c("message", "power", "superiority", "conf_int", "pvalue", "local_size_data")
  for (arg_name in logical_args) {
    val <- get(arg_name)
    if (is.character(val) && toupper(val) %in% c("Y", "N")) {
      stop(sprintf(
        "`%s` must be TRUE/FALSE, not \"%s\". (Note: this package uses R logicals, not the SAS Y/N convention.)",
        arg_name, val
      ))
    }
    if (!is.logical(val) || length(val) != 1L || is.na(val)) {
      stop(sprintf("`%s` must be a single TRUE/FALSE value (got %s).", arg_name, val))
    }
  }

  ## Cross argument consistency -----------------------------------------------

  if (!power && (!missing(power_at_pi1) || !missing(power_at_pi2))) {
    warning("`power_at_pi1`/`power_at_pi2` were supplied but `power = FALSE`; they will be ignored.")
  }

  # START FUNCTION ============================================================

  # Data inputs:

  # u -- observed a (# of successes for column 1)
  # m -- observed column sum 1
  # v -- observed b (# of successes for column 2)
  # n -- observed column sum 2

  # Total successes:
  t <- u+v

  # and set:
  z <- c(u, t)

  ## Initialise results objects -----------------------------------------------

  RESULT_pval <- NULL
  RESULT_estimate <- NULL
  RESULT_conf_int <- NULL

  plot_data <- NULL
  power_modified <- NULL

  #DATANAME <- deparse(substitute(data))
  DATANAME <- paste0("u = ", u, ", v = ", v, ", m = ", m, ", n = ", n)
  METHOD <- paste0("Non-Conservative Size-\u03b1 Modified Fisher's Exact Test")

  ## Construct test frame -----------------------------------------------------

  df <- construct_test_frame(.odds_ratio = odds_ratio, .m = m, .n = n,
                             .alpha = alpha, .precision = precision,
                             .message = message)

  # If power is needed, calculate:

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

    power_modified <- power_modified(p, gamma0, or, .m = m, .n = n, .df = df_power,
                             .alpha = alpha, .precision = precision,
                             .superiority = superiority)*100

  }

  ## Store OR point estimate --------------------------------------------------

  # Input values, with +0.5 adjustment if 0:

  a <- max(u, 0.5)
  b <- max(m-u, 0.5)
  c <- max(v, 0.5)
  d <- max(n-v, 0.5)

  or0 <- (a*d)/(b*c)
  RESULT_estimate <- or0

  z_alpha <- -stats::qnorm(alpha/2)

  upper0 <- exp(log(or0) + z_alpha*sqrt(1/a+1/b+1/c+1/d))
  lower0 <- exp(log(or0) - z_alpha*sqrt(1/a+1/b+1/c+1/d))

  ## MFET - Calculate OR confidence intervals ---------------------------------

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

        answer <- .accept(z, .odds_ratio = or, .m = m, .n = n, .df = df2,
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

        answer <- .accept(z, .odds_ratio = or, .m = m, .n = n, .df = df2,
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

  ## MFET - find optimal gamma ------------------------------------------------

  # Find optimal gamma0 for size data and for output later:

  opt_gamma0 <- optimise_gamma0(.odds_ratio = odds_ratio, .m = m, .n = n,
                                .alpha = alpha, .precision = precision,
                                .method = method, .maze = maze,
                                .zoom_iter = zoom_iter)

  ## MFET - store size data ---------------------------------------------------

  if(local_size_data){

    plot_data <- data.frame(
      "pi1" = seq(0, 1, by = 1/100)
    )

    for(row in 1:101){

      point <- plot_data$pi1[[row]]

      plot_data$size[[row]] <- local_size_modified(point, .gamma0 = opt_gamma0,
                                               .odds_ratio = odds_ratio,
                                               .m = m, .n = n, .df = df,
                                               .alpha = alpha,
                                               .precision = precision)*100

    }

    plot_data$size <- as.numeric(plot_data$size)
    plot_data$method <- method

  }

  ## MFET - Calculate p-value -------------------------------------------------

  if(pvalue){

    pval_lower <- 0
    pval_upper <- 1
    a0 <- 0.5
    crit <- abs(pval_upper-pval_lower)

    # p-value for H0 :

    while(crit > precision){

      df3 <- construct_test_frame(.odds_ratio = odds_ratio, .m = m, .n = n,
                                  .alpha = a0, .precision = precision)

      reject <- 1-.accept(z, .odds_ratio = odds_ratio, .m = m, .n = n,
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
  if (conf_int) attr(RESULT_conf_int, "conf.level") <- 1 - alpha

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
    power = if(power) power_modified,
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
