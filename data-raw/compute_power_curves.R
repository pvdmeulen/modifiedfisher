# WRAPPER
# Compute power curves for all three tests over a grid of pi2 values

compute_power_curves <- function(m, n, pi1 = 0.5, alpha = 0.05, or = 1,
                                 precision = 1e-3,
                                 pi2_grid  = seq(0.51, 0.99, by = 0.01),
                                 method    = "zoom", maze = 10, zoom_iter = 6,
                                 include_randomised = FALSE) {

  message(
    paste0(
      "Building test frame and finding optimal gamma0 (OR = ",
      or,
      ")")
  )

  df <- modifiedfisher:::construct_test_frame(
    .odds_ratio = or, .m = m, .n = n, .alpha = alpha, .precision = precision)

  gamma0_opt <- modifiedfisher:::optimise_gamma0(
    .odds_ratio = or, .m = m, .n = n, .alpha = alpha, .precision = precision,
    .method = method, .maze = maze, .zoom_iter = zoom_iter)

  results <- data.frame(pi2 = pi2_grid)

  message("Computing power: modified FET...")
  results$modified <- sapply(pi2_grid, function(pi2)
    modifiedfisher:::power_modified(
      c(pi1, pi2), .gamma0 = gamma0_opt, .odds_ratio = or,
      .m = m, .n = n, .df = df, .alpha = alpha, .precision = precision,
      .superiority = FALSE) * 100)

  message("Computing power: Woolf asymptotic...")
  results$woolf <- sapply(pi2_grid, function(pi2)
    modifiedfisher:::power_asymptotic(c(pi1, pi2), .m = m, .n = n, .alpha = alpha) * 100)

  message("Computing power: SAS Proc FREQ exact...")
  results$sas_freq <- sapply(pi2_grid, function(pi2)
    modifiedfisher:::power_probability(c(pi1, pi2), .m = m, .n = n, .alpha = alpha) * 100)

  message("Computing power: conservative FET...")
  results$conservative <- sapply(pi2_grid, function(pi2)
    modifiedfisher:::power_conservative(c(pi1, pi2), .m = m, .n = n, .df = df,
                       .alpha = alpha, .precision = precision) * 100)

  if (include_randomised) {
    message("Computing power: randomised (UMPU) FET...")
    results$randomised <- sapply(pi2_grid, function(pi2)
      modifiedfisher:::power_randomised(c(pi1, pi2), .m = m, .n = n, .df = df,
                       .alpha = alpha, .precision = precision) * 100)
  }

  attr(results, "pi1")   <- pi1
  attr(results, "m")     <- m
  attr(results, "n")     <- n
  attr(results, "alpha") <- alpha

  return(results)

}
