# WRAPPER:
# Compute size curves for all three tests over a grid of p0 values

compute_size_curves <- function(m, n, alpha = 0.05, or = 1, precision = 1e-3,
                                p0_grid = seq(0.01, 0.99, by = 0.01),
                                method = "zoom", maze = 10, zoom_iter = 6,
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

  message(paste0("gamma0 = ", round(gamma0_opt, 6),
                 "  |  grid points: ", length(p0_grid)))

  results <- data.frame(p0 = p0_grid)

  message("Computing: modified FET...")
  results$modified <- sapply(p0_grid, function(p)
    modifiedfisher:::local_size_modified(
      p, .gamma0 = gamma0_opt, .odds_ratio = or,
      .m = m, .n = n, .df = df, .alpha = alpha, .precision = precision) * 100)

  message("Computing: Woolf asymptotic...")
  results$woolf <- sapply(p0_grid, function(p)
    modifiedfisher:::local_size_asymptotic(p, .odds_ratio = or, .m = m, .n = n,
                                      .alpha = alpha) * 100)

  message("Computing: SAS Proc FREQ exact...")
  results$sas_freq <- sapply(p0_grid, function(p)
    modifiedfisher:::local_size_probability(p, .odds_ratio = or, .m = m, .n = n,
                                         .alpha = alpha) * 100)

  message("Computing: conservative FET...")
  results$conservative <- sapply(p0_grid, function(p)
    modifiedfisher:::local_size_modified(
      p, .gamma0 = 1, .odds_ratio = or,
      .m = m, .n = n, .df = df, .alpha = alpha, .precision = precision) * 100)

  if (include_randomised) {
    message("Computing: randomised FET...")
    results$randomised <- sapply(p0_grid, function(p)
      modifiedfisher:::local_size_randomised(
        p, .odds_ratio = or, .m = m, .n = n,
        .df = df, .alpha = alpha, .precision = precision) * 100)
  }

  attr(results, "gamma0_opt") <- gamma0_opt
  attr(results, "m")          <- m
  attr(results, "n")          <- n
  attr(results, "alpha")      <- alpha

  return(results)

}
