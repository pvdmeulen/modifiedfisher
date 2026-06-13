# Test - power is non-decreasing ----------------------------------------------

test_that("mfet_power is monotonically non-decreasing in pi2 for fixed pi1", {

  m <- 12
  n <- 11
  alpha <- 0.05
  odds_ratio <- 1

  df <- modifiedfisher:::construct_test_frame(.odds_ratio = odds_ratio, .m = m, .n = n,
                                              .alpha = alpha, .precision = 1e-3)

  g0 <- modifiedfisher:::optimise_gamma0(.odds_ratio = odds_ratio, .m = m, .n = n,
                                         .alpha = alpha, .precision = 1e-3,
                                         .method = "zoom", .maze = 19, .zoom_iter = 6)

  pi1 <- 0.2
  pi2_grid <- seq(0.2, 0.9, by = 0.05)

  powers <- sapply(pi2_grid, function(pi2) {

    modifiedfisher:::power_mfet(p = c(pi1, pi2), .gamma0 = g0, .odds_ratio = odds_ratio,
                                .m = m, .n = n, .df = df, .alpha = alpha, .precision = 1e-3,
                                .superiority = FALSE)

  })

  # Power should never decrease as pi2 moves away from pi1
  expect_true(all(diff(powers) >= -1e-8))

  # And power at pi2 = pi1 should be close to alpha (size), not power
  power_at_null <- modifiedfisher:::power_mfet(p = c(pi1, pi1), .gamma0 = g0,
                                               .odds_ratio = odds_ratio,
                                               .m = m, .n = n, .df = df,
                                               .alpha = alpha, .precision = 1e-3,
                                               .superiority = FALSE)
  expect_lte(power_at_null, alpha + 1e-6)
})

# Test - power increases with sample size -------------------------------------

test_that("power_mfet increases with sample size for a fixed effect size", {
  alpha <- 0.05
  odds_ratio <- 1
  pi1 <- 0.2
  pi2 <- 0.5

  sample_sizes <- list(c(m = 10, n = 10), c(m = 20, n = 20), c(m = 40, n = 40))

  powers <- sapply(sample_sizes, function(sz) {
    m <- sz[["m"]]; n <- sz[["n"]]
    df <- modifiedfisher:::construct_test_frame(.odds_ratio = odds_ratio, .m = m, .n = n,
                                                .alpha = alpha, .precision = 1e-3)

    g0 <- modifiedfisher:::optimise_gamma0(.odds_ratio = odds_ratio, .m = m, .n = n,
                                           .alpha = alpha, .precision = 1e-3,
                                           .method = "zoom", .maze = 19, .zoom_iter = 6)

    modifiedfisher:::power_mfet(p = c(pi1, pi2), .gamma0 = g0, .odds_ratio = odds_ratio,
                                .m = m, .n = n, .df = df, .alpha = alpha, .precision = 1e-3,
                                .superiority = FALSE)

  })

  expect_true(all(diff(powers) >= -1e-8))

})
