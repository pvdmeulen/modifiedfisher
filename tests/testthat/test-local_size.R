test_that("MFET size never exceeds alpha across a grid of nuisance parameters", {
  df <- construct_test_frame(.odds_ratio = 1, .m = 12, .n = 11,
                             .alpha = 0.05, .precision = 1e-3)

  g0 <- optimise_gamma0(.odds_ratio = 1, .m = 12, .n = 11,
                        .alpha = 0.05, .precision = 1e-3,
                        .method = "zoom", .maze = 19, .zoom_iter = 6)

  sizes <- sapply(seq(0.01, 0.99, by = 0.02), function(p0) {
    local_size_modified(p0, .gamma0 = g0, .odds_ratio = 1, .m = 12, .n = 11, .df = df,
                    .alpha = 0.05, .precision = 1e-3)
  })

  expect_true(all(sizes <= 0.05 + 1e-6))
})

test_that("sas_probability_pvalue size does not exceed alpha for OR != 1 (regression for null-mismatch bug)", {
  sizes <- sapply(seq(0.05, 0.95, by = 0.05), function(p0) {
    local_size_probability(p0, .odds_ratio = 2, .m = 12, .n = 11, .alpha = 0.05)
  })
  expect_true(all(sizes <= 0.05 + 1e-6))
})

test_that("Woolf test centred on theta0 controls size for OR != 1 (regression for null-mismatch bug)", {
  sizes <- sapply(seq(0.05, 0.95, by = 0.05), function(p0) {
    local_size_asymptotic(p0, .odds_ratio = 2, .m = 65, .n = 71, .alpha = 0.05)
  })
  # Woolf doesn't strictly control size, but it shouldn't be wildly above alpha for large n
  expect_true(all(sizes <= 0.10))
})
