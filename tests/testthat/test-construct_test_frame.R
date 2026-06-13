# Absolute-tolerance helper. testthat edition 3 interprets expect_equal()'s
# `tolerance` as relative, which is too tight for the small probability values
# in Table 1: a 1e-3 relative tolerance on gamma = 0.178 permits only ~1.8e-4
# absolute difference, smaller than the paper's own 3-decimal rounding (+/- 5e-4).
# Table 1 values are reported to 3 decimals, so an absolute tolerance of 1e-3
# safely covers both the paper rounding and the function's numerical precision.

expect_within <- function(object, expected, abs_tol = 1e-3) {
  expect_lt(abs(object - expected), abs_tol)
}

# Test - Check that our gammas are the same as Table 1 ------------------------

test_that("construct_test_frame matches Table 1 of van der Meulen et al. (2021)", {
  df <- modifiedfisher:::construct_test_frame(.odds_ratio = 1, .m = 6, .n = 4,
                                              .alpha = 0.05, .precision = 1e-3)

  # Critical values are integers: exact equality is appropriate
  expect_equal(df$c1[df$t == 5], 2)
  expect_equal(df$c2[df$t == 6], 5)

  # Randomisation probabilities: compare against Table 1 (3 dp) with an
  # absolute tolerance, not expect_equal()'s relative tolerance.
  expect_within(df$gamma1[df$t == 5], 0.005)
  expect_within(df$gamma2[df$t == 6], 0.178)

  # t = 0 edge case (gamma1 = gamma2 = alpha/2 exactly)
  expect_within(df$gamma1[df$t == 0], 0.025)
  expect_within(df$gamma2[df$t == 0], 0.025)

})

# Test - Boundary cases match closed-form results -----------------------------

test_that("boundary cases of t match closed-form results", {

  df <- modifiedfisher:::construct_test_frame(.odds_ratio = 1, .m = 6, .n = 4,
                                              .alpha = 0.05, .precision = 1e-3)

  # t = 0: gamma1 = gamma2 = alpha/2, c1 = c2 = 0
  expect_equal(df$c1[df$t == 0], 0)
  expect_equal(df$c2[df$t == 0], 0)
  expect_within(df$gamma1[df$t == 0], 0.025)

  # t = n+m: gamma1 = gamma2 = alpha/2, c1 = c2 = m
  last <- nrow(df)
  expect_equal(df$c1[last], 6)
  expect_within(df$gamma1[last], 0.025)

  # t = 1: gamma1 = gamma2 = alpha, c1 = 0, c2 = 1
  expect_equal(df$c1[df$t == 1], 0)
  expect_equal(df$c2[df$t == 1], 1)
  expect_within(df$gamma1[df$t == 1], 0.05)

})

# Output is stable on rerun ---------------------------------------------------

test_that("construct_test_frame output is stable (snapshot)", {

  skip_on_cran()

  df <- modifiedfisher:::construct_test_frame(.odds_ratio = 1, .m = 20, .n = 20,
                                              .alpha = 0.05, .precision = 1e-3)

  expect_snapshot(df)

})
