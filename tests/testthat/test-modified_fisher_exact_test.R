# Test - MFET invalid argument rejection --------------------------------------

test_that("modified_fisher_exact_test validates inputs", {
  # counts
  expect_error(modified_fisher_exact_test(u = -1, m = 10, v = 3, n = 10,
                                          odds_ratio = 1), "non-negative")
  expect_error(modified_fisher_exact_test(u = 11, m = 10, v = 3, n = 10,
                                          odds_ratio = 1), "cannot exceed")
  expect_error(modified_fisher_exact_test(u = 2.5, m = 10, v = 3, n = 10,
                                          odds_ratio = 1), "integer")
  expect_error(modified_fisher_exact_test(u = 5, m = 0, v = 3, n = 10,
                                          odds_ratio = 1), "at least 1")

  # odds_ratio
  expect_error(modified_fisher_exact_test(u = 5, m = 10, v = 3, n = 10, odds_ratio = 0), "positive")
  expect_error(modified_fisher_exact_test(u = 5, m = 10, v = 3, n = 10, odds_ratio = -1), "positive")
  expect_error(modified_fisher_exact_test(u = 5, m = 10, v = 3, n = 10, odds_ratio = Inf), "finite")

  # alpha
  expect_error(modified_fisher_exact_test(u = 5, m = 10, v = 3, n = 10,
                                          odds_ratio = 1, alpha = 0), "between 0 and 1")

  expect_error(modified_fisher_exact_test(u = 5, m = 10, v = 3, n = 10,
                                          odds_ratio = 1, alpha = 1), "between 0 and 1")

  # method
  expect_error(modified_fisher_exact_test(u = 5, m = 10, v = 3, n = 10,
                                          odds_ratio = 1, method = "newton"))

  # logical / SAS-style flags
  expect_error(modified_fisher_exact_test(u = 5, m = 10, v = 3, n = 10,
                                          odds_ratio = 1, power = "Y"), "TRUE/FALSE")

  # power_at_pi range
  expect_error(modified_fisher_exact_test(u = 5, m = 10, v = 3, n = 10,
                                          odds_ratio = 1, power_at_pi1 = 1.5), "between 0 and 1")
})

# Test - MFET pval and confidence interval agree ------------------------------

test_that("p-value and CI agree: theta0=1 rejected iff 1 is outside the CI", {

  res <- modified_fisher_exact_test(u = 13, m = 41, v = 6, n = 47,
                                    odds_ratio = 1)

  rejects_null <- res$p.value < 0.05
  one_outside_ci <- !(1 >= res$conf.int[1] && 1 <= res$conf.int[2])

  expect_equal(rejects_null, one_outside_ci)

})

# Test - conf_int = FALSE works -----------------------------------------------

test_that("conf_int = FALSE returns without error (regression test)", {
  res <- modified_fisher_exact_test(u = 5, m = 12, v = 7, n = 11,
                                    odds_ratio = 1,
                                    power = FALSE, conf_int = FALSE)
  expect_null(res$conf.int)
  expect_true(!is.null(res$p.value))
})

# Test - MFET symmetry and edge cases -----------------------------------------

test_that("swapping groups inverts the odds ratio", {
  res1 <- modified_fisher_exact_test(u = 13, m = 41, v = 6, n = 47,
                                     odds_ratio = 1)
  res2 <- modified_fisher_exact_test(u = 6, m = 47, v = 13, n = 41,
                                     odds_ratio = 1)

  expect_equal(unname(res1$estimate), 1 / unname(res2$estimate), tolerance = 1e-3)
  expect_equal(res1$p.value, res2$p.value, tolerance = 1e-3)
})

test_that("handles u = 0 or v = 0 (zero cells)", {
  expect_no_error(modified_fisher_exact_test(u = 0, m = 10, v = 3, n = 10,
                                             odds_ratio = 1))
  expect_no_error(modified_fisher_exact_test(u = 10, m = 10, v = 0, n = 10,
                                             odds_ratio = 1))
})

test_that("handles t = 0 and t = m+n (no successes / all successes)", {
  expect_no_error(modified_fisher_exact_test(u = 0, m = 5, v = 0, n = 5,
                                             odds_ratio = 1))
  expect_no_error(modified_fisher_exact_test(u = 5, m = 5, v = 5, n = 5,
                                             odds_ratio = 1))
})

# Test - MFET results match paper results from Table 2 ------------------------

# Validation against Table 2 of original paper.

# ERRATA IN THE GALLEY PROOF (verified via closed-form Woolf confidence limits
# and base R's fisher.test(), which require no package code):
#
#   - Example 1 (5/12 vs 7/11): all three LOWER confidence limits in the table
#     lost a leading zero in typesetting. Woolf's closed-form lower limit is
#     0.076 (printed as 0.760); fisher.test gives 0.055 (printed as 0.550);
#     the modified FE lower limit is 0.0716 (printed as 0.716). Note the
#     printed MFET row is internally impossible: a test-based CI of
#     (0.716, 2.210) would exclude its own point estimate of 0.408.
#
#   - Example 4: the data header "72/128 versus 58/142" is a typo for
#     "71/128 versus 58/142". With u = 71, Woolf's closed form gives a point
#     estimate of 1.113 and 95% CI of 1.804 / 2.925, matching the printed row
#     digit-for-digit, and fisher.test matches the printed Proc FREQ row. With
#     u = 72 as printed, Woolf gives 1.148 / 1.862 / 3.020, which matches
#     neither.
#
#   - Examples 2 and 3: printed values are correct (closed-form Woolf rows
#     reproduce exactly).

# Tolerances: the package's bisection precision (default 1e-3 on the OR scale)
# means agreement with the paper's 3-4 significant figures is expected to
# ~1e-2 on CI limits and ~1e-3 on p-values.

# Note - tolerances are absolute via expect_within() below. Rationale:

# The p-value bisection terminates on a bracket of width `precision` (default
# 1e-3) and returns its midpoint, so package p-values carry +/- 5e-4; the paper
# rounds to 3 decimals (+/- 5e-4). Combined: 1.5e-3.

# CI limits are bisected to within `precision` on the OR scale (+/- 1e-3); the
# paper rounds to 3-4 significant figures (up to +/- 5e-3 for values like
# 10.25). Combined: 1e-2 is a safe absolute bound.

# Do not switch these to expect_equal(tolerance =): testthat's tolerance is
# RELATIVE, so small p-values would need misleadingly large values (e.g. a
# 3e-4 absolute difference on p = 0.034 is a 0.9% relative difference).

  expect_within <- function(object, expected, abs_tol) {
    expect_lt(abs(object - expected), abs_tol)
  }

test_that("Results for (5/12 vs 7/11) match Table 2 Example 1", {

  #skip_on_cran()

  res <- modified_fisher_exact_test(u = 5, m = 12, v = 7, n = 11,
                                    odds_ratio = 1)

  expect_within(unname(res$estimate), 0.408, abs_tol = 1e-3)
  expect_within(res$p.value, 0.321, abs_tol = 1.5e-3)

  # Paper prints 0.716; correct value is 0.0716 (leading-zero typo, see header)
  expect_within(res$conf.int[1], 0.0716, abs_tol = 1e-2)
  expect_within(res$conf.int[2], 2.210, abs_tol = 1e-2)

})

test_that("Results for (13/41 vs 6/47) match Table 2 Example 2", {

  skip_on_cran()

  res <- modified_fisher_exact_test(u = 13, m = 41, v = 6, n = 47,
                                    odds_ratio = 1)

  expect_within(unname(res$estimate), 3.173, abs_tol = 1e-3)
  expect_within(res$p.value, 0.034, abs_tol = 1.5e-3)
  expect_within(res$conf.int[1], 1.082, abs_tol = 1e-2)
  # 10.25 is printed to 2 dp, so paper rounding alone contributes +/- 5e-3
  expect_within(res$conf.int[2], 10.25, abs_tol = 2e-2)

})

test_that("Results for (37/65 vs 28/71) match Table 2 Example 3", {

  skip_on_cran()

  res <- modified_fisher_exact_test(u = 37, m = 65, v = 28, n = 71,
                                    odds_ratio = 1)

  expect_within(unname(res$estimate), 2.029, abs_tol = 1e-3)
  expect_within(res$p.value, 0.0439, abs_tol = 1.5e-3)
  expect_within(res$conf.int[1], 1.023, abs_tol = 1e-2)
  expect_within(res$conf.int[2], 4.083, abs_tol = 1e-2)

})

test_that("Results for (71/128 vs 58/142) match Table 2 Example 4", {

  skip_on_cran()

  # Example 4: use the data the paper actually analysed (u=71, not 72)
  res <- modified_fisher_exact_test(u = 71, m = 128, v = 58, n = 142,
                                    odds_ratio = 1)

  expect_within(unname(res$estimate), 1.804, abs_tol = 1e-3)
  expect_within(res$p.value, 0.0175, abs_tol = 1.5e-3)
  expect_within(res$conf.int[1], 1.108, abs_tol = 1e-2)
  expect_within(res$conf.int[2], 2.936, abs_tol = 1e-2)

})
