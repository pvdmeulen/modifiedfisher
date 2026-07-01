# Tests for comparison functions ----------------------------------------------

# pvalue_probability
# local_size_randomised
# power_asymptotic
# power_conservative
# power_probability
# power_randomised

# These functions are used internally by the vignettes and pkgdown site, but are
# not essential to the modified_fisher_exact_test() function itself.

## Shared setup ---------------------------------------------------------------

m_s <- 12L; n_s <- 11L; alpha_s <- 0.05

df_s <- construct_test_frame(.odds_ratio = 1, .m = m_s, .n = n_s,
                             .alpha = alpha_s, .precision = 1e-3)

g0_s <- optimise_gamma0(.odds_ratio = 1, .m = m_s, .n = n_s,
                        .alpha = alpha_s, .precision = 1e-3,
                        .method = "zoom", .maze = 10, .zoom_iter = 6)

## pvalue_probability ---------------------------------------------------------

test_that("pvalue_probability returns value in [0, 1]", {
  pval <- pvalue_probability(s = 5, t = 12, m = 12, n = 11, odds_ratio = 1)
  expect_gte(pval, 0)
  expect_lte(pval, 1)
})

test_that("pvalue_probability is 1 when the support has only one point", {
  # t = 0: only one possible table (u = 0, v = 0)
  pval <- pvalue_probability(s = 0, t = 0, m = 6, n = 4, odds_ratio = 1)
  expect_equal(pval, 1)
})

test_that("pvalue_probability for OR != 1 still lies in [0, 1]", {
  pval <- pvalue_probability(s = 5, t = 12, m = 12, n = 11, odds_ratio = 2)
  expect_gte(pval, 0)
  expect_lte(pval, 1)
})

## local_size_randomised ------------------------------------------------------

test_that("local_size_randomised equals alpha at every nuisance value (UMPU property)", {
  # The fully randomised UMPU test has exact unconditional size alpha for all p0.
  sizes <- sapply(seq(0.1, 0.9, by = 0.1), function(p0) {
    local_size_randomised(nuisance = p0, .m = m_s, .n = n_s, .df = df_s,
                          .odds_ratio = 1, .alpha = alpha_s, .precision = 1e-3)
  })
  expect_true(all(abs(sizes - alpha_s) < 1e-8))
})

test_that("local_size_randomised >= local_size_modified (randomised is more powerful)", {
  g0 <- g0_s
  sizes_rand <- sapply(seq(0.1, 0.9, by = 0.2), function(p0) {
    local_size_randomised(nuisance = p0, .m = m_s, .n = n_s, .df = df_s,
                          .odds_ratio = 1, .alpha = alpha_s, .precision = 1e-3)
  })
  sizes_mod <- sapply(seq(0.1, 0.9, by = 0.2), function(p0) {
    local_size_modified(p0, .gamma0 = g0, .odds_ratio = 1,
                        .m = m_s, .n = n_s, .df = df_s,
                        .alpha = alpha_s, .precision = 1e-3)
  })
  expect_true(all(sizes_rand >= sizes_mod - 1e-10))
})

## power_asymptotic -----------------------------------------------------------

test_that("power_asymptotic at null (pi1 = pi2) does not exceed alpha by much", {
  # Woolf's test is not exact, but for n = 20 the size should be near alpha
  p_null <- power_asymptotic(p = c(0.3, 0.3), .m = 20, .n = 20, .alpha = alpha_s)
  expect_lte(p_null, alpha_s + 0.02)
})

test_that("power_asymptotic is non-decreasing in pi2 for fixed pi1", {
  powers <- sapply(seq(0.2, 0.9, by = 0.1), function(pi2) {
    power_asymptotic(p = c(0.2, pi2), .m = 20, .n = 20, .alpha = alpha_s)
  })
  expect_true(all(diff(powers) >= -1e-8))
})

test_that("power_asymptotic with superiority = TRUE is <= superiority = FALSE", {
  p1 <- power_asymptotic(p = c(0.2, 0.6), .m = m_s, .n = n_s,
                         .alpha = alpha_s, .superiority = FALSE)
  p2 <- power_asymptotic(p = c(0.2, 0.6), .m = m_s, .n = n_s,
                         .alpha = alpha_s, .superiority = TRUE)
  expect_lte(p2, p1 + 1e-10)
})

## power_conservative ---------------------------------------------------------

test_that("power_conservative does not exceed power_modified", {
  # Conservative test (gamma0 = 1) rejects a subset of what the MFET rejects.
  p_cons <- power_conservative(p = c(0.3, 0.6), .m = m_s, .n = n_s,
                               .df = df_s, .alpha = alpha_s, .precision = 1e-3)
  p_mod  <- power_modified(p = c(0.3, 0.6), .gamma0 = g0_s, .odds_ratio = 1,
                           .m = m_s, .n = n_s, .df = df_s,
                           .alpha = alpha_s, .precision = 1e-3,
                           .superiority = FALSE)
  expect_lte(p_cons, p_mod + 1e-10)
})

test_that("power_conservative with superiority = TRUE is <= superiority = FALSE", {
  p1 <- power_conservative(p = c(0.2, 0.6), .m = m_s, .n = n_s,
                           .df = df_s, .alpha = alpha_s, .precision = 1e-3,
                           .superiority = FALSE)
  p2 <- power_conservative(p = c(0.2, 0.6), .m = m_s, .n = n_s,
                           .df = df_s, .alpha = alpha_s, .precision = 1e-3,
                           .superiority = TRUE)
  expect_lte(p2, p1 + 1e-10)
})

## power_probability ----------------------------------------------------------

test_that("power_probability at null does not greatly exceed alpha", {
  p_null <- power_probability(p = c(0.3, 0.3), .m = 20, .n = 20, .alpha = alpha_s)
  expect_lte(p_null, alpha_s + 0.01)
})

test_that("power_probability is non-decreasing in pi2 for fixed pi1", {
  powers <- sapply(seq(0.2, 0.9, by = 0.1), function(pi2) {
    power_probability(p = c(0.2, pi2), .m = m_s, .n = n_s, .alpha = alpha_s)
  })
  expect_true(all(diff(powers) >= -1e-8))
})

test_that("power_probability with superiority = TRUE is <= superiority = FALSE", {
  p1 <- power_probability(p = c(0.2, 0.6), .m = m_s, .n = n_s,
                          .alpha = alpha_s, .superiority = FALSE)
  p2 <- power_probability(p = c(0.2, 0.6), .m = m_s, .n = n_s,
                          .alpha = alpha_s, .superiority = TRUE)
  expect_lte(p2, p1 + 1e-10)
})

## power_randomised -----------------------------------------------------------

test_that("power_randomised is exactly alpha at null (pi1 = pi2)", {
  p_null <- power_randomised(p = c(0.3, 0.3), .m = m_s, .n = n_s,
                             .df = df_s, .alpha = alpha_s, .precision = 1e-3)
  expect_equal(p_null, alpha_s, tolerance = 1e-8)
})

test_that("power_randomised >= power_modified (UMPU dominates MFET)", {
  p_rand <- power_randomised(p = c(0.3, 0.6), .m = m_s, .n = n_s,
                             .df = df_s, .alpha = alpha_s, .precision = 1e-3)
  p_mod  <- power_modified(p = c(0.3, 0.6), .gamma0 = g0_s, .odds_ratio = 1,
                           .m = m_s, .n = n_s, .df = df_s,
                           .alpha = alpha_s, .precision = 1e-3,
                           .superiority = FALSE)
  expect_gte(p_rand, p_mod - 1e-10)
})

test_that("power_randomised with superiority = TRUE is <= superiority = FALSE", {
  p1 <- power_randomised(p = c(0.2, 0.6), .m = m_s, .n = n_s,
                         .df = df_s, .alpha = alpha_s, .precision = 1e-3,
                         .superiority = FALSE)
  p2 <- power_randomised(p = c(0.2, 0.6), .m = m_s, .n = n_s,
                         .df = df_s, .alpha = alpha_s, .precision = 1e-3,
                         .superiority = TRUE)
  expect_lte(p2, p1 + 1e-10)
})
