test_that("calc_exp_value reduces to hypergeometric mean when OR = 1", {

  # E(S|T=t) for the central hypergeometric has a known closed form: m*t/(m+n)
  e_s <- modifiedfisher:::.calc_expected_value(.odds_ratio = 1, .m = 6, .n = 4,
                                         .t = 5, .precision = 1e-6)

  expect_equal(e_s, 6 * 5 / 10, tolerance = 1e-4)

})
