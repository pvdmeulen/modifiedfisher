test_that("boundary cases of t match closed-form results", {

  # Check that boundary cases match the closed form results in paper:

  df <- modifiedfisher:::construct_test_frame(.odds_ratio = 1, .m = 6, .n = 4,
                                              .alpha = 0.05, .precision = 1e-3)

  # t = 0: gamma1 = gamma2 = alpha/2, c1 = c2 = 0
  expect_equal(df$c1[df$t == 0], 0)
  expect_equal(df$c2[df$t == 0], 0)
  expect_equal(df$gamma1[df$t == 0], 0.025)

  # t = n+m: gamma1 = gamma2 = alpha/2, c1 = c2 = m
  last <- nrow(df)
  expect_equal(df$c1[last], 6)
  expect_equal(df$gamma1[last], 0.025)

  # t = 1: gamma1 = gamma2 = alpha, c1=0, c2=1
  expect_equal(df$c1[df$t == 1], 0)
  expect_equal(df$c2[df$t == 1], 1)
  expect_equal(df$gamma1[df$t == 1], 0.05)

})


test_that(".find_gamma12 returns gamma in [0,1]^2 or (-1,-1)", {

  g <- modifiedfisher:::.find_gamma12(c1 = 2, c2 = 4, .m = 6, .n = 4, .t = 5,
                                     .odds_ratio = 1, .alpha = 0.05, .precision = 1e-3)
  expect_length(g, 2)
  expect_true(all(g >= -1 & g <= 1))

})
