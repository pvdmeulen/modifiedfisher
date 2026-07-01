# .build_rejection_matrix ------------------------------------------------------

test_that(".build_rejection_matrix agrees with .modified_reject() for all cells", {

  df_small <- modifiedfisher:::construct_test_frame(.odds_ratio = 1, .m = 6, .n = 4,
                                                    .alpha = 0.05, .precision = 1e-3)
  g0 <- 0.03

  R <- modifiedfisher:::.build_rejection_matrix(df_small, g0, 6, 4)

  expect_equal(dim(R), c(7L, 5L))

  for (u in 0:6) {
    for (v in 0:4) {

      expect_equal(
        R[u + 1L, v + 1L],
        as.double(modifiedfisher:::.modified_reject(c(u, u + v), df_small, g0)),
        info = sprintf("u = %d, v = %d", u, v)
      )

    }
  }
})
