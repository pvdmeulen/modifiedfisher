# .modified_reject() indexes its test frame by t + 1 (row t+1 holds the critical
# values and randomisation probabilities for total T = t). The stub frames
# below must therefore have a row for every t up to the value being tested,
# otherwise `.df$c1[[t+1]]` is out of bounds. The helper builds a full-length
# frame and places the values of interest at the target t.

make_stub_frame <- function(t_target, c1, c2, gamma1, gamma2, t_max = 10) {

  df <- data.frame(
    t      = 0:t_max,
    c1     = 0L,
    c2     = 0L,
    gamma1 = 0,
    gamma2 = 0
  )

  df$c1[t_target + 1]     <- c1
  df$c2[t_target + 1]     <- c2
  df$gamma1[t_target + 1] <- gamma1
  df$gamma2[t_target + 1] <- gamma2

  return(df)

}

test_that(".modified_reject covers all five branches", {

  df <- make_stub_frame(t_target = 5, c1 = 2, c2 = 4, gamma1 = 0.3, gamma2 = 0.7)

  # strictly interior -> accept
  expect_equal(modifiedfisher:::.modified_reject(c(3, 5), df, .gamma0 = 0.5), 0)
  # strictly exterior -> reject
  expect_equal(modifiedfisher:::.modified_reject(c(0, 5), df, .gamma0 = 0.5), 1)
  expect_equal(modifiedfisher:::.modified_reject(c(6, 5), df, .gamma0 = 0.5), 1)
  # at c1, gamma1 = 0.3 < gamma0 = 0.5 -> accept
  expect_equal(modifiedfisher:::.modified_reject(c(2, 5), df, .gamma0 = 0.5), 0)
  # at c2, gamma2 = 0.7 > gamma0 = 0.5 -> reject
  expect_equal(modifiedfisher:::.modified_reject(c(4, 5), df, .gamma0 = 0.5), 1)

})

test_that(".modified_reject handles c1 == c2", {

  df <- make_stub_frame(t_target = 5, c1 = 3, c2 = 3, gamma1 = 0.2, gamma2 = 0.4)

  # gamma1 + gamma2 = 0.6 > gamma0 = 0.5 -> reject
  expect_equal(modifiedfisher:::.modified_reject(c(3, 5), df, .gamma0 = 0.5), 1)
  # gamma1 + gamma2 = 0.6 < gamma0 = 0.7 -> accept
  expect_equal(modifiedfisher:::.modified_reject(c(3, 5), df, .gamma0 = 0.7), 0)

})
