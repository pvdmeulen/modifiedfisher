# Woolf's asymptotic 95% CI for the odds ratio of a 2x2 table.
# Cells: a = successes group 1, b = failures group 1,
#        c = successes group 2, d = failures group 2.

woolf_ci <- function(a, b, c, d, conf_level = 0.95) {
  or <- (a * d) / (b * c)
  se <- sqrt(1/a + 1/b + 1/c + 1/d)
  z  <- qnorm(1 - (1 - conf_level) / 2)
  c(lower = exp(log(or) - z * se),
    est   = or,
    upper = exp(log(or) + z * se))
}

cat("===========================================================\n")
cat(" Table 2, Example 1:  5/12 vs 7/11\n")
cat(" Claim: the printed lower CI limits dropped a leading zero\n")
cat("===========================================================\n")

# 5/12 vs 7/11  ->  a = 5, b = 7, c = 7, d = 4
w1 <- woolf_ci(5, 7, 7, 4)
f1 <- fisher.test(matrix(c(5, 7, 7, 4), nrow = 2, byrow = TRUE))

cat(sprintf("Woolf lower limit:      %.3f   (paper prints 0.760)\n", w1["lower"]))
cat(sprintf("Proc FREQ lower limit:  %.3f   (paper prints 0.550)\n", f1$conf.int[1]))
cat("\nThe printed modified-test interval (0.716, 2.210) is also impossible:\n")
cat(sprintf("  point estimate = %.3f, which lies OUTSIDE (0.716, 2.210).\n",
            w1["est"]))
cat("  A test-based confidence interval must contain its point estimate,\n")
cat("  so 0.716 should read 0.0716 (the leading zero was dropped).\n\n")

cat("===========================================================\n")
cat(" Table 2, Example 4:  printed '72/128 vs 58/142'\n")
cat(" Claim: the printed figures were computed with u = 71\n")
cat("===========================================================\n")

# Compare u = 72 (as printed) against u = 71 (hypothesised true value).
# 72/128 -> a = 72, b = 56;  71/128 -> a = 71, b = 57;  58/142 -> c = 58, d = 84
w4_72 <- woolf_ci(72, 56, 58, 84)
w4_71 <- woolf_ci(71, 57, 58, 84)
f4_71 <- fisher.test(matrix(c(71, 57, 58, 84), nrow = 2, byrow = TRUE))

cat("Paper Table 2 prints, for Example 4:\n")
cat("  estimate 1.804,  Woolf CI (1.113, 2.925),  Proc FREQ CI (1.082, 3.018)\n\n")

cat(sprintf("With u = 72 (as printed):  est %.3f,  Woolf CI (%.3f, %.3f)\n",
            w4_72["est"], w4_72["lower"], w4_72["upper"]))
cat(sprintf("With u = 71 (hypothesis):  est %.3f,  Woolf CI (%.3f, %.3f)\n",
            w4_71["est"], w4_71["lower"], w4_71["upper"]))
cat(sprintf("With u = 71, fisher.test:  Proc FREQ CI (%.3f, %.3f),  p = %.4f\n",
            f4_71$conf.int[1], f4_71$conf.int[2], f4_71$p.value))
cat("\nThe u = 71 figures match the printed table; the u = 72 figures do not.\n")
cat("So the data label '72/128' should read '71/128'.\n")
