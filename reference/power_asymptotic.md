# Power of Woolf's asymptotic test

Computes the unconditional power of Woolf's asymptotic Wald test for
\\H_0\\: OR = 1, at response rates \\(\pi_1, \pi_2)\\. Rejects \\H_0\\
when \\\|\log \hat{\theta}\| / SE(\log \hat{\theta}) \> z\_{\alpha/2}\\,
where \\\hat{\theta} = u(n-v) / ((m-u)v)\\. Uses the Haldane correction
(replacing zero cells with 0.5) for tables where one or more cells are
zero.

## Usage

``` r
power_asymptotic(p, .m, .n, .alpha, .superiority = FALSE)
```

## Arguments

- p:

  Length-2 numeric vector \\(\pi_1, \pi_2)\\: success probability in
  group 1 (\\\pi_1\\) and group 2 (\\\pi_2\\) at which power is
  evaluated.

- .m:

  Number of trials in group 1.

- .n:

  Number of trials in group 2.

- .alpha:

  Nominal significance level \\\alpha\\. No default.

- .superiority:

  Logical. If `TRUE`, power is computed only over tables where the
  observed rate in group 2 exceeds that in group 1. Defaults to `FALSE`.

## Value

A single numeric value: the power of Woolf's asymptotic test at the
response rates \\(\pi_1, \pi_2)\\, in \\\[0, 1\]\\.

## See also

Other power:
[`power_conservative()`](https://pvdmeulen.github.io/modifiedfisher/reference/power_conservative.md),
[`power_modified()`](https://pvdmeulen.github.io/modifiedfisher/reference/power_modified.md),
[`power_probability()`](https://pvdmeulen.github.io/modifiedfisher/reference/power_probability.md),
[`power_randomised()`](https://pvdmeulen.github.io/modifiedfisher/reference/power_randomised.md)

## Examples

``` r
power_asymptotic(p = c(0.2, 0.6), .m = 6, .n = 4, .alpha = 0.05)
#> [1] 0.08497562
```
