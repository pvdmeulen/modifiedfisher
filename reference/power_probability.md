# Power of the SAS Proc FREQ exact test

Computes the unconditional power of SAS Proc FREQ's exact test for
\\H_0\\: OR = 1, at response rates \\(\pi_1, \pi_2)\\. For each possible
table \\(u, v)\\, the rejection decision uses
[`pvalue_probability()`](pvalue_probability.md): the sum of central
hypergeometric probabilities given \\T = u + v\\ that are no greater
than the observed probability, compared against \\\alpha\\. This mirrors
the rejection rule used in
[`local_size_probability()`](local_size_probability.md).

## Usage

``` r
power_probability(p, .m, .n, .alpha, .superiority = FALSE)
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

A single numeric value: the power of the SAS Proc FREQ exact test at the
response rates \\(\pi_1, \pi_2)\\, in \\\[0, 1\]\\.

## See also

Other power: [`power_asymptotic()`](power_asymptotic.md),
[`power_conservative()`](power_conservative.md),
[`power_modified()`](power_modified.md),
[`power_randomised()`](power_randomised.md)

## Examples

``` r
power_probability(p = c(0.2, 0.6), .m = 6, .n = 4, .alpha = 0.05)
#> [1] 0.1755824
```
