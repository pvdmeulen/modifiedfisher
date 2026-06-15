# P-value from the SAS Proc FREQ exact test

Computes the two-sided Fisher exact p-value for \\H_0\\: OR =
`odds_ratio` by summing all probabilities in the conditional
distribution of \\U\\ given \\T = t\\ that are no greater than the
observed probability \\P\_{\theta_0}(U = s \mid T = t)\\. The
conditional distribution is the Fisher non-central hypergeometric with
the specified odds ratio, computed via
[`BiasedUrn::dFNCHypergeo()`](https://rdrr.io/pkg/BiasedUrn/man/BiasedUrn-2-Univariate.html).
When `odds_ratio = 1` this reduces to the standard central
hypergeometric, matching the p-value reported by SAS Proc FREQ for
\\H_0\\: OR = 1.

## Usage

``` r
pvalue_probability(s, t, m, n, odds_ratio = 1)
```

## Arguments

- s:

  Number of successes observed in group 1. No default.

- t:

  Total number of successes across both groups (\\t = s + \\ successes
  in group 2); the conditioning variable in the hypergeometric
  distribution. No default.

- m:

  Number of trials in group 1. No default.

- n:

  Number of trials in group 2. No default.

- odds_ratio:

  Null hypothesis odds ratio \\\theta_0\\. Defaults to 1.

## Value

A single numeric value: the two-sided conditional exact p-value, in
\\\[0, 1\]\\.

## See also

[`local_size_probability()`](https://pvdmeulen.github.io/modifiedfisher/reference/local_size_probability.md)
for the local size of the SAS Proc FREQ exact test;
[`power_probability()`](https://pvdmeulen.github.io/modifiedfisher/reference/power_probability.md)
for the power of the SAS Proc FREQ exact test;
[`modified_fisher_exact_test()`](https://pvdmeulen.github.io/modifiedfisher/reference/modified_fisher_exact_test.md)
for the main user-facing function.

## Examples

``` r
# Two-sided exact p-value for 5/12 vs 7/11 under H0: OR = 1
# (here t = u + v = 5 + 7 = 12):
pvalue_probability(s = 5, t = 12, m = 12, n = 11, odds_ratio = 1)
#> [1] 0.4136492
```
