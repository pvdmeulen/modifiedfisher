# Construct testing frame and randomisation values

Constructs a data frame of critical values (\\c_1\\, \\c_2\\) and
randomisation probabilities (\\\gamma_1\\, \\\gamma_2\\) for every
possible value of the total \\T = 0, \ldots, m + n\\, given the null
odds ratio, significance level \\\alpha\\, and precision. Starting from
the \\\alpha/2\\ quantiles of the Fisher non-central hypergeometric
distribution, a spiral search over \\(c_1, c_2)\\ is used whenever the
initial solution for \\(\gamma_1, \gamma_2)\\ falls outside \\\[0,
1\]\\.

## Usage

``` r
construct_test_frame(.odds_ratio, .m, .n, .alpha, .precision, .message = FALSE)
```

## Arguments

- .odds_ratio:

  The null hypothesis odds ratio \\\theta_0\\. No default.

- .m:

  Number of trials in group 1.

- .n:

  Number of trials in group 2.

- .alpha:

  Nominal significance level \\\alpha\\. No default.

- .precision:

  Numerical precision for quantile calculations and
  [`BiasedUrn::dFNCHypergeo()`](https://rdrr.io/pkg/BiasedUrn/man/BiasedUrn-2-Univariate.html).
  No default.

- .message:

  A logical. Defaults to `FALSE`. Setting this to `TRUE` will print
  progress messages; useful for debugging.

## Value

A data frame with `m + n + 1` rows, one per possible total \\T = 0,
\ldots, m + n\\, and columns `t` (the total), `c1` and `c2` (lower and
upper critical values), `d1` and `d2` (the \\\alpha/2\\ quantiles used
as starting points), and `gamma1` and `gamma2` (the randomisation
probabilities at `c1` and `c2`).

## See also

[`modified_fisher_exact_test()`](modified_fisher_exact_test.md) for the
main user-facing function; [`optimise_gamma0()`](optimise_gamma0.md)
which uses this frame to find the optimal gamma0;
[`size_modified()`](size_modified.md) for the resulting test size.

Other modified: [`local_size_modified()`](local_size_modified.md),
[`modified_fisher_exact_test()`](modified_fisher_exact_test.md),
[`optimise_gamma0()`](optimise_gamma0.md),
[`power_modified()`](power_modified.md),
[`size_modified()`](size_modified.md)

## Examples

``` r
# Critical values and randomisation probabilities for m = 6, n = 4
# (reproduces Table 1 of van der Meulen et al., 2021):
construct_test_frame(.odds_ratio = 1, .m = 6, .n = 4,
                     .alpha = 0.05, .precision = 1e-3)
#>     t c1 d1    gamma1 c2 d2    gamma2
#> 1   0  0  0 0.0250000  0  1 0.0250000
#> 2   1  0  0 0.0500000  1  2 0.0500000
#> 3   2  0  0 0.1500000  2  3 0.0900000
#> 4   3  0  1 0.6000000  3  4 0.1800000
#> 5   4  1  2 0.1777778  4  4 0.3488889
#> 6   5  2  2 0.0050000  4  5 0.0050000
#> 7   6  2  3 0.3488889  5  6 0.1777778
#> 8   7  3  4 0.1800000  6  6 0.6000000
#> 9   8  4  5 0.0900000  6  6 0.1500000
#> 10  9  5  6 0.0500000  6  6 0.0500000
#> 11 10  6  6 0.0250000  6  6 0.0250000
```
