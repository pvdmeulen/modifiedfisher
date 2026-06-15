# Local size of the SAS Proc FREQ exact test

Computes the unconditional rejection probability (local size) of the SAS
Proc FREQ exact test for \\H_0\\: OR = \\\theta_0\\ at nuisance
parameter \\p_0\\, where \\p_1 = p_0 / (p_0 + \theta_0 (1 - p_0))\\. For
each possible table \\(i, j)\\, rejection is determined by
[`pvalue_probability()`](pvalue_probability.md): the sum of Fisher
non-central hypergeometric probabilities (given \\T = i + j\\ and OR =
\\\theta_0\\) that are no greater than the observed probability,
compared against \\\alpha\\.

## Usage

``` r
local_size_probability(nuisance, .odds_ratio, .m, .n, .alpha)
```

## Arguments

- nuisance:

  The nuisance parameter \\p_0\\: success probability in group 1 under
  \\H_0\\. Must be in \\\[0, 1\]\\.

- .odds_ratio:

  The null hypothesis odds ratio \\\theta_0\\. No default.

- .m:

  Number of trials in group 1.

- .n:

  Number of trials in group 2.

- .alpha:

  Nominal significance level \\\alpha\\. No default.

## Value

A single numeric value: the local size of the SAS Proc FREQ exact test
at the given nuisance parameter, in \\\[0, 1\]\\.

## See also

Other size: [`local_size_asymptotic()`](local_size_asymptotic.md),
[`local_size_modified()`](local_size_modified.md),
[`local_size_randomised()`](local_size_randomised.md),
[`size_modified()`](size_modified.md)

## Examples

``` r
local_size_probability(nuisance = 0.5, .odds_ratio = 1, .m = 6, .n = 4,
                    .alpha = 0.05)
#> [1] 0.02148438
```
