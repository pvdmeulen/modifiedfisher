# Local size of Woolf's asymptotic test

Computes the unconditional rejection probability (local size) of Woolf's
asymptotic Wald test for \\H_0\\: OR = \\\theta_0\\ at nuisance
parameter \\p_0\\, where \\p_1 = p_0 / (p_0 + \theta_0 (1 - p_0))\\.
Rejects \\H_0\\ when \\\|\log \hat{\theta} - \log \theta_0\| / SE(\log
\hat{\theta}) \> z\_{\alpha/2}\\, where \\\hat{\theta} = u(n-v) /
((m-u)v)\\. Uses the Haldane correction (replacing zero cells with 0.5)
to handle tables where one or more cells are zero, matching the SAS
macro.

## Usage

``` r
local_size_asymptotic(nuisance, .odds_ratio, .m, .n, .alpha)
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

A single numeric value: the local size of Woolf's asymptotic test at the
given nuisance parameter, in \\\[0, 1\]\\.

## See also

Other size: [`local_size_modified()`](local_size_modified.md),
[`local_size_probability()`](local_size_probability.md),
[`local_size_randomised()`](local_size_randomised.md),
[`size_modified()`](size_modified.md)

## Examples

``` r
local_size_asymptotic(nuisance = 0.5, .odds_ratio = 1, .m = 6, .n = 4,
                 .alpha = 0.05)
#> [1] 0.01367188
```
