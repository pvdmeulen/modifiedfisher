# The non-conservative, size-\\\alpha\\ modified Fisher exact test

Computes the non-conservative, size-\\\alpha\\ modified Fisher exact
test for comparing two proportions \\u/m\\ and \\v/n\\, testing \\H_0\\:
OR = `odds_ratio`. Returns an `htest` object containing test-based
two-sided p-values and \\(1 - \alpha)\\ confidence intervals for the
odds ratio, which agree by construction. The test maximises the actual
size over the nuisance parameter \\p_0\\ subject to it remaining no
greater than \\\alpha\\, making it less conservative than the standard
Fisher exact test while strictly controlling the Type I error rate.

## Usage

``` r
modified_fisher_exact_test(
  u,
  m,
  v,
  n,
  odds_ratio,
  alpha = 0.05,
  precision = 0.001,
  message = FALSE,
  maze = 10,
  method = "zoom",
  zoom_iter = 6,
  conf_int = TRUE,
  pvalue = TRUE,
  local_size_data = FALSE,
  power = TRUE,
  superiority = FALSE,
  power_at_pi1 = 0.5,
  power_at_pi2 = 0.75
)
```

## Arguments

- u:

  Number of successes observed in group 1 (out of `m` trials). No
  default.

- m:

  Number of trials in group 1. No default.

- v:

  Number of successes observed in group 2 (out of `n` trials). No
  default.

- n:

  Number of trials in group 2. No default.

- odds_ratio:

  The null hypothesis odds ratio \\\theta_0\\. No default.

- alpha:

  Nominal significance level \\\alpha\\. Defaults to 0.05.

- precision:

  Numerical precision for p-values, confidence limits, and size
  calculations. Defaults to 1e-03.

- message:

  Logical. If `TRUE`, prints progress messages during execution.
  Defaults to `FALSE`.

- maze:

  Number of grid points evaluated at each zoom iteration. Defaults to
  10.

- method:

  Numerical method for maximising the local size over the nuisance
  parameter. Either `"zoom"` (default) or `"trust"` (trust-region method
  from the trust package).

- zoom_iter:

  Number of zoom iterations. Defaults to 6.

- conf_int:

  Logical. If `FALSE`, skips the \\(1 - \alpha)\\ two-sided confidence
  interval. Defaults to `TRUE`.

- pvalue:

  Logical. If `FALSE`, skips the two-sided test-based p-value. Defaults
  to `TRUE`.

- local_size_data:

  Logical. If `TRUE`, attaches a `local.size.data` data frame to the
  results for plotting the size of the test as a function of the
  nuisance parameter. Defaults to `FALSE`.

- power:

  Logical. If `FALSE`, skips the power calculation. Defaults to `TRUE`.

- superiority:

  Logical. If `TRUE`, power is computed only over tables where the
  observed rate in group 2 exceeds that in group 1. Defaults to `FALSE`.

- power_at_pi1:

  Success probability in group 1 at which to evaluate power. Defaults to
  0.5.

- power_at_pi2:

  Success probability in group 2 at which to evaluate power. Defaults to
  0.75.

## Value

An object of class `htest`: a list with components `p.value` (the
two-sided test-based p-value, if `pvalue = TRUE`), `estimate` (the
sample odds ratio), `conf.int` (the two-sided \\(1 - \alpha)\\
confidence interval for the odds ratio, if `conf_int = TRUE`),
`null.value` (the null odds ratio), `alternative`, `method`, and
`data.name`. It additionally carries `support.data` (the test frame from
[`construct_test_frame()`](construct_test_frame.md)), `gamma0` (the
optimal threshold), `power` (the power at `power_at_pi1` and
`power_at_pi2`, if `power = TRUE`), `local.size.data` (a data frame for
plotting the size as a function of the nuisance parameter, if
`local_size_data = TRUE`), and `fn_args` (the arguments used).

## See also

[`construct_test_frame()`](construct_test_frame.md) for the underlying
test frame of critical values and randomisation probabilities;
[`optimise_gamma0()`](optimise_gamma0.md) for the gamma0 optimisation;
[`size_modified()`](size_modified.md) for the size of the test maximised
over the nuisance parameter;
[`local_size_modified()`](local_size_modified.md) for the local size at
a fixed nuisance parameter value;
[`power_modified()`](power_modified.md) for the power of the test.

Other modified: [`construct_test_frame()`](construct_test_frame.md),
[`local_size_modified()`](local_size_modified.md),
[`optimise_gamma0()`](optimise_gamma0.md),
[`power_modified()`](power_modified.md),
[`size_modified()`](size_modified.md)

## Examples

``` r
# Example 1 of van der Meulen, Raymond & van der Meulen (2021):
# 5/12 successes versus 7/11, testing H0: OR = 1 at alpha = 0.05.
modified_fisher_exact_test(u = 5, m = 12, v = 7, n = 11, odds_ratio = 1)
#> 
#>  Non-Conservative Size-α Modified Fisher's Exact Test
#> 
#> data:  u = 5, v = 7, m = 12, n = 11
#> p-value = 0.3208
#> alternative hypothesis: true odds ratio is not equal to 1
#> 95 percent confidence interval:
#>  0.072034 2.209549
#> sample estimates:
#> odds ratio 
#>  0.4081633 
#> 

# \donttest{
# The size of the test as a function of the nuisance parameter can be
# returned for diagnostic plotting (slower, as the size is evaluated on a
# grid of 101 nuisance-parameter values):
res <- modified_fisher_exact_test(u = 5, m = 12, v = 7, n = 11,
                                  odds_ratio = 1, local_size_data = TRUE)
head(res$local.size.data)
#>    pi1         size method
#> 1 0.00 0.0000000000   zoom
#> 2 0.01 0.0006921852   zoom
#> 3 0.02 0.0092790921   zoom
#> 4 0.03 0.0393074943   zoom
#> 5 0.04 0.1038314178   zoom
#> 6 0.05 0.2116501340   zoom
# }
```
