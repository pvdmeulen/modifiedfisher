# Find the optimal \\\gamma_0\\

Finds the optimal \\\gamma_0\\ threshold: the largest value of
\\\gamma_0\\ such that the actual size of the MFET (the maximum of the
local size over the nuisance parameter \\p_0\\) does not exceed
\\\alpha\\. Uses a bisection search over the \\2(m + n + 1)\\ sorted
randomisation probability values from the test frame.

## Usage

``` r
optimise_gamma0(
  .odds_ratio,
  .m,
  .n,
  .alpha,
  .precision,
  .method,
  .maze,
  .zoom_iter
)
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

  Numerical precision. No default.

- .method:

  Numerical method for maximising the local size: `"zoom"` or `"trust"`.

- .maze:

  Number of grid points evaluated at each zoom iteration.

- .zoom_iter:

  Number of zoom iterations (`"zoom"` method only).

## Value

A single numeric value: the optimal threshold \\\gamma_0\\.

## See also

[`modified_fisher_exact_test()`](https://pvdmeulen.github.io/modifiedfisher/reference/modified_fisher_exact_test.md)
for the main user-facing function;
[`construct_test_frame()`](https://pvdmeulen.github.io/modifiedfisher/reference/construct_test_frame.md)
for the test frame passed to this function;
[`size_modified()`](https://pvdmeulen.github.io/modifiedfisher/reference/size_modified.md)
for the size function being maximised.

Other modified:
[`construct_test_frame()`](https://pvdmeulen.github.io/modifiedfisher/reference/construct_test_frame.md),
[`local_size_modified()`](https://pvdmeulen.github.io/modifiedfisher/reference/local_size_modified.md),
[`modified_fisher_exact_test()`](https://pvdmeulen.github.io/modifiedfisher/reference/modified_fisher_exact_test.md),
[`power_modified()`](https://pvdmeulen.github.io/modifiedfisher/reference/power_modified.md),
[`size_modified()`](https://pvdmeulen.github.io/modifiedfisher/reference/size_modified.md)

## Examples

``` r
optimise_gamma0(.odds_ratio = 1, .m = 6, .n = 4, .alpha = 0.05,
                .precision = 1e-3, .method = "zoom", .maze = 10,
                .zoom_iter = 6)
#> [1] 0.3488889
```
