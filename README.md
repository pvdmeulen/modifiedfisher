
<!-- README.md is generated from README.Rmd. Please edit that file -->

# modifiedfisher: Non-Conservative Size-α Modified Fisher Exact Test

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/pvdmeulen/modifiedfisher/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pvdmeulen/modifiedfisher/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `modifiedfisher` package implements the non-randomised,
non-conservative, size-α modified Fisher exact test (MFET) for comparing
two binomial proportions, as introduced by van der Meulen (2008) and
extended with agreeing test-based p-values, confidence intervals, and
power in van der Meulen, Raymond, and van der Meulen (2021).

The classical Fisher exact test is the uniformly most powerful unbiased
(UMPU) test for two independent binomial proportions, but requires
randomisation at the critical values to achieve size α exactly — making
it impractical without a natural p-value or confidence interval. The
standard non-randomised alternatives are unsatisfactory: rejecting only
at strictly extreme values is needlessly conservative (as in the SAS
Proc FREQ exact test), while the mid p-value and Woolf’s asymptotic test
do not strictly control the Type I error rate. The MFET resolves this by
additionally rejecting when the randomisation probability at a critical
value exceeds a data-driven threshold γ₀, chosen to bring the actual
size as close to α as possible without exceeding it. P-values and
confidence intervals are both derived directly from this test, so they
agree by construction.

## Installation

``` r
# Install from GitHub:
devtools::install_github("pvdmeulen/modifiedfisher")
```

## Quick start

The main function is `modified_fisher_exact_test()`. It takes the number
of successes (`u`, `v`) and number of trials (`m`, `n`) in each group,
and the odds ratio under the null hypothesis (`odds_ratio`):

``` r
library(modifiedfisher)

# Example 2 from Table 2 in van der Meulen et al. (2021):
# 13/41 vs 6/47, testing H0: OR = 1
result <- modified_fisher_exact_test(u = 13, m = 41, v = 6, n = 47,
                                     odds_ratio = 1)
result
#> 
#>  Non-Conservative Size-α Modified Fisher's Exact Test
#> 
#> data:  u = 13, v = 6, m = 41, n = 47
#> p-value = 0.03369
#> alternative hypothesis: true odds ratio is not equal to 1
#> 95 percent confidence interval:
#>   1.082336 10.249348
#> sample estimates:
#> odds ratio 
#>   3.172619
```

The returned object is of class `htest`, so standard components are
accessible directly:

``` r
result$estimate   # odds ratio estimate
#> odds ratio 
#>   3.172619
result$p.value    # two-sided p-value
#> [1] 0.03369141
result$conf.int   # 95% confidence interval
#> [1]  1.082336 10.249348
#> attr(,"conf.level")
#> [1] 0.95
```

These match Table 2 of the paper (OR = 3.173, p = 0.034, 95% CI
1.082–10.25).

## Comparing with other tests

The table below reproduces Table 2 from the paper, comparing the MFET
against Woolf’s asymptotic test and the SAS Proc FREQ exact test across
four worked examples. The disagreement between the SAS Proc FREQ p-value
and confidence interval in Example 2 (p \< 0.05 but CI includes OR = 1)
illustrates the consistency advantage of the MFET’s test-based approach.

| Example          | Test          | OR    | p-value   | 95% CI            |
|------------------|---------------|-------|-----------|-------------------|
| 5/12 vs 7/11     | Modified FE   | 0.408 | 0.321     | 0.716 – 2.210     |
|                  | Woolf         | 0.408 | 0.296     | 0.760 – 2.193     |
|                  | SAS Proc FREQ | 0.408 | 0.414     | 0.550 – 2.871     |
| 13/41 vs 6/47    | Modified FE   | 3.173 | **0.034** | **1.082 – 10.25** |
|                  | Woolf         | 3.173 | 0.036     | 1.077 – 9.34      |
|                  | SAS Proc FREQ | 3.173 | 0.039     | 0.969 – 11.31     |
| 37/65 vs 28/71   | Modified FE   | 2.029 | 0.0439    | 1.023 – 4.083     |
|                  | Woolf         | 2.029 | 0.0425    | 1.024 – 4.021     |
|                  | SAS Proc FREQ | 2.029 | 0.0584    | 0.970 – 4.258     |
| 72/128 vs 58/142 | Modified FE   | 1.804 | 0.0175    | 1.108 – 2.936     |
|                  | Woolf         | 1.804 | 0.0167    | 1.113 – 2.925     |
|                  | SAS Proc FREQ | 1.804 | 0.0203    | 1.082 – 3.018     |

## Diagnosing size control

Setting `local_size_data = TRUE` attaches a data frame of the local
size, evaluated across the nuisance parameter space, to the returned
object. This is the diagnostic plot recommended in the paper for
confirming that the size optimisation has not been trapped at a local
maximum:

``` r
result_diag <- modified_fisher_exact_test(u = 13, m = 41, v = 6, n = 47,
                                          odds_ratio = 1,
                                          local_size_data = TRUE)

alpha <- 0.05

head(result_diag$local.size.data)
#> # A tibble: 6 × 3
#>     pi1   size method
#>   <dbl>  <dbl> <chr> 
#> 1  0    0      zoom  
#> 2  0.01 0.0543 zoom  
#> 3  0.02 0.462  zoom  
#> 4  0.03 1.25   zoom  
#> 5  0.04 2.16   zoom  
#> 6  0.05 2.96   zoom
```

<picture>
<source media="(prefers-color-scheme: dark)" srcset="man/figures/README-dark_plot_data-1.svg">
<source media="(prefers-color-scheme: light)" srcset="man/figures/README-light_plot_data-1.svg">
<img alt="Plot of the local size of the modified Fisher exact test as a function of the nuisance parameter pi1, with a horizontal reference line at alpha = 0.05. The size reaches but does not exceed the reference line over a broad central region.">
</picture>

The size should reach but not exceed the dashed α reference line,
confirming non-conservative size control across the full range of the
nuisance parameter.

## Learn more

The package vignette covers additional features in detail:

- Testing a null hypothesis other than OR = 1
- Controlling numerical precision and choosing the optimisation method
- Power calculations, including for superiority testing
- Accessing the test frame (critical values and randomisation
  probabilities)
- Reproducing the size and power comparison plots from the paper using
  the exported `local_size_*()` and `power_*()` functions

``` r
vignette("modifiedfisher")
```

## References

van der Meulen EA, Raymond K, van der Meulen PJ (2021). Consistent
Confidence Limits, P Values, and Power of the Non-Conservative, Size-α
Modified Fisher Exact Test. *Journal of Biostatistics and Biometric
Applications*, 6(1):102. <https://doi.org/10.13140/RG.2.2.34661.73448>

van der Meulen EA (2008). A Nonrandomized, Nonconservative Version of
the Fisher Exact Test. *Communications in Statistics — Theory and
Methods*, 37:699–708.
