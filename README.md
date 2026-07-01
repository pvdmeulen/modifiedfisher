
<!-- README.md is generated from README.Rmd. Please edit that file -->

# modifiedfisher <a href="https://pvdmeulen.github.io/modifiedfisher/"><img src="man/figures/logo.png" align="right" height="139" alt="modifiedfisher website" /></a>

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/pvdmeulen/modifiedfisher/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pvdmeulen/modifiedfisher/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/pvdmeulen/modifiedfisher/branch/main/graph/badge.svg?token=%3Ctoken%3E)](https://codecov.io/gh/pvdmeulen/modifiedfisher)
<!-- badges: end -->

## The Non-Conservative Size-α Modified Fisher Exact Test

The `modifiedfisher` package implements the non-randomised,
non-conservative, size-α modified Fisher exact test for comparing two
binomial proportions, as introduced by van der Meulen (2008) and
extended with agreeing test-based p-values, confidence intervals, and
power in van der Meulen, Raymond, and van der Meulen (2021).

The classical Fisher exact test is the uniformly most powerful unbiased
(UMPU) test for two independent binomial proportions, but requires
randomisation at the critical values to achieve size α exactly, which
makes it impractical without a natural p-value or confidence interval.
The standard non-randomised alternatives are unsatisfactory: rejecting
only at strictly extreme values is needlessly conservative (as in the
SAS Proc FREQ exact test or R’s `fisher.test()`), while the mid p-value
and Woolf’s asymptotic test do not strictly control the Type I error
rate.

The modified test resolves this by additionally rejecting when the
randomisation probability at a critical value exceeds a data-driven
threshold γ₀, chosen to bring the actual size as close to α as possible
without exceeding it. P-values and confidence intervals are both derived
directly from this test, so they agree by construction.

## Installation

`modifiedfisher` is not yet on CRAN. Install the current development
version from GitHub using `devtools`:

``` r
# install.packages("devtools")
devtools::install_github("pvdmeulen/modifiedfisher")
```

A CRAN submission is planned once the package clears the remaining
`R CMD check` notes. In regulated or locked-down environments where
non-CRAN packages cannot be installed, contact the maintainer for a
source tarball.

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

These match Table 2 of the paper (OR = 3.173, p = 0.034, 95% CI 1.082 to
10.25).

## Comparing with other tests

The table below reproduces Table 2 Example 2 from the paper, comparing
the modified test against Woolf’s asymptotic test and the SAS Proc FREQ
exact test for the same worked example as above. The disagreement
between the SAS Proc FREQ p-value and confidence interval (p \< 0.05 but
CI includes OR = 1) illustrates the consistency advantage of the
modified test’s test-based approach.

| Test          | OR    | p-value   | 95% CI             |
|---------------|-------|-----------|--------------------|
| Modified FE   | 3.173 | **0.034** | **1.082 to 10.25** |
| Woolf         | 3.173 | 0.036     | 1.077 to 9.34      |
| SAS Proc FREQ | 3.173 | 0.039     | 0.969 to 11.31     |

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

Full documentation is available on the package website:
<https://pvdmeulen.github.io/modifiedfisher/>. The articles there cover:

- [Get
  started](https://pvdmeulen.github.io/modifiedfisher/articles/modifiedfisher.html):
  an introduction to `modified_fisher_exact_test()`, its arguments,
  power, and diagnostics.
- [Background and
  comparison](https://pvdmeulen.github.io/modifiedfisher/articles/background-and-comparison.html):
  why the test exists, and a comparison with `fisher.test()` and
  `exact2x2::fisher.exact()` (both similar to SAS PROC Freq).
- [Overview of
  algorithm](https://pvdmeulen.github.io/modifiedfisher/articles/overview-of-algorithm.html):
  a step by step overview of how the test is built in R.
- [Reproducing the paper’s
  figures](https://pvdmeulen.github.io/modifiedfisher/articles/reproducing-paper-figures.html):
  rebuilding the size and power curves, and adapting them to your own
  design.

Broadly the same material is in the package vignette, which covers:

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
