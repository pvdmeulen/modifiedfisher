# modifiedfisher

A non-conservative, size-α modified Fisher exact test for comparing two
proportions. Its p-value, confidence interval, and power all agree with
one another by construction, and it controls the Type I error rate
strictly while recovering most of the power that the ordinary exact test
gives up to conservativeness. It implements the test of van der Meulen,
Raymond and van der Meulen (2021), building on van der Meulen (2008).

## Installation

``` r

devtools::install_github("pvdmeulen/modifiedfisher")
```

## Usage

The main function takes the two success counts (`u`, `v`), the two group
sizes (`m`, `n`), and a null odds ratio, and returns a familiar `htest`
object.

``` r

# Example 1 of the paper: 5/12 successes versus 7/11, H0: OR = 1, alpha = 0.05.
res <- modified_fisher_exact_test(u = 5, m = 12, v = 7, n = 11, odds_ratio = 1)

res$estimate    # sample odds ratio
res$p.value     # two-sided, test-based p-value
res$conf.int    # agreeing (1 - alpha) confidence interval
```

## Learn more

- [Introduction to
  modifiedfisher](https://pvdmeulen.github.io/modifiedfisher/articles/modifiedfisher.md):
  an introduction to
  [`modified_fisher_exact_test()`](https://pvdmeulen.github.io/modifiedfisher/reference/modified_fisher_exact_test.md),
  its arguments, power, and diagnostics.
- [Background and
  comparison](https://pvdmeulen.github.io/modifiedfisher/articles/background-and-comparison.md):
  why the test exists, and a comparison with
  [`fisher.test()`](https://rdrr.io/r/stats/fisher.test.html) and
  [`exact2x2::fisher.exact()`](https://rdrr.io/pkg/exact2x2/man/exact2x2.html)
  (both similar to SAS PROC Freq).
- [Overview of
  algorithm](https://pvdmeulen.github.io/modifiedfisher/articles/overview-of-algorithm.md):
  a step by step overview of how the test is built in R.
- [Reproducing the paper’s
  figures](https://pvdmeulen.github.io/modifiedfisher/articles/reproducing-paper-figures.md):
  rebuilding the size and power curves, and adapting them to your own
  design.

## Reference

van der Meulen EA, Raymond K, van der Meulen PJ (2021). *Consistent
Confidence Limits, P Values, and Power of the Non-Conservative, Size-α
Modified Fisher Exact Test.* Journal of Biostatistics and Biometric
Applications 6(1):102.

van der Meulen EA (2008). A Nonrandomized, Nonconservative Version of
the Fisher Exact Test. *Communications in Statistics - Theory and
Methods*, 37:699-708.
