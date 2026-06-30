# Overview of algorithm

This article traces how the modified Fisher exact test is built in R,
from the test frame through to the agreeing p-value and confidence
interval. It assumes the statistical background, which is in the
[Background and
comparison](https://pvdmeulen.github.io/modifiedfisher/articles/background-and-comparison.md)
article; the derivations and proofs are in van der Meulen, Raymond and
van der Meulen (2021).

The notation: \\u\\ successes from \\m\\ trials in Group A, \\v\\ from
\\n\\ in Group B, with total \\T = u + v\\.

|                  | Group A: Success | Group A: Failure | Total  |
|:-----------------|:----------------:|:----------------:|:------:|
| Group B: Success |        ..        |        ..        | \\v\\  |
| Group B: Failure |        ..        |        ..        | \\v'\\ |
| Total            |      \\u\\       |      \\u'\\      | \\t\\  |

## The basic premise

The modified test makes a single change to the randomised (UMPU) Fisher
test: instead of always rejecting at a critical value (randomised) or
never rejecting there (conservative), it uses one global threshold
\\\gamma_0\\ and rejects at a boundary only if that boundary’s
randomisation probability \\\gamma\\ exceeds \\\gamma_0\\. The
randomisation is taken away, and most of the boundary rejection mass is
recovered.

The implementation chooses \\\gamma_0\\ for a given null odds ratio,
then inverts the resulting test for the p-value and confidence interval.
The entry point is
[`modified_fisher_exact_test()`](https://pvdmeulen.github.io/modifiedfisher/reference/modified_fisher_exact_test.md):

``` r

modified_fisher_exact_test(u = 5, m = 12, v = 7, n = 11, odds_ratio = 1)
```

The steps below show what the above function does internally. Helper
functions not exported to the namespace are prefixed by `.`, these can
still be accessed through `modifiedfisher:::` (e.g. for troubleshooting
purposes, or to visualise the steps with your own example as in this
article).

## Step 1: Build the test frame of critical values and randomisations

For a given null odds ratio \\\theta_0\\ and each total \\T = 0, \ldots,
m+n\\, find the critical values \\c_1, c_2\\ and admissible
randomisation probabilities \\\gamma_1, \gamma_2\\. This is
[`construct_test_frame()`](https://pvdmeulen.github.io/modifiedfisher/reference/construct_test_frame.md),
returning one row per \\T\\ with columns `c1`, `c2`, `gamma1`, `gamma2`
(and the starting quantiles `d1`, `d2`).

With the critical values fixed, the two requirements (conditional size
\\\alpha\\ and unbiasedness) form a \\2\times2\\ linear system solved
for \\(\gamma_1, \gamma_2)\\ by `.find_gamma12()`.
[`construct_test_frame()`](https://pvdmeulen.github.io/modifiedfisher/reference/construct_test_frame.md)
searches for the critical values giving an admissible solution (both
\\\gamma\\’s in \\\[0, 1\]\\), starting at the \\\alpha/2\\ quantiles of
the non-central hypergeometric distribution and spiralling outward in
larger squares until one is found for every \\T\\. At \\T = 0\\ and \\T
= m+n\\ there is no room to randomise, so \\\gamma_1 = \gamma_2 =
\alpha/2\\. The frame holds \\2(m+n+1)\\ randomisation probabilities in
total, which matters in Step 4.

``` r

# Example 1 from the paper (m = 12, n = 11):
construct_test_frame(.odds_ratio = 1, .m = 12, .n = 11, .alpha = 0.05, .precision = 1e-3)
```

The first ten rows of this dataframe are given by:

    #>    t c1 d1    gamma1 c2 d2    gamma2
    #> 1  0  0  0 0.0250000  0  1 0.0250000
    #> 2  1  0  0 0.0500000  1  2 0.0500000
    #> 3  2  0  0 0.1100000  2  3 0.1000000
    #> 4  3  0  0 0.2566667  3  4 0.2100000
    #> 5  4  0  1 0.6416667  4  4 0.4666667
    #> 6  5  1  1 0.1081481  4  5 0.0000337
    #> 7  6  1  2 0.3630208  5  6 0.1892519
    #> 8  7  1  2 0.9953571  6  6 0.5526948
    #> 9  8  2  2 0.2761364  6  7 0.0543831
    #> 10 9  2  3 0.8060101  7  8 0.3582323

## Step 2: Define the rejection rule for a candidate \\\gamma_0\\

Given the frame and a candidate \\\gamma_0\\, the verdict on any table
\\(u, T)\\ is deterministic:

- \\u\\ outside \\\[c_1, c_2\]\\: reject;
- \\u\\ strictly inside \\(c_1, c_2)\\: accept;
- \\u = c_1\\: reject iff \\\gamma_1 \> \gamma_0\\;
- \\u = c_2\\: reject iff \\\gamma_2 \> \gamma_0\\;
- \\c_1 = c_2\\ (a single boundary point): reject iff \\\gamma_1 +
  \gamma_2 \> \gamma_0\\.

This is `.modified_reject()`. Since the verdict depends only on the
frame and \\\gamma_0\\, not on the unknown success probability, the
whole rejection region is built once into an \\(m+1)\times(n+1)\\ 0/1
matrix (`.build_rejection_matrix()`), so later steps avoid recomputing
it.

## Step 3: Measure the (worst case) size over the nuisance parameter \\\pi_1\\

A \\\gamma_0\\ is acceptable only if the test’s true (unconditional)
size stays below \\\alpha\\. Conditional on \\T\\ the test is
size-\\\alpha\\ by construction, but **unconditionally** the rejection
probability depends on the unknown success rate \\\pi_1\\. The real size
is therefore a function of \\\pi_1\\, and we need to guard against the
‘worst case’.

For a fixed \\\pi_1\\, the rejection probability is a sum over all
tables weighted by their two binomial probabilities, computed as the
bilinear form \\p_u^\top R\\ p_v\\ with \\R\\ the rejection matrix from
Step 2
([`local_size_modified()`](https://pvdmeulen.github.io/modifiedfisher/reference/local_size_modified.md)).
The size is the maximum over \\\pi_1 \in (0, 1)\\
([`size_modified()`](https://pvdmeulen.github.io/modifiedfisher/reference/size_modified.md)).

Maximisation uses either `"zoom"` (grid, then re-grid finer around the
peak, repeat; the default, and also the SAS macro’s approach) or
`"trust"` (a trust-region optimiser with the analytic gradient
`.local_size_gradient_modified()`). The size curve can be multi-peaked,
so the safest check is to plot it via `local_size_data = TRUE`.

``` r

testframe <- construct_test_frame(.odds_ratio = 1, .m = 12, .n = 11, .alpha = 0.05, .precision = 1e-3)

# Worst-case size at a candidate gamma0 of 0.05:
size_modified(.c = 0.05, .odds_ratio = 1, .m = 12, .n = 11, .df = testframe,
              .alpha = 0.05, .precision = 1e-3, .method = "zoom",
              .maze = 10, .zoom_iter = 6)
```

## Step 4: Find the optimal \\\gamma_0\\ by bisection

We still want the largest worst-case size to be less than \\\alpha\\.
The size only changes when \\\gamma_0\\ crosses one of the \\2(m+n+1)\\
actual \\\gamma\\ values from Step 1, so rather than search a continuum,
we sort those values and search among them.

Raising \\\gamma_0\\ can only *remove* boundaries, so size is monotone
in the sorted index, which makes a bisection appropriate. This is done
in
[`optimise_gamma0()`](https://pvdmeulen.github.io/modifiedfisher/reference/optimise_gamma0.md):
it sorts the pooled `gamma1`/`gamma2` vector and bisects to the largest
threshold whose size does not exceed \\\alpha\\. The optimum is not a
single number but a half-open interval between two adjacent sorted
\\\gamma\\’s and any value in it gives the same test.

``` r

optimise_gamma0(.odds_ratio = 1, .m = 12, .n = 11, .alpha = 0.05,
                .precision = 1e-3, .method = "zoom", .maze = 10, .zoom_iter = 6)
```

## Step 5: Invert the test for the p-value and confidence interval

Steps 1 to 4 decide, at some given null odds ratio \\\theta_0\\, whether
the observed table is rejected (and the `.accept()` helper function ties
it together: build the frame, optimise \\\gamma_0\\, apply the rule).
Confidence intervals and p-values are obtained by inverting that
decision:

- **Confidence interval**: all \\\theta_0\\ not rejected at level
  \\\alpha\\. Its edges are found by bisection, starting from Woolf’s
  asymptotic limits and moving inward to the first rejection, to the
  requested `precision`.
- **P-value**: the smallest \\\alpha\\ at which the table is rejected
  for the given null, found by bisecting on \\\alpha\\.

Both invert the same test, so they agree by construction.

The main
[`modified_fisher_exact_test()`](https://pvdmeulen.github.io/modifiedfisher/reference/modified_fisher_exact_test.md)
function runs the following pipeline:

[`construct_test_frame()`](https://pvdmeulen.github.io/modifiedfisher/reference/construct_test_frame.md)
→ `.modified_reject()` / `.build_rejection_matrix()` →
[`local_size_modified()`](https://pvdmeulen.github.io/modifiedfisher/reference/local_size_modified.md)
/
[`size_modified()`](https://pvdmeulen.github.io/modifiedfisher/reference/size_modified.md)
→
[`optimise_gamma0()`](https://pvdmeulen.github.io/modifiedfisher/reference/optimise_gamma0.md)
→ `.accept()`

The function returns an `htest` with `$p.value`,`$estimate`,
`$conf.int`, the optimal `$gamma0`, the test frame, and optionally the
size-versus-\\\pi_1\\ dataframe:

``` r

result <- modified_fisher_exact_test(u = 5, m = 12, v = 7, n = 11, odds_ratio = 1, local_size_data = TRUE)

result$estimate
result$p.value
result$conf.int
result$local.size.data
```

## Reference

van der Meulen EA, Raymond K, van der Meulen PJ (2021). *Consistent
Confidence Limits, P Values, and Power of the Non-Conservative, Size-α
Modified Fisher Exact Test.* Journal of Biostatistics and Biometric
Applications 6(1):102.
