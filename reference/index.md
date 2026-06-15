# Package index

## Modified Fisher exact test

The main user-facing test.

- [`modified_fisher_exact_test()`](https://pvdmeulen.github.io/modifiedfisher/reference/modified_fisher_exact_test.md)
  : The non-conservative, size-\\\alpha\\ modified Fisher exact test

## Building blocks

Internal steps of the algorithm, exported for advanced use.

- [`construct_test_frame()`](https://pvdmeulen.github.io/modifiedfisher/reference/construct_test_frame.md)
  : Construct testing frame and randomisation values
- [`optimise_gamma0()`](https://pvdmeulen.github.io/modifiedfisher/reference/optimise_gamma0.md)
  : Find the optimal \\\gamma_0\\
- [`local_size_modified()`](https://pvdmeulen.github.io/modifiedfisher/reference/local_size_modified.md)
  : Local size of the MFET
- [`size_modified()`](https://pvdmeulen.github.io/modifiedfisher/reference/size_modified.md)
  : Actual size of the MFET for a given \\\gamma_0\\
- [`power_modified()`](https://pvdmeulen.github.io/modifiedfisher/reference/power_modified.md)
  : Power of the MFET

## Comparison tests

Size and power of the alternative tests, for benchmarking.

- [`local_size_asymptotic()`](https://pvdmeulen.github.io/modifiedfisher/reference/local_size_asymptotic.md)
  : Local size of Woolf's asymptotic test
- [`power_asymptotic()`](https://pvdmeulen.github.io/modifiedfisher/reference/power_asymptotic.md)
  : Power of Woolf's asymptotic test
- [`local_size_probability()`](https://pvdmeulen.github.io/modifiedfisher/reference/local_size_probability.md)
  : Local size of the SAS Proc FREQ exact test
- [`power_probability()`](https://pvdmeulen.github.io/modifiedfisher/reference/power_probability.md)
  : Power of the SAS Proc FREQ exact test
- [`pvalue_probability()`](https://pvdmeulen.github.io/modifiedfisher/reference/pvalue_probability.md)
  : P-value from the SAS Proc FREQ exact test
- [`local_size_randomised()`](https://pvdmeulen.github.io/modifiedfisher/reference/local_size_randomised.md)
  : Local size of the randomised (UMPU) Fisher exact test
- [`power_randomised()`](https://pvdmeulen.github.io/modifiedfisher/reference/power_randomised.md)
  : Power of the randomised (UMPU) Fisher exact test
- [`power_conservative()`](https://pvdmeulen.github.io/modifiedfisher/reference/power_conservative.md)
  : Power of the conservative Fisher exact test
