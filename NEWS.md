# modifiedfisher 0.0.2

This release is a major overhaul of the package structure, documentation, and
test coverage relative to v0.0.1. The underlying statistical methodology is
unchanged.

## New exports

v0.0.1 exported only `modified_fisher_exact_test()`. The following diagnostic
and comparison functions are now exported too:

- `construct_test_frame()` — builds the full table of critical values and
  randomisation probabilities for each possible total T.
- `optimise_gamma0()` — finds the optimal γ₀ threshold by bisection over the
  sorted randomisation probabilities.
- `size_mfet()` — the local size of the modified Fisher exact test, maximised
  over the nuisance parameter.
- `local_size_mfet()`, `local_size_woolf()`, `local_size_procfreq()`,
  `local_size_randomised()` — local size at a fixed nuisance parameter, for
  each of the four tests compared in van der Meulen et al. (2021).
- `power_mfet()`, `power_woolf()`, `power_procfreq()`, `power_randomised()`,
  `power_conservative()` — power for each test at specified response rates.
- `pvalue_procfreq()` — the SAS Proc FREQ-style conditional exact p-value.

## Renamed functions and files

Internal functions have been renamed for consistency, and file names now match
function names throughout. The main changes:

| Old name | New name |
|---|---|
| `calc_exp_value` | `.calc_expected_value` |
| `find_gamma` | `.find_gamma12` |
| `local_size` | `local_size_mfet` |
| `local_size_gradient` | `.local_size_gradient_mfet` |
| `local_power` / `power_mfet` | `power_mfet` |
| `mod_fe_test` | `.mfet_reject` |
| `mod_fe_size` / `mfet_size` | `size_mfet` |
| `sas_procfreq_pvalue` | `pvalue_procfreq` |

Pure internal helpers are now dot-prefixed (`.calc_expected_value`,
`.find_gamma12`, `.mfet_reject`, `.accept`, `.build_rejection_matrix`,
`.local_size_gradient_mfet`) and carry `@noRd`, so they no longer appear in the
help index.

## Documentation

- Updated the `Title` and `Description` fields in `DESCRIPTION`.
- All exported functions now have full roxygen documentation, with `@param`,
  `@return`, `@examples`, `@family`, and `@seealso` cross-links.
- Added `knitr` and `rmarkdown` to `Suggests`, along with a vignette
  (`vignettes/modifiedfisher.Rmd`).
- The PDF reference manual now builds cleanly via `R CMD Rd2pdf`.

## Tests

Test coverage has grown from no test files to 7 (451 lines), covering
`.calc_expected_value`, `.find_gamma12`, `.mfet_reject`, `construct_test_frame`,
`local_size_*`, `modified_fisher_exact_test`, and `power_mfet`. The four worked
examples from Table 2 of van der Meulen et al. (2021) are included as regression
tests in `test-modified_fisher_exact_test.R`.

### A note on Table 2 in the source paper

While checking the package against those Table 2 examples, I came across two
typos in the published version. The package produces the correct numbers — it's
the printed table that's off. Both were straightforward to confirm using
closed-form Woolf limits and base R's `fisher.test()`, neither of which relies
on any package code.

- **Example 1 (5/12 vs 7/11): the lower CI limits dropped a leading zero.** The
  modified test's lower limit is 0.0716, not the printed "0.716" (and likewise
  Woolf's is 0.076 not 0.760, and Proc FREQ's is 0.055 not 0.550). The printed
  interval gives itself away, too: (0.716, 2.210) doesn't contain its own point
  estimate of 0.408, which a test-based interval can't do.
- **Example 4: printed as "72/128 vs 58/142", but the numbers come from
  u = 71.** Running it with u = 71 lines everything up — estimate 1.804,
  p-value 0.0175, CI (1.108, 2.936). With u = 72, nothing matches.
- Examples 2 and 3 are fine as printed.

The tests use the corrected values, with comments noting why, so they don't get
"fixed" back to the printed ones later.
