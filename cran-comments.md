## Submission

This is a first submission of modifiedfisher (version 0.0.3).

The package implements the non-randomised, non-conservative, size-alpha
modified Fisher exact test, with agreeing test-based p-values, confidence
intervals, and power.

## Test environments

<!-- Replace with the environments you actually ran, and the date. -->

- Local: macOS 26.5.1, R 4.5.2, R CMD check --as-cran
- win-builder: R-devel and R-release (devtools::check_win_devel(), check_win_release())
- macOS builder: <result> (mac.r-project.org/macbuilder)
- GitHub Actions: ubuntu-latest, windows-latest, macOS-latest (R-release and R-devel)

## R CMD check results

<!-- Replace with the actual output. The expected result is below. -->

0 errors | 0 warnings | 1 note

The note is the standard one for a first submission:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Peter van der Meulen <peter.vd.meulen@icloud.com>'
  New submission

The DESCRIPTION references the underlying paper using the <doi:...> form. The
DOI is assigned by ResearchGate and resolves to the article.

The package uses British spelling, declared via the Language: en-GB field, so
spellings such as 'randomised' are handled by the en-GB dictionary. Remaining
non-dictionary terms, such as the surname 'Woolf', are listed in inst/WORDLIST.

## Downstream dependencies

There are no downstream dependencies, as this is a first submission.
