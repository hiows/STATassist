## Test environments

* local Windows 11, R 4.6.1
* R CMD check --as-cran --no-manual

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes for CRAN

* Heavy ML / plotting / embedding tests call `skip_on_cran()` so check time
  stays manageable; the full suite still runs under `devtools::test()`.
* Longer cross-validation and embedding examples are wrapped in `\donttest`.
