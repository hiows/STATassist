## Test environments

* local: Windows 11 x64, R 4.6.1, `R CMD check --as-cran --no-manual`
<!-- Add once the builders have answered:
* win-builder: R-release, R-devel
* macOS builder: R release (arm64)
-->

## R CMD check results

0 errors | 0 warnings | 1 note

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Wonseok Oh <hiows97@gmail.com>'
New submission
```

## Notes

This is a new submission.

Examples that fit a model by cross-validation, tune a hyperparameter or compute
a t-SNE or UMAP embedding are wrapped in `\donttest`. A faster form of each is
still run, so every exported function has an executable example.

Test files that fit models, draw plots or embed points call `skip_on_cran()` to
keep check time down. The full suite runs locally with `devtools::test()`.
