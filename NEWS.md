# STATassist 0.2.0

## New features

* `compare_two_groups()`, `compare_multiple_groups()` and `compare_one_sample()`
  gain `input_scale`. With `input_scale = "log2"` each observation is raised
  back through `2^x` before the centres of `$effect` are taken, so a fold change
  computed from already-logged data is a ratio of centres on the original scale
  rather than a ratio of logarithms, which can carry the wrong sign. The tests
  still run on the values as they were given, which is the reason for logging
  them in the first place.

* `fc_mean` now defaults to `"geom"` under `input_scale = "log2"` and stays at
  `"arith"` otherwise. On the log2 scale the geometric mean is the one centre
  that reduces `log2fc` to `mean(x) - mean(y)`. An explicit `fc_mean` always
  wins.

* `compare_one_sample()` reads `mu` on the same scale as the data. The
  comparison value under `input_scale = "log2"` is `2^mu`, which is always
  positive, so the `mu = 0` case that leaves `fold_change` undefined on the raw
  scale cannot arise there.

* New `simulate_two_groups()` generates two-group log2 expression data with a
  chosen number of features moved up and down on purpose and returns the planted
  answer with it. Its `args` element is named after the arguments of
  `compare_two_groups()`, so the comparison is one `do.call()` away, and
  features that were not planted have a true fold change of exactly zero, so a
  false positive rate can be computed by definition rather than estimated.

* `draw_butterfly_hist()` gains `type`, which draws the histogram bars alone
  (`"freq"`), a kernel density estimate alone (`"dens"`), or the two overlaid
  (`"both"`), along with `dens_adjust`, `dens_lwd`, `dens_col` and `dens_alpha`.
  A density and a bar can only be read against one axis when the bar is a
  density too, so anything other than `"freq"` moves `scale` to `"density"` and
  rejects a conflicting `scale` that was asked for explicitly. When a density is
  drawn the return value carries `group_densities`, one `"density"` object per
  group level.

* The result contract is at `0.2.1`: `parameters$input_scale` records the scale
  the data was read on. Nothing was removed or renamed, so a consumer written
  against `0.2.0` still reads every result.

## Breaking changes

* The first argument of `estimate_significance()` is now `comparison_result`
  rather than `res`, which is what its error messages have to name. Calls that
  pass the comparison positionally, as every example did, are unaffected; a call
  that named `res = ` has to be updated.

* `draw_grouped_boxplot()` draws the light theme unless asked otherwise:
  `dark` now defaults to `FALSE`, as it already did everywhere else in the
  package. Pass `dark = TRUE` for the previous default.

## Bug fixes

* `draw_volcano_plot()` no longer fails when `anno_feats = TRUE` and no feature
  clears both cutoffs. `graphics::text()` treats a zero-length `labels` as an
  error rather than as nothing to draw, so the plot is now drawn with a
  `message()` in place of the labels.

# STATassist 0.1.0

* The result contract widened from `0.1.0` to `0.2.0` with two additions and no
  changes: `posthoc`, one table per test holding one row per feature and pair of
  levels, and `diagnostics`. Both are present in every result, empty where the
  scenario has nothing to put in them, so all three scenarios are read the same
  way.

* New `compare_multiple_groups()` pairs each omnibus test with the pairwise
  procedure that shares its assumptions: one-way ANOVA with Tukey HSD, Welch's
  ANOVA with Games-Howell, Yuen's trimmed-mean ANOVA with pairwise Yuen,
  Kruskal-Wallis with Dunn's test. Repeated conditions swap in RM-ANOVA, which
  reports Mauchly's sphericity test and both epsilon corrections on the same
  row, and Friedman with Conover's comparisons. The pairwise stage runs only for
  features whose omnibus test cleared `posthoc_alpha`; a feature that did not
  qualify is absent from the table rather than present with `NA`.

* New `compare_one_sample()` runs the one-sample t-test, the signed-rank test
  with a Hodges-Lehmann pseudo-median, and a score test with a Wilson interval
  for binary features, against a hypothesised value in `design$mu`.

* New `diagnose_distribution()` and `screen_outliers()`. Each assumption is
  checked twice, by tests that fail differently, and every comparison carries
  the same checks as `$diagnostics` unless `diagnose = FALSE`. A failed check
  never changes which tests run. `screen_outliers()` flags observations and does
  not remove them.

* New `plot()` method for `sa_comparison`, with an estimate, a pairwise-contrast
  and a p-value view. It reads only the columns the contract guarantees, so one
  method covers all three scenarios, and `type = "auto"` picks the first view
  the chosen table can support.

* `estimate_significance()` and `print.sa_comparison()` no longer assume a
  two-group result, so they work on every scenario without being told which one
  produced the object.

* Added a `testthat` layer that uses base R equivalents and hand-computed
  expectations only, including regression tests pinning the v0.0.1 behaviour.

# STATassist 0.0.1

* First release. The `sa_comparison` contract fixes one shape for every
  comparison: each test table carries `features`, `n_used`, `pval`, `pval_adj`,
  `lower_conf` and `upper_conf`, and the effect estimates live in one table that
  all of the tests share.

* `compare_two_groups()` runs six methods, three at a time: Welch's t-test,
  Mann-Whitney U and Brunner-Munzel for independent groups, paired t, Wilcoxon
  signed-rank and Yuen's paired trimmed-mean test for repeated ones, with the
  fold change computed alongside them.

* `summarize_descriptive_stats()`, `estimate_significance()`,
  `draw_grouped_boxplot()`, `draw_volcano_plot()` and `draw_butterfly_hist()`.

* Dependencies are limited to base R.
