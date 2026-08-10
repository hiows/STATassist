# STATassist 0.7.0

This adds the fourth result contract, the one for a function that searches rather
than fits. A model is handed its predictors and answers about them; a selection
is handed candidates and answers which of them to keep, which is a question no
existing slot had a shape for. Two functions answer it and disagree about the
price: one holds rows out and keeps the subset that scored best on them, the other
keeps the model whose likelihood is worth what its parameters cost.

## New features

* New `perform_rfe()` runs a recursive feature elimination: the candidates are
  ranked, the weakest is dropped, and what is left is scored, until one predictor
  is standing. The elimination is inside the resampling rather than before it, so
  the ranking is recomputed in each fold and no subset size is scored on the rows
  that chose it. There is no `cv` argument for that reason. An elimination with
  nothing held out has no score to choose a size by, so it would not be a shorter
  version of this function but a different and wrong one.

* New `perform_stepwise()` runs a stepwise search by AIC or BIC: the model is refitted
  with each term taken out or put back, the move that lowers the criterion most is
  taken, and the search stops when no single move lowers it any further.
  `direction` chooses which moves are allowed, and `"both"` is the one that
  reconsiders a term it has already dropped. It has no `cv` and no `seed` either,
  for a different reason than `perform_rfe()` has none: a criterion is a penalised
  likelihood computed on the rows the model was fitted to, so nothing is held out,
  nothing is resampled and nothing is random.

* The two criteria are one search at two prices, 2 per parameter against `log(n)`,
  so past seven observations BIC charges more and keeps fewer predictors.
  `profile` reports both at every step whichever one is searching, so a path
  chosen by AIC can be read against what BIC would have said about the same
  models. They are `AIC()`'s and `BIC()`'s own numbers, the ones a model's
  `fit_stats` reports, rather than the `extractAIC()` scale `step()` searches on;
  the two differ by a constant across models fitted to the same rows, so the path
  is ordered identically and only the printed values differ.

* `ranking$estimate` is one number for both groups of candidates: what the
  criterion would be with that predictor left out of the selected model, minus
  what it is with it in. A predictor the search kept is worth the rise that
  dropping it would cause and so is positive; one the search left out would raise
  the criterion by being added and so is negative. The sign is the search's own
  decision about it and the size is by how much, which puts the selection at the
  top of the table with no second sort.

* A search that walks back to the intercept is an error rather than a result. That
  no candidate pays for itself at this charge is an answer, but not one the
  contract can carry, since `selected` would be empty and `ranking` would have
  nothing to be read against. The message says which charge was levied and, for a
  BIC search, that AIC levies 2 instead.

* New `sa_selection`, the fourth result contract. `candidates` takes the place
  `features` holds in a comparison, `terms` in a model and `points` in a
  reduction, and it is in the order the search ranked rather than the order the
  columns arrived. Two tables hang off it, `ranking` with one row per candidate
  and `profile` with one row per model the search compared, because "which
  predictors" and "how many" are two answers and only the second can be read from
  a table of models. What `profile` repeats depends on how the search moved: one
  row per subset size for an elimination, one per step of the path for a stepwise
  search, with the `n_vars` column and the single `chosen` row in common.
  `resampling` is what tells the two apart at a glance, since a search that holds
  nothing out leaves it `NULL`. The eleven slots are `analysis`, `candidates`,
  `design`, `parameters`, `selected`, `ranking`, `profile`, `resampling`,
  `engine`, `fit` and `metadata`.

* `$selected` is a set of column names and nothing else, so it goes straight back
  into `predictors =` of any `fit_*()` call. That is what the ranking is built
  for, and it is why a linear or logistic search folds its coefficients back onto
  the columns they came from: `lm()` and `glm()` see a `k`-level factor as
  `k - 1` dummy columns, and none of the three is something a later fit would
  accept. A factor is ranked by the largest statistic among its levels, so it is
  kept as long as one level is worth keeping, and eliminated as a column.
  `perform_stepwise()` has none of that work to do: `step()` moves whole terms of the
  formula, and a factor is one term however many dummy columns it becomes.

* The ranking is the absolute t or Wald statistic rather than the coefficient,
  which is where this departs from `caret::lmFuncs`. A coefficient is an effect
  per unit of its predictor, so ranking by its size ranks by the units the
  predictors happened to be measured in: the same model with a predictor in
  grams rather than in kilograms eliminates in a different order. The statistic
  has divided the units out, and it is what `caret::lrFuncs` already ranks a
  logistic regression by, so the two models now rank on one scale. A forest ranks
  by the permutation importance `fit_rf()` reports as `estimate`, so its ranking
  here and its importance table there are the same measure.

* `model` names what is fitted inside the search — `"linear"`, `"logistic"` or
  `"rf"` — and the outcome has to agree with it. A disagreement is an error
  naming the model that would have fitted rather than a silently different
  analysis. The forest inside the search grows at `randomForest()`'s own `mtry`
  for each subset size rather than at one value throughout, since a fixed `mtry`
  exceeds the predictor count at the small end of the profile, which is where the
  whole question is.

* `outcome_lv` gains a companion, `control_label`, in `perform_rfe()` and in the
  helper every model function reads its levels through. It names the reference
  class on its own, for the usual case where the sort has it backwards and the
  other level needs no saying, and it defaults to `outcome_lv[1]`, so a call that
  names one of the two names the reference either way. Naming both and
  disagreeing is an error rather than a precedence rule: either argument is a
  complete answer, so a call that contradicts itself has no reading more likely
  than the other.

* `control_label` now reaches every function whose result has a direction to
  report: `fit_logistic_regression()`, `fit_elastic_net()`, `fit_rf()` and
  `fit_svm()` beside `perform_rfe()`, and `compare_two_groups()` and
  `compare_multiple_groups()` beside them. It defaults to the reference the call
  already named, so no existing call changes.

* The two families read it differently, and the difference is what the argument
  is for. A model's `outcome_lv` holds the two classes and nothing else, so
  `fit_logistic_regression(outcome_lv = c("a", "b"), control_label = "b")` is a
  call that contradicts itself and is an error. A comparison's `group_lv` also
  carries the display order of every level, which `control_label` says nothing
  about, so `compare_two_groups(group_lv = c("a", "b"), control_label = "b")`
  moves `b` to the front and reverses every difference and ratio. In
  `compare_multiple_groups()` the same move carries the post-hoc stage with it:
  the reference is the denominator of the fold change and the subtracted side of
  every contrast at once, and the remaining levels keep contrasting each other
  in the order they were given.

* Naming it is also enough to say that a numeric column of zeroes and ones is
  two classes rather than two numbers, in the four models as in `perform_rfe()`,
  and the message that announces the guess now says so.

* `$fit` is the `rfe` object and carries no `sa_fit` class, unlike the `$fit` of
  a model. That class exists to route `coef()` and `summary()` to a
  `$finalModel`, and a search has none. Everything else in the object is a
  scalar, a character vector, a named list or a data.frame, so dropping that one
  slot leaves an object that writes out as JSON, which is the rule all four
  contracts keep.


# STATassist 0.6.0

This release adds the two result contracts that have no feature axis. Everything
before it answered a question about features; a model answers about terms and a
dimension reduction about points, and both needed a shape of their own rather
than a feature table with the wrong column names.

## New features

* New `split_data()` partitions rows into a training and a test half, and it is
  the first function in the package that computes no statistic. The split is
  where a training set first learns something it must not know, so both routes
  are closed by arguments rather than by advice: `stratified` keeps the balance
  of the whole data on both sides, and `id` sends every row of one sampling unit
  to the same side. `caret::createDataPartition()` has no notion of a sampling
  unit, so the rows are folded into units, each unit takes the stratum its rows
  agree on, the partition is drawn over units and then unfolded back to row
  indices. A unit whose rows disagree about `stratified` is an error and not a
  majority vote, since a majority produces a split that looks stratified and is
  not.

* With `id`, `p_train` is a proportion of units rather than of rows, and the two
  differ whenever units differ in size, so the row proportion actually reached is
  reported in `parameters$achieved_p`. Both `stratified` and `id` take either a
  column name or a vector of their own, and `design` records which it was, so a
  result says what it was split on. The returned shape does not depend on
  `times`: `datasets` is a list of length one when one split was asked for, and
  `train_rows` / `test_rows` carry the original row numbers.

* New `sa_model`, the second result contract, and `fit_linear_regression()` and
  `fit_logistic_regression()` on top of it. Every table in the comparison
  contract repeats `features` in one order; a model has one outcome and a set of
  terms, and the terms are not the columns that were handed in, since one factor
  predictor becomes several. `terms` takes that place and `coefficients$terms`
  repeats its order. The eleven slots are `analysis`, `terms`, `design`,
  `parameters`, `coefficients`, `fit_stats`, `performance`, `resampling`,
  `engine`, `fit` and `metadata`.

* `$fit` is the one exception to the rule that no engine object is stored. A
  model that cannot be handed to `predict()` is not a model, and there is no
  plain list in R that stands in for one. The other ten slots keep the old rule,
  so dropping `$fit` leaves an object that goes out as JSON.

* `$fit` carries an `sa_fit` class **in front of** `caret`'s `train` rather than
  in place of it, so every method `caret` defines is still inherited and
  `inherits(fit$fit, "train")` is still true. It exists because `caret` defines
  no `coef()` for its own class, so `coef()` on a fitted model returned `NULL`,
  which is not an answer a fitted model may give. `predict.sa_fit()` also accepts
  `type = "response"`, which `caret::predict.train()` rejects in its first line
  even though the value is already there as the second column of `type = "prob"`.

* Cross-validation scores these two models and does not choose them. The final
  fit uses every usable row either way, so `cv = TRUE` and `cv = FALSE` produce
  identical `coefficients` and differ only in `performance` and `resampling`.
  `parameters` records the values that were used rather than the ones that were
  asked for: LOOCV reads neither a fold count nor a repeat count, so both are
  `NA`, and `cv = FALSE` leaves `cv_method` `NA` as well.

* Listwise deletion happens once, before the folds are drawn, rather than inside
  each fold. Left to the engine it happens per fold, so the folds are scored on
  different subsets and their numbers stop being comparable. `design$n_dropped`
  says how many rows went. A predictor with one value cannot contribute and is
  excluded with a message and a record in `design$dropped_predictors`; a term
  that could not be estimated because another predictor already spans it keeps
  its row with `NA` inference, because dropping it makes the table quietly
  shorter than the model.

* Intervals are built from the summary table by hand rather than by
  `stats::confint()`, for two reasons. A rank-deficient fit has a coefficient
  vector and a covariance matrix of different lengths and fails, and the default
  `confint()` for a `glm` is a profile-likelihood interval, which is a different
  quantity from the Wald standard error and z in the same row. A linear model now
  gets a t interval and a logistic one a Wald interval, and both always agree
  with the standard error beside them.

* New `simulate_regression()`, `simulate_classification()` and
  `make_block_cor()`, so a coefficient table can be scored the way
  `simulate_two_groups()` made a verdict scorable. A predictor that was not
  planted has a coefficient of exactly zero, which makes a false positive a count
  rather than an estimate. `args` is named after the arguments of the `fit_*()`
  functions and `split_args` after those of `split_data()`, because a single list
  would kill one of the two `do.call()` calls with an unused argument.

* The answer comes back on both axes: `truth` has one row per predictor and
  `truth_term` one row per term, aligned with `coefficients` by position. The
  term names are computed from what `lm()` would paste together rather than read
  back off a fit. `max_cor_signal` records why `make_block_cor()` exists at all:
  a null predictor correlated with a planted one is pulled off zero by the data,
  no number of rows repairs it, and the lookup explains what would otherwise be
  filed as a false positive.

* `make_block_cor()` is exported and now validates. Overlapping blocks are an
  error rather than a later block quietly winning, and the matrix is checked for
  positive definiteness, since symmetric with a unit diagonal and every value in
  `[-1, 1]` still describes data that cannot exist. The correlated draws use
  `chol()` rather than a new dependency, which is stricter as a side effect: an
  eigen decomposition would pass a matrix that is not positive definite.

* New `fit_elastic_net()`, covering LASSO, ridge and elastic net as three corners
  of one model rather than three models, and both outcome types from the column
  it is given. This is the first model where resampling **chooses**, so
  `parameters$alpha` and `parameters$lambda` are the values that won and the grid
  is the rows of `performance`, with `n_candidates` recording how many pairs were
  scored. `print.sa_model()` no longer assumes the first row of `performance` is
  the best one; it finds the row `caret` reports as `bestTune`.

* A penalized coefficient table has no `stderr`, `pval`, `statistic`, `df` or
  interval columns, and they are **absent rather than `NA`**. A penalized
  estimate is deliberately biased and the usual standard error assumes an
  unbiased one, so there is no honest number for them; a column that is entirely
  `NA` reads as a table with its values missing, which is a different claim from
  a question this model does not answer. `selected` takes their place, the
  inference columns became an all-or-nothing group that `sa_new_model()`
  enforces, and `is.null(fit$coefficients$pval)` is how a consumer tells the two
  kinds of table apart.

* Factor predictors are coded here rather than left to the engine. `glmnet` takes
  a numeric matrix and what `caret` does in front of it is `as.matrix()`, which
  turns a three-level factor into one numeric predictor of integer codes and
  fits without complaint — a model that assumes an order and a spacing nobody
  stated. `sa_design_matrix()` builds the same coding `lm()` would.

* New `predict()` and `coef()` methods for `sa_model`. `predict()` belongs on the
  result rather than on `$fit`, because an engine object knows only the column
  names it was handed: `glmnet` and `kernlab` were given a design matrix and read
  it by position, so a frame whose numeric columns are in a different order is
  multiplied by the wrong coefficients with no error at all, and a factor
  predictor disappears from it entirely. Only the result knows which columns were
  predictors and what the levels of a factor were, which is also what makes one
  line work for all five models.

* One prediction comes back per row of `newdata`, and a row that cannot be read
  is `NA` rather than absent, so the answer stays aligned with the input. This
  has to be stated rather than left to the engine because of the penalized
  models: a coefficient of exactly zero is never read in a sparse product, so
  without the rule a missing value would break some rows and not others depending
  on which predictor it landed in. A level missing from `newdata` is not an
  error — its dummy column is zero, since the levels come from the new
  `design$predictor_lv` rather than from the new rows — while a level the
  training data never saw is an error naming both the column and the level.

* New `fit_rf()`, the first model with no coefficients. A forest holds hundreds
  of trees and their splits rather than one effect per predictor, so `estimate`
  is permutation importance (`%IncMSE` for regression,
  `MeanDecreaseAccuracy` for classification) and `impurity` carries the other
  measure the same fit reports, because the two disagree in a way worth seeing.
  Importance is not divided by the between-tree standard deviation the way
  `randomForest::importance()` does by default, since that ratio is referred to
  no distribution. A negative importance is an answer and not a missing value: a
  predictor that carries nothing can do worse than its own permutation.

* The table is sorted by importance, since that is the order worth reading first
  in a model with no coefficients, and the order the columns arrived in is still
  in `design$predictors`. `fit_stats` is out-of-bag rather than in-sample and
  says so in every name (`oob_rmse`, `oob_accuracy` and the rest), because a
  third of the rows are out of bag for each tree and the forest has already
  predicted them from trees that never saw them.

* `mtry` is the only tuned argument, and `NULL` resolves to the single
  conventional value rather than to a grid, so `fit_rf(data, outcome,
  cv = FALSE)` is a complete call. An `mtry` above the number of predictors is
  refused by name, because `randomForest()` would warn and silently fit a
  different forest from the one reported.

* New `fit_svm()`, the second model with no coefficients, for the opposite
  reason: a forest has too many numbers per predictor to report one, and a radial
  kernel machine has none, holding support vectors and weights that are points in
  the data rather than directions in the predictor space. `estimate` is
  permutation importance measured in the metric the resampling tuned on, so the
  table and `performance` are in the same unit, and it is measured on the fitted
  rows because a machine sees every row at once and has no out-of-bag half. That
  is weaker than the forest's and the documentation says so.

* `sigma = NULL` reads the kernel width from the data as the median of
  `kernlab::sigest()` rather than tuning it blind, the predictors are centred and
  scaled so that width is on the standardised scale, and `seed` now fixes all
  three random steps a machine takes: the permutations, `sigest()`'s sampling of
  row pairs, and the Platt scaling behind class probabilities.

* New `sa_reduction`, the third result contract, with `perform_pca()`,
  `perform_tsne()` and `perform_umap()`. `points` takes the place of `features`
  and `terms`, and it is `points` rather than `samples` because this is the one
  row axis in the package the caller chooses: `embedding_scale` decides which
  margin is embedded and `design$point_type` reports the answer. `variance` and
  `loadings` are present exactly when the method is PCA, and `sa_new_reduction()`
  refuses an object carrying one without the other.

* These are three functions rather than one call with a `methods` argument. The
  comparison functions can put four tests in one result because they answer the
  same question on one feature axis; these three answer in coordinates that share
  no scale, so nothing but `points` could be joined between them, and
  `perplexity`, `n_neighbors` and `metric` each belong to exactly one of them.
  What makes them comparable is the input instead: all three read `data` through
  one function, so the same rows drop for the same reason.

* `perform_pca()` and `perform_tsne()` see literally the same matrix, which takes
  turning off two of `Rtsne`'s defaults: `normalize = TRUE` would overwrite the
  `center` and `scale` that were asked for, and `pca = TRUE, initial_dims = 50`
  would show t-SNE a rotation rather than the matrix. Both are recorded in
  `engine$overridden`. `perform_umap()` standardises nothing by default, because
  `metric` is its own argument and `"cosine"` or `"pearson"` already compares the
  shape of a row rather than its size.

* `embedding_scale` does different work in the three, which is the other half of
  why one function could not hold them. PCA is a singular value decomposition, so
  one fit answers both margins and the matrix is never turned around:
  `embedding_scale = "features"` rescales the rotation from unit length to
  variance-weighted length and puts the sample scale in `$loadings`, while
  `$variance` and `$fit` are untouched, so the axis labels are the same either
  way. t-SNE and UMAP have no such pair and embed the rows they were given.

* Transposing the input by hand is documented at length in all three, because it
  is the one mistake here that produces a plot instead of an error.
  `prcomp()` always centres and scales the columns it is handed, so
  `perform_pca(t(data))` standardises samples rather than features; the shapes
  match, the picture reads, and the answer is to a third question.

* Both neighbourhood sizes are derived from the engine's own limits when they are
  `NULL`, and from the number of **points** rather than the number of samples. A
  value that breaks the limit is an error naming it, and so is a derived value
  that does not hold; below about sixteen points the derived value comes with a
  message, since a neighbourhood is most of what these two methods are.

## Dependencies

* The package no longer depends on base R alone. `caret` came in with
  `split_data()`, on the grounds that this package builds machine learning on top
  of an established framework rather than reimplementing it, and `glmnet`,
  `randomForest`, `kernlab`, `Rtsne` and `umap` came in behind the models and
  embeddings that call them through it. The comparison, diagnostic and
  visualisation functions still use base R only, so nothing in the contract that
  existed before this release acquired a dependency.

## Documentation

* README rewritten around one set of scenarios, in three parts: comparison,
  modelling, and dimension reduction. All twenty figures are drawn by
  `tools/render_readme_figures.R` from the calls the text quotes, on one seed and
  one correlation structure, so a figure cannot disagree with a number beside it.
  The two-group example moves to thirty features with eight planted each way, and
  the same `make_block_cor()` blocks run through all five models, which lets the
  correlated null predictor be followed from section to section instead of
  explained once.

# STATassist 0.5.0

## New features

* New `draw_heatmap()` draws one cell per feature and sample from the same wide
  input the comparison functions take, with the sample groups as a coloured
  strip above the columns and a dendrogram on each axis that was clustered. The
  input is transposed, so features run down the rows the way an expression
  heatmap is read, and `feats` selects which of them to draw and in what order.
  Nothing was added to `Imports`: the cells, the group strip and the two trees
  are `stats::heatmap()`, and the colour key it does not draw is overlaid beside
  them afterwards — an upright bar with its numbers to its right and the group
  legend, under a bold title, to the right of those, all starting level with the
  top row of cells. The strip they go in is as wide as they measure rather than
  a fixed share of the device, and it sits against the feature labels rather
  than at the edge of it. The bar is sized in units of text, so it stays a bar
  of about an inch rather than growing with the device.

* The clustering is done before `stats::heatmap()` is called and the trees are
  handed to it already built, so the tree that comes back on the result is the
  tree the plot shows. It also keeps the scaling ahead of the clustering, which
  is what lets a feature with no variance be centred instead of divided by zero.
  At least two features and two samples are needed, which is what
  `stats::heatmap()` can draw a grid and a tree from.

* Features are z-scored across the samples by default. One colour scale is
  shared by every cell, and features are not measured on a common scale, so
  without it a single high-abundance feature takes the whole range and the rest
  of the plot is left white. `scale` also offers `"sample"` and `"none"`. The
  plot does not name which one ran, but the numbers beside the colour key are
  what it produced, and `matrix` on the result is the scaled data.

* The clustering comes back on the result, as the `hclust` objects, the two
  permutations, and the scaled matrix in the order it was drawn. What the
  picture claims about the data can therefore be checked instead of eyeballed.

* Missing values are drawn as grey cells rather than dropped, and a feature with
  no variance is centred rather than divided by zero. An axis whose distances
  are not all defined, which happens when a pair of features shares no observed
  sample, keeps its input order and says so instead of failing.

# STATassist 0.4.0

## Breaking changes

* `group_lv[1]` is now the reference in `compare_two_groups()` as it already
  was in `compare_multiple_groups()`. Differences read
  `group_lv[2] - group_lv[1]`, ratios read `group_lv[2] / group_lv[1]`, and
  `alternative = "greater"` asks whether `group_lv[2]` is the larger group, so
  `mean_diff`, `hl_shift`, `trim_diff`, `fold_change`, `log2fc` and the
  confidence intervals all come back with the opposite sign, and
  `relative_effect` is mirrored around 0.5. The first level used to be the
  numerator with two groups and the reference with three or more, which meant
  the same control handed to both functions in the same position pointed two
  different ways. One rule now covers both. Code that named the level of
  interest first has to swap the two labels to get the numbers it got before.

* `simulate_two_groups()` defaults to `group_lv = c("control", "case")` and
  plants its effects in the second level. A planted increase still comes back
  as a positive `log2fc`, since the control is now the reference. The draws
  consume the random stream in the order they always did, so a given `seed`
  still produces the same numbers; only the row order of `args$data` and
  `args$group` follows the new level order.

## Improvements

* `draw_grouped_boxplot()` and `draw_butterfly_hist()` need no new argument to
  put a control on the left. Both already draw the levels in the order
  `group_lv` gives them, so the reference-first rule places it there, and a
  boxplot and a comparison of the same call now agree on which group is which
  side of the difference.

* The per-level tables in `$diagnostics` of a two-group comparison are built
  reference first, matching `design$group_lv`, rather than in the order the
  tests consume the two samples.

# STATassist 0.3.0

## New features

* New `draw_forest_plot()` exports the drawing that was reachable only through
  `plot()` on a comparison result. The name says what is drawn, so the plot can
  be found the way every other visualisation in the package is found, and it
  takes the same arguments the method took: `test`, `type`, `sort_by`, `dark`
  and the label and colour controls.

* `plot()` on a `sa_comparison` now calls `draw_forest_plot()` and carries no
  logic of its own. Existing code is unaffected: the two entry points draw the
  same rows and resolve `type = "auto"` to the same view.

* `draw_forest_plot()` gains `feats`, which names the features to draw and the
  order to draw them in, from the top of the plot down. It replaces the earlier
  `feature`, which selected one feature in the post-hoc view only. The
  post-hoc view now takes several features at once and labels each contrast
  with the feature it belongs to; a feature whose omnibus test never qualified
  it for the pairwise stage is reported in a `message()` rather than drawn as
  an empty row. The selection is applied before `type = "auto"` resolves, so
  the view is chosen from the rows that will actually be drawn.

* `draw_forest_plot()` gains `use_adjusted`, matching `draw_volcano_plot()`.
  With `use_adjusted = FALSE` the plot reads `pval` instead of `pval_adj`, and
  the colouring, the sorting, the p-value view and the legend follow, so the
  plot always names the p-value it used.

* `draw_forest_plot()` gains `xlim`, which fixes the x axis range instead of
  deriving it from the values being drawn, so two plots can be read against
  each other. The interval bounds are clamped to the range given, the way they
  already were to the derived one.

* The legend of `draw_forest_plot()` moved into a narrow panel of its own on
  the right, the arrangement `draw_grouped_boxplot()` already used, where it
  can no longer cover the rows it describes. Its entries read `Significant` and
  `Not significant` under a title carrying the threshold, and the new
  `cex.legend` sizes it independently of the axis annotation. A long row label
  widens the left margin only as far as half the panel; past that the labels
  are shrunk to fit rather than squeezing the plot.

* New `simulate_multiple_groups()` generates log2-scale abundance data for one
  control group and any number of treatment groups, and returns the effects
  planted in it. Its `args` element is named after the arguments of
  `compare_multiple_groups()`, so the comparison is one `do.call()` away, and
  `paired = TRUE` builds a repeated design with `id` in place of independent
  groups.

* The sizes are given per group. `n_control` is the control group and `n_treat`
  is one size per treatment group, so `n_treat = c(50, 40, 30)` is three
  treatment groups of different sizes and its length is the only statement of
  how many there are. A single size is spread over the groups when `group_lv`
  says how many there are, and labels and sizes that count differently are an
  error rather than a guess. Under `paired = TRUE` every group holds the same
  subjects, so unequal sizes are rejected.

* `n_up` and `n_down` default to a fraction of `n_feats` rather than to fixed
  counts, so asking for fewer features plants fewer effects instead of failing.
  At the default `n_feats = 100` they are the 15 and 15 the other defaults were
  tuned against.

* The answer comes back in three tables, because a multi-group comparison has
  two stages that fail separately. `truth` has one row per feature and scores
  the omnibus tests and `$effect`; `truth_group` has one row per feature and
  level, holding the delta, the centre, the standard deviation and the group
  size behind each one; and `truth_contrast` has one row per feature and pair
  of levels, in the row order and the `group1 - group2` direction the post-hoc
  tables use, so scoring the pairwise stage is a `merge()`.

* Each planted feature is given one of three effect shapes, in proportions set
  by `pattern_mix`: `"all"` moves every treatment group alike, `"gradient"`
  moves them in a dose-response ramp, and `"single"` moves one of them and
  leaves the rest at exactly zero. They are recovered at visibly different
  rates by the same omnibus test, which is the point of planting more than one.

* A comparison result carries a new `pairwise` slot, holding the numbers of
  `posthoc` rearranged so that one contrast can be read on its own. It is keyed
  first by test and then by contrast label, so with four levels
  `res$pairwise$anova_test` is a list of six tables and
  `res$pairwise$anova_test[["treat_1 - control"]]` is that one comparison
  across every feature. Empty exactly where `posthoc` is empty, which is the
  two-group and the one-sample scenario and any run with `posthoc = FALSE`.

* Those tables are rectangular where `posthoc` is ragged. Each holds every
  feature, in the row order the rest of the object uses, so contrasts line up
  against each other and against the omnibus tables by position; a feature that
  did not clear `posthoc_alpha` is present with its inference columns `NA`
  rather than absent. `posthoc` keeps the old shape, where an absent row means
  the question was never asked.

* They add `fold_change` and `log2fc`, which no post-hoc procedure reports.
  The ratio divides `group1` by `group2`, the direction `estimate` already
  reads in, so the two agree in sign within a row. It is a different quantity
  from `effect$log2fc`, which is the most extreme level over the reference, and
  it is filled even for a feature that never entered the pairwise stage, since
  dividing two centres does not require a test to have been run.

* `estimate_significance()` gains `by`. The default, `by = "omnibus"`, is what
  it has always done. `by = "contrast"` reads the pairwise stage instead and
  returns one verdict table per contrast, each carrying that contrast's own
  `log2fc` and p-value and enough attributes to be handed straight to
  `draw_volcano_plot()`. The adjustment axis differs between the two: the
  omnibus reading corrects across features, the contrast reading reuses the
  `pval_adj` the pairwise stage produced across the contrasts within a feature,
  unless `adj_type` names a method, which corrects across the features of that
  one contrast.

* The x axis of a volcano plot names what its `log2fc` compares when the answer
  is not a pair the caller chose. A multi-group omnibus verdict carries the one
  fold change its `effect` table holds, the level furthest from the reference
  rather than any named level, and the axis now reads
  `log2 FC (most extreme level vs setosa)`. A two-group, one-sample or
  single-contrast verdict does compare two fixed centres and keeps the plain
  `log2 FC`. The new `xlab` argument overrides either, which `...` could not do
  while the label was fixed inside the call to `plot()`.

* `print()` on the verdict summarises the test, the two cutoffs, the adjustment
  and how many features cleared them, with one line per contrast under
  `by = "contrast"`, instead of dumping the table and its attributes.

## Bug fixes

* The x axis of `draw_forest_plot()` now runs the width of the panel. It was
  drawn by `plot()` and so stopped at the outermost tick, which read as a
  broken axis whenever an interval end reached past the last labelled value.

* `draw_forest_plot()` and `draw_grouped_boxplot()` put back only the graphical
  parameters they set, instead of restoring a whole
  `par(no.readonly = TRUE)` snapshot. That snapshot carries `fin`, `pin` and
  `mai`, which are absolute sizes, so restoring it pinned the figure to the
  size the plot happened to be drawn at and the next plot on a resized device
  came out small in a corner of it. A panel grid the caller set up with
  `par(mfrow = )` still survives the call.

## Breaking changes

* A post-hoc contrast now subtracts the reference rather than being subtracted
  from it. `group1` is the later of the two levels of `group_lv`, so a contrast
  reads `treat_1 - control` where it read `control - treat_1`, and `estimate`,
  `statistic`, `lower_conf` and `upper_conf` come back with the opposite sign.
  The old arrangement put `group_lv[1]` first, which is right for two groups,
  where the first level is the numerator, but with three or more the first level
  is the reference: a feature the treatment raised then had a positive
  `effect$log2fc` and a negative post-hoc `estimate`, so a volcano plot and a
  forest plot of one comparison disagreed on which way it moved. They now agree,
  and the labels match the `"later-earlier"` rows of `stats::TukeyHSD()`. The
  contrasts themselves, their row order, `pval`, `pval_adj`, `stderr`, `df` and
  every omnibus and `effect` column are unchanged. Code that indexed
  `res$pairwise[[test]]` or the `by = "contrast"` verdict by contrast label has
  to name the label the other way round.

* A comparison result no longer carries `schema_version`, and neither does the
  `sa_diagnosis` object of `diagnose_distribution()`. It named the version of
  the object's layout for a consumer in another language, so it was the first
  thing a user read for something no R code ever looked at.
  `Configuration/schema/sa_result.schema.json` remains the normative definition
  of the layout.

* `posthoc` and `pairwise` are present only where a pairwise stage actually
  ran. A two-group or one-sample comparison leaves both out, as does
  `compare_multiple_groups(posthoc = FALSE)`, rather than carrying an empty map
  for a question that does not arise. Reading an absent slot as
  `res$posthoc[[test]]` still gives `NULL`, and iterating
  `names(res$posthoc)` still yields nothing, so the usual ways of reading them
  are unaffected; code that asserted the slot exists has to check for it
  instead. A stage that ran and qualified no feature still reports a present,
  zero-row table, which is a different fact from a stage that never ran.

* `estimate_significance()` returns an `sa_significance` object rather than the
  verdict table itself. It holds two elements: `analysis_type`, the `analysis`
  of the comparison the verdict was read from, and `significance`, the table
  that used to be returned, or the list of one table per contrast under
  `by = "contrast"`. The scenario name is beside the table because `log2fc` does
  not answer the same question in all three scenarios, and a consumer that only
  ever saw the table had to hope the caller remembered which one produced it.
  `sig$is_signif` becomes `sig$significance$is_signif`, and
  `by_pair[["virginica - setosa"]]` becomes
  `by_pair$significance[["virginica - setosa"]]`. The cutoffs, the test name and
  the adjustment stay on the table as attributes, so `draw_volcano_plot()` reads
  them from wherever the table came from, and it accepts the new object as well
  as a bare table. Handing it the whole `by = "contrast"` reading now errors and
  names the element to plot instead.

* The first argument of `draw_volcano_plot()` is now `significance_result`
  rather than `x`, matching the `comparison_result` of
  `estimate_significance()`. The error messages name the argument, so they now
  say what to pass rather than repeating a letter. `draw_forest_plot()` takes
  `comparison_result` for the same reason; the `plot()` method keeps `x`,
  which is the name the generic fixes. Positional calls are unaffected; a call
  that named `x =` has to be renamed.


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
