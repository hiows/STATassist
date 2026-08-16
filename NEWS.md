# STATassist (development version)

# STATassist 1.0.0

Scoring the models. The five fitting functions returned a model and the README
scored it by hand, with a `pROC` snippet for the curve and a scatter written out
each time. `predict.sa_model()` already abstracts the five engines, so a layer
above it can score any of them side by side, and holding several models against
one baseline is where the metrics that need a pair — DeLong's test, the IDI and
the NRI — become available.

And the other half of the unsupervised family. The three reductions answered
where each point sits once the features have been squeezed into two dimensions,
which left the neighbouring question — which points belong together — with no
function to ask it. The four `cluster_*()` functions ask it, and they reuse the
reductions' row axis and input reading rather than starting a fourth of their
own, so a clustering and an embedding of the same frame are about the same rows.
`draw_dim_reduction_plot()` is what that shared axis was for: it paints the
labels one of them found onto the coordinates the other placed, beside the
grouping that was known before either ran.

And the pair of features, which nothing had asked about. Every function so far
reduced one feature at a time or every feature at once; the question of how two
of them move together had no entry point, and a correlation matrix computed by
hand had no p-value beside it and no way of being drawn.

And the summary that had a table but no picture of its own.
`summarize_descriptive_stats()` has reduced every feature and group to a row of
locations and spreads since the first release, and seeing one of those columns
across a set of groups meant reshaping the table and calling `barplot()` by hand,
with the interval either side of each mean worked out again on the way.

## New features

* New `summarize_association_stats()` reduces every pair of features to the
  association between them. Pearson, Spearman and Kendall come back side by side
  on the same pairs, the way a comparison reports a parametric, a rank-based and
  a robust test side by side, so that a linear coefficient and a monotonic one
  disagreeing is a fact about the data rather than about which call was made.

  Each method gets a slot of four features-by-features matrices: `corr` from
  `stats::cor()`, `pvalue` from `stats::cor.test()`, `adj_pvalue` adjusted across
  the pairs, and `n`, the observations the pair shared. Beside them `design`
  records the call.

  The adjustment runs over the pairs that produced a p-value rather than over
  every cell of the triangle. A pair `stats::cor.test()` refused, having fewer
  than three shared observations or no variance on one side, is not a test that
  was performed, and counting it would shrink the other pairs for a comparison
  that never happened. Such a pair comes back `NA` rather than aborting the
  screen. The diagonal is a property of the matrix and not an estimate: `corr` is
  1 there and both p-values are `NA`.

  It is not an `sa_comparison`. This is a screen, the pairwise counterpart of
  `summarize_descriptive_stats()`, and it returns a plain list for the same
  reason that one returns a plain data.frame.

* New `draw_corrplot()` draws that matrix, and adds three decisions to
  `draw_heatmap()` rather than a second drawing engine. Nothing is standardised,
  a coefficient already being on a common scale. The colours are fixed at -1 to
  1, so the same colour means the same strength from one plot to the next. And
  the two axes hold the same features, so they are clustered once and permuted
  together; a symmetric matrix clustered twice can come back with its rows and
  columns in different orders, and the diagonal then wanders off the diagonal.

  Given the p-values, the pairs that did not clear `sig_level` are drawn as blank
  cells. The blanking happens after the clustering: a cell removed first would
  change the tree, and the picture would no longer be the matrix the reader is
  being shown. The diagonal is never blanked, a feature not being tested against
  itself, which is why the diagonal carries no p-value in the first place.

  The distance is `1 - cor()`, the same one `cluster_hclust()` and
  `draw_heatmap(dist_method = "correlation")` mean, so a corrplot and a heatmap
  of the same features group them the same way.

* New `draw_grouped_barplot()` draws a column of `summarize_descriptive_stats()`
  as clusters of bars, one cluster per feature and one bar per group level. It is
  the summary counterpart of `draw_grouped_boxplot()`: a box shows the
  distribution a group's observations have, a bar shows one number standing for
  them, which is what a figure wants when the point is about a location rather
  than a spread. The heights are that function's own column rather than a second
  pass over the data, so a bar and a row of the table are the same number.

  `errorbar` is read under `mainbar` rather than independently of it, because
  only two of the summary columns are locations an interval either side says
  anything about. A mean takes `"se"`, `"sd"` and Student's `"ci"` at
  `conf_level`. A median takes `"ci"` alone, the notch interval
  `median +/- 1.58 * IQR / sqrt(n)`, which is the interval
  `draw_grouped_boxplot()` already returns as `median_confidence_stats`, so a bar
  and the notch of the box beside it are the same width on the same data. Every
  other height — a count, a spread, a shape — has no second quantity to be
  uncertain by and takes `"none"`. Asking for one of the others is an error
  rather than a silently dropped interval, the alternative being a figure that
  answers a question other than the one it was asked.

  A bar is read against the zero it stands on, so a derived `ylim` always holds
  zero and puts its headroom on the side the bars run to; heights that go both
  ways, as `"skewness"` does, get the baseline drawn as well. A height the
  summary could not compute leaves its bar blank rather than shifting the ones
  beside it, which is what makes a group too small for a shape estimate visible
  instead of absent.

* New `evaluate_regression_models()` scores one or more fitted regressions on
  the same held-out rows. `$metrics` reports `cor`, `r_squared`, `rmse`, `mae`,
  `bias` and the calibration line, and `$comparisons` reports each model less
  the baseline. `cor` and `r_squared` are both there because they answer
  different questions: `r_squared` is `1 - SSE/SST` on the scored rows, so it is
  the variance the predictions actually removed and can be negative, while
  `cor^2` is what that would be after rescaling by a line fitted to those same
  rows. The gap between them is what `calib_slope` and `calib_intercept` report.

  There is no p-value beside the differences. A difference of held-out errors
  has no null this package is in a position to state, and the two numbers it is
  made of are in `$metrics`.

* New `evaluate_classification_models()` does the same for a two-class outcome.
  `$metrics` reports `auc` with a DeLong interval, `brier`, and `accuracy`,
  `sensitivity` and `specificity` at a stated `threshold`; `$curves` holds the
  ROC operating points as a rectangular table rather than a curve object.

  `$comparisons` asks three questions of each model against the baseline rather
  than picking one. `delta_auc` asks whether the ranking improved and is tested
  by DeLong's paired test. `idi` asks how much further apart the two classes'
  predicted probabilities moved, which an AUC is blind to, since reordering
  nothing leaves it untouched. `nri` asks how often a probability moved the
  right way, counting direction only, so a model that helps many rows slightly
  and hurts a few badly scores well here and badly on the IDI.

  The class order is the fits'. `outcome_lv` and `control_label` are checked
  against the order the models were fitted with rather than used to change it,
  since a fitted classification predicts the probability of its own second level
  and cannot be re-pointed after the fact.

* Both return the new `sa_performance` contract, a plain list of `analysis`,
  `models`, `design`, `parameters`, `predictions`, `metrics`, `comparisons`,
  `curves` and `metadata`. No engine object is stored: the calibration line is
  two numbers rather than the `lm` that produced them, so the whole object
  writes out as JSON.

  Every model is scored on the rows **all** of them could predict.
  `predict.sa_model()` answers `NA` for a row incomplete across that model's
  predictors, so models fitted on different predictor sets otherwise come back
  with different rows filled in, and DeLong's test, the IDI and the NRI are
  paired statistics with no meaning across different rows. One message reports
  how many rows went and why.

* New `draw_prediction_plot()` and `draw_roc_curve()`, with
  `plot.sa_performance()` dispatching to whichever the result calls for. The
  scatter draws the identity line and the calibration line from the two numbers
  `$metrics` holds rather than fitting again, so the picture and the table
  cannot drift apart. Its `type` is `"panel"` past one model, because two clouds
  of points on shared axes make a third that belongs to neither;
  `type = "overlay"` with `points = FALSE` compares the calibration lines alone.

* The AUC, ROC, DeLong, IDI, NRI and Brier kernels are written out in the
  package rather than taken from `pROC`. `pROC` covers the first three and none
  of the last three, and DeLong's test has no counterpart in `scikit-learn` or
  `scipy` for the Python port to call, so depending on it would buy the easy
  half and leave the hard half to be written twice against two sets of defaults.
  `pROC` is in `Suggests` as the oracle the tests check these against. Runtime
  dependencies are unchanged.

* New `cluster_hclust()`, `cluster_kmeans()`, `cluster_dbscan()` and
  `cluster_snn()`, the second half of the unsupervised family. The reductions ask
  where each point sits; these ask which points belong together. There are four
  because they disagree about what a cluster is, and the disagreement is the
  information: the first two are told how many groups to find and will always
  find that many, while the other two are told how dense a group has to be and
  derive the count from that, so they can return two clusters, or nine, or none,
  and they can refuse to place a point.

  `cluster_hclust()` returns the tree as well as the cut, so
  `cutree($fit, k = 5)` asks for another cut without rebuilding the distance
  matrix. It is also the only one of the four with a choice of distance, because
  it is the only one where the choice is real: k-means minimises squared
  Euclidean distance by construction and the two density methods measure a radius
  in it.

  `cluster_kmeans()` runs `nstart = 25` rather than `stats::kmeans()`'s own 1,
  which makes the answer a search rather than a draw; the override is reported in
  `engine`. `seed` makes a run repeatable and restores the random stream
  afterwards.

* All four return the new `sa_cluster` contract, a plain list of `analysis`,
  `points`, `design`, `parameters`, `assignments`, `clusters`, `engine`, `fit`
  and `metadata`. The row axis is `points`, the same one `sa_reduction` uses, and
  the four read their input through the same helper the three reductions do — so
  the rows are the same rows and an assignment can be painted straight onto a set
  of scores. `cluster_scale` sits where `embedding_scale` does and chooses the
  same way between clustering the samples and clustering the features.

  Noise is cluster `0`, which is `dbscan`'s own convention. `clusters` never has
  a row for it, so `nrow(clusters)` is the number of groups found and
  `design$n_noise` is what did not join one; a partitioning method reports
  `n_noise = 0` and both shapes read down the same columns. Every point being
  noise is a possible answer rather than an error, and it says so.

  `assignments$silhouette` is the one number that compares across all four, since
  coordinates from two methods share no scale but a ratio of distances does. It
  is Rousseeuw's definition written out in the package rather than taken from
  `cluster`: a singleton scores 0, noise scores `NA` and takes no part in anyone
  else's arithmetic, and a single cluster scores `NA` throughout.

* `cluster_dbscan()` derives both of its arguments when they are not given, and
  says what it derived. `min_pts` follows the textbook `d + 1`, floored at 4 and
  capped at half the points, since this package's usual input is wider than it is
  tall and a threshold above half can only ever return one cluster. `eps` is a
  radius in the units of the data and so can have no constant default; it is
  taken as the 95th percentile of the k-distance curve, which reads as "assume
  about one point in twenty is noise". `parameters$eps_source` records whether
  the value was supplied or derived.

  `cluster_snn()` is shared nearest neighbour clustering: each point keeps its
  `k` nearest neighbours and two points are linked when their neighbour lists
  overlap. Its `eps` counts shared neighbours where `cluster_dbscan()`'s is a
  distance. Both engines call the argument `eps`, so renaming one here would put
  this documentation at odds with theirs; the two say so instead.

* `dbscan` is a new dependency, the engine behind those two. `draw_heatmap()` now
  measures its dendrogram distances through the same function `cluster_hclust()`
  does, so the tree it draws and the tree that gets cut are one tree.

* New `draw_dim_reduction_plot()`, which `plot()` on an `sa_reduction` now calls.
  Two coordinates of a reduction as a scatter, with a clustering of the same
  points taking the colours and a grouping that was known already taking the
  shapes.

  Two channels rather than one, because the question a clustering usually raises
  is whether it recovered a grouping nobody told it about, and one channel would
  answer that by hiding it: either the clustering or the grouping would have to
  be left out, or the two crossed into a single set of labels whose count is the
  product of theirs. On two channels it is read directly — one colour per shape
  is a clustering that found the groups, and a shape split across colours is a
  group the data does not see as one thing. Give one and it takes the colours on
  its own; give neither and the points are drawn in the foreground colour.

  The two objects have to be about the same points, which is checked rather than
  assumed: both contracts carry `points` and both read their input through the
  same function, so a mismatch means two different frames were reduced and
  clustered and is refused instead of lined up by position. Points that joined no
  cluster are grey rather than given a colour beside the others, since a point
  left out is the absence of a cluster and not one more of them, and a PCA's axes
  carry the share of the variance they explain while an embedding's do not.

  `group_lv` here orders the levels and selects no rows, which is the opposite of
  what it does in `draw_grouped_boxplot()`. Dropping a row at this point would
  erase a point the reduction had already placed, so a level `group` uses and
  `group_lv` leaves out is an error.

* `draw_dim_reduction_plot()` gains `cluster_lv`, one label per cluster for the
  legend. `NULL` keeps `#1`, `#2`, and so on; noise is still listed as
  `noise (n)` when present.

* `draw_dim_reduction_plot()`'s `col` now names group levels when only `group`
  is given, and its new `pch` argument names group shapes or, with no `group`,
  the points themselves. Defaults are unchanged: shapes from the built-in
  sequence and `16` on a plain scatter.

* Clustering on a reduction's `$scores` now reads sample names from its `points`
  column when the table has no row names, so a clustering of `PC1` and `PC2`
  lines up with the reduction it is drawn on.

* `make_block_cor()` gains `against`, which splits a block in two: the indices in
  `features` and the indices in `against` each correlate at `cor` among
  themselves and at `-cor` with the other side. A block that shares one value has
  a floor on how negative that value can be — `-1/(k - 1)` for `k` predictors, so
  three of them could not disagree past -0.5 and four past -0.333 — which left a
  strong negative correlation unavailable at the sizes worth simulating. A split
  block is one factor with a sign per predictor rather than a demand that every
  pair disagree, so it holds for any `cor` below 1 whatever its size.

  That floor, and the ceiling of 1 beside it, are now checked block by block
  before the assembled matrix is factorised. A `cor` no block of that size could
  hold is reported with the bound it would have to clear instead of as a matrix
  that is not positive definite, and the message for the assembled matrix carries
  its smallest eigenvalue and can now say that every block held on its own, since
  by then only a `default_cor` the blocks cannot sit inside is left to reach it.

* `draw_corrplot()` now names `cex.anno`, `cex.axis`, `cex.main` and
  `cex.legend` on its signature. They were already forwarded through `...` to
  [draw_heatmap()], but were easy to miss in `?draw_corrplot`.

## Breaking changes

* `draw_corrplot()`'s first argument is now `cor_matrix` rather than `x`, since
  that is what it holds whether the caller hands in a matrix or the result of
  [summarize_association_stats()].

* `summarize_association_stats()` no longer returns `cor_matrix`. It was the
  `corr` of the first entry of `methods`, which is already available as
  `res[[res$design$methods[1]]]$corr`, and keeping both names for one matrix
  made it look as though `cor_matrix` always meant Pearson when it followed
  whatever came first in `methods`. [draw_corrplot()] still takes the result and
  reads the slot named by `method`, defaulting to the first.

## Bug fixes

* `make_block_cor()` no longer drops most of a block specification without a
  word. `list(features = 1:3, cor = 0.9, features = 4:6, cor = -0.4)` is one
  block with a repeated name rather than the two blocks it reads as, and `$`
  returns the first of a repeated name, so the matrix came back holding the first
  pair of values and nothing else. Repeated names, names a block has no use for,
  and elements with no name at all are now errors, and the message shows what the
  several-blocks form looks like.

# STATassist 0.8.0

The interaction plot, and the slot it needed. A factorial result reported that a
term was significant without keeping anything a reader could see the term in:
`$effect` reduces the grid to two cells and a term test reduces it to a p-value,
and neither says which way the lines run. The kernel was already computing the
cell means and throwing them away after the post-hoc stage, so the picture cost a
slot rather than a second pass over the data.

## New features

* `compare_factorial_groups()` now returns `$cells`, one row per feature and cell
  of the crossed grid, holding the cell mean, the within-cell `sd`, the count and
  a pooled `se`, plus one column per factor named after the factor and holding
  the level. `mean` is the arithmetic mean the crossed model was fitted on rather
  than a centre on the `fc_mean` scale, so a plot of this table and the F tests
  beside it describe one fit.

  `se` is `sqrt(ms_error / n)`, pooled over the whole model rather than taken
  within the cell. That is what makes the variance of any marginal mean
  recoverable from this column alone: over a set of cells it is
  `sum((1 / length(S))^2 * se^2)`, the expression the Tukey stage scales its
  contrasts by, so an error bar drawn from this and a contrast in `$posthoc`
  cannot disagree. Because the factor columns are named after the factors, the
  six fixed names are reserved and a factor may no longer be called `features`,
  `cell`, `n`, `mean`, `sd` or `se`.

* New `draw_interaction_plot()` joins the cell means of one factor across the
  levels of another, one line per level of the tracing factor, which is the
  picture an interaction term is a test of. Factors beyond the two being drawn
  are averaged away unweighted, the same marginal mean the post-hoc stage
  contrasts, so the gap between two points of one line is the `estimate` its
  contrast reports.

  `type` divides the three ways a design of more than two factors can be shown.
  `"pairwise"` draws one pair and spends its panels on features, so several read
  side by side. `"matrix"` draws every pair at once in an upper triangle, the row
  tracing and the column on the x axis. `"facet"` keeps the levels of a third
  factor in panels of their own instead of averaging them away, which is what
  shows a two-factor interaction that itself depends on a third; the last two
  spend their panels on the factors and so draw one feature at a time.
  `type = "auto"` reads the arguments: naming `facet` asks for the facet view,
  naming `x` or `trace` for the pairwise one, and naming none of them gives the
  pairwise view for two factors and the matrix for three or more.

  `errorbar` is off by default. A bar here is the standard error of the mean
  being drawn rather than of the difference between two of them, and the
  difference is what an interaction is about, so two overlapping bars do not
  settle the term test.

## Bug fixes

* `compare_factorial_groups()` no longer shortens its list of fits when a feature
  could not be fitted. `fits[[i]] <- NULL` deletes the element rather than
  emptying it, and a following success happened to restore the length, which is
  why this only surfaced where the last feature failed: `$terms` then lost its
  statistics columns altogether rather than filling them with `NA`.

The categorical family. Two columns of labels are a contingency table, not a
feature-wise comparison, so the result is an `sa_categorical` rather than an
`sa_comparison` and `estimate_significance()` refuses it. The mosaic is what
reads it, in the place a volcano plot holds for the numeric scenarios, and
`estimate_categorical_significance()` is the verdict beside it.

What the contract turns on is that "expected" is not a property of a table. It
is a property of a table and a claim about it, and the three designs here make
three different claims. So the result names its own claim and everything read
off it -- the residuals, the diagnostics, the shading of the mosaic -- is read
under that one.

## New features

* New `compare_categorical_groups()` crosses two categorical variables, or
  reads the columns as repeated binary conditions on the same row. An
  independent design returns the chi-square test of independence beside
  Fisher's exact test; a matched design of two conditions returns McNemar's
  test and three or more return Cochran's Q. There is no argument naming a
  test. `$association` carries Cramer's V and the contingency coefficient on
  every independent table, phi and the odds ratio on a 2 x 2, the three paired
  measures on a matched 2 x 2, and Kendall's W on Cochran's Q. `control_label`
  points the reference the way it does for a crossed design, and that is the
  direction of the odds ratio.

* `design$null` names the hypothesis the result is about, one of
  `"independence"`, `"symmetry"` and `"marginal_homogeneity"`, and
  `$cells$expected` is read under it. A cell is expected at the product of its
  margins under independence, at the average of it and its transpose under
  symmetry, and at the pooled response rate under marginal homogeneity. The
  matched case is the one where this is load bearing: McNemar's test asks about
  symmetry, so the diagonal is expected at exactly what it holds and carries no
  residual, and the squared Pearson residuals of the two discordant cells sum to
  `(b - c)^2 / (b + c)`, McNemar's uncorrected statistic exactly. The cell table
  and the p-value beside it are one piece of arithmetic read two ways.

* Fisher's exact test returns an `NA` p-value, with a message, where its
  enumeration cannot be walked. The algorithm visits every table with the
  observed margins and a large r x c table has more of them than its workspace
  holds; that is a limit of the enumeration rather than a fault in the data, and
  failing the whole call over it would take the chi-square result down with it.
  `simulate_p_value = TRUE` is what answers there.

* `$cells$std_residual` is `NA` under symmetry. Its variance correction is
  derived for a table held against its own margins and has no counterpart there,
  so the column is empty rather than holding a number that looks referable to a
  standard normal and is not.

* `as.table()` on a result gives the contingency table the tests were run on,
  built from `$cells` on request rather than stored beside it. Every slot of the
  result is therefore a scalar, a character vector, a named list or a data.frame,
  which is what the `sa_result` contract asks for and what lets the whole object
  be written out as JSON. `Configuration/schema/sa_categorical.schema.json`
  records the shape.

* `max_levels`, default 20, refuses a variable taking more distinct values than
  that, naming the variable and its count. `as.character()` turns a continuous
  measurement into one label per observation, and a table with a cell per
  observation is not something a test of association has anything to say about.
  It is checked against the levels actually used, so naming three levels of a
  fifty-valued column in `category_lv` is a way through.

* New `simulate_categorical_groups()` plants the association the comparison
  estimates. `assoc = 0` is the product of the margins exactly, so it is null in
  the strict sense, and equal transition probabilities are the same for a matched
  design. A matched design plants the paired odds ratio as a ratio of transition
  probabilities, and `truth_cell` carries `p_symmetric` and
  `expected_symmetry_n` beside the independent share, so the planted table can be
  scored against the null its own tests use. Both designs key `truth_cell` on
  `c("row_level", "col_level")`, so it merges with `$cells` with neither side
  renamed.

* New `estimate_categorical_significance()` reduces a contingency table to
  verdicts, which the whole-table reading of the result cannot be reduced to: an
  association is not signed, so there is nothing to put on an effect axis beside
  the p-value. A cell is different. It was expected at some count and observed at
  another, so `observed / expected` says how far it moved and which way, and it
  is defined whatever the shape of the table. That ratio is `lift`, the quantity
  `simulate_categorical_groups()` already plants under the same name, and
  `log2_lift` is the effect axis. The p-value is the cell's own standardized
  residual referred to a standard normal, and the cells of one table are one
  family, so `adj_type` adjusts across them — the first adjustment rather than a
  replacement for one, there being no `pval_adj` column on a categorical result.

  So the same `log2_lift_cutoff` means the same thing on a 2 x 2, a 2 x 3 and
  anything larger. It defaults to 1, deliberately the number
  `estimate_significance()` uses, though a doubling is a stricter demand of a
  cell than of a fold change.

  `by = "table"` is the second reading: one row, an association measure beside an
  omnibus p-value and no adjusted column, one table being one question.
  `measure = "auto"` takes the odds ratio on an independent 2 x 2, Cramer's V on
  anything larger, the paired odds ratio under symmetry and Kendall's W under
  marginal homogeneity. `effect_cutoff` is `NULL` by default, so the p-value
  alone decides: the conventional thresholds for Cramer's V are conventions
  rather than facts about the measure, and a default is where a convention is
  hardest to notice.

  A matched pair of conditions has no cell reading and says so. Its standardized
  residuals are `NA` throughout, for the reason the column is empty there at all,
  so there is no p-value axis to read; `by = "table"` answers instead. Three or
  more matched conditions are read off a condition-by-response table whose
  arithmetic is that of independence, so the cell reading works there.

  The verdict is an `sa_categorical_significance` and deliberately not an
  `sa_significance`: its rows are cells rather than features, so
  `draw_volcano_plot()` refusing it is the point of the separate class. It keys
  on `c("row_level", "col_level")`, which merges with `$cells` and with
  `simulate_categorical_groups()$truth_cell` with neither side renamed, so the
  estimated lift can be scored against the planted one.

* New `draw_mosaic_plot()` splits the x axis by the first variable's marginal
  shares and each strip by the second variable's conditional shares, so the area
  of a tile is the cell's share of the table. `plot()` on an `sa_categorical`
  calls it.

  It reads a `compare_categorical_groups()` result and nothing else, the way
  `draw_interaction_plot()` reads a factorial one. The shading and the expected
  lines are statements about the null the result was tested against, which a
  bare table does not carry, so the levels that take part and the reference they
  are read against are settled by the comparison rather than restated here.

  The shading reads the residual of `design$null`, so a matched design is a
  picture of departure from symmetry rather than from a hypothesis nothing in the
  result has a p-value for, and `residual = "standardized"` is refused there. A
  dashed segment marks each strip where that null would have cut it, which under
  independence is the same height in every strip and under symmetry is not. The
  level names on the y axis sit on the reference strip, the first one, since no
  single set of positions can label strips that are cut at different heights.
  `anno_cells` writes as much of the count and the conditional percentage as the
  tile has room for, measured against the label rather than against a fixed
  fraction of the plot, and the residual key sizes itself from its own text.
  `gap` is the width of one gap rather than of all of them together, capped so
  that the gaps never take more than two fifths of an axis however many levels
  ask for one.

The factorial family, built from the answer key inwards. The simulator came
first, which is the reverse of how every other simulator in the package was
written: the others read the names in their `args` off a `compare_*()` or
`fit_*()` that already existed, and here there was nothing to read them off. That
order is what makes the comparison checkable on arrival. The data generator was
shown correct against `aov()` before a line of the analysis depended on it being
so, and the analysis is now scored against the same answer key term by term.

## New features

* New `compare_factorial_groups()` analyses a crossed design as one model. Two
  factors are a two-way ANOVA, three a three-way one and more than three a
  factorial one, and those are three names for the same fully crossed linear
  model rather than three procedures. Which name applies follows from
  `length(factor_lv)`, so there is no argument for choosing it; it is reported in
  `design$anova_type` and in the label of the test. This is the one comparison in
  the package that fits a single model instead of running competing tests side by
  side, because what multiplies in a factorial design is not the procedures but
  the questions.

* The result grows an axis for those questions. `$terms` holds one row per
  feature and per model term, every main effect and every interaction of every
  order, with `df`, `ss`, `ms`, `f_stat`, `eta_sq`, `partial_eta_sq`,
  `log2_effect` and a `pval_adj` adjusted **within each term** across features,
  since a term is one family and pooling two terms would correct each for the
  other's multiplicity. Its `terms` and `term_order` columns are the ones
  `truth_term` carries, so a result and an answer key merge on
  `c("features", "terms")` with neither side renamed.

* `log2_effect` is the size of a term, with a sign. `eta_sq` and
  `partial_eta_sq` say how much variance a term accounts for and neither has a
  direction, so neither can be an effect axis. This one is the largest ANOVA
  component of the term, taken by decomposing `log2()` of the same cell centres
  `$effect` is built from, which is the decomposition
  `simulate_factorial_groups()` records in `truth_term$max_abs_delta`: the
  measurement and the answer key are one quantity, and over eight seeds the two
  agree at r = 0.78, per term between 0.71 and 0.90. Read it as a deviation from
  what the rest of the model predicts rather than as a fold change — a two-level
  factor whose levels differ by one log2 unit contributes -0.5 and +0.5 — and
  note that a cutoff meant for a fold change is stricter here than it looks.

* `$tests$anova_test` keeps the whole-model F test at one row per feature, which
  is the one-way ANOVA that treats the cells as groups: a crossed model is the
  cell means model in another basis, so the two leave the same residuals and the
  term sums of squares add up to the cell sum of squares. Keeping that table in
  its usual shape is what lets `estimate_significance()`, `print()`,
  `draw_forest_plot()` and `draw_volcano_plot()` read a factorial result without
  knowing that it is one.

* Reading a factorial result **as** one is the other half of that. The default
  `estimate_significance()` reading is still `"omnibus"`, the whole-model F
  paired with the most extreme cell against the reference, as in every other
  scenario. A third reading, `by = "term"`, returns one verdict table per model
  term. Each term table is named after the term and ordered as `$terms` lists
  them, the same shape `by = "contrast"` already returns for pairwise contrasts.
  The p-value is that term's own and the effect size is its `log2_effect`. There
  is no choice of adjustment axis to make: both branches adjust across the
  features of one term, which is the family `$terms$pval_adj` was built over.

* An omnibus verdict of a factorial or multi-group comparison now carries
  `extreme_cell` or `extreme_level` beside `log2fc`, naming which cell or level
  was furthest from the reference on the log2 scale. [draw_volcano_plot()] labels
  the x axis accordingly — most extreme cell vs the reference cell for a crossed
  design, most extreme level vs the reference level for a multi-group one.

* `control_label` reaches the crossed design, in `compare_factorial_groups()` and
  in `draw_grouped_boxplot()`. A crossed design has one reference per factor
  rather than one in total, so it takes one level name per factor it points, as
  `list(treatment = "control", sex = "male")` or as
  `c(treatment = "control")` — the two shapes say the same thing. The named level
  moves to the front of its own factor, the factors it says nothing about are
  left alone, and the reference cell is where all of them land. It does more here
  than it does for a single factor, where `group_lv` is required and naming a
  reference can only re-point one: `factor_lv` is optional, and levels taken from
  the data arrive sorted, so `control_label` is the way to say which cell the
  fold changes divide by without listing every level of every factor. Moving the
  reference leaves the ANOVA where it was — a sum-to-zero coding spans the same
  space whichever level leads a factor — and moves the cell labels, `$effect`,
  the volcano's x axis and the direction of every post-hoc contrast. Both
  functions take it right after `factor_lv`, so a call that passed a later
  argument by position has to name it now.

* `draw_volcano_plot()` draws a term reading as a figure of panels, one per term,
  rather than asking which table to name. A crossed design decomposes into a
  fixed and small set of terms, three for two factors, and seeing which of them a
  feature responded to is a comparison between panels; the number of pairwise
  contrasts follows from the level counts and is arbitrary, so a contrast list is
  still sent back to name one. The panels share one rule and one pair of axes,
  `terms` chooses them — by default the two main effects of the first two factors
  and their interaction, with the terms left out named in a `message()` — and
  `panel_nrow` arranges them. The x axis of a term panel is labelled as an effect
  rather than as a fold change, because a component is not a ratio of two
  centres. Its first two arguments are new, so a call that passed `use_adjusted`
  by position has to name it now.

* On the `crossover` shape, whose main effects are exactly zero while its cells
  plainly differ, the term panels put the finding where it belongs: at
  `abs(log2fc) >= 0.5` and `adj_pvalue <= 0.05` the interaction panel calls half
  of them and the two main effect panels call 0.016 of their rows, against
  nothing planted there. The whole-model verdict calls the same features changed
  and cannot say more than that.

* `terms` is a new slot of the `sa_comparison` contract, and the first table in it
  that is not one row per feature. `sa_new_comparison()` takes `terms = NULL` and
  drops the slot when it is not given, the way `posthoc` and `pairwise` are
  dropped, so the three existing comparison functions return exactly what they
  returned before. `sa_result.schema.json` grows a `termTable` definition to
  match. There is no `pairwise` on a factorial result: its contrasts are indexed
  by factor and stratum rather than by a single label, so
  `estimate_significance(by = "contrast")` reports that the result has no
  pairwise stage, which is the true statement.

* `print()` on a factorial result leads with the design rather than the cells.
  Eight cell labels strung together with `vs` is not what a reader of a crossed
  design wants first; `factors : treatment (4) x sex (2)`, the ANOVA that ran, and
  then one line per term is.

* `ss_type` chooses what the term tests are built from, `"III"` by default,
  `"II"` or `"I"`. The three agree on a balanced design and part company when the
  cells are unequal. Type III makes a main effect a statement about the levels
  rather than about how many observations landed in each, which is what agrees
  with the unweighted marginal means the post-hoc stage compares and with the
  decomposition the simulator plants. Type I is what `stats::aov()` reports, and
  it is there so that these numbers can be checked against an outside
  implementation on unbalanced data.

* The pairwise stage answers both questions a factorial design has.
  A marginal contrast compares two levels with the other factors averaged away, a
  simple effect compares them inside one combination of the other factors, and
  `stratum` tells them apart, `NA` for the first. Both are Tukey-Kramer on the
  mean square error of the whole model rather than a second analysis of the same
  data, and the marginal means are unweighted, so a level is not pulled towards
  whichever combination was sampled most heavily. `posthoc_scope` chooses which
  of the two kinds runs.

* Which contrasts run is decided term by term rather than by the whole-model
  test: a factor's marginal contrasts need that factor's main effect to clear
  `posthoc_alpha`, and a simple effect needs the interaction with the factors held
  fixed to clear it. Gating on the whole model instead would compare the levels
  of a factor the model says has no effect, on the strength of a different factor
  that has one. There is no `posthoc_p_adjust`, since Tukey's p-values are
  already family-wise within a block.

* `within` is accepted and refused. A repeated-measurement factor needs a second
  error stratum, which is not implemented, so a non-empty value is an error that
  names what is missing rather than a design analysed as though the repeated
  measurements were independent. Leaving the argument out altogether would make
  `do.call()` on a within-subject simulation fail with `unused argument`, which
  says nothing.

* Scored against the simulator over eight seeds and 480 features, the term tests
  find a planted main effect about four times in five and a planted interaction
  about three times in four, calling an unplanted term about once in twenty. On
  the `crossover` shape, whose main effects are exactly zero while its cells
  plainly differ, the main effects are called at 0.06 and the interaction at
  0.94: the case a design read one factor at a time cannot see is the case this
  function is for.

* New `simulate_factorial_groups()` crosses any number of factors and returns the
  effects planted in them. `factor_lv` names the factors and their levels, so its
  length is how many are crossed and there is no separate argument for the count.
  The first factor is the primary one and the ones after it are the factors the
  treatment may or may not depend on. Two is the fewest there can be: one factor
  is `simulate_multiple_groups()`, and the error says so.

* `within` names the factors measured on the same subjects. Naming none gives a
  factorial ANOVA, naming all of them a fully repeated design and naming some a
  mixed one, so the three are one function rather than three. A subject sits in
  one combination of the between factors and is seen under every combination of
  the within ones, which leaves the within-subject rectangle complete. `args`
  then carries `id` and `within` the way `simulate_multiple_groups()` carries
  `id` and `paired`.

* The answer has a layer the one-factor simulators do not: `truth_term`, one row
  per feature and per model term, every main effect and every interaction of
  every order, with `is_effect` saying whether it was planted and `is_within`
  saying which error stratum tests it. Crossing factors turns the answer into a
  statement per term rather than per level, and this is the table that scores an
  ANOVA table row by row. Beside it are `truth` per feature, `truth_cell` per
  feature and cell, and `truth_contrast` per feature and pair of levels.

* The effect is planted in five shapes, chosen so that the terms come apart
  rather than move together. `main_only` moves the primary factor alone;
  `additive` moves it and a partner with parallel profiles, which makes an
  interaction call a false positive by construction; `interaction` makes the
  treatment effect depend on the partner; `nuisance_only` moves the partner and
  leaves the treatment at zero; and `crossover` is pure interaction, with **both
  main effects exactly zero** while the cells plainly differ. That last one is
  the case a design read one factor at a time cannot see, and it is exact rather
  than approximate because the effect is built in the space an ANOVA decomposes
  into and then shifted to put the reference cell at zero, which touches the
  grand mean and no term.

* `pattern_mix` is the same argument, and the same three profiles, that
  `simulate_multiple_groups()` has. It is on a different axis from `term_mix`
  here: one says which terms move, the other what the profile along a factor
  looks like, and the two are crossed at random while both sets of counts stay
  exactly what the weights ask for.

* An `aov()` on the defaults finds the treatment main effect about four times in
  five, which is the rate `simulate_multiple_groups()` was tuned to, and calls a
  term that was not planted about once in twenty. `interaction_scale` is set so
  that the `interaction` shape's interaction is found about three times in five,
  against about nine in ten for `crossover`: between trivially recovered and
  indistinguishable from noise, which is where a default belongs.

* `truth_contrast` carries both pairwise questions a factorial design has. A
  `stratum` of `NA` is the marginal contrast, averaged over the other factors,
  which is what a main effect is a statement about; anything else names the
  combination the contrast was taken inside, which is the simple effect. Under
  `crossover` the first is exactly zero everywhere and the second is not, which
  is the pair of statements the two kinds of row exist to keep apart.

* `args` is named for `compare_factorial_groups()`, so a simulated between-subject
  design is one `do.call()` away from being analysed and its `$terms` merges with
  `truth_term` without either side being renamed.

* `draw_grouped_boxplot()` draws a crossed design. `factors` and `factor_lv` are
  the arguments `compare_factorial_groups()` takes, and they are read by the same
  `sa_fact_layout()`, so the boxes are the cells the model fits rather than a
  second reading of the same data. Naming both `group` and `factors` is an error,
  since they are two ways of saying what the boxes are.

* A crossed design is drawn one panel per feature, with the factors after the
  first along the x axis and the **primary factor as the coloured boxes** inside
  each of those clusters, so that the two factors sit side by side in one panel.
  An interaction is then local: a treatment effect that reverses between the
  sexes is a pattern of colours that visibly flips a couple of centimetres away.
  Splitting the two factors between the legend and the panels, which is what the
  transpose does, means reading a colour profile in one panel against the same
  profile in another.

* `panel_by = "factor"` asks for that transpose: one panel per combination of
  the factors after the first, with the features along the x axis, which is what
  to ask for when the question is about the features rather than the crossing.
  The boxes and the returned statistics are the same either way; only their
  grouping into panels and clusters differs.

* `panel_nrow` now defaults to `NULL`, "let the arrangement decide": one row for
  panels over the factors, and a grid as near square as the panel count allows
  for panels over the features, where a row of ten of them cannot be read. A
  number given explicitly still wins.

* The columns of the returned statistics are the cell labels, the same strings
  `compare_factorial_groups()` puts in `design$group_lv` and
  `simulate_factorial_groups()` in `truth_cell`, so what the picture shows can be
  read against what the analysis found and against what was planted. A cell
  holding no observation keeps its column and its position, filled with `NA` and
  drawn blank, and is named in a message: an empty cell that shifted its
  neighbours along would relabel every box after it.

* Whether the panels share a y axis range follows from what a panel holds.
  Panels over the factors hold the same features and share one range taken from
  every value drawn, since panels of the same quantity scaled to their own
  values cannot be read against each other. Panels over the features hold
  different quantities, each on its own baseline, so a shared range would
  flatten every one of them: each keeps the range `graphics::boxplot()` gives it
  and carries its own axis. An explicit `ylim` is shared by every panel in
  either arrangement, and a single panel is left to `graphics::boxplot()` as
  before, so a single-factor plot is drawn and returned exactly as it was.

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
