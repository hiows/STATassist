#' STATassist: Run Every Applicable Statistical Test at Once
#'
#' STATassist is organised around comparison scenarios rather than around
#' individual tests. One function covers one situation and runs everything that
#' applies to it, returning feature-wise result tables with a shared column
#' layout. The guiding idea is that a group comparison should never rest on a
#' single test: a parametric, a rank-based and a robust procedure are reported
#' side by side so that disagreement between them becomes visible instead of
#' hidden, and the choice of what to report stays with the analyst.
#'
#' Everything is vectorised over features, so the same call serves one
#' measurement or several thousand.
#'
#' @section Comparison:
#' One function per situation, chosen by how many groups there are:
#'
#' \describe{
#'   \item{[compare_one_sample()]}{One sample against a hypothesised value: a
#'     one-sample t-test, a Wilcoxon signed-rank test and a proportion test.}
#'   \item{[compare_two_groups()]}{Exactly two levels: Welch's or paired t,
#'     Wilcoxon rank-sum or signed-rank, and Brunner-Munzel or Yuen's trimmed
#'     mean test for dependent samples. The fold change comes with them.}
#'   \item{[compare_multiple_groups()]}{Three or more levels, independent or
#'     repeated: one-way, Welch's, trimmed mean and Kruskal-Wallis omnibus
#'     tests, or repeated measures ANOVA and Friedman, each followed by the
#'     post-hoc procedure that shares its assumptions.}
#' }
#'
#' Direction is the order of `group_lv` in both of the grouped ones, and
#' `control_label` is the second way of stating it: the level it names moves to
#' the front and the rest keep the order they were given, so a fold change and
#' every post-hoc contrast can be turned around without the levels being retyped.
#'
#' @section Assumption diagnostics:
#' [diagnose_distribution()] reports the normality tests, the homogeneity of
#' variance tests and [screen_outliers()] together. The same checks are attached
#' to every comparison result as `$diagnostics`, computed on the observations
#' that were actually tested, so an assumption is never silently ignored. A
#' failed check never changes which tests run: it changes which member of the
#' reported family deserves the most weight, and that judgement stays with the
#' analyst.
#'
#' @section Significance and visualisation:
#' [estimate_significance()] reduces a comparison to one row per feature holding
#' both axes of a volcano plot and a multiplicity adjusted verdict, and
#' [draw_volcano_plot()] draws it. Because the significance table is derived from
#' the comparison object rather than assembled from loose vectors, the effect size
#' and the p-value beside it always describe the same observations and the same
#' direction. The verdict comes back as `$significance` beside the
#' `$analysis_type` it was read from, so what the `log2fc` column compares is
#' never in doubt.
#'
#' ```
#' res <- compare_two_groups(data, feats, group, group_lv)
#' sig <- estimate_significance(res, test = "t_test")
#' draw_volcano_plot(sig)
#' ```
#'
#' [draw_forest_plot()] reads the comparison object itself, drawing each
#' estimate beside its confidence interval, or the pairwise contrasts of a
#' multi-group result. It is what [plot()] on a `sa_comparison` calls.
#'
#' [draw_grouped_boxplot()] draws a grouped boxplot for the same wide-format
#' input and optionally returns the box summary statistics and median
#' confidence intervals behind the plot. [draw_butterfly_hist()] puts the two
#' group distributions of a single feature back to back on shared breaks.
#' [draw_heatmap()] takes that same input the other way round, drawing every
#' feature and every sample at once as one cell each, with the sample groups
#' annotated above the columns and a dendrogram on each axis it clustered.
#'
#' @section Descriptive summary:
#' [summarize_descriptive_stats()] reduces each feature to one row of sample
#' size, central tendency, dispersion, quartiles, outlier fences and
#' distribution shape, split by group level when a grouping vector is supplied.
#'
#' @section Result contract:
#' Every comparison returns a `sa_comparison` object, a plain named list of
#' scalars, character vectors and data.frames with an S3 class on top that only
#' supplies [print()] and [plot()]. No fitted model or other R-only object is
#' stored anywhere, so the result can be written out as JSON and rebuilt
#' elsewhere. Every test kernel is likewise a plain function of numeric vectors.
#' Both choices exist to keep a future Python implementation a transcription
#' rather than a redesign.
#'
#' Every result carries `effect` and `tests`, and one attaches `diagnostics`
#' whenever `diagnose = TRUE`. A comparison of three or more levels adds
#' `posthoc`, one table per test holding one row per feature and pair of levels,
#' and `pairwise`, the same numbers rearranged into one rectangular table per
#' contrast. Those two slots are there only when a pairwise stage actually ran,
#' so a two-group or one-sample result does not offer an empty version of them.
#'
#' @section Practising on a known answer:
#' [simulate_two_groups()] generates log2-scale expression data with a chosen
#' number of features moved up and down on purpose, and hands back the planted
#' answer with it. Its `args` element is named after the arguments of
#' [compare_two_groups()], so a comparison is one `do.call()` away and its
#' verdict can be scored against what was actually there.
#'
#' ```
#' sim <- simulate_two_groups(seed = 1)
#' res <- do.call(compare_two_groups, sim$args)
#' sig <- estimate_significance(res)$significance
#' table(planted = sim$truth$direction != "none", called = sig$is_signif)
#' ```
#'
#' [simulate_multiple_groups()] does the same for one control group and any
#' number of treatment groups, independent or measured on the same subjects.
#' A multi-group comparison has two stages that fail separately, so the answer
#' comes back in three tables rather than one: `truth` scores the omnibus tests
#' feature by feature, `truth_group` holds the effect planted in each level, and
#' `truth_contrast` scores the pairwise stage in the row order and the direction
#' the post-hoc tables use. Each planted feature also carries the shape of its
#' effect, since a change every treatment group shares and a change in one group
#' alone are found at different rates by the same test.
#'
#' ```
#' sim <- simulate_multiple_groups(seed = 1)
#' res <- do.call(compare_multiple_groups, sim$args)
#' ph <- merge(res$posthoc$anova_test, sim$truth_contrast,
#'             by = c("features", "contrast"))
#' table(differs = ph$is_diff, called = ph$pval_adj <= 0.05)
#' ```
#'
#' @section Machine learning:
#' [split_data()] draws the train/test partition every later step is fitted
#' inside. It stratifies on whatever the model will predict, so both halves
#' carry the balance of the whole data set, and it splits over sampling units
#' when `id` is given, so repeated measurements of one subject cannot appear on
#' both sides. This family is built on `caret` rather than written from scratch.
#'
#' [fit_linear_regression()] and [fit_logistic_regression()] fit what it hands
#' back. Each takes the wide data frame the comparison functions take, one column
#' of it as the outcome, and returns a `sa_model`: a coefficient table, the
#' goodness of fit of the model as a whole, and how the same procedure scored on
#' rows it had not seen inside each resampling fold. Cross-validation decides
#' whether the model is scored, not how it is fitted, so the coefficients do not
#' depend on the scheme. `outcome_lv` fixes the direction of a classification the
#' way `group_lv` fixes it for a comparison: the first level is the reference, so
#' the coefficients describe the odds of the second. `control_label` names that
#' reference on its own, here and in [fit_elastic_net()], [fit_rf()],
#' [fit_svm()] and [perform_rfe()], and naming it is also the one way to say that
#' a numeric column of zeroes and ones is two classes rather than two numbers.
#' Where it is read differently from the comparisons is in disagreement: an
#' `outcome_lv` holds the two classes and nothing else, so naming the other one
#' as the reference is an error, while the same argument re-points a `group_lv`
#' that carries the display order of every level besides.
#'
#' ```
#' sp <- split_data(data, stratified = "outcome", seed = 1)
#' fit <- fit_logistic_regression(sp$datasets[[1]]$train_data,
#'                                outcome = "outcome")
#' predict(fit, newdata = sp$datasets[[1]]$test_data, type = "response")
#' ```
#'
#' [fit_elastic_net()] is the penalized version of both of them, and the first
#' model here that the resampling picks rather than merely scores. `penalty` names
#' which corner of it to fit — a lasso, a ridge or the mixture — and the outcome
#' decides whether it is a regression or a two-class classification, so one
#' function covers what the other two cover between them. The penalty is what a
#' coefficient table cannot have both ways: the estimates are shrunk, some of them
#' to exactly zero, so the fit selects predictors as well as estimating them and
#' `selected` says which survived, but there is no standard error to test against,
#' and the columns that would report one are absent from its table rather than
#' present and empty. `is.null(fit$coefficients$pval)` is therefore the whole test
#' for which kind of table is in hand. Both are in the same `sa_model` shape, so
#' the choice between them is about what is being asked rather than about how the
#' answer will be read.
#'
#' [fit_rf()] grows a random forest over the same two outcome types, and is the
#' first model here that answers with no coefficients at all. A forest holds splits
#' rather than one number per predictor, so its table reports how much each
#' predictor was worth to it — the loss when its values are shuffled among the rows
#' each tree had not seen — and says which predictors carried the fit rather than
#' which way they pushed it. Two more things follow from the trees. A factor is
#' split on its levels directly, so the terms are the predictors themselves rather
#' than `k - 1` of them per factor, and every tree leaves about a third of the rows
#' out of its bootstrap sample, so `fit_stats` is an out-of-bag score rather than
#' the in-sample one the other models report.
#'
#' [fit_svm()] fits a support vector machine with a radial kernel over the same two
#' outcome types, and answers with importance for a different reason than the
#' forest does. A forest has too many numbers per predictor to report one; a kernel
#' machine has none at all, since the surface it finds lives in a space the kernel
#' never builds. What it can still be asked is what it would lose without a given
#' term, so `estimate` is the rise in error when that term's values are shuffled
#' among the rows — measured on the rows it was fitted to, because a machine has no
#' out-of-bag rows to measure it on. Both of its arguments are tuned rather than
#' one: `C` is what a margin violation costs and `sigma` how far one support vector
#' reaches, and `sigma = NULL` reads a starting width off the distances the data
#' holds.
#'
#' [predict()] on the result is what reads new rows, and it takes the data frame
#' they arrive in rather than anything prepared for the engine: the columns are
#' matched by name, a factor is coded at the levels the fit recorded, and a row
#' with a missing predictor comes back `NA` so that the answer stays aligned with
#' what was asked. That is the one call that works the same way for every model in
#' the family, since a penalized fit and a machine were given a design matrix and
#' cannot be handed the frame again.
#'
#' A `sa_model` keeps one element the rest of the package refuses to hold: `fit`,
#' the engine object itself, because a model that cannot be handed to [predict()]
#' is not much of a model. It answers to [coef()] and [summary()] as well, which
#' read through to the `lm`, `glm` or `glmnet` inside it, so `coef(fit$fit)` is
#' the named vector those models give. `coef(fit)` is the coefficient table, since
#' that is what the result holds and answering twice with the same vector would
#' waste one of the two calls. Every other element is a scalar, a character vector,
#' a named list or a data.frame, so dropping that one leaves an object that still
#' writes out as JSON.
#'
#' [simulate_regression()] and [simulate_classification()] are the known answer
#' for this family, and they carry two argument lists rather than one because two
#' functions consume them: `args` is named after the fitting function and
#' `split_args` after [split_data()]. The answer comes back on the term axis as
#' well as the predictor axis, since a categorical predictor is several
#' coefficients and one that takes a single value is none. Everything a model has
#' to survive can be asked for on purpose: predictors that correlate, missing
#' cells, and repeated measurements of one subject.
#'
#' ```
#' sim <- simulate_regression(seed = 1)
#' fit <- do.call(fit_linear_regression, c(sim$args, cv = FALSE))
#' scored <- merge(fit$coefficients, sim$truth_term, by = "terms")
#' table(planted = scored$beta != 0, called = scored$pval <= 0.05)
#' ```
#'
#' [make_block_cor()] builds the `cor_mat` they take. A null predictor correlated
#' with a planted one is the case a coefficient table gets wrong however many rows
#' it is given, and `truth$max_cor_signal` is what accounts for it afterwards.
#'
#' @section Feature selection:
#' Two functions here search rather than fit. A model is handed its predictors and
#' answers about them; a selection is handed candidates and answers which of them
#' to keep. Both return `sa_selection`, the fourth row axis in the package.
#' `candidates` takes the place `features` holds in a comparison, `terms` in a
#' model and `points` in a reduction, and it is in the order the search ranked
#' rather than the order the columns arrived. Two tables hang off it: `ranking`,
#' one row per candidate, and `profile`, one row per model the search compared.
#' `$selected` is a set of column names and nothing else, which is what lets either
#' result go straight into `predictors =` of a `fit_*()` call.
#'
#' [perform_rfe()] ranks candidates, drops the weakest, and scores what is left,
#' over and over until one predictor is standing. Every subset size is scored
#' inside the resampling rather than after it, so the ranking is recomputed in each
#' fold and no size is scored on the rows that chose it. There is no `cv` argument
#' for that reason: an elimination with nothing held out has no score to choose a
#' size by. Its `profile` is one row per subset size.
#'
#' ```
#' sel <- perform_rfe(train_data, outcome = "outcome", model = "rf")
#' fit_rf(train_data, outcome = "outcome", predictors = sel$selected)
#' ```
#'
#' That takes work inside an elimination that fits a linear or logistic model,
#' since those see a factor as `k - 1` dummy columns rather than as the column that
#' was passed in, and it is why the ranking is folded back onto the input columns
#' before anything is eliminated. It is also why the ranking is the absolute t or
#' Wald statistic rather than the coefficient: a coefficient is an effect per unit
#' of its predictor, so ranking by its size ranks by the units the predictors
#' happened to be measured in.
#'
#' [perform_stepwise()] asks the same question and pays for the answer with a
#' penalised likelihood instead of a resampled score. The model is refitted with one
#' term taken out or put back at a time, the move that lowers AIC or BIC the most is
#' taken, and the search stops when no single move lowers it further. AIC levies 2
#' per parameter and BIC levies `log(n)`, so past seven observations BIC charges
#' more and keeps fewer predictors. Nothing is held out, so nothing is resampled,
#' so nothing is random, and there is neither a `cv` argument nor a `seed`.
#' `profile` holds one row per step of the path, with both criteria at every step,
#' and `resampling` is `NULL`, which is the slot that tells the two searches apart
#' at a glance. `step()` moves whole terms, so a factor is one candidate however
#' many dummy columns it becomes and nothing has to be folded back onto the input.
#'
#' ```
#' sel <- perform_stepwise(train_data, outcome = "outcome", criterion = "BIC")
#' fit_linear_regression(train_data, outcome = "outcome",
#'                       predictors = sel$selected)
#' ```
#'
#' What an information criterion is not is a validation. It is computed on the rows
#' the model was fitted to, and the model was kept because it scored best on them,
#' so the p-values a later fit reports for the selected predictors on the same
#' rows read smaller than they are. The honest score is on the test half of
#' [split_data()], which the search never saw.
#'
#' @section Unsupervised learning:
#' [perform_pca()], [perform_tsne()] and [perform_umap()] are the first functions
#' here with no outcome at all. Nothing is predicted and nothing is scored, so what
#' they answer is where each sample sits once hundreds of features have been
#' squeezed into two dimensions. There are three of them because they disagree, and
#' the disagreement is the information: a principal component analysis is a
#' rotation, so it is readable as "which features moved this sample" but can only
#' find straight structure, while t-SNE and UMAP find structure that curves and
#' cannot say which feature made it. A cluster all three find is a different fact
#' from one that appears in a single embedding.
#'
#' All three return `sa_reduction`, a third row axis of its own: `points` takes the
#' place `features` holds in a comparison and `terms` holds in a model, and `scores`
#' repeats it in that order, which is what lets two of the three be read side by
#' side. [perform_pca()] and [perform_tsne()] standardise the features by default so
#' that the two pictures differ by method rather than by preprocessing, and the
#' `Rtsne` defaults that would have broken that are overridden on purpose and
#' reported in `engine`. [perform_umap()] leaves the values as they arrived, since
#' its `metric` is where the scale question is usually answered instead.
#'
#' ```
#' res <- perform_pca(data)
#' plot(res$scores[c("PC1", "PC2")],
#'      xlab = paste0("PC1 (", round(res$variance$prop_var[1], 2), "%)"))
#' ```
#'
#' The coordinates are `$scores` in all three, and the engine object is `$fit`. The
#' input is one row per sample as everywhere else in the package, and which margin
#' becomes the points is `embedding_scale`, with `design$point_type` reporting it.
#' Transposing the input by hand instead standardises samples rather than features,
#' which produces a picture that looks right and answers a third question; the three
#' say so in their own documentation because the mistake is easy to make and hard to
#' see.
#'
#' @keywords internal
#' @importFrom caret createDataPartition rfe rfeControl train trainControl
# `glmnet` is reached by name rather than by call: `caret::train()` is asked for
# `method = "glmnet"` and loads it itself, so nothing in this package calls the
# function that is imported here. The import is what declares the dependency and
# what keeps `coef()` and `predict()` on the fitted object dispatching to
# `glmnet`'s own methods.
#' @importFrom glmnet glmnet
# `randomForest` is imported for the same reason and in the same way: `caret` is
# asked for `method = "rf"` and loads it, and the import is what declares the
# dependency and keeps `predict()` and `importance()` on the fitted object
# dispatching to `randomForest`'s own methods.
#' @importFrom randomForest randomForest
# `Rtsne` and `umap` are the two engines this package calls itself rather than
# through `caret`, which has no method for either. Both are reached by `::` in
# `perform_tsne()` and `perform_umap()` as well; the imports are here so that the
# dependency is declared in one place with the other four.
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#' @importFrom grDevices colorRampPalette hcl.colors
#' @importFrom graphics abline axis boxplot grconvertX grconvertY grid hist
#' @importFrom graphics layout legend par
#' @importFrom graphics plot.default plot.new plot.window points rect segments
#' @importFrom graphics strwidth text
#' @importFrom stats AIC BIC as.dendrogram as.dist bartlett.test binomial coef
#' @importFrom stats complete.cases cor cov dist dnorm friedman.test hclust
#' @importFrom stats heatmap kruskal.test ks.test mad median p.adjust
#' @importFrom stats p.adjust.methods pchisq pf pnorm prop.test pt ptukey qnorm
#' @importFrom stats qt qtukey quantile rnorm runif sd setNames shapiro.test
#' @importFrom stats t.test var wilcox.test
#' @importFrom utils combn head packageVersion tail
"_PACKAGE"
