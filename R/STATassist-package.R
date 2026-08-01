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
#' direction.
#'
#' ```
#' res <- compare_two_groups(data, feats, group, group_lv)
#' sig <- estimate_significance(res, test = "t_test")
#' draw_volcano_plot(sig)
#' ```
#'
#' [draw_grouped_boxplot()] draws a grouped boxplot for the same wide-format
#' input and optionally returns the box summary statistics and median
#' confidence intervals behind the plot. [draw_butterfly_hist()] puts the two
#' group distributions of a single feature back to back on shared breaks.
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
#' The contract is at version `0.2.0`. Alongside the `effect` and `tests` tables
#' of `0.1.0` it carries `posthoc`, one table per test holding one row per
#' feature and pair of levels, and `diagnostics`. Both are present in every
#' result, empty where the scenario has nothing to put in them, so a consumer
#' reads all three scenarios the same way.
#'
#' @keywords internal
#' @importFrom grDevices hcl.colors
#' @importFrom graphics abline axis boxplot grid hist layout legend par
#' @importFrom graphics plot.default plot.new points rect segments text
#' @importFrom stats bartlett.test complete.cases cov dnorm friedman.test
#' @importFrom stats kruskal.test ks.test mad median p.adjust p.adjust.methods
#' @importFrom stats pchisq pf pnorm prop.test pt ptukey qnorm qt qtukey
#' @importFrom stats quantile sd setNames shapiro.test t.test var wilcox.test
#' @importFrom utils combn packageVersion
"_PACKAGE"
