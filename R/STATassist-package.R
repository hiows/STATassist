#' STATassist: Standardised Statistical Comparison Workflows
#'
#' STATassist runs several complementary statistical tests in one call and
#' returns feature-wise result tables with a shared column layout. The guiding
#' idea is that a group comparison should never rest on a single test: a
#' parametric, a rank-based and a robust procedure are reported side by side so
#' that disagreement between them becomes visible instead of hidden.
#'
#' @section Two-group comparison:
#' [compare_two_groups()] compares exactly two group levels across any number
#' of numeric features and returns three tables:
#'
#' \describe{
#'   \item{`t_test`}{Welch's t-test (independent) or paired t-test.}
#'   \item{`wilcox_test`}{Wilcoxon rank-sum test (independent) or signed-rank
#'     test (paired), with the Hodges-Lehmann shift estimate.}
#'   \item{`robust_test`}{Brunner-Munzel test (independent) or Yuen's trimmed
#'     mean test for dependent samples (paired).}
#' }
#'
#' @section Descriptive summary:
#' [summarize_descriptive_stats()] reduces each feature to one row of sample
#' size, central tendency, dispersion, quartiles, outlier fences and
#' distribution shape, split by group level when a grouping vector is supplied.
#'
#' @section Effect size and significance:
#' [calculate_fold_change()] reduces each feature to a ratio of group means, and
#' [evaluate_significance()] pairs those ratios with p-values from
#' [compare_two_groups()] to produce one table holding both axes of a volcano
#' plot together with a multiplicity adjusted verdict:
#'
#' ```
#' fc  <- calculate_fold_change(data, feats, group, group_lv)
#' pv  <- compare_two_groups(data, feats, group, group_lv)$test_results$t_test$pval
#' sig <- evaluate_significance(feats, fc, pv)
#' draw_volcano_plot(sig$features, sig$log2fc, sig$adj_pvalue)
#' ```
#'
#' These three functions are joined by feature order alone, so each one checks
#' that its per-feature vectors line up, and reorders them by name when the names
#' are available.
#'
#' @section Visualisation:
#' [draw_grouped_boxplot()] draws a grouped boxplot for the same wide-format
#' input and optionally returns the box summary statistics and median
#' confidence intervals behind the plot. [draw_butterfly_hist()] puts the two
#' group distributions of a single feature back to back on shared breaks.
#' [draw_volcano_plot()] plots the output of [evaluate_significance()].
#'
#' @keywords internal
#' @importFrom colorspace qualitative_hcl
#' @importFrom graphics abline axis boxplot grid hist layout legend par
#' @importFrom graphics plot.default plot.new points rect text
#' @importFrom lawstat brunner.munzel.test
#' @importFrom stats mad p.adjust p.adjust.methods pt qt quantile sd setNames
#' @importFrom stats t.test var wilcox.test
#' @importFrom tidyr all_of pivot_longer
#' @importFrom WRS2 yuend
"_PACKAGE"
