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
#' [compare_two_groups()] compares exactly two group levels across any number
#' of numeric features. It returns the fold change between the groups plus three
#' test tables:
#'
#' \describe{
#'   \item{`t_test`}{Welch's t-test (independent) or paired t-test.}
#'   \item{`wilcox_test`}{Wilcoxon rank-sum test (independent) or signed-rank
#'     test (paired), with the Hodges-Lehmann shift estimate.}
#'   \item{`robust_test`}{Brunner-Munzel test (independent) or Yuen's trimmed
#'     mean test for dependent samples (paired).}
#' }
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
#' [compare_two_groups()] returns a `sa_comparison` object, a plain named list of
#' scalars, character vectors and data.frames with an S3 class on top that only
#' supplies [print()]. No fitted model or other R-only object is stored anywhere,
#' so the result can be written out as JSON and rebuilt elsewhere. Every test
#' kernel is likewise a plain function of numeric vectors. Both choices exist to
#' keep a future Python implementation a transcription rather than a redesign.
#'
#' @keywords internal
#' @importFrom grDevices hcl.colors
#' @importFrom graphics abline axis boxplot grid hist layout legend par
#' @importFrom graphics plot.default plot.new points rect text
#' @importFrom stats dnorm mad p.adjust p.adjust.methods pt qnorm qt quantile
#' @importFrom stats sd setNames t.test var wilcox.test
#' @importFrom utils packageVersion
"_PACKAGE"
