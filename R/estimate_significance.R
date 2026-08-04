#' Reduce a comparison to one significance verdict per feature
#'
#' Puts the two axes of a volcano plot side by side and flags the features that
#' clear both cutoffs. The effect size comes from the comparison's `effect`
#' table and the p-value from whichever test is named, so a single comparison can
#' be read out through the parametric, the rank-based or the robust lens without
#' recomputing anything.
#'
#' @param comparison_result A comparison result, as returned by
#'   [compare_two_groups()], [compare_multiple_groups()] or
#'   [compare_one_sample()].
#' @param test Which test in `comparison_result` supplies the p-values. One of
#'   `names(comparison_result$tests)`, so `"t_test"`, `"wilcox_test"` or
#'   `"robust_test"` for a two-group comparison. Defaults to the first test the
#'   scenario ran, which is the parametric one in every scenario.
#' @param log2fc_cutoff Minimum `abs(log2fc)` required to call a feature
#'   significant. The default of 1 is a two-fold change.
#' @param pval_cutoff Largest `adj_pvalue` allowed for a feature to be called
#'   significant.
#' @param adj_type Multiplicity adjustment. `NULL`, the default, reuses the
#'   `pval_adj` column that [compare_two_groups()] already computed. Naming a
#'   method from [stats::p.adjust.methods] re-adjusts the raw p-values instead,
#'   which is the only way to get a different adjustment without rerunning the
#'   comparison. Passing a method here does not adjust twice: the input is always
#'   the unadjusted `pval` column.
#'
#' @return A data.frame with one row per feature and the columns `features`,
#'   `log2fc`, `pvalue`, `adj_pvalue` and `is_signif`. The cutoffs, the test name
#'   and the adjustment actually used are attached as attributes, which is where
#'   [draw_volcano_plot()] picks them up so that the plotted guides cannot
#'   disagree with the verdict.
#'
#' @details
#' `is_signif` combines `abs(log2fc) >= log2fc_cutoff` with
#' `adj_pvalue <= pval_cutoff`, and is therefore judged on the adjusted
#' p-values. Pass `adj_type = "none"` to test the raw ones.
#'
#' The two ways a fold change can fall outside the domain of `log2()` are not
#' equivalent. A fold change of exactly zero gives `log2fc = -Inf`, an
#' infinitely large decrease, which clears any magnitude cutoff. A fold change
#' whose two group centres have opposite signs gives `log2fc = NaN`, which makes
#' `is_signif` `NA`: such a feature is undecided rather than decided against.
#' [compare_two_groups()] reports both kinds when it builds the `effect` table.
#' Note that `subset()` and `[` drop `NA` rows silently, so filter with
#' `which(x$is_signif)` if the count matters.
#'
#' @seealso [compare_two_groups()] for the input and [draw_volcano_plot()] for
#'   the output.
#'
#' @examples
#' iris2 <- iris[iris$Species != "setosa", ]
#' res <- compare_two_groups(
#'   data     = iris2,
#'   feats    = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
#'   group    = iris2$Species,
#'   group_lv = c("virginica", "versicolor")
#' )
#'
#' estimate_significance(res)
#'
#' ## The same comparison judged on the rank-based and the robust test
#' estimate_significance(res, test = "wilcox_test", log2fc_cutoff = 0.1)
#' estimate_significance(res, test = "robust_test", log2fc_cutoff = 0.1)
#'
#' ## A stricter fold change cutoff leaves nothing significant here, since the
#' ## two species are well under two-fold apart on every measurement.
#' estimate_significance(res, log2fc_cutoff = 1.5)
#'
#' ## Bonferroni instead of the Benjamini-Hochberg the comparison used
#' estimate_significance(res, adj_type = "bonferroni", log2fc_cutoff = 0.1)
#'
#' @export
estimate_significance <- function(comparison_result,
                                  test = names(comparison_result$tests)[1],
                                  log2fc_cutoff = 1,
                                  pval_cutoff = 0.05,
                                  adj_type = NULL) {

  sa_check_scalar_num(log2fc_cutoff, "log2fc_cutoff", 0)
  sa_check_scalar_num(pval_cutoff, "pval_cutoff", 0, 1, lower_open = TRUE)
  if (!is.null(adj_type)) {
    sa_check_p_adjust(adj_type, "adj_type")
  }

  tbl <- sa_pick_test(comparison_result, test, arg = "comparison_result")
  pvalue <- tbl$pval
  sa_check_pvalues(pvalue)

  # Re-adjusting `pval_adj` would apply the correction twice. The raw column is
  # always the input, so `adj_type` replaces the comparison's choice rather than
  # compounding it.
  if (is.null(adj_type)) {
    adj_pvalue <- tbl$pval_adj
    adj_used <- comparison_result$parameters$p_adjust
  } else {
    adj_pvalue <- stats::p.adjust(pvalue, method = adj_type)
    adj_used <- adj_type
  }

  out <- data.frame(
    features   = comparison_result$features,
    log2fc     = comparison_result$effect$log2fc,
    pvalue     = pvalue,
    adj_pvalue = adj_pvalue,
    stringsAsFactors = FALSE
  )
  out$is_signif <- abs(out$log2fc) >= log2fc_cutoff &
    out$adj_pvalue <= pval_cutoff
  rownames(out) <- NULL

  structure(out,
            analysis      = comparison_result$analysis,
            group_lv      = comparison_result$design$group_lv,
            test          = test,
            test_label    = comparison_result$test_info[[test]]$label,
            adj_type      = adj_used,
            log2fc_cutoff = log2fc_cutoff,
            pval_cutoff   = pval_cutoff)
}
