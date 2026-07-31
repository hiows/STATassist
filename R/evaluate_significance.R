#' Combine fold changes and p-values into a significance table
#'
#' Puts the two axes of a volcano plot side by side, adjusts the p-values for
#' multiplicity and flags the features that clear both cutoffs.
#'
#' @param feats Character vector of feature names. One output row per entry, in
#'   this order.
#' @param fold_change Numeric vector of fold changes, one per feature, such as
#'   the output of [calculate_fold_change()]. Reported as `log2fc`.
#' @param pvalue Numeric vector of unadjusted p-values, one per feature.
#' @param adj_type Multiplicity adjustment passed to [stats::p.adjust()], one of
#'   `stats::p.adjust.methods`. Use `"none"` to leave the p-values alone.
#' @param log2fc_cutoff Minimum `abs(log2fc)` required to call a feature
#'   significant. The default of 1 is a two-fold change.
#' @param pval_cutoff Largest `adj_pvalue` allowed for a feature to be called
#'   significant.
#'
#' @return A data.frame with one row per feature and the columns `features`,
#'   `log2fc`, `pvalue`, `adj_pvalue` and `is_signif`.
#'
#' @details
#' `fold_change` and `pvalue` must have one value per feature. When either
#' carries names, as [calculate_fold_change()] output does, the names must match
#' `feats` and are used to reorder the values; otherwise position is trusted.
#'
#' `is_signif` combines `abs(log2fc) >= log2fc_cutoff` with
#' `adj_pvalue <= pval_cutoff`, and is therefore judged on the adjusted
#' p-values. Pass `adj_type = "none"` to test the raw ones.
#'
#' The two ways a fold change can fall outside the domain of `log2()` are not
#' equivalent. A fold change of exactly zero gives `log2fc = -Inf`, an
#' infinitely large decrease, which clears any magnitude cutoff. A negative fold
#' change gives `log2fc = NaN`, which makes `is_signif` `NA`: such a feature is
#' undecided rather than decided against. Both are listed in a `message()`.
#' Note that `subset()` and `[` drop `NA` rows silently, so filter with
#' `which(x$is_signif)` if the count matters.
#'
#' @seealso [calculate_fold_change()] for the `fold_change` argument,
#'   [compare_two_groups()] for `pvalue`, and [draw_volcano_plot()] to plot the
#'   result.
#'
#' @examples
#' iris2 <- iris[iris$Species != "setosa", ]
#' feats <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
#' lv <- c("versicolor", "virginica")
#'
#' fc <- calculate_fold_change(iris2, feats, iris2$Species, lv,
#'                             case_label = "virginica")
#' pv <- compare_two_groups(iris2, feats, iris2$Species,
#'                          lv)$test_results$t_test$pval
#'
#' evaluate_significance(feats, fc, pv)
#'
#' ## A stricter fold change cutoff leaves nothing significant here, since the
#' ## two species are well under two-fold apart on every measurement.
#' evaluate_significance(feats, fc, pv, log2fc_cutoff = 1.5)
#'
#' ## Bonferroni instead of the default Benjamini-Hochberg
#' evaluate_significance(feats, fc, pv, adj_type = "bonferroni",
#'                       log2fc_cutoff = 0.1)
#'
#' @export
evaluate_significance <- function(feats,
                                  fold_change,
                                  pvalue,
                                  adj_type = "BH",
                                  log2fc_cutoff = 1,
                                  pval_cutoff = 0.05) {

  # The choices used to be spelled out as c("BH", "none", "bonferonni"), where
  # match.arg() accepted the misspelled third one and p.adjust() then rejected
  # it. Validating against p.adjust.methods removes the trap and opens up the
  # remaining methods.
  if (!is.character(adj_type) || length(adj_type) != 1L ||
      !adj_type %in% stats::p.adjust.methods) {
    stop("`adj_type` must be one of: ",
         paste(stats::p.adjust.methods, collapse = ", "), ".", call. = FALSE)
  }
  sa_check_scalar_num(log2fc_cutoff, "log2fc_cutoff", 0)
  sa_check_scalar_num(pval_cutoff, "pval_cutoff", 0, 1, lower_open = TRUE)

  sa_check_feat_names(feats)
  fold_change <- sa_align_to_feats(fold_change, feats, "fold_change")
  pvalue <- sa_align_to_feats(pvalue, feats, "pvalue")
  sa_check_pvalues(pvalue)

  # log2() of a negative number warns once per call; the features responsible
  # are more useful than the warning, so they are collected and reported here.
  log2fc <- suppressWarnings(log2(fold_change))
  zero_fc <- feats[!is.na(fold_change) & fold_change == 0]
  if (length(zero_fc) > 0L) {
    message("Fold change is zero for ", length(zero_fc),
            " feature(s), so `log2fc` is -Inf and clears any ",
            "`log2fc_cutoff`: ", paste(zero_fc, collapse = ", "), ".")
  }
  neg_fc <- feats[!is.na(fold_change) & fold_change < 0]
  if (length(neg_fc) > 0L) {
    message("Fold change is negative for ", length(neg_fc),
            " feature(s), so `log2fc` is NaN and `is_signif` is NA: ",
            paste(neg_fc, collapse = ", "), ".")
  }

  out <- data.frame(
    features   = feats,
    log2fc     = log2fc,
    pvalue     = pvalue,
    adj_pvalue = stats::p.adjust(pvalue, method = adj_type),
    stringsAsFactors = FALSE
  )
  out$is_signif <- abs(out$log2fc) >= log2fc_cutoff &
    out$adj_pvalue <= pval_cutoff

  rownames(out) <- NULL
  out
}
