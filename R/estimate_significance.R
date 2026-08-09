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
#' @param by Which p-value the verdict is read from. `"omnibus"`, the default,
#'   uses the named test's own table, which is the only thing a two-group or a
#'   one-sample comparison has. `"contrast"` uses the pairwise stage of a
#'   multi-group comparison instead and returns one verdict table per contrast.
#'
#' @return An `sa_significance` object, a list of two elements:
#'   \describe{
#'     \item{`analysis_type`}{The `analysis` of the comparison this verdict was
#'       read from, so a table that has been passed around still says which
#'       scenario produced it.}
#'     \item{`significance`}{With `by = "omnibus"`, a data.frame with one row per
#'       feature and the columns `features`, `log2fc`, `pvalue`, `adj_pvalue` and
#'       `is_signif`. With `by = "contrast"`, a list of those same data.frames,
#'       one per pairwise contrast and named after it, in the order
#'       `comparison_result$pairwise[[test]]` fixes.}
#'   }
#'
#'   The cutoffs, the test name and the adjustment actually used are attached to
#'   each data.frame as attributes, which is where [draw_volcano_plot()] picks
#'   them up so that the plotted guides cannot disagree with the verdict. A
#'   contrast table carries `contrast`, `group1` and `group2` on top of those, so
#'   a single element of `significance` can be handed straight to
#'   [draw_volcano_plot()].
#'
#' @details
#' `is_signif` combines `abs(log2fc) >= log2fc_cutoff` with
#' `adj_pvalue <= pval_cutoff`, and is therefore judged on the adjusted
#' p-values. Pass `adj_type = "none"` to test the raw ones.
#'
#' The two ways of reading a multi-group comparison answer different questions
#' and use different numbers. `by = "omnibus"` asks whether a feature differs
#' across the levels at all, and pairs that with the one `log2fc` the `effect`
#' table carries, the most extreme level against the reference.
#' `by = "contrast"` asks the same question of one pair of levels at a time, and
#' each table carries the `log2fc` of that pair, in the direction its `contrast`
#' label reads. Both divide by the reference, so the two readings agree on which
#' way a feature moved.
#'
#' The adjustment axis differs too. Under `by = "contrast"` with
#' `adj_type = NULL` the `pval_adj` of the pairwise stage is reused, which was
#' adjusted across the contrasts within each feature by `posthoc_p_adjust`, or
#' not at all for Tukey's HSD and Games-Howell, whose p-values are already
#' family-wise. Naming a method instead adjusts across the features within each
#' contrast, which is the axis `by = "omnibus"` always works on.
#'
#' A feature whose omnibus test did not clear `posthoc_alpha` was never compared
#' pairwise, so its `pvalue` is `NA` in every contrast table. Its `log2fc` is
#' still reported, since a ratio of group centres does not depend on a test
#' having been run, and `is_signif` follows the same three-valued rule it does
#' everywhere else: `NA`, undecided, unless the magnitude cutoff already rules
#' the feature out, which makes it `FALSE` whatever the p-value would have been.
#'
#' The two ways a fold change can fall outside the domain of `log2()` are not
#' equivalent. A fold change of exactly zero gives `log2fc = -Inf`, an
#' infinitely large decrease, which clears any magnitude cutoff. A fold change
#' whose two group centres have opposite signs gives `log2fc = NaN`, which makes
#' `is_signif` `NA`: such a feature is undecided rather than decided against.
#' [compare_two_groups()] reports both kinds when it builds the `effect` table.
#' Note that `subset()` and `[` drop `NA` rows silently, so filter with
#' `which(x$significance$is_signif)` if the count matters.
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
#' ## One verdict table per pairwise contrast of a multi-group comparison
#' multi <- compare_multiple_groups(iris, c("Sepal.Length", "Petal.Length"),
#'                                  iris$Species, levels(iris$Species))
#' by_pair <- estimate_significance(multi, by = "contrast",
#'                                  log2fc_cutoff = 0.1)
#' names(by_pair$significance)
#' by_pair$significance[["virginica - setosa"]]
#'
#' @export
estimate_significance <- function(comparison_result,
                                  test = names(comparison_result$tests)[1],
                                  log2fc_cutoff = 1,
                                  pval_cutoff = 0.05,
                                  adj_type = NULL,
                                  by = c("omnibus", "contrast")) {

  by <- match.arg(by)
  sa_check_scalar_num(log2fc_cutoff, "log2fc_cutoff", 0)
  sa_check_scalar_num(pval_cutoff, "pval_cutoff", 0, 1, lower_open = TRUE)
  if (!is.null(adj_type)) {
    sa_check_p_adjust(adj_type, "adj_type")
  }

  tbl <- sa_pick_test(comparison_result, test, arg = "comparison_result")

  if (by == "contrast") {
    return(sa_new_significance(
      comparison_result$analysis,
      sa_significance_by_contrast(comparison_result, test, adj_type,
                                  log2fc_cutoff, pval_cutoff)
    ))
  }

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

  out <- sa_significance_table(comparison_result$features,
                               comparison_result$effect$log2fc,
                               pvalue, adj_pvalue,
                               log2fc_cutoff, pval_cutoff)

  out <- do.call(structure,
                 c(list(out),
                   sa_significance_attrs(comparison_result, test, adj_used,
                                         log2fc_cutoff, pval_cutoff)))

  sa_new_significance(comparison_result$analysis, out)
}


#' Wrap a verdict table in the object the user sees
#'
#' The scenario name sits beside the table rather than only in its attributes,
#' since a consumer reading the verdict has to know which question the `log2fc`
#' column answers, and an attribute is easy to lose and easy to miss.
#'
#' @keywords internal
#' @noRd
sa_new_significance <- function(analysis_type, significance) {
  structure(list(analysis_type = analysis_type,
                 significance  = significance),
            class = c("sa_significance", "sa_result"))
}


#' The attributes a verdict table describes itself with
#'
#' [draw_volcano_plot()] reads the cutoffs back off here, so a table that has
#' been passed around still knows which comparison and which thresholds produced
#' it.
#'
#' @keywords internal
#' @noRd
sa_significance_attrs <- function(comparison_result, test, adj_used,
                                  log2fc_cutoff, pval_cutoff) {
  list(analysis      = comparison_result$analysis,
       group_lv      = comparison_result$design$group_lv,
       test          = test,
       test_label    = comparison_result$test_info[[test]]$label,
       adj_type      = adj_used,
       log2fc_cutoff = log2fc_cutoff,
       pval_cutoff   = pval_cutoff)
}


#' The verdict table both reading of a comparison produce
#'
#' Shared so that the two axes and the rule combining them cannot drift apart
#' between the omnibus and the contrast paths.
#'
#' @keywords internal
#' @noRd
sa_significance_table <- function(features, log2fc, pvalue, adj_pvalue,
                                  log2fc_cutoff, pval_cutoff) {
  out <- data.frame(
    features   = features,
    log2fc     = log2fc,
    pvalue     = pvalue,
    adj_pvalue = adj_pvalue,
    stringsAsFactors = FALSE
  )
  out$is_signif <- abs(out$log2fc) >= log2fc_cutoff &
    out$adj_pvalue <= pval_cutoff
  rownames(out) <- NULL
  out
}


#' One verdict table per pairwise contrast
#'
#' @keywords internal
#' @noRd
sa_significance_by_contrast <- function(comparison_result, test, adj_type,
                                        log2fc_cutoff, pval_cutoff) {
  pairwise <- comparison_result$pairwise[[test]]
  if (is.null(pairwise) || length(pairwise) == 0L) {
    stop("`by = \"contrast\"` needs a pairwise stage, and ",
         "`comparison_result$pairwise$", test, "` is absent. Only ",
         "compare_multiple_groups() runs one, and only when ",
         "`posthoc = TRUE`.", call. = FALSE)
  }

  lapply(pairwise, function(tbl) {
    pvalue <- tbl$pval
    sa_check_pvalues(pvalue)
    if (is.null(adj_type)) {
      adj_pvalue <- tbl$pval_adj
      adj_used <- comparison_result$parameters$posthoc_p_adjust
    } else {
      # Across the features of this one contrast, which is a different family
      # from the one the pairwise stage corrected over.
      adj_pvalue <- stats::p.adjust(pvalue, method = adj_type)
      adj_used <- adj_type
    }

    out <- sa_significance_table(tbl$features, tbl$log2fc, pvalue, adj_pvalue,
                                 log2fc_cutoff, pval_cutoff)

    do.call(structure,
            c(list(out),
              sa_significance_attrs(comparison_result, test, adj_used,
                                    log2fc_cutoff, pval_cutoff),
              list(contrast = tbl$contrast[1],
                   group1   = tbl$group1[1],
                   group2   = tbl$group2[1])))
  })
}


#' Print a significance verdict
#'
#' Summarises the rule that was applied and how many features cleared it, rather
#' than printing the table itself. Reach into `x$significance` for that.
#'
#' @param x An `sa_significance` object, as returned by
#'   [estimate_significance()].
#' @param ... Ignored, present for consistency with [print()].
#'
#' @return `x` invisibly.
#'
#' @examples
#' iris2 <- iris[iris$Species != "setosa", ]
#' res <- compare_two_groups(iris2, c("Sepal.Length", "Petal.Length"),
#'                           iris2$Species, c("virginica", "versicolor"))
#' estimate_significance(res, log2fc_cutoff = 0.1)
#'
#' @export
print.sa_significance <- function(x, ...) {
  # A contrast reading holds one table per contrast, and every one of them was
  # judged by the same rule, so the header is read off whichever comes first.
  tables <- if (is.data.frame(x$significance)) {
    list(x$significance)
  } else {
    x$significance
  }
  head_tbl <- tables[[1]]

  cat("<", class(x)[1], "> ", x$analysis_type, "\n", sep = "")
  cat("  test     : ", attr(head_tbl, "test"), "  (",
      attr(head_tbl, "test_label"), ")\n", sep = "")
  # Tukey's HSD and Games-Howell report family-wise p-values, so the pairwise
  # stage adjusted nothing and the attribute is absent rather than "none".
  adj_used <- attr(head_tbl, "adj_type")
  if (is.null(adj_used)) adj_used <- "none"
  cat("  cutoffs  : abs(log2fc) >= ", attr(head_tbl, "log2fc_cutoff"),
      ", adj_pvalue <= ", attr(head_tbl, "pval_cutoff"),
      "  (", adj_used, ")\n", sep = "")

  count <- function(tbl) {
    n_signif <- sum(tbl$is_signif %in% TRUE)
    n_undecided <- sum(is.na(tbl$is_signif))
    paste0(n_signif, " of ", nrow(tbl), " significant",
           if (n_undecided > 0L) paste0("  (", n_undecided, " undecided)"))
  }

  if (is.data.frame(x$significance)) {
    cat("  verdict  : ", count(x$significance), "\n", sep = "")
  } else {
    cat("\n  $significance, one table per contrast\n")
    width <- max(nchar(names(tables)))
    for (nm in names(tables)) {
      cat("    ", formatC(nm, width = -width), "  ", count(tables[[nm]]), "\n",
          sep = "")
    }
  }

  invisible(x)
}
