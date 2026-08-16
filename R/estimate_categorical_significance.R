#' Reduce a contingency table to significance verdicts
#'
#' The categorical counterpart of [estimate_significance()], and the reason it is
#' a second function rather than a branch of that one is that the two axes a
#' verdict is made of sit at a different granularity here. A comparison asks its
#' question once per feature; a contingency table is asked about as a whole, and
#' the place where a signed effect and a p-value exist side by side is one level
#' down, at the cell.
#'
#' @section Why a cell has both axes and a table has one:
#' [compare_categorical_groups()] says that an association is not signed, which is
#' true of the **table**: `cramers_v` reports how far the table sits from its null
#' and not in which direction, because past a 2 x 2 there is no single direction
#' to name. A **cell** is different. It was expected at some count and observed at
#' another, and `observed / expected` says both how far it moved and which way.
#'
#' That ratio is `lift`, the same quantity [simulate_categorical_groups()] plants
#' and reports in `truth_cell$lift`, and `log2(lift)` is what the effect axis of
#' this function is. It is defined on every table, 2 x 2, 2 x 3 or larger, so one
#' cutoff means the same thing whatever shape the table is.
#'
#' The p-value axis is `std_residual`, which is built to be referred to a standard
#' normal, so `2 * pnorm(-abs(std_residual))` is the cell's own two-sided
#' p-value. This is the standard post-hoc reading of a contingency table.
#'
#' The two axes are not the same information. `lift` is a ratio and does not
#' change when the table is observed on twice as many rows; `std_residual` grows
#' with the square root of the count. So a cell can be far from what was expected
#' and poorly evidenced, or close to it and firmly established, exactly as
#' `log2fc` and a p-value come apart in the numeric scenarios.
#'
#' A cell axis also restores a multiplicity axis. The table has as many cells as
#' it has, they are one family, and `adj_type` adjusts across them. That is why
#' this function computes an adjustment rather than reusing one: `sa_categorical`
#' carries no `pval_adj` column, there being nothing to adjust across at the level
#' the tests are reported at.
#'
#' @section The two readings:
#' \describe{
#'   \item{`by = "cell"`}{One row per cell of the table, with `log2_lift` on the
#'     effect axis and the adjusted p-value of the cell's standardized residual on
#'     the other. This is the reading that generalises, so it is the default.}
#'   \item{`by = "table"`}{One row, the whole-table verdict: an association
#'     measure out of `$association` beside the p-value of one of the tests. There
#'     is no multiplicity here, one table being one question, so no `adj_pvalue`
#'     column is carried.}
#' }
#'
#' Each reading ignores the arguments the other one reads, and says so rather than
#' letting a setting that changes nothing pass unremarked.
#'
#' @section Which measure a table reading reports:
#' `measure = "auto"` reads the design, since which measures exist at all depends
#' on it and on the size of the table.
#'
#' | design | table | measure |
#' |---|---|---|
#' | independent | 2 x 2 | `odds_ratio` |
#' | independent | larger | `cramers_v` |
#' | matched, two conditions | 2 x 2 | `odds_ratio_paired` |
#' | matched, three or more | condition by response | `kendalls_w` |
#'
#' Any other row of `$association` can be named instead.
#'
#' `effect_cutoff` is `NULL` by default, and then the verdict is the p-value
#' alone. The conventional thresholds for `cramers_v` are not defaults here
#' because they are conventions rather than facts about the measure, and a
#' default is the one place a convention is hardest to notice. Naming a number
#' reads it on the measure's own scale: `odds_ratio` and `odds_ratio_paired` are
#' ratios centred at 1, so the cutoff is a fold either way and has to be at least
#' 1, while every other measure is centred at zero and the cutoff is compared
#' against the magnitude.
#'
#' @param categorical_comparison_result A categorical comparison result, as
#'   returned by [compare_categorical_groups()]. A numeric comparison is refused
#'   and pointed at [estimate_significance()].
#' @param by `"cell"` for one verdict per cell of the table, `"table"` for one
#'   verdict for the table as a whole. See "The two readings".
#' @param test Which test supplies the p-value of a table reading. One of
#'   `names(categorical_comparison_result$tests)`, so `"chisq_test"` or
#'   `"fisher_test"` for an independent design. Defaults to the first test the
#'   design ran. Read only by `by = "table"`: a cell reading takes its p-value
#'   from the cell's own standardized residual, which no test reports.
#' @param log2_lift_cutoff Minimum `abs(log2_lift)` required to call a cell
#'   significant. Read only by `by = "cell"`. The default of 1 is a cell holding
#'   twice, or half, what its null expected, and is deliberately the same number
#'   [estimate_significance()] uses. It is a strict demand of a contingency
#'   table, where departures well under a doubling are ordinary; halve it to ask
#'   of a cell what the default asks of a fold change.
#' @param pval_cutoff Largest p-value allowed for a verdict of significant. Read
#'   by both readings, against `adj_pvalue` under `by = "cell"` and against the
#'   test's own `pvalue` under `by = "table"`.
#' @param adj_type Multiplicity adjustment across the cells, one of
#'   [stats::p.adjust.methods]. Read only by `by = "cell"`. Unlike the `adj_type`
#'   of [estimate_significance()] this is the first adjustment rather than a
#'   replacement for one, since a categorical result carries no adjusted column.
#'   `"none"` tests the raw p-values.
#' @param measure Which row of `$association` a table reading puts on its effect
#'   axis, or `"auto"` to take the one the design defines. Read only by
#'   `by = "table"`.
#' @param effect_cutoff Magnitude `measure` has to reach for a table reading to
#'   call the table significant, or `NULL` to judge on the p-value alone. Read
#'   only by `by = "table"`.
#'
#' @return An `sa_categorical_significance` object, a list of two elements:
#'   \describe{
#'     \item{`analysis_type`}{`"categorical_comparison"`.}
#'     \item{`significance`}{With `by = "cell"`, a data.frame with one row per
#'       cell and the columns `row_level`, `col_level`, `observed`, `expected`,
#'       `lift`, `log2_lift`, `std_residual`, `pvalue`, `adj_pvalue` and
#'       `is_signif`. With `by = "table"`, a one-row data.frame with `measure`,
#'       `estimate`, `lower_conf`, `upper_conf`, `pvalue` and `is_signif`.}
#'   }
#'
#'   The cutoffs, the reading, the null hypothesis and whichever of the test name
#'   and the adjustment applies are attached to the data.frame as attributes, so
#'   a table that has been passed around still says which rule produced it.
#'
#'   A cell table keys on `c("row_level", "col_level")`, the key
#'   `categorical_comparison_result$cells` and
#'   `simulate_categorical_groups()$truth_cell` also use, so a verdict merges with
#'   either without renaming. Note that `truth_cell` carries a `lift` column of
#'   its own, the planted one, so a merge of the two distinguishes them as
#'   `lift.x` and `lift.y`.
#'
#'   This is deliberately not an `sa_significance`. Its columns are cells rather
#'   than features and `log2_lift` rather than `log2fc`, so [draw_volcano_plot()]
#'   refusing it is the point of the separate class rather than an omission.
#'   [draw_mosaic_plot()] is what draws this scenario.
#'
#' @details
#' `is_signif` is three-valued and follows the rule the numeric scenarios use.
#' `NA` is undecided rather than decided against, which is what an `NA`
#' `std_residual` or an undefined `lift` leaves a cell as.
#'
#' A cell holding no observation at all has a `lift` of exactly zero and so a
#' `log2_lift` of `-Inf`, an infinitely large shortfall, which clears any
#' magnitude cutoff. A cell whose expected count is zero has no ratio to take and
#' is `NA` in both columns. The two are different findings and are reported
#' differently, the same way [estimate_significance()] separates a fold change of
#' zero from one that cannot be formed.
#'
#' @section Why a matched pair of conditions has no cell reading:
#' A matched two-condition design is tested for symmetry, and the variance
#' correction the standardized residual divides by is derived for a table held
#' against its own margins. `$cells$std_residual` is therefore `NA` throughout
#' such a result, so there is no p-value axis to read and the request is refused
#' rather than answered with a different quantity. `by = "table"` reads it, and
#' [draw_mosaic_plot()] refuses `residual = "standardized"` there for the same
#' reason.
#'
#' Three or more matched conditions are a different case. Their null is marginal
#' homogeneity, read off a condition-by-response table whose arithmetic is that of
#' independence, so the standardized residual exists and the cell reading works.
#'
#' @seealso [compare_categorical_groups()] for the input,
#'   [draw_mosaic_plot()] to draw the same residuals this reads, and
#'   [estimate_significance()] for the numeric scenarios.
#'
#' @references
#' Haberman, S. J. (1973). The analysis of residuals in cross-classified tables.
#' *Biometrics*, 29(1), 205-220.
#'
#' Beasley, T. M. and Schumacker, R. E. (1995). Multiple regression approach to
#' analyzing contingency tables: Post hoc and planned comparison procedures.
#' *The Journal of Experimental Education*, 64(1), 79-93.
#'
#' Agresti, A. (2002). *Categorical Data Analysis*, 2nd ed. Wiley.
#'
#' @examples
#' smoking <- data.frame(
#'   smoker = rep(c("y", "n"), each = 60),
#'   grade  = c(rep(c("high", "mid", "low"), c(10, 20, 30)),
#'              rep(c("high", "mid", "low"), c(30, 20, 10)))
#' )
#' res <- compare_categorical_groups(smoking)
#'
#' ## One verdict per cell. The default cutoff is a doubling either way, which on
#' ## this table only the two corners reach.
#' sig <- estimate_categorical_significance(res)
#' sig
#' sig$significance
#'
#' ## `lift` is what the effect axis is built from, and it is a ratio: the cell
#' ## holding half of what independence expected is at 0.5.
#' sig$significance[c("row_level", "col_level", "lift", "log2_lift")]
#'
#' ## A gentler magnitude cutoff asks less of a cell than a doubling.
#' estimate_categorical_significance(res, log2_lift_cutoff = 0.5)
#'
#' ## The whole-table verdict instead. A 2 x 3 table has no odds ratio, so the
#' ## measure `"auto"` reports is Cramer's V.
#' estimate_categorical_significance(res, by = "table")$significance
#'
#' ## Two levels of `grade` make a 2 x 2 table, which is where an odds ratio
#' ## exists, and `effect_cutoff` is then read as a fold either way.
#' two_by_two <- compare_categorical_groups(
#'   smoking,
#'   category_lv = list(smoker = c("n", "y"), grade = c("low", "high"))
#' )
#' estimate_categorical_significance(two_by_two, by = "table",
#'                                   effect_cutoff = 2)$significance
#'
#' ## The cell axis against the association that was planted on it.
#' sim <- simulate_categorical_groups(n_samples = 400, assoc = 0.4, seed = 1)
#' fit <- do.call(compare_categorical_groups, sim$args)
#' scored <- merge(estimate_categorical_significance(fit)$significance,
#'                 sim$truth_cell, by = c("row_level", "col_level"))
#' scored[c("row_level", "col_level", "lift.x", "lift.y")]
#'
#' @export
estimate_categorical_significance <- function(
    categorical_comparison_result,
    by = c("cell", "table"),
    test = names(categorical_comparison_result$tests)[1],
    log2_lift_cutoff = 1,
    pval_cutoff = 0.05,
    adj_type = "BH",
    measure = "auto",
    effect_cutoff = NULL) {

  by <- match.arg(by)

  # Before anything reaches for a slot, so that the wrong object is told what it
  # is rather than failing somewhere inside on a slot it does not have. The
  # default of `test` reads one, and being a promise it is not forced until then.
  res <- categorical_comparison_result
  if (!inherits(res, "sa_categorical")) {
    if (inherits(res, "sa_comparison")) {
      stop("`categorical_comparison_result` is a numeric comparison result. ",
           "estimate_significance() is what reads one; this function reads a ",
           "contingency table.", call. = FALSE)
    }
    stop("`categorical_comparison_result` must be a categorical comparison ",
         "result, as returned by compare_categorical_groups().", call. = FALSE)
  }

  sa_check_scalar_num(log2_lift_cutoff, "log2_lift_cutoff", 0)
  sa_check_scalar_num(pval_cutoff, "pval_cutoff", 0, 1, lower_open = TRUE)
  sa_check_p_adjust(adj_type, "adj_type")

  sa_warn_unread_args(by, missing(test), missing(measure),
                      missing(effect_cutoff), missing(log2_lift_cutoff),
                      missing(adj_type))

  if (identical(by, "cell")) {
    if (identical(res$design$null, "symmetry")) {
      stop("`by = \"cell\"` needs a p-value per cell, and this result was ",
           "tested for symmetry, where `$cells$std_residual` is NA throughout: ",
           "the variance correction it divides by is derived for a table held ",
           "against its own margins and has no counterpart there. Use ",
           "`by = \"table\"`, whose verdict reads McNemar's p-value and the ",
           "paired odds ratio.", call. = FALSE)
    }
    out <- sa_categorical_significance_cells(res, adj_type, log2_lift_cutoff,
                                             pval_cutoff)
    out <- do.call(structure,
                   c(list(out),
                     sa_categorical_significance_attrs(res, by, pval_cutoff),
                     list(log2_lift_cutoff = log2_lift_cutoff,
                          adj_type         = adj_type)))
    return(sa_new_categorical_significance(res$analysis, out))
  }

  tbl <- sa_pick_categorical_test(res, test)
  measure <- sa_resolve_measure(res, measure)
  sa_check_effect_cutoff(effect_cutoff, measure)

  out <- sa_categorical_significance_row(res, tbl, measure, effect_cutoff,
                                         pval_cutoff)
  out <- do.call(structure,
                 c(list(out),
                   sa_categorical_significance_attrs(res, by, pval_cutoff),
                   list(test          = test,
                        test_label    = res$test_info[[test]]$label,
                        measure       = measure,
                        effect_cutoff = effect_cutoff)))
  sa_new_categorical_significance(res$analysis, out)
}


#' One verdict per cell of the table
#'
#' `lift` and the standardized residual are two readings of the same departure,
#' and they are kept apart because they answer different questions: the first is a
#' ratio and does not move with the sample size, the second is the departure
#' divided by what its own sampling variability would have been.
#'
#' @param res The categorical comparison result.
#' @param adj_type,log2_lift_cutoff,pval_cutoff The arguments as received.
#'
#' @keywords internal
#' @noRd
sa_categorical_significance_cells <- function(res, adj_type, log2_lift_cutoff,
                                              pval_cutoff) {
  cells <- res$cells

  # An expected count of zero leaves no ratio to take, which is not a lift of
  # zero: `sa_finite_or_na()` is what keeps the two apart. A lift of exactly zero
  # survives it, and its log2 is -Inf, an infinitely large shortfall.
  lift <- sa_finite_or_na(cells$observed / cells$expected)
  log2_lift <- log2(lift)

  # The standardized residual is built to be referred to a standard normal, so
  # this is the cell's own two-sided p-value and not an approximation of one.
  pvalue <- 2 * stats::pnorm(-abs(cells$std_residual))
  sa_check_pvalues(pvalue)

  # The cells of one table are one family. This is the first adjustment rather
  # than a replacement for one, since `sa_categorical` carries no adjusted
  # column: at the level its tests are reported at there is nothing to adjust
  # across.
  adj_pvalue <- stats::p.adjust(pvalue, method = adj_type)

  out <- data.frame(
    row_level    = cells$row_level,
    col_level    = cells$col_level,
    observed     = cells$observed,
    expected     = cells$expected,
    lift         = lift,
    log2_lift    = log2_lift,
    std_residual = cells$std_residual,
    pvalue       = pvalue,
    adj_pvalue   = adj_pvalue,
    stringsAsFactors = FALSE
  )
  out$is_signif <- abs(out$log2_lift) >= log2_lift_cutoff &
    out$adj_pvalue <= pval_cutoff
  rownames(out) <- NULL
  out
}


#' The one-row verdict a table reading produces
#'
#' @param res The categorical comparison result.
#' @param tbl The one-row table of the test that was named.
#' @param measure,effect_cutoff,pval_cutoff The arguments, `measure` resolved.
#'
#' @keywords internal
#' @noRd
sa_categorical_significance_row <- function(res, tbl, measure, effect_cutoff,
                                            pval_cutoff) {
  row <- res$association[res$association$measure == measure, , drop = FALSE]
  pvalue <- tbl$pval
  sa_check_pvalues(pvalue)

  out <- data.frame(
    measure    = measure,
    estimate   = row$estimate,
    lower_conf = row$lower_conf,
    upper_conf = row$upper_conf,
    pvalue     = pvalue,
    stringsAsFactors = FALSE
  )
  out$is_signif <- sa_assoc_clears(measure, row$estimate, effect_cutoff) &
    out$pvalue <= pval_cutoff
  rownames(out) <- NULL
  out
}


#' The measure a design defines, or the one that was named
#'
#' @param res The categorical comparison result.
#' @param measure `"auto"` or a row of `$association`.
#'
#' @keywords internal
#' @noRd
sa_resolve_measure <- function(res, measure) {
  if (!is.character(measure) || length(measure) != 1L || is.na(measure)) {
    stop("`measure` must be a single measure name, or \"auto\".", call. = FALSE)
  }

  if (identical(measure, "auto")) {
    dims <- res$design$dim
    return(switch(
      res$design$null,
      symmetry             = "odds_ratio_paired",
      marginal_homogeneity = "kendalls_w",
      if (length(dims) == 2L && all(dims == 2L)) "odds_ratio" else "cramers_v"
    ))
  }

  if (!measure %in% res$association$measure) {
    stop("`measure` must name one of the measures this design defines: ",
         paste(res$association$measure, collapse = ", "), ". Got ", measure,
         ". Which measures exist depends on the design and on the size of the ",
         "table, so a measure absent here is one this result has no value for ",
         "rather than one that was left out.", call. = FALSE)
  }
  measure
}


#' Which scale an association measure lives on
#'
#' Two of them are ratios centred at 1 and the rest are centred at zero, which is
#' the whole of what a cutoff has to know. The zero-centred ones are compared on
#' their magnitude, which is a no-op for the measures that cannot be negative and
#' the correct reading for the ones that can: `cohens_g` and
#' `risk_difference_paired` are below zero when the later condition lowers the
#' response, and a departure downwards is as large as the same one upwards.
#'
#' @keywords internal
#' @noRd
sa_assoc_scale <- function(measure) {
  if (measure %in% c("odds_ratio", "odds_ratio_paired")) "ratio" else "magnitude"
}


#' Whether an estimate reaches the cutoff on its own scale
#'
#' @keywords internal
#' @noRd
sa_assoc_clears <- function(measure, estimate, cutoff) {
  if (is.null(cutoff)) {
    return(TRUE)
  }
  if (identical(sa_assoc_scale(measure), "ratio")) {
    return(estimate >= cutoff | estimate <= 1 / cutoff)
  }
  abs(estimate) >= cutoff
}


#' Refuse a cutoff that cannot mean what it says on this measure
#'
#' A ratio cutoff below 1 admits every table rather than a stricter set of them,
#' since `estimate >= c | estimate <= 1 / c` covers the whole line once `c` drops
#' under 1. That is a silently empty demand, so it is an error instead.
#'
#' @keywords internal
#' @noRd
sa_check_effect_cutoff <- function(cutoff, measure) {
  if (is.null(cutoff)) {
    return(invisible(NULL))
  }
  sa_check_scalar_num(cutoff, "effect_cutoff", 0, lower_open = TRUE)
  if (identical(sa_assoc_scale(measure), "ratio") && cutoff < 1) {
    stop("`effect_cutoff` is read on the scale of `", measure,
         "`, a ratio centred at 1, so it is a fold either way and has to be ",
         "at least 1. A cutoff of ", cutoff, " asks for an estimate above ",
         cutoff, " or below ", format(1 / cutoff),
         ", which every value meets.", call. = FALSE)
  }
  invisible(NULL)
}


#' The test table a table reading names
#'
#' `sa_pick_test()` is the counterpart for the comparison scenarios and requires
#' an `sa_comparison`, which this result deliberately is not.
#'
#' @keywords internal
#' @noRd
sa_pick_categorical_test <- function(res, test) {
  if (!is.character(test) || length(test) != 1L || is.na(test)) {
    stop("`test` must be a single test name.", call. = FALSE)
  }
  if (!test %in% names(res$tests)) {
    stop("`test` must name one of the tests in ",
         "`categorical_comparison_result`: ",
         paste(names(res$tests), collapse = ", "), ". Got ", test, ".",
         call. = FALSE)
  }
  res$tests[[test]]
}


#' Say so when a setting the chosen reading does not read was supplied
#'
#' The two readings take their p-value from different places, so about half of the
#' arguments do nothing under either one. A setting that changes nothing is worth
#' a sentence rather than silence, which is the choice
#' `compare_categorical_groups()` makes about `exact` and `simulate_p_value`.
#'
#' @keywords internal
#' @noRd
sa_warn_unread_args <- function(by, test_missing, measure_missing,
                                cutoff_missing, lift_missing, adj_missing) {
  supplied <- if (identical(by, "cell")) {
    c("test", "measure", "effect_cutoff")[
      !c(test_missing, measure_missing, cutoff_missing)]
  } else {
    c("log2_lift_cutoff", "adj_type")[!c(lift_missing, adj_missing)]
  }
  if (length(supplied) == 0L) {
    return(invisible(NULL))
  }

  warning(paste0("`", supplied, "`", collapse = " and "),
          if (length(supplied) > 1L) " are" else " is",
          " not read by `by = \"", by, "\"` and ",
          if (length(supplied) > 1L) "were" else "was", " ignored. ",
          if (identical(by, "cell")) {
            paste0("A cell reading takes its p-value from the cell's own ",
                   "standardized residual, which no test and no association ",
                   "measure reports.")
          } else {
            paste0("A table reading has one p-value and no family to adjust ",
                   "across, and its effect axis is `measure` rather than a ",
                   "lift.")
          }, call. = FALSE)
  invisible(NULL)
}


#' Wrap a categorical verdict in the object the user sees
#'
#' @keywords internal
#' @noRd
sa_new_categorical_significance <- function(analysis_type, significance) {
  structure(list(analysis_type = analysis_type,
                 significance  = significance),
            class = c("sa_categorical_significance", "sa_result"))
}


#' The attributes both readings describe themselves with
#'
#' @keywords internal
#' @noRd
sa_categorical_significance_attrs <- function(res, by, pval_cutoff) {
  # `table_dim` rather than `dim`, which on a data.frame is the number of rows
  # and columns of the verdict itself and cannot be borrowed for anything else.
  list(analysis    = res$analysis,
       null        = res$design$null,
       by          = by,
       table_dim   = res$design$dim,
       pval_cutoff = pval_cutoff)
}


#' Print a categorical significance verdict
#'
#' Summarises the rule that was applied and what cleared it, rather than printing
#' the table itself. Reach into `x$significance` for that.
#'
#' @param x An `sa_categorical_significance` object, as returned by
#'   [estimate_categorical_significance()].
#' @param ... Ignored, present for consistency with [print()].
#'
#' @return `x` invisibly.
#'
#' @examples
#' res <- compare_categorical_groups(
#'   data.frame(smoker = rep(c("y", "n"), each = 30),
#'              grade  = c(rep(c("high", "low"), c(22, 8)),
#'                         rep(c("high", "low"), c(9, 21))))
#' )
#' estimate_categorical_significance(res)
#' estimate_categorical_significance(res, by = "table")
#'
#' @export
print.sa_categorical_significance <- function(x, ...) {
  tbl <- x$significance
  by <- attr(tbl, "by")
  dims <- attr(tbl, "table_dim")

  cat("<", class(x)[1], "> ", x$analysis_type, "\n", sep = "")
  cat("  reading  : ", by, "  (", paste(dims, collapse = " x "), " table)\n",
      sep = "")
  cat("  null     : ", sa_null_label(attr(tbl, "null")), "\n", sep = "")

  if (identical(by, "cell")) {
    cat("  cutoffs  : abs(log2_lift) >= ", attr(tbl, "log2_lift_cutoff"),
        ", adj_pvalue <= ", attr(tbl, "pval_cutoff"),
        "  (", attr(tbl, "adj_type"), ")\n", sep = "")

    n_signif <- sum(tbl$is_signif %in% TRUE)
    n_undecided <- sum(is.na(tbl$is_signif))
    cat("  verdict  : ", n_signif, " of ", nrow(tbl), " cell(s) significant",
        if (n_undecided > 0L) paste0("  (", n_undecided, " undecided)"),
        "\n", sep = "")

    hits <- which(tbl$is_signif %in% TRUE)
    if (length(hits) > 0L) {
      cat("\n  cells\n")
      labels <- paste(tbl$row_level[hits], tbl$col_level[hits], sep = " : ")
      width <- max(nchar(labels))
      for (i in seq_along(hits)) {
        at <- hits[i]
        cat("    ", formatC(labels[i], width = -width),
            "  lift = ", sa_fmt_est(tbl$lift[at]),
            ", adj_pvalue = ",
            format.pval(tbl$adj_pvalue[at], digits = 3, eps = 1e-16),
            "\n", sep = "")
      }
    }

    return(invisible(x))
  }

  cat("  test     : ", attr(tbl, "test"), "  (", attr(tbl, "test_label"),
      ")\n", sep = "")
  cutoff <- attr(tbl, "effect_cutoff")
  cat("  cutoffs  : ",
      if (is.null(cutoff)) {
        ""
      } else if (identical(sa_assoc_scale(tbl$measure), "ratio")) {
        paste0(tbl$measure, " >= ", cutoff, " or <= ", format(1 / cutoff), ", ")
      } else {
        paste0("abs(", tbl$measure, ") >= ", cutoff, ", ")
      },
      "pvalue <= ", attr(tbl, "pval_cutoff"), "\n", sep = "")
  cat("  verdict  : ", tbl$measure, " = ", sa_fmt_est(tbl$estimate),
      if (!is.na(tbl$lower_conf)) {
        paste0("  [", sa_fmt_est(tbl$lower_conf), ", ",
               sa_fmt_est(tbl$upper_conf), "]")
      },
      "  (", if (isTRUE(tbl$is_signif)) {
        "significant"
      } else if (isFALSE(tbl$is_signif)) {
        "not significant"
      } else {
        "undecided"
      }, ")\n", sep = "")

  invisible(x)
}
