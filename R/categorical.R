# The result contract for the categorical scenario, and the one place that says
# why it is a contract of its own rather than another `sa_comparison`.
#
# Every comparison scenario answers one question per numeric feature, which is
# what lets `sa_new_comparison()` hold each of its tables to `identical(df$features,
# features)` and what lets `estimate_significance()` read one `log2fc` per row.
# A contingency table has no such axis. The question is asked once, of the table
# as a whole, and at that level there is no signed effect to put beside the
# p-value: an association is not signed, and the odds ratio that is signed is
# only defined when both variables are binary.
#
# So this scenario is deliberately not an `sa_comparison`, for the same reason
# `diagnose_distribution()` is not one. What it keeps is the vocabulary --
# `tests` beside `test_info`, `design`, `parameters`, `metadata` -- so that a
# reader who knows one result object can read this one.
#
# One level down there is an axis, and `estimate_categorical_significance()` is
# what reads it. A cell was expected at some count and observed at another, so
# `observed / expected` says how far it moved and which way, and it is defined
# whatever the shape of the table; the standardized residual beside it is built
# to be referred to a standard normal, so the cell has a p-value of its own. That
# is a verdict axis, and it is why the absence of one here is a statement about
# the table rather than about the scenario. `$cells` is what both of them are
# computed from, so nothing is stored for them.
#
# What this contract adds to that vocabulary is `design$null`. A contingency
# table can be held against more than one null hypothesis, the tests in this
# scenario hold it against three, and `expected` means a different number under
# each. Naming which one produced the cell table is what keeps the residuals, the
# diagnostics and the shading of the mosaic talking about the same hypothesis the
# p-value beside them is about.
#
# Everything here is a scalar, a character vector, a named list or a data.frame.
# No matrix, no table and no engine object is stored, which is what lets the
# object be written straight out as JSON and rebuilt in another language;
# `as.table()` is how the table itself is read back out.


#' The null hypotheses a cell table can be built against
#'
#' Three, one per design, and every one of them is a statement the object also
#' carries a p-value for.
#'
#' \describe{
#'   \item{`"independence"`}{Two variables cross-classified on one sample. The
#'     expected count of a cell is the product of its margins over the total,
#'     which is what the chi-square test of independence and Fisher's exact test
#'     are both about.}
#'   \item{`"symmetry"`}{One thing measured twice, crossed against itself. The
#'     expected count of a cell is the average of it and its transpose, so the
#'     diagonal is expected to be exactly what it is and only the discordant
#'     cells carry a residual. That is McNemar's test.}
#'   \item{`"marginal_homogeneity"`}{Three or more repeated conditions summarised
#'     as a condition-by-response table. Every condition is expected to show the
#'     pooled response rate, which on that table is arithmetically the same
#'     formula as independence and is a different claim about the world. That is
#'     Cochran's Q.}
#' }
#'
#' @keywords internal
#' @noRd
sa_categorical_nulls <- function() {
  c("independence", "symmetry", "marginal_homogeneity")
}


#' A one-line reading of what a null hypothesis says
#'
#' Used by `print.sa_categorical()` and by the mosaic key, so the two cannot
#' describe the same shading differently.
#'
#' @keywords internal
#' @noRd
sa_null_label <- function(null) {
  switch(
    null,
    independence = paste0(
      "independence -- a cell is expected at the product of its margins"
    ),
    symmetry = paste0(
      "symmetry -- a cell is expected at the average of it and its transpose"
    ),
    marginal_homogeneity = paste0(
      "marginal homogeneity -- every condition is expected at the pooled rate"
    ),
    null
  )
}


#' Column names the cell table must carry
#'
#' The canonical form of the table. A matrix would say the same thing more
#' compactly, but it is the one shape the result contract does not take: a
#' dimnamed matrix survives a trip through JSON in some writers and loses its
#' labels in others, and the labels are the whole content of a contingency
#' table. One row per cell keeps them beside the number they belong to.
#'
#' `expected`, `residual` and `std_residual` are all read under `design$null` and
#' mean nothing without it. `residual` is the Pearson residual, the quantity that
#' squares and sums to the statistic of the test that null belongs to, so it says
#' how that statistic was made up. `std_residual` is the standardized (adjusted)
#' residual, which is referred to a standard normal and so says which cells are
#' individually surprising; its variance correction is derived under independence
#' and has no counterpart under symmetry, so it is `NA` there rather than a
#' number that looks comparable and is not.
#'
#' @keywords internal
#' @noRd
sa_categorical_cell_columns <- function() {
  c("row_level", "col_level", "observed", "expected", "residual",
    "std_residual", "prop_total", "prop_row", "prop_col")
}


#' Column names every test table must carry
#'
#' One row rather than one row per feature, and no `pval_adj`. There is a single
#' table and therefore a single question, so there is no family to adjust across;
#' a column holding a copy of `pval` under another name would suggest otherwise.
#'
#' `statistic` and `df` are `NA` for a test that has neither, which Fisher's
#' exact test does not: it conditions on the margins and reads a probability off
#' the hypergeometric distribution rather than referring a statistic to a null
#' one. The columns exist for every test, as in the comparison contract; being
#' finite is not required.
#'
#' @keywords internal
#' @noRd
sa_categorical_test_columns <- function() {
  c("n_used", "statistic", "df", "pval", "lower_conf", "upper_conf")
}


#' Column names the association table must carry
#'
#' One row per measure rather than one column per measure, because which measures
#' exist depends on the design and on the size of the table, and a wide table
#' would carry a column of `NA` for every measure the design does not define.
#'
#' @keywords internal
#' @noRd
sa_association_columns <- function() {
  c("measure", "estimate", "lower_conf", "upper_conf")
}


#' Assemble a categorical comparison result
#'
#' The checks here guard the contract rather than the user's input, so they fire
#' only on a mistake inside the package and say so.
#'
#' @param analysis Scenario identifier, `"categorical_comparison"`.
#' @param variables The variable names the analysis was made on, in the order
#'   `category_lv` fixed.
#' @param design Named list describing the layout. Must carry `null`, one of
#'   `sa_categorical_nulls()`, and `dim`, the shape of the table the cells came
#'   from.
#' @param parameters Named list of the analysis choices the caller made.
#' @param cells data.frame of one row per cell of the table.
#' @param tests Named list of one one-row data.frame per test.
#' @param test_info Named list with one entry per element of `tests`.
#' @param association data.frame of one row per association measure.
#' @param diagnostics The approximation check, or `NULL` when it was not
#'   requested.
#' @param subclass Extra class prepended ahead of `sa_categorical`.
#'
#' @keywords internal
#' @noRd
sa_new_categorical <- function(analysis,
                               variables,
                               design,
                               parameters,
                               cells,
                               tests,
                               test_info,
                               association,
                               diagnostics = NULL,
                               subclass = character(0)) {

  if (!is.list(design) || is.null(design$null) ||
        !isTRUE(design$null %in% sa_categorical_nulls())) {
    stop("internal error: `design$null` must name one of ",
         paste(sa_categorical_nulls(), collapse = ", "), ", and holds ",
         paste(format(design$null), collapse = " "), ".", call. = FALSE)
  }
  if (length(design$dim) != 2L || any(design$dim < 1L)) {
    stop("internal error: `design$dim` must be the two dimensions of the ",
         "table.", call. = FALSE)
  }

  if (!is.list(tests) || length(tests) == 0L || is.null(names(tests))) {
    stop("internal error: `tests` must be a non-empty named list.",
         call. = FALSE)
  }
  if (!setequal(names(tests), names(test_info))) {
    stop("internal error: `tests` and `test_info` name different tests: ",
         paste(names(tests), collapse = ", "), " vs ",
         paste(names(test_info), collapse = ", "), ".", call. = FALSE)
  }

  for (nm in names(tests)) {
    what <- paste0("`tests$", nm, "`")
    if (!is.data.frame(tests[[nm]])) {
      stop("internal error: ", what, " must be a data.frame.", call. = FALSE)
    }
    if (nrow(tests[[nm]]) != 1L) {
      stop("internal error: ", what, " must hold exactly one row, one table ",
           "being one question, but holds ", nrow(tests[[nm]]), ".",
           call. = FALSE)
    }
    absent <- setdiff(sa_categorical_test_columns(), names(tests[[nm]]))
    if (length(absent) > 0L) {
      stop("internal error: ", what, " is missing contract column(s): ",
           paste(absent, collapse = ", "), ".", call. = FALSE)
    }
  }

  if (!is.data.frame(cells)) {
    stop("internal error: `cells` must be a data.frame.", call. = FALSE)
  }
  absent <- setdiff(sa_categorical_cell_columns(), names(cells))
  if (length(absent) > 0L) {
    stop("internal error: `cells` is missing contract column(s): ",
         paste(absent, collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(cells) != prod(design$dim)) {
    stop("internal error: `cells` holds ", nrow(cells), " row(s) for a ",
         paste(design$dim, collapse = " x "), " table.", call. = FALSE)
  }

  if (!is.data.frame(association)) {
    stop("internal error: `association` must be a data.frame.", call. = FALSE)
  }
  absent <- setdiff(sa_association_columns(), names(association))
  if (length(absent) > 0L) {
    stop("internal error: `association` is missing contract column(s): ",
         paste(absent, collapse = ", "), ".", call. = FALSE)
  }

  slots <- list(
    analysis    = analysis,
    variables   = variables,
    design      = design,
    parameters  = parameters,
    cells       = cells,
    tests       = tests,
    test_info   = test_info,
    association = association,
    diagnostics = diagnostics,
    metadata    = sa_metadata()
  )

  structure(slots, class = c(subclass, "sa_categorical", "sa_result"))
}


#' The contingency table a categorical comparison was run on
#'
#' `$cells` is the canonical form, one row per cell, because that is the shape
#' that survives being written out as JSON with its labels attached. A table is
#' the shape to read it in, so it is built on request rather than stored twice
#' and left to drift.
#'
#' The rows dropped for a missing value or for a level outside `category_lv` are
#' already gone, so this is the table the tests were run on rather than a fresh
#' [table()] of the input.
#'
#' @param x A categorical comparison, as returned by
#'   [compare_categorical_groups()].
#' @param ... Ignored, present for consistency with [as.table()].
#'
#' @return A two-dimensional [table()] of counts, its dimnames named after the
#'   two axes of the design.
#'
#' @examples
#' res <- compare_categorical_groups(
#'   data.frame(smoker = rep(c("y", "n"), each = 30),
#'              grade  = c(rep(c("high", "low"), c(22, 8)),
#'                         rep(c("high", "low"), c(9, 21))))
#' )
#' as.table(res)
#'
#' ## Which is the table the tests were run on, so it feeds `stats` directly.
#' stats::chisq.test(as.table(res))
#'
#' @export
as.table.sa_categorical <- function(x, ...) {
  sa_categorical_table(x$cells, x$design$row_var, x$design$col_var)
}


#' Fold a cell table back into a contingency table
#'
#' Indexed by level name rather than by position, so it does not depend on the
#' order the cells happen to be in.
#'
#' @keywords internal
#' @noRd
sa_categorical_table <- function(cells, row_var = "row", col_var = "column") {
  row_lv <- unique(cells$row_level)
  col_lv <- unique(cells$col_level)

  out <- matrix(0, nrow = length(row_lv), ncol = length(col_lv),
                dimnames = stats::setNames(list(row_lv, col_lv),
                                           c(row_var, col_var)))
  at <- cbind(match(cells$row_level, row_lv), match(cells$col_level, col_lv))
  out[at] <- cells$observed
  # Counts are whole numbers, and an integer table is what `stats::chisq.test()`
  # and the rest of the engines are handed everywhere else.
  if (all(out == round(out))) {
    storage.mode(out) <- "integer"
  }
  as.table(out)
}


#' Print a categorical comparison
#'
#' Summarises the table, the null hypothesis it was held against, the tests run
#' on it and the association measures, rather than printing the cell table.
#' `as.table()` gives the table as counts, `$cells` the residuals and
#' `$association` the measures with their intervals.
#'
#' @param x A categorical comparison, as returned by
#'   [compare_categorical_groups()].
#' @param alpha Threshold the reported verdict on each test is read at.
#' @param ... Ignored, present for consistency with [print()].
#'
#' @return `x` invisibly.
#'
#' @examples
#' compare_categorical_groups(
#'   data.frame(smoker = c(rep("y", 30), rep("n", 30)),
#'              grade  = c(rep(c("high", "low"), c(22, 8)),
#'                         rep(c("high", "low"), c(9, 21))))
#' )
#'
#' @export
print.sa_categorical <- function(x, alpha = 0.05, ...) {
  sa_check_scalar_num(alpha, "alpha", 0, 1, lower_open = TRUE)

  design <- x$design
  params <- x$parameters

  cat("<", class(x)[1], "> ", x$analysis, "\n", sep = "")
  cat("  table    : ", design$row_var, " (", design$dim[1], ") x ",
      design$col_var, " (", design$dim[2], ")",
      "  (", prod(design$dim), " cells, ",
      if (isTRUE(design$paired)) {
        paste0("matched by ", design$pairing)
      } else {
        "independent"
      }, ")\n", sep = "")
  cat("  null     : ", sa_null_label(design$null), "\n", sep = "")
  cat("  observed : ", design$n_used, " row(s)\n", sep = "")
  cat("  settings : conf_level = ", params$conf_level,
      ", correct = ", params$correct,
      if (isTRUE(params$simulate_p_value)) {
        paste0(", simulated on ", params$n_resamples, " resample(s)")
      }, "\n", sep = "")

  cat("\n  tests\n")
  width <- max(nchar(names(x$tests)))
  for (nm in names(x$tests)) {
    tbl <- x$tests[[nm]]
    verdict <- if (is.na(tbl$pval)) {
      "not computed"
    } else if (tbl$pval <= alpha) {
      paste0("null rejected at ", alpha)
    } else {
      paste0("null retained at ", alpha)
    }
    cat("    $", formatC(nm, width = -width), "  pval = ",
        format.pval(tbl$pval, digits = 3, eps = 1e-16),
        "  (", verdict, ")\n", sep = "")
    cat("    ", strrep(" ", width + 2), x$test_info[[nm]]$label, "\n", sep = "")
  }

  cat("\n  association\n")
  width <- max(nchar(x$association$measure))
  for (i in seq_len(nrow(x$association))) {
    row <- x$association[i, ]
    cat("    ", formatC(row$measure, width = -width), "  ",
        sa_fmt_est(row$estimate),
        if (!is.na(row$lower_conf)) {
          paste0("  [", sa_fmt_est(row$lower_conf), ", ",
                 sa_fmt_est(row$upper_conf), "]")
        }, "\n", sep = "")
  }

  if (!is.null(x$diagnostics)) {
    cat("\n  $diagnostics attached, rule ", x$diagnostics$rule, ": ",
        if (isTRUE(x$diagnostics$approx_ok)) "met" else "not met",
        "\n", sep = "")
  }
  if (isTRUE(design$n_dropped > 0L)) {
    cat("\n  dropped  : ", design$n_dropped,
        " row(s) outside `category_lv`\n", sep = "")
  }
  if (isTRUE(design$n_incomplete > 0L)) {
    cat("  dropped  : ", design$n_incomplete,
        " row(s) missing a value the table needs\n", sep = "")
  }

  invisible(x)
}


#' Format an estimate for printing
#'
#' @keywords internal
#' @noRd
sa_fmt_est <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }
  format(x, digits = 3, trim = TRUE)
}
