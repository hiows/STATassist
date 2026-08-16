# The result contract shared by every comparison scenario. `compare_two_groups()`
# fills these slots and `estimate_significance()` reads them back, so a second
# scenario function reusing the same slots inherits the whole downstream path
# without either side knowing which tests were actually run.
#
# Everything in the object is a scalar, a character vector, a named list or a
# data.frame. No engine object is stored anywhere, which is what lets the object
# be written straight out as JSON and rebuilt in another language.

#' Reproducibility metadata attached to every result
#'
#' @keywords internal
#' @noRd
sa_metadata <- function() {
  list(
    package_version = as.character(utils::packageVersion("STATassist")),
    r_version       = R.version.string,
    platform        = R.version$platform,
    # ISO-8601 with an explicit offset, so the value survives a trip through
    # JSON into a language with a different default timezone.
    timestamp       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
}


#' Column names every test table in a comparison result must carry
#'
#' The per-test statistics differ, but these are the columns a consumer is
#' allowed to rely on regardless of which test it is looking at.
#'
#' @keywords internal
#' @noRd
sa_test_table_columns <- function() {
  c("features", "n_used", "pval", "pval_adj", "lower_conf", "upper_conf")
}


#' Column names every post-hoc table must carry
#'
#' A post-hoc table holds one row per feature and pair of levels rather than one
#' row per feature, which is the whole reason it lives in its own slot instead of
#' alongside the omnibus tables. `contrast` is the readable pair label and
#' `group1` / `group2` are the two levels it is made of, in the direction
#' `estimate` reads as `group1 - group2`.
#'
#' @keywords internal
#' @noRd
sa_posthoc_table_columns <- function() {
  c("features", "contrast", "group1", "group2", "n1", "n2", "estimate",
    "stderr", "statistic", "df", "pval", "pval_adj", "lower_conf",
    "upper_conf")
}


#' The post-hoc columns that carry a number rather than a label
#'
#' @keywords internal
#' @noRd
sa_posthoc_stat_columns <- function() {
  setdiff(sa_posthoc_table_columns(),
          c("features", "contrast", "group1", "group2"))
}


#' Column names a term table must carry
#'
#' A factorial analysis answers on two axes. Whether a feature responds to the
#' design at all is one question per feature and lives in `$tests`; which part of
#' the design it responds to is one question per feature and model term and lives
#' here, because a table of that shape cannot satisfy the one-row-per-feature
#' alignment every other table in the object is held to.
#'
#' `terms` is the term label in [stats::terms()] form, `a` for a main effect and
#' `a:b` for an interaction, and `term_order` is how many factors it is over. The
#' pair is what a simulator's `truth_term` is keyed on, so the two tables merge
#' without either side being renamed.
#'
#' `log2_effect` is here for the same reason `log2fc` is required of an effect
#' table: it is the only signed size on this axis, so it is the only column a
#' term can be given an effect axis from. `estimate_significance(by = "term")`
#' reads it by name.
#'
#' @keywords internal
#' @noRd
sa_term_table_columns <- function() {
  c("features", "terms", "term_order", "n_used", "df", "ss", "f_stat",
    "log2_effect", "pval", "pval_adj")
}


#' Column names a cell table must carry
#'
#' The cell means of a factorial comparison, one row per feature and cell of the
#' crossed grid. `$effect` reduces the same numbers to one signed size per
#' feature and `$tests` to one p-value, and neither can be read backwards into
#' the pattern that produced it, which is what an interaction plot draws.
#'
#' `cell` is the dot-joined level label `$design$group_lv` lists, and the table
#' also carries one column per factor, named after the factor and holding the
#' level name, so that a subset can be taken on a factor without parsing the
#' label. Those columns are why these six names are reserved: a factor may not
#' be called any of them.
#'
#' `mean` is the arithmetic cell mean the crossed model was fitted on, not a
#' centre on the `fc_mean` scale, so that a plot of this table and the F tests
#' beside it describe the same fit. `se` is `sqrt(ms_error / n)`, pooled over the
#' whole model rather than computed within the cell, which is what makes the
#' variance of any marginal mean recoverable from it: over a set of cells S it is
#' `sum((1 / length(S))^2 * se^2)`, the same expression the Tukey stage scales
#' its contrasts by. `sd` is the within-cell spread, which is not that and is
#' kept for the reader rather than for the arithmetic.
#'
#' @keywords internal
#' @noRd
sa_cell_table_columns <- function() {
  c("features", "cell", "n", "mean", "sd", "se")
}


#' Column names every pairwise table must carry
#'
#' The same numbers as a post-hoc table, rearranged into one rectangular table
#' per contrast so that a single contrast can be read on its own. Two columns
#' are added: `fold_change` and `log2fc`, the ratio of the two group centres,
#' which no post-hoc procedure reports because it does not depend on which test
#' was run.
#'
#' Unlike a post-hoc table, a pairwise table holds every feature of the
#' comparison, in the order `features` fixes. A feature whose omnibus test did
#' not clear `posthoc_alpha` is therefore present, with its ratio columns filled
#' and every inference column `NA`.
#'
#' @keywords internal
#' @noRd
sa_pairwise_table_columns <- function() {
  c("features", "contrast", "group1", "group2", "fold_change", "log2fc",
    sa_posthoc_stat_columns())
}


#' Assemble a comparison result object
#'
#' The checks here guard the contract rather than the user's input, so they fire
#' only on a mistake inside the package and say so.
#'
#' @param analysis Scenario identifier, such as `"two_group_comparison"`.
#' @param features Feature names, in the row order every table uses.
#' @param design Named list describing the data layout: group levels, whether
#'   the samples are paired, how pairs were formed, what was dropped.
#' @param parameters Named list of the analysis choices the caller made.
#' @param effect data.frame of the effect estimates shared by all tests, one row
#'   per feature.
#' @param tests Named list of one data.frame per test, one row per feature.
#' @param test_info Named list with one entry per element of `tests`, describing
#'   which test was actually run.
#' @param terms data.frame of one row per feature and model term, or `NULL` for a
#'   scenario whose model has a single term, in which case the slot is left out of
#'   the result.
#' @param cells data.frame of one row per feature and cell of a crossed grid, or
#'   `NULL` for a scenario that has no grid to report, in which case the slot is
#'   left out of the result.
#' @param posthoc Named list of post-hoc tables, named after the test each one
#'   follows, or `list()` when the scenario has no post-hoc stage, in which case
#'   the slot is left out of the result. Rows are feature by pair rather than one
#'   per feature.
#' @param pairwise Named list holding, per test, a list of one table per
#'   contrast, or `list()` when the scenario has no post-hoc stage. The same
#'   numbers as `posthoc` seen one contrast at a time, with every feature
#'   present.
#' @param diagnostics Assumption checks attached to this analysis, or `NULL`
#'   when they were not requested.
#' @param subclass Extra class prepended ahead of `sa_comparison`.
#'
#' @keywords internal
#' @noRd
sa_new_comparison <- function(analysis,
                              features,
                              design,
                              parameters,
                              effect,
                              tests,
                              test_info,
                              terms = NULL,
                              cells = NULL,
                              posthoc = list(),
                              pairwise = list(),
                              diagnostics = NULL,
                              subclass = character(0)) {

  if (!is.list(tests) || length(tests) == 0L || is.null(names(tests))) {
    stop("internal error: `tests` must be a non-empty named list.",
         call. = FALSE)
  }
  if (!setequal(names(tests), names(test_info))) {
    stop("internal error: `tests` and `test_info` name different tests: ",
         paste(names(tests), collapse = ", "), " vs ",
         paste(names(test_info), collapse = ", "), ".", call. = FALSE)
  }
  if (!is.list(posthoc)) {
    stop("internal error: `posthoc` must be a list.", call. = FALSE)
  }
  if (length(posthoc) > 0L && !all(names(posthoc) %in% names(tests))) {
    stop("internal error: `posthoc` names a test that was not run: ",
         paste(setdiff(names(posthoc), names(tests)), collapse = ", "), ".",
         call. = FALSE)
  }
  if (!is.list(pairwise)) {
    stop("internal error: `pairwise` must be a list.", call. = FALSE)
  }
  if (length(pairwise) > 0L && !all(names(pairwise) %in% names(tests))) {
    stop("internal error: `pairwise` names a test that was not run: ",
         paste(setdiff(names(pairwise), names(tests)), collapse = ", "), ".",
         call. = FALSE)
  }

  check_table <- function(df, what) {
    if (!is.data.frame(df)) {
      stop("internal error: ", what, " must be a data.frame.", call. = FALSE)
    }
    if (!identical(df$features, features)) {
      stop("internal error: ", what, " is not aligned with `features`.",
           call. = FALSE)
    }
  }

  check_table(effect, "`effect`")
  for (nm in names(tests)) {
    check_table(tests[[nm]], paste0("`tests$", nm, "`"))
    absent <- setdiff(sa_test_table_columns(), names(tests[[nm]]))
    if (length(absent) > 0L) {
      stop("internal error: `tests$", nm, "` is missing contract column(s): ",
           paste(absent, collapse = ", "), ".", call. = FALSE)
    }
  }

  # A term table holds one row per feature and term, so it is checked the way a
  # post-hoc table is: the features it names have to be features of this
  # comparison, but there is no position to align them by.
  if (!is.null(terms)) {
    if (!is.data.frame(terms)) {
      stop("internal error: `terms` must be a data.frame.", call. = FALSE)
    }
    absent <- setdiff(sa_term_table_columns(), names(terms))
    if (length(absent) > 0L) {
      stop("internal error: `terms` is missing contract column(s): ",
           paste(absent, collapse = ", "), ".", call. = FALSE)
    }
    unknown <- setdiff(terms$features, features)
    if (length(unknown) > 0L) {
      stop("internal error: `terms` holds feature(s) absent from the ",
           "comparison: ", paste(unique(unknown), collapse = ", "), ".",
           call. = FALSE)
    }
  }

  # A cell table is rectangular but not one row per feature, so it is checked on
  # membership like a term table and then on the block size the grid fixes.
  if (!is.null(cells)) {
    if (!is.data.frame(cells)) {
      stop("internal error: `cells` must be a data.frame.", call. = FALSE)
    }
    absent <- setdiff(sa_cell_table_columns(), names(cells))
    if (length(absent) > 0L) {
      stop("internal error: `cells` is missing contract column(s): ",
           paste(absent, collapse = ", "), ".", call. = FALSE)
    }
    if (!identical(unique(cells$features), features)) {
      stop("internal error: `cells` does not hold every feature of the ",
           "comparison once, in order.", call. = FALSE)
    }
  }

  # A post-hoc table is checked on membership rather than on identity, because
  # it holds several rows per feature and may legitimately skip a feature whose
  # omnibus test was not significant.
  for (nm in names(posthoc)) {
    what <- paste0("`posthoc$", nm, "`")
    if (!is.data.frame(posthoc[[nm]])) {
      stop("internal error: ", what, " must be a data.frame.", call. = FALSE)
    }
    absent <- setdiff(sa_posthoc_table_columns(), names(posthoc[[nm]]))
    if (length(absent) > 0L) {
      stop("internal error: ", what, " is missing contract column(s): ",
           paste(absent, collapse = ", "), ".", call. = FALSE)
    }
    unknown <- setdiff(posthoc[[nm]]$features, features)
    if (length(unknown) > 0L) {
      stop("internal error: ", what, " holds feature(s) absent from the ",
           "comparison: ", paste(unique(unknown), collapse = ", "), ".",
           call. = FALSE)
    }
  }

  # A pairwise table, unlike a post-hoc one, is rectangular, so it is held to
  # the same alignment every other table in the object is held to.
  for (nm in names(pairwise)) {
    if (!is.list(pairwise[[nm]]) ||
          !is.character(names(pairwise[[nm]])) ||
          anyNA(names(pairwise[[nm]]))) {
      stop("internal error: `pairwise$", nm, "` must be a list named by ",
           "contrast.", call. = FALSE)
    }
    for (ct in names(pairwise[[nm]])) {
      what <- paste0("`pairwise$", nm, "[[\"", ct, "\"]]`")
      check_table(pairwise[[nm]][[ct]], what)
      absent <- setdiff(sa_pairwise_table_columns(),
                        names(pairwise[[nm]][[ct]]))
      if (length(absent) > 0L) {
        stop("internal error: ", what, " is missing contract column(s): ",
             paste(absent, collapse = ", "), ".", call. = FALSE)
      }
    }
  }

  slots <- list(
    analysis    = analysis,
    features    = features,
    design      = design,
    parameters  = parameters,
    effect      = effect,
    tests       = tests,
    terms       = terms,
    cells       = cells,
    posthoc     = posthoc,
    pairwise    = pairwise,
    test_info   = test_info,
    diagnostics = diagnostics,
    metadata    = sa_metadata()
  )

  # A scenario whose model has one term says everything it has to say about that
  # term in `$tests`, so the slot is dropped rather than left empty, the same way
  # the three below are. `list()` keeps a `NULL` element rather than discarding
  # it, which is why this is written out. `cells` goes the same way: a scenario
  # with no crossed grid has no cell to report a mean for, and one level of a
  # single factor is a group rather than a cell.
  #
  # A scenario with no pairwise stage has nothing to say in these two slots, and
  # an empty map reads as a result that is missing something rather than as one
  # for which the question does not arise. They are absent instead, so a reader
  # asks `is.null(res$posthoc)`, which is also what `res$posthoc[[test]]`
  # already answers for a test that has no post-hoc stage.
  if (is.null(terms)) slots$terms <- NULL
  if (is.null(cells)) slots$cells <- NULL
  if (length(posthoc) == 0L) slots$posthoc <- NULL
  if (length(pairwise) == 0L) slots$pairwise <- NULL

  structure(slots, class = c(subclass, "sa_comparison", "sa_result"))
}


#' Pull one test table out of a comparison result
#'
#' Shared by every function that lets the user name a test, so the error message
#' listing the valid choices is written once.
#'
#' @keywords internal
#' @noRd
sa_pick_test <- function(res, test, arg) {
  if (!inherits(res, "sa_comparison")) {
    stop("`", arg, "` must be a comparison result, as returned by ",
         "compare_two_groups().", call. = FALSE)
  }
  if (!is.character(test) || length(test) != 1L || is.na(test)) {
    stop("`test` must be a single test name.", call. = FALSE)
  }
  if (!test %in% names(res$tests)) {
    stop("`test` must name one of the tests in `", arg, "`: ",
         paste(names(res$tests), collapse = ", "), ". Got ", test, ".",
         call. = FALSE)
  }
  res$tests[[test]]
}


#' Print a comparison result
#'
#' Summarises which tests were run and how many features each one called
#' significant, rather than printing the tables themselves. Reach into
#' `x$tests` for those, and into `x$posthoc` for the pairwise stage.
#'
#' @param x A comparison result, as returned by [compare_two_groups()],
#'   [compare_multiple_groups()] or [compare_one_sample()].
#' @param alpha Threshold applied to `pval_adj` for the significant counts.
#' @param ... Ignored, present for consistency with [print()].
#'
#' @return `x` invisibly.
#'
#' @examples
#' iris2 <- iris[iris$Species != "setosa", ]
#' compare_two_groups(iris2, c("Sepal.Length", "Petal.Length"),
#'                    iris2$Species, c("versicolor", "virginica"))
#'
#' @export
print.sa_comparison <- function(x, alpha = 0.05, ...) {
  sa_check_scalar_num(alpha, "alpha", 0, 1, lower_open = TRUE)

  design <- x$design
  params <- x$parameters

  cat("<", class(x)[1], "> ", x$analysis, "\n", sep = "")
  # A one-sample comparison has a hypothesised value where the others have group
  # levels, so the header line is chosen from what the design actually holds
  # rather than from the analysis name.
  if (!is.null(design$factor_lv)) {
    # A crossed design's group levels are its cells, and eight of them listed as
    # "a1.b1 vs a1.b2 vs ..." says less than the two factors they came from.
    cat("  factors  : ",
        paste(paste0(names(design$factor_lv), " (",
                     lengths(design$factor_lv), ")"), collapse = " x "),
        "  (", length(design$group_lv), " cells, independent)\n", sep = "")
    cat("  anova    : ", sub("_", "-", design$anova_type),
        ", Type ", params$ss_type, " sums of squares\n", sep = "")
  } else if (is.null(design$group_lv)) {
    cat("  mu       : ", design$mu, "\n", sep = "")
  } else {
    cat("  groups   : ", paste(design$group_lv, collapse = " vs "),
        if (isTRUE(design$paired)) {
          paste0("  (paired by ", design$pairing, ")")
        } else {
          "  (independent)"
        }, "\n", sep = "")
  }
  cat("  features : ", length(x$features), "\n", sep = "")
  cat("  settings : alternative = ", params$alternative,
      ", conf_level = ", params$conf_level,
      ", p_adjust = ", params$p_adjust, "\n", sep = "")

  cat("\n  tests\n")
  width <- max(nchar(names(x$tests)))
  for (nm in names(x$tests)) {
    tbl <- x$tests[[nm]]
    n_signif <- sum(!is.na(tbl$pval_adj) & tbl$pval_adj <= alpha)
    n_failed <- sum(is.na(tbl$pval))
    cat("    $", formatC(nm, width = -width), "  ",
        n_signif, " of ", nrow(tbl), " at pval_adj <= ", alpha,
        if (n_failed > 0L) paste0("  (", n_failed, " not computed)"),
        "\n", sep = "")
    cat("    ", strrep(" ", width + 2), x$test_info[[nm]]$label, "\n", sep = "")
    ph <- x$posthoc[[nm]]
    if (!is.null(ph) && nrow(ph) > 0L) {
      n_pairs <- sum(!is.na(ph$pval_adj) & ph$pval_adj <= alpha)
      cat("    ", strrep(" ", width + 2), "post-hoc: ", n_pairs, " of ",
          nrow(ph), " contrast(s) over ", length(unique(ph$features)),
          " feature(s), ", x$test_info[[nm]]$posthoc_label, "\n", sep = "")
    }
  }

  # The whole-model row above says that a feature responds to the design. Which
  # part of the design it responds to is the question a crossed model was fitted
  # to answer, so it is summarised here rather than left for the reader to reach
  # into `$terms` for.
  if (!is.null(x$terms)) {
    cat("\n  terms\n")
    labels <- unique(x$terms$terms)
    width <- max(nchar(labels))
    for (nm in labels) {
      rows <- x$terms[x$terms$terms == nm, , drop = FALSE]
      n_signif <- sum(!is.na(rows$pval_adj) & rows$pval_adj <= alpha)
      cat("    ", formatC(nm, width = -width), "  ", n_signif, " of ",
          nrow(rows), " at pval_adj <= ", alpha, "\n", sep = "")
    }
  }

  if (!is.null(x$diagnostics)) {
    cat("\n  $diagnostics attached\n")
  }
  if (length(design$unmatched_ids) > 0L) {
    cat("\n  dropped  : ", length(design$unmatched_ids),
        " unpaired id(s)\n", sep = "")
  }
  if (isTRUE(design$n_dropped > 0L)) {
    cat("  dropped  : ", design$n_dropped, " row(s) outside `group_lv`\n",
        sep = "")
  }

  invisible(x)
}
