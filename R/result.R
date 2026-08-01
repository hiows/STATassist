# The result contract shared by every comparison scenario. `compare_two_groups()`
# fills these slots and `estimate_significance()` reads them back, so a second
# scenario function reusing the same slots inherits the whole downstream path
# without either side knowing which tests were actually run.
#
# Everything in the object is a scalar, a character vector, a named list or a
# data.frame. No engine object is stored anywhere, which is what lets the object
# be written straight out as JSON and rebuilt in another language.

#' Version of the comparison result contract
#'
#' Independent of the package version. It moves only when the shape of the
#' object changes: a new field is a minor bump, a removed or reinterpreted field
#' is a major one. Consumers in other languages read this to decide whether they
#' understand the object they were handed.
#'
#' @keywords internal
#' @noRd
sa_schema_version <- function() {
  "0.2.0"
}


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
#' @param posthoc Named list of post-hoc tables, named after the test each one
#'   follows, or `list()` when the scenario has no post-hoc stage. Rows are
#'   feature by pair rather than one per feature.
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
                              posthoc = list(),
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
  # A list with no names serialises to a JSON array, and the contract says
  # `posthoc` is a map from test name to table in every result. Naming the empty
  # case keeps it an object there, so a consumer that iterates over the map does
  # not have to special-case the scenarios that have no pairwise stage.
  if (length(posthoc) == 0L) {
    posthoc <- structure(list(), names = character(0))
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

  structure(
    list(
      schema_version = sa_schema_version(),
      analysis       = analysis,
      features       = features,
      design         = design,
      parameters     = parameters,
      effect         = effect,
      tests          = tests,
      posthoc        = posthoc,
      test_info      = test_info,
      diagnostics    = diagnostics,
      metadata       = sa_metadata()
    ),
    class = c(subclass, "sa_comparison", "sa_result")
  )
}


#' Pull one test table out of a comparison result
#'
#' Shared by every function that lets the user name a test, so the error message
#' listing the valid choices is written once.
#'
#' @keywords internal
#' @noRd
sa_pick_test <- function(res, test, arg = "res") {
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
  if (is.null(design$group_lv)) {
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
