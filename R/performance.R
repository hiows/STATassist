# The result contract shared by the two evaluation scenarios. `sa_model` is
# organised around terms, because a model has an outcome and a set of them. An
# evaluation has neither: it has a set of *models*, scored on rows none of them
# was fitted on, so `models` takes the place `terms` holds over there and every
# table in the object repeats it in the same order.
#
# The rule `sa_model` breaks is kept here. Nothing in this object is an engine
# object: the calibration line of a regression is two numbers rather than the
# `lm` that produced them, and a ROC curve is a table of thresholds rather than
# whatever computed it. So everything is a scalar, a character vector, a named
# list or a data.frame, and the whole object writes out as JSON. An evaluation
# has nothing to predict with, which is what made `$fit` unavoidable for a model
# and makes it unnecessary here.


#' Column names the prediction table carries
#'
#' Long rather than wide, so that a model name cannot collide with a column name
#' and so that the table survives being written out with its labels attached.
#' `row` is the position in the `newdata` that was scored, kept because the rows
#' that survive are an intersection and are therefore not `seq_len(n)`.
#'
#' @keywords internal
#' @noRd
sa_prediction_table_columns <- function() {
  c("model", "row", "observed", "predicted")
}


#' Column names the metric table of a regression evaluation carries
#'
#' `cor` and `r_squared` are both here because they are different questions and
#' agree only for predictions that need no calibration. `r_squared` is
#' `1 - SSE/SST` on the held-out rows, so it is the fraction of the variance the
#' predictions actually removed and can be negative; `cor^2` is what it would be
#' if the predictions were first rescaled by a line fitted to these same rows.
#' The gap between the two is exactly what `calib_slope` and `calib_intercept`
#' report, which is why all four travel together.
#'
#' `rmse` and `mae` carry the names `fit_*()` already uses in `$performance`, so
#' that a resampled number and a held-out number read in the same unit under the
#' same name.
#'
#' @keywords internal
#' @noRd
sa_regression_metric_columns <- function() {
  c("model", "n_used", "cor", "r_squared", "rmse", "mae", "bias",
    "calib_slope", "calib_intercept")
}


#' Column names the metric table of a classification evaluation carries
#'
#' `auc` and `brier` need no threshold. The three after them do, and it is the
#' one `parameters$threshold` records, so the table says what it was measured
#' at rather than implying a cut that was never named. Their names are the ones
#' `fit_*()` uses in `fit_stats`, and they are read against `outcome_lv[2]` the
#' same way.
#'
#' @keywords internal
#' @noRd
sa_classification_metric_columns <- function() {
  c("model", "n_used", "n_events", "auc", "auc_lower_conf", "auc_upper_conf",
    "brier", "accuracy", "sensitivity", "specificity")
}


#' Column names the comparison table of a regression evaluation carries
#'
#' Every column is `new - baseline`, so a positive `delta_cor` and a negative
#' `delta_rmse` both say the new model did better. There is no p-value here on
#' purpose: a difference of held-out errors has no null this function is in a
#' position to state, and the two numbers it is made of are in `$metrics`.
#'
#' @keywords internal
#' @noRd
sa_regression_comparison_columns <- function() {
  c("model", "delta_cor", "delta_r_squared", "delta_rmse", "delta_mae")
}


#' Column names the comparison table of a classification evaluation carries
#'
#' Three questions about the same pair of models, which is why all three are
#' here rather than one being picked. `delta_auc` asks whether the ranking
#' improved, and DeLong's test is paired because both models ranked the same
#' rows. `idi` asks how much further apart the predicted probabilities of the
#' two classes moved, and `nri` how often a probability moved the right way,
#' neither of which a difference of AUCs is sensitive to.
#'
#' @keywords internal
#' @noRd
sa_classification_comparison_columns <- function() {
  c("model",
    "delta_auc", "delta_auc_lower_conf", "delta_auc_upper_conf",
    "delta_auc_pval",
    "idi", "idi_lower_conf", "idi_upper_conf", "idi_pval",
    "nri", "nri_event", "nri_nonevent", "nri_lower_conf", "nri_upper_conf",
    "nri_pval")
}


#' Column names the ROC curve table carries
#'
#' One row per model and threshold, rather than one curve object per model. A
#' curve is a set of points and this is the set, which both draws directly and
#' writes out as JSON, neither of which a `pROC` object does.
#'
#' The first `threshold` of each curve is `Inf`, the cut above every prediction,
#' where nothing is called an event. It is the one cell of the whole object that
#' JSON has no number for and comes back as `null`; the coordinates it labels
#' survive, and a reader that only draws the curve loses nothing. No finite
#' threshold can stand in for it, since a curve reaching both corners needs one
#' cut outside the range of the predictions at one end or the other.
#'
#' @keywords internal
#' @noRd
sa_roc_curve_columns <- function() {
  c("model", "threshold", "sensitivity", "specificity")
}


#' Assemble an evaluation result object
#'
#' The checks here guard the contract rather than the user's input, so they fire
#' only on a mistake inside the package and say so.
#'
#' @param analysis `"regression_performance"` or `"classification_performance"`.
#' @param models Model names, in the row order every table uses. The baseline is
#'   the first of them.
#' @param design Named list describing the rows that were scored: the `outcome`
#'   label, its type, the `baseline` model's name, and the row counts.
#' @param parameters Named list of the scoring choices.
#' @param predictions data.frame of one row per model and scored row.
#' @param metrics data.frame of one row per model.
#' @param comparisons data.frame of one row per model other than the baseline,
#'   or `NULL` when nothing was compared against it, in which case the slot is
#'   left out of the result.
#' @param curves data.frame of one row per model and threshold, or `NULL` for a
#'   regression, in which case the slot is left out of the result.
#'
#' @keywords internal
#' @noRd
sa_new_performance <- function(analysis,
                               models,
                               design,
                               parameters,
                               predictions,
                               metrics,
                               comparisons = NULL,
                               curves = NULL) {

  if (!is.character(models) || length(models) == 0L || anyNA(models)) {
    stop("internal error: `models` must be a non-empty character vector.",
         call. = FALSE)
  }
  if (anyDuplicated(models) > 0L) {
    stop("internal error: `models` must be unique.", call. = FALSE)
  }
  if (!analysis %in% c("regression_performance",
                       "classification_performance")) {
    stop("internal error: unknown analysis `", analysis, "`.", call. = FALSE)
  }
  if (!identical(design$baseline, models[1])) {
    stop("internal error: `design$baseline` is not the first of `models`.",
         call. = FALSE)
  }

  metric_columns <- if (identical(analysis, "regression_performance")) {
    sa_regression_metric_columns()
  } else {
    sa_classification_metric_columns()
  }
  comparison_columns <- if (identical(analysis, "regression_performance")) {
    sa_regression_comparison_columns()
  } else {
    sa_classification_comparison_columns()
  }

  check_columns <- function(df, what, want) {
    absent <- setdiff(want, names(df))
    if (length(absent) > 0L) {
      stop("internal error: ", what, " is missing contract column(s): ",
           paste(absent, collapse = ", "), ".", call. = FALSE)
    }
  }

  if (!is.data.frame(metrics)) {
    stop("internal error: `metrics` must be a data.frame.", call. = FALSE)
  }
  check_columns(metrics, "`metrics`", metric_columns)
  if (!identical(metrics$model, models)) {
    stop("internal error: `metrics` is not aligned with `models`.",
         call. = FALSE)
  }

  if (!is.data.frame(predictions)) {
    stop("internal error: `predictions` must be a data.frame.", call. = FALSE)
  }
  check_columns(predictions, "`predictions`", sa_prediction_table_columns())
  # Every model was scored on the same rows, which is the whole point of the
  # intersection the caller took, so the table holds each model once and in
  # order rather than merely holding known names.
  if (!identical(unique(predictions$model), models)) {
    stop("internal error: `predictions` does not hold every model once, in ",
         "order.", call. = FALSE)
  }

  # A comparison is against the baseline, so the baseline has no row of its own
  # and a name that is not a model of this evaluation has nowhere to come from.
  if (!is.null(comparisons)) {
    if (!is.data.frame(comparisons)) {
      stop("internal error: `comparisons` must be a data.frame or NULL.",
           call. = FALSE)
    }
    check_columns(comparisons, "`comparisons`", comparison_columns)
    if (!identical(comparisons$model, models[-1])) {
      stop("internal error: `comparisons` must hold every model other than ",
           "the baseline, once and in order.", call. = FALSE)
    }
  }

  if (!is.null(curves)) {
    if (identical(analysis, "regression_performance")) {
      stop("internal error: a regression evaluation has no ROC curve.",
           call. = FALSE)
    }
    if (!is.data.frame(curves)) {
      stop("internal error: `curves` must be a data.frame or NULL.",
           call. = FALSE)
    }
    check_columns(curves, "`curves`", sa_roc_curve_columns())
    unknown <- setdiff(curves$model, models)
    if (length(unknown) > 0L) {
      stop("internal error: `curves` holds model(s) absent from the ",
           "evaluation: ", paste(unique(unknown), collapse = ", "), ".",
           call. = FALSE)
    }
  }

  slots <- list(
    analysis    = analysis,
    models      = models,
    design      = design,
    parameters  = parameters,
    predictions = predictions,
    metrics     = metrics,
    comparisons = comparisons,
    curves      = curves,
    metadata    = sa_metadata()
  )

  # An evaluation of one model compared it to nothing, and a regression has no
  # curve to report. An empty table in either slot would read as a result that
  # lost its values rather than as one for which the question does not arise, so
  # the slot goes instead and `is.null(res$comparisons)` is the test.
  if (is.null(comparisons)) slots$comparisons <- NULL
  if (is.null(curves)) slots$curves <- NULL

  structure(slots, class = c("sa_performance", "sa_result"))
}


#' Refuse anything that is not the evaluation this plot draws
#'
#' The two pictures are not interchangeable, and the object that carries one
#' also carries the other's slots' worth of names, so a classification handed to
#' the scatter would otherwise fail on a missing column rather than on the
#' mistake that was made.
#'
#' @keywords internal
#' @noRd
sa_performance_input <- function(x, want, arg, other) {
  if (!inherits(x, "sa_performance")) {
    if (inherits(x, "sa_model")) {
      stop("`", arg, "` is a fitted model rather than an evaluation of one. ",
           "Score it on held-out rows with evaluate_regression_models() or ",
           "evaluate_classification_models() first.", call. = FALSE)
    }
    stop("`", arg, "` must be an evaluation result, as returned by ",
         "evaluate_regression_models() or evaluate_classification_models().",
         call. = FALSE)
  }
  if (!identical(x$analysis, want)) {
    stop("`", arg, "` is a ", x$analysis, " result. Use ", other,
         " for that one.", call. = FALSE)
  }
  invisible(x)
}


#' Choose which models to draw, and in what order
#'
#' `NULL` draws all of them in the order the evaluation holds, which puts the
#' baseline first.
#'
#' @keywords internal
#' @noRd
sa_performance_models <- function(x, models) {
  if (is.null(models)) {
    return(x$models)
  }
  if (!is.character(models) || length(models) == 0L || anyNA(models)) {
    stop("`models` must be a non-empty character vector of model names, or ",
         "NULL for every model in the order the evaluation holds them.",
         call. = FALSE)
  }
  unknown <- setdiff(models, x$models)
  if (length(unknown) > 0L) {
    stop("`models` names model(s) the evaluation does not hold: ",
         paste(unknown, collapse = ", "), ". Available: ",
         paste(x$models, collapse = ", "), ".", call. = FALSE)
  }
  dup <- unique(models[duplicated(models)])
  if (length(dup) > 0L) {
    stop("`models` contains duplicated names: ",
         paste(dup, collapse = ", "), ".", call. = FALSE)
  }
  models
}


#' One colour per drawn model
#'
#' @keywords internal
#' @noRd
sa_performance_colours <- function(n, col) {
  if (is.null(col)) {
    return(grDevices::hcl.colors(max(n, 2L), "Dark 2")[seq_len(n)])
  }
  if (length(col) != 1L && length(col) != n) {
    stop("`col` must hold one colour, or one per drawn model (", n, ").",
         call. = FALSE)
  }
  rep(col, length.out = n)
}


#' Print one indented entry with its continuation hanging under the value
#'
#' `sa_cat_field()` pads its label to a fixed nine characters, which is what
#' lines the top-level fields of a result up with each other. A model name is
#' not that: it is as wide as the widest name in this evaluation and no wider,
#' and it sits one level further in.
#'
#' @keywords internal
#' @noRd
sa_cat_entry <- function(label, text, width = 78, indent = 4L) {
  head_text <- paste0(strrep(" ", indent), label, "  ")
  lines <- strwrap(text, width = width - nchar(head_text))
  cat(head_text, lines[1], "\n", sep = "")
  for (line in lines[-1]) {
    cat(strrep(" ", nchar(head_text)), line, "\n", sep = "")
  }
  invisible(NULL)
}


#' Format an estimate with the interval and p-value that belong to it
#'
#' The name of the quantity is the caller's, since it is what the line is
#' aligned on, so this begins at the equals sign.
#'
#' @keywords internal
#' @noRd
sa_fmt_inference <- function(estimate, lower, upper, pval) {
  paste0("= ", sa_fmt_num(estimate, 3),
         "  [", sa_fmt_num(lower, 3), ", ", sa_fmt_num(upper, 3), "]",
         "  p = ", sa_fmt_num(pval, 3))
}


#' Print an evaluation result
#'
#' Summarises what was scored on what, and how each model did, rather than
#' printing every table. The per-model scores are in `x$metrics`, what each
#' model did against the baseline is in `x$comparisons`, and the predictions
#' themselves are in `x$predictions`.
#'
#' @param x An evaluation, as returned by [evaluate_regression_models()] or
#'   [evaluate_classification_models()].
#' @param n Maximum number of models to show. The rest are counted.
#' @param ... Ignored, present for consistency with [print()].
#'
#' @return `x` invisibly.
#'
#' @examples
#' ## Two models of the same outcome, scored on rows neither was fitted on.
#' train <- mtcars[1:24, ]
#' test <- mtcars[25:32, ]
#' full <- fit_linear_regression(train, outcome = "mpg",
#'                               predictors = c("wt", "hp", "disp"),
#'                               cv = FALSE)
#' small <- fit_linear_regression(train, outcome = "mpg",
#'                                predictors = "wt", cv = FALSE)
#' evaluate_regression_models(full, list(weight_only = small), newdata = test)
#'
#' @export
print.sa_performance <- function(x, n = 10L, ...) {
  n <- sa_check_count(n, "n", 0)
  design <- x$design
  classification <- identical(x$analysis, "classification_performance")

  cat("<sa_performance> ", x$analysis, "\n", sep = "")
  cat("  outcome  : ", design$outcome, "  (", design$outcome_type, ")\n",
      sep = "")
  if (!is.null(design$outcome_lv)) {
    cat("             scoring the probability of ", design$outcome_lv[2],
        " against ", design$outcome_lv[1], ", ", design$n_events, " of ",
        design$n_used, " row(s)\n", sep = "")
  }
  cat("  rows     : ", design$n_used, " scored",
      if (design$n_dropped > 0L) {
        paste0("  (", design$n_dropped, " row(s) dropped)")
      },
      "\n", sep = "")
  cat("  models   : ", length(x$models), ", baseline = ", design$baseline,
      "\n", sep = "")
  if (classification) {
    cat("  threshold: ", x$parameters$threshold,
        "  (accuracy, sensitivity and specificity only)\n", sep = "")
  }

  shown <- utils::head(x$metrics, n)
  width <- if (nrow(shown) > 0L) max(nchar(shown$model)) else 0L
  cat("\n  metrics\n")
  for (i in seq_len(nrow(shown))) {
    row <- shown[i, ]
    text <- if (classification) {
      paste0("auc = ", sa_fmt_num(row$auc, 3),
             "  [", sa_fmt_num(row$auc_lower_conf, 3), ", ",
             sa_fmt_num(row$auc_upper_conf, 3), "]",
             ", brier = ", sa_fmt_num(row$brier, 3),
             ", accuracy = ", sa_fmt_num(row$accuracy, 3))
    } else {
      paste0("cor = ", sa_fmt_num(row$cor, 3),
             ", r_squared = ", sa_fmt_num(row$r_squared, 3),
             ", rmse = ", sa_fmt_num(row$rmse, 3),
             ", mae = ", sa_fmt_num(row$mae, 3))
    }
    sa_cat_entry(formatC(row$model, width = -width), text)
  }
  if (nrow(x$metrics) > nrow(shown)) {
    cat("    ... and ", nrow(x$metrics) - nrow(shown),
        " more model(s) in $metrics\n", sep = "")
  }

  if (!is.null(x$comparisons)) {
    cat("\n  comparisons  (against ", design$baseline, ")\n", sep = "")
    shown <- utils::head(x$comparisons, n)
    width <- if (nrow(shown) > 0L) max(nchar(shown$model)) else 0L
    for (i in seq_len(nrow(shown))) {
      row <- shown[i, ]
      if (classification) {
        # Three estimates with an interval and a p-value apiece do not fit on a
        # line beside a model name, and folding them onto one would put the
        # name of the model further from the numbers than the numbers are from
        # each other.
        cat("    ", row$model, "\n", sep = "")
        sa_cat_entry(formatC("delta_auc", width = -9),
                     sa_fmt_inference(row$delta_auc,
                                      row$delta_auc_lower_conf,
                                      row$delta_auc_upper_conf,
                                      row$delta_auc_pval),
                     indent = 6L)
        sa_cat_entry(formatC("IDI", width = -9),
                     sa_fmt_inference(row$idi, row$idi_lower_conf,
                                      row$idi_upper_conf, row$idi_pval),
                     indent = 6L)
        sa_cat_entry(formatC("NRI", width = -9),
                     sa_fmt_inference(row$nri, row$nri_lower_conf,
                                      row$nri_upper_conf, row$nri_pval),
                     indent = 6L)
      } else {
        sa_cat_entry(
          formatC(row$model, width = -width),
          paste0("delta_cor = ", sa_fmt_num(row$delta_cor, 3),
                 ", delta_r_squared = ", sa_fmt_num(row$delta_r_squared, 3),
                 ", delta_rmse = ", sa_fmt_num(row$delta_rmse, 3),
                 ", delta_mae = ", sa_fmt_num(row$delta_mae, 3))
        )
      }
    }
    if (nrow(x$comparisons) > nrow(shown)) {
      cat("    ... and ", nrow(x$comparisons) - nrow(shown),
          " more model(s) in $comparisons\n", sep = "")
    }
  }

  invisible(x)
}


#' Draw an evaluation result
#'
#' Dispatches on what was evaluated, since the two scenarios have different
#' pictures rather than one picture with a switch: a regression is drawn against
#' the outcome it predicted and a classification against the two classes it
#' ranked. Call [draw_prediction_plot()] or [draw_roc_curve()] directly to reach
#' their own arguments by name.
#'
#' @param x An evaluation, as returned by [evaluate_regression_models()] or
#'   [evaluate_classification_models()].
#' @param ... Passed to whichever of the two is called.
#'
#' @return Whatever the function called returns, invisibly.
#'
#' @seealso [draw_prediction_plot()] and [draw_roc_curve()].
#'
#' @examples
#' train <- mtcars[1:24, ]
#' test <- mtcars[25:32, ]
#' fit <- fit_linear_regression(train, outcome = "mpg",
#'                              predictors = c("wt", "hp"), cv = FALSE)
#' plot(evaluate_regression_models(fit, newdata = test))
#'
#' @export
plot.sa_performance <- function(x, ...) {
  if (identical(x$analysis, "classification_performance")) {
    draw_roc_curve(x, ...)
  } else {
    draw_prediction_plot(x, ...)
  }
}
