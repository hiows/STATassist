# The result contract shared by every model fitting function, and the counterpart
# of `result.R` for the machine learning family. A comparison result is organised
# around a feature axis: every table repeats `features` in the same order. A model
# has no feature axis. It has one outcome and a set of terms, and the terms are
# not the columns that were passed in, since a factor predictor becomes several of
# them. `terms` is therefore what `features` is over there: the row order every
# table in the object follows.
#
# One slot breaks the rule the comparison contract keeps. `fit` holds the
# `caret::train` object itself, because a model that cannot be handed to
# `predict()` is not much of a model, and nothing else in R can stand in for it.
# It carries an `sa_fit` class in front of `caret`'s own so that `coef()` and
# `summary()` reach the fit inside it, which is what the two methods at the foot
# of this file do. `coef()` on the result itself is a different question with a
# different answer: the table, since that is what this object holds and the named
# vector is one call away in `$fit`. Everything else is a scalar, a character
# vector, a named list or a data.frame, so dropping that one slot leaves an object
# that still writes out as JSON.

#' Column names every coefficient table must carry
#'
#' What every model answers is which terms it has and what it estimated for each,
#' so those two are the contract and the rest is per-model. `terms` is also the row
#' order every other table in the object follows.
#'
#' @keywords internal
#' @noRd
sa_model_coef_columns <- function() {
  c("terms", "estimate")
}


#' The inference columns, which a model either fills or does not have
#'
#' `statistic` is a t value for a linear model and a Wald z for a logistic one,
#' which is why the column is not named after either. A penalized fit has neither:
#' its estimates are deliberately biased and the standard error assumes an
#' unbiased one, so there is nothing to put in these columns and they are absent
#' from its table rather than present and `NA`.
#'
#' Absent rather than `NA` because the two say different things. An `NA` in a
#' table that has the column is a question that was asked and came back
#' unanswerable, which is what an aliased term in a linear model is. Six columns
#' of `NA` down their whole length read as a table that lost its values, not as a
#' model that never had them, and `is.null(fit$coefficients$pval)` tells a
#' consumer which kind of table it is holding in one line.
#'
#' They come as a group. A table with a `pval` and no `stderr` would be a third
#' kind of table that nothing in the package produces on purpose, so
#' `sa_new_model()` refuses it.
#'
#' @keywords internal
#' @noRd
sa_model_inference_columns <- function() {
  c("stderr", "statistic", "df", "pval", "lower_conf", "upper_conf")
}


#' Assemble a fitted model result object
#'
#' The checks here guard the contract rather than the user's input, so they fire
#' only on a mistake inside the package and say so.
#'
#' @param analysis Model identifier, such as `"linear_regression"`.
#' @param terms Coefficient term names, in the row order every table uses.
#' @param design Named list describing the data the model saw: the outcome, its
#'   type and levels, the predictors, and how many rows were usable.
#' @param parameters Named list of the fitting choices, with the resampling
#'   arguments as they were actually used rather than as they were passed.
#' @param coefficients data.frame of one row per term.
#' @param fit_stats Named list of goodness-of-fit scalars for the model as a
#'   whole, which are per-model rather than per-term and so do not fit the table.
#' @param performance data.frame of the resampled performance, one row per
#'   hyperparameter combination, or `NULL` when nothing was resampled.
#' @param resampling data.frame of one row per resample, or `NULL`.
#' @param engine Named list naming what actually fitted the model.
#' @param fit The engine object, kept so that the model can predict.
#'
#' @keywords internal
#' @noRd
sa_new_model <- function(analysis,
                         terms,
                         design,
                         parameters,
                         coefficients,
                         fit_stats,
                         performance = NULL,
                         resampling = NULL,
                         engine,
                         fit) {

  if (!is.character(terms) || length(terms) == 0L) {
    stop("internal error: `terms` must be a non-empty character vector.",
         call. = FALSE)
  }
  if (!is.data.frame(coefficients)) {
    stop("internal error: `coefficients` must be a data.frame.", call. = FALSE)
  }
  if (!identical(coefficients$terms, terms)) {
    stop("internal error: `coefficients` is not aligned with `terms`.",
         call. = FALSE)
  }
  absent <- setdiff(sa_model_coef_columns(), names(coefficients))
  if (length(absent) > 0L) {
    stop("internal error: `coefficients` is missing contract column(s): ",
         paste(absent, collapse = ", "), ".", call. = FALSE)
  }
  inference <- sa_model_inference_columns()
  present <- intersect(inference, names(coefficients))
  if (length(present) > 0L && length(present) < length(inference)) {
    stop("internal error: `coefficients` carries some inference column(s) and ",
         "not others, so it is neither a table with inference nor one without: ",
         paste(present, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.list(fit_stats) || is.null(names(fit_stats))) {
    stop("internal error: `fit_stats` must be a named list.", call. = FALSE)
  }
  for (nm in c("package", "method", "label", "metrics")) {
    if (is.null(engine[[nm]])) {
      stop("internal error: `engine` is missing `", nm, "`.", call. = FALSE)
    }
  }
  if (!is.null(performance) && !is.data.frame(performance)) {
    stop("internal error: `performance` must be a data.frame or NULL.",
         call. = FALSE)
  }
  if (!is.null(resampling) && !is.data.frame(resampling)) {
    stop("internal error: `resampling` must be a data.frame or NULL.",
         call. = FALSE)
  }

  structure(
    list(
      analysis     = analysis,
      terms        = terms,
      design       = design,
      parameters   = parameters,
      coefficients = coefficients,
      fit_stats    = fit_stats,
      performance  = performance,
      resampling   = resampling,
      engine       = engine,
      fit          = fit,
      metadata     = sa_metadata()
    ),
    class = c("sa_model", "sa_result")
  )
}


#' Format a number for printing without deciding how big it is
#'
#' A coefficient can be a rate per unit of body weight or a log odds per
#' millimetre, so a fixed number of decimals is either noise or nothing at all.
#' Significant digits keep both readable. Each value is formatted on its own and
#' trimmed, since the padding a vector would share is the caller's business: only
#' the coefficient table lines up in columns, and it says so with its own widths.
#'
#' @keywords internal
#' @noRd
sa_fmt_num <- function(x, digits = 4) {
  vapply(x, function(v) {
    if (is.null(v) || length(v) == 0L || is.na(v)) {
      "NA"
    } else {
      format(signif(v, digits), trim = TRUE)
    }
  }, character(1), USE.NAMES = FALSE)
}


#' Print a labelled field, wrapping the continuation under the value
#'
#' The goodness-of-fit line holds ten scalars and the resampled line holds three
#' with their standard deviations, so both outrun one console line. Wrapping the
#' label together with the value would collapse the padding that lines the labels
#' up, since [strwrap()] normalises whitespace, so the two are wrapped apart.
#'
#' @keywords internal
#' @noRd
sa_cat_field <- function(label, text, width = 78) {
  head_text <- paste0("  ", formatC(label, width = -9), ": ")
  lines <- strwrap(text, width = width - nchar(head_text))
  cat(head_text, lines[1], "\n", sep = "")
  for (line in lines[-1]) {
    cat(strrep(" ", nchar(head_text)), line, "\n", sep = "")
  }
  invisible(NULL)
}


#' The row of the performance table the model was actually fitted at
#'
#' `performance` holds one row per hyperparameter combination, in the order they
#' were scored rather than in the order they placed, so the chosen combination is
#' the one `parameters` names and not the first row. An unpenalized model has a
#' single combination and a single row, which is why this used to be that first
#' row: the two agreed as long as nothing was tuned.
#'
#' The chosen values are read from the result's own `parameters` rather than from
#' `caret`'s `bestTune`, so that this asks nothing of `$fit`. That slot is the one
#' part of the object that cannot be serialised, and a result printed without it
#' would otherwise fall through to the first row while the line above it still
#' named the chosen penalty.
#'
#' Falling back to the first row rather than failing is deliberate. This is only
#' ever used to summarise for printing, and a table that cannot be matched is not
#' a reason to refuse to print the rest of the object.
#'
#' @param performance The `performance` table, known to have at least one row.
#' @param chosen The result's `parameters`, matched to `performance` by the names
#'   the two share, which are the tuned ones. A model that tuned nothing shares
#'   no name and has one row anyway.
#'
#' @return One row of `performance`.
#'
#' @keywords internal
#' @noRd
sa_performance_row <- function(performance, chosen) {
  if (is.null(chosen) || nrow(performance) == 1L) {
    return(performance[1L, , drop = FALSE])
  }
  at <- rep(TRUE, nrow(performance))
  for (nm in intersect(names(chosen), names(performance))) {
    at <- at & performance[[nm]] == chosen[[nm]]
  }
  performance[if (any(at)) which(at)[1L] else 1L, , drop = FALSE]
}


#' Print a fitted model
#'
#' Summarises what was fitted to what, and how it did, rather than printing every
#' table. The coefficient table is in `x$coefficients`, which is also what
#' [coef()] answers with; the resampled folds are in `x$resampling`, and the
#' engine object `predict()` and [summary()] take is `x$fit`.
#'
#' @param x A fitted model, as returned by [fit_linear_regression()],
#'   [fit_logistic_regression()], [fit_elastic_net()], [fit_rf()] or
#'   [fit_svm()].
#' @param n Maximum number of coefficient rows to show. The rest are counted.
#' @param ... Ignored, present for consistency with [print()].
#'
#' @return `x` invisibly.
#'
#' @examples
#' fit_linear_regression(mtcars, outcome = "mpg",
#'                       predictors = c("wt", "hp"), cv = FALSE)
#'
#' @export
print.sa_model <- function(x, n = 10L, ...) {
  n <- sa_check_count(n, "n", 0)
  design <- x$design
  params <- x$parameters

  cat("<sa_model> ", x$analysis, "\n", sep = "")
  cat("  outcome  : ", design$outcome, "  (", design$outcome_type, ")\n",
      sep = "")
  if (!is.null(design$outcome_lv)) {
    # Only a model with an odds ratio in its table is modelling odds. A forest
    # models the same class in the same direction and reports no odds at all, so
    # the line names the class without naming a scale it never used.
    cat("             modelling ",
        if (!is.null(x$coefficients$odds_ratio)) "the odds of " else "",
        design$outcome_lv[2], " against ", design$outcome_lv[1], ", ",
        design$n_events, " of ", design$n_used, " row(s)\n", sep = "")
  }
  cat("  rows     : ", design$n_used, " used",
      if (design$n_dropped > 0L) {
        paste0("  (", design$n_dropped, " incomplete row(s) dropped)")
      },
      "\n", sep = "")
  cat("  terms    : ", length(x$terms), " over ", length(design$predictors),
      " predictor(s)\n", sep = "")
  cat("  settings : ",
      if (isTRUE(params$cv)) {
        paste0(params$cv_method,
               if (!is.na(params$n_fold)) {
                 paste0(", ", params$n_fold, " fold(s)")
               },
               if (!is.na(params$n_repeat)) {
                 paste0(" x ", params$n_repeat, " repeat(s)")
               })
      } else {
        "no resampling"
      },
      # A model that reports no interval has no confidence level to report
      # either, so the label is left out rather than printed against nothing.
      if (!is.null(params$conf_level)) {
        paste0(", conf_level = ", params$conf_level)
      },
      "\n", sep = "")
  if (!is.null(params$penalty)) {
    cat("  penalty  : ", params$penalty,
        ", alpha = ", sa_fmt_num(params$alpha, 3),
        ", lambda = ", sa_fmt_num(params$lambda, 3),
        if (isTRUE(params$cv)) {
          paste0("  (chosen from ", params$n_candidates, " candidate(s))")
        },
        "\n", sep = "")
  }
  if (!is.null(params$ntree)) {
    # Only a grid of more than one `mtry` was chosen between. A forest's default
    # `mtry` is a rule of thumb rather than a search, so the usual line would
    # report one candidate as having won something.
    cat("  forest   : ", params$ntree, " tree(s), mtry = ", params$mtry,
        ", nodesize = ", params$nodesize,
        if (isTRUE(params$cv) && params$n_candidates > 1L) {
          paste0("  (mtry chosen from ", params$n_candidates, " candidate(s))")
        },
        "\n", sep = "")
  }
  if (!is.null(params$kernel)) {
    cat("  kernel   : ", params$kernel,
        ", C = ", sa_fmt_num(params$C, 3),
        ", sigma = ", sa_fmt_num(params$sigma, 3),
        if (isTRUE(params$cv) && params$n_candidates > 1L) {
          paste0("  (chosen from ", params$n_candidates, " candidate(s))")
        },
        "\n", sep = "")
  }

  # A penalized fit has no standard error and so no interval and no p-value, and
  # its table does not carry the columns at all. What it answers instead is
  # whether the penalty kept the term.
  inference <- !is.null(x$coefficients$pval)
  # A forest and a machine have no coefficients at all. What they answer with is
  # how much each term was worth to them, and the heading says so rather than
  # calling those numbers something they are not. Neither column that would stand
  # in for an estimate is there, which is what tells such a table from a penalized
  # one: `selected` is a statement about a coefficient that a permuted loss has no
  # counterpart for. The word in brackets names which measure the printed column
  # is, which matters for a forest, whose table holds an impurity-based one too.
  importance <- !inference && is.null(x$coefficients$selected)
  cat("\n  ",
      if (importance) "importance  (permutation)" else "coefficients",
      "\n", sep = "")
  shown <- utils::head(x$coefficients, n)
  width <- if (nrow(shown) > 0L) max(nchar(shown$terms)) else 0L
  for (i in seq_len(nrow(shown))) {
    row <- shown[i, ]
    cat("    ", formatC(row$terms, width = -width), "  ",
        formatC(sa_fmt_num(row$estimate), width = 10),
        if (inference) {
          paste0("  [", sa_fmt_num(row$lower_conf, 3), ", ",
                 sa_fmt_num(row$upper_conf, 3), "]  p = ",
                 sa_fmt_num(row$pval, 3))
        } else if (!is.null(row$selected)) {
          paste0("  ", if (isTRUE(row$selected)) "selected" else "dropped")
        },
        "\n", sep = "")
  }
  if (nrow(x$coefficients) > nrow(shown)) {
    cat("    ... and ", nrow(x$coefficients) - nrow(shown),
        " more term(s) in $coefficients\n", sep = "")
  }

  cat("\n")
  sa_cat_field("fit", paste0(names(x$fit_stats), " = ",
                             sa_fmt_num(unlist(x$fit_stats), 3),
                             collapse = ", "))

  if (!is.null(x$performance)) {
    best <- sa_performance_row(x$performance, params)
    scored <- vapply(x$engine$metrics, function(m) {
      sd_col <- paste0(m, "SD")
      paste0(m, " = ", sa_fmt_num(best[[m]], 3),
             if (!is.null(best[[sd_col]])) {
               paste0(" (SD ", sa_fmt_num(best[[sd_col]], 2), ")")
             })
    }, character(1))
    sa_cat_field("resample",
                 paste0(paste(scored, collapse = ", "),
                        if (!is.null(x$resampling)) {
                          paste0(" over ", nrow(x$resampling), " resample(s)")
                        }))
  }
  if (length(design$dropped_predictors) > 0L) {
    cat("  dropped  : ", paste(design$dropped_predictors, collapse = ", "),
        " (single valued)\n", sep = "")
  }

  invisible(x)
}


#' Coefficients of a fitted model
#'
#' The coefficient table, `object$coefficients` itself: one row per term, in the
#' order of `object$terms`, carrying everything the model estimated about each of
#' them rather than the estimate alone.
#'
#' The `lm`-shaped answer, a named numeric vector, is `coef(x$fit)`. Two objects
#' are reached by two calls, so they answer in their own terms rather than both
#' giving the same thing: the engine object answers as the engine does, and the
#' result object answers with the table it was assembled to hold. Ask `$fit` when
#' a vector is what is wanted to index or multiply.
#'
#' What the table holds depends on the model, and the columns say which kind it
#' is. An unpenalized fit carries the standard error, the statistic, the p-value
#' and the confidence limits, and a classification the odds ratio; a penalized one
#' has no standard error to report and carries `selected` instead, so
#' `is.null(coef(x)$pval)` is the test for a fit that cannot be asked for
#' inference. A forest and a support vector machine have no coefficient of any
#' kind, and their table is the importance of each term rather than an effect per
#' unit of it; see [fit_rf()] and [fit_svm()]. Every term keeps its row either way — a term a penalty dropped with
#' an `estimate` of exactly 0, one the engine could not estimate with `NA` — since
#' a table shorter than the model would be a table that had lost terms rather than
#' one describing terms that went to zero.
#'
#' @param object A fitted model, as returned by [fit_linear_regression()],
#'   [fit_logistic_regression()], [fit_elastic_net()], [fit_rf()] or
#'   [fit_svm()].
#' @param ... Ignored, present for consistency with [coef()]. `stats::coef()`'s
#'   own `complete` argument is among what is ignored: it drops the `NA` cells of
#'   a vector, which on a table would mean indexing a data.frame with a logical
#'   matrix, and the rows are what say which terms the model has.
#'
#' @return The `coefficients` element of `object`, a data.frame of one row per
#'   term.
#'
#' @seealso [fit_linear_regression()], [fit_elastic_net()], [fit_rf()],
#'   [fit_svm()], and [coef.sa_fit()] for the same question asked of the engine
#'   object in `$fit`.
#'
#' @examples
#' fit <- fit_linear_regression(mtcars, outcome = "mpg",
#'                              predictors = c("wt", "hp"), cv = FALSE)
#' coef(fit)
#'
#' ## The named vector is one object further in.
#' coef(fit$fit)
#'
#' ## A penalized fit answers with what it has: no inference, but `selected`.
#' pen <- fit_elastic_net(mtcars, outcome = "mpg", penalty = "lasso",
#'                        lambda = 0.5, cv = FALSE)
#' coef(pen)
#'
#' @export
coef.sa_model <- function(object, ...) {
  object$coefficients
}


#' Predict from a fitted model on rows it was not fitted to
#'
#' Takes the data frame the fit took — the test half of a [split_data()] result,
#' say — and answers one prediction per row of it. The columns are read by name
#' and coded the way the fit coded them, so the rows to predict can be handed over
#' exactly as they came, outcome column and all.
#'
#' This is the method to use rather than `predict(object$fit, newdata = )`. The
#' engine object knows the names of the columns it was given and nothing about
#' where they came from, which is enough for [fit_linear_regression()],
#' [fit_logistic_regression()] and [fit_rf()], whose engine was handed the
#' predictor frame itself, and not enough for [fit_elastic_net()] or [fit_svm()],
#' whose engines were handed a design matrix. `glmnet` and `kernlab` read that
#' matrix by position, and `caret` prepares `newdata` for it by keeping the columns
#' whose names it recognises, in the order `newdata` happens to hold them. A factor
#' predictor is therefore dropped whole, since `x_cat` is not one of the names —
#' the model has `x_catmid` and `x_cathigh` — and a set of numeric predictors in
#' another order is silently matched to the wrong coefficients, or measured along
#' the wrong axis of the kernel. Here the terms are rebuilt from the levels the fit
#' recorded and put in the model's own order by name, which is what makes one call
#' work for every model in the family.
#'
#' @details
#' Columns the model never saw are ignored, so the outcome column and anything
#' else the frame carries can stay. A predictor that is absent is an error naming
#' it, and so is a factor level the fit never saw, since neither has a coefficient
#' to be predicted with.
#'
#' There is one prediction per row of `newdata` whatever the row holds, and a row
#' with a missing value among the predictors gets `NA`. That is the rule the fit
#' already follows in reverse: those are the rows `design$n_dropped` counted, and
#' they cannot be predicted for the same reason they could not be fitted. Saying
#' so with an `NA` in place keeps the answer aligned with `newdata`, which a
#' shorter vector would not be. It also has to be said explicitly, because a
#' penalized fit would otherwise predict some of them: a coefficient of exactly
#' zero drops out of the sparse product before the missing value it multiplies is
#' ever read, so whether a row could be predicted would depend on which predictor
#' the hole fell in.
#'
#' @param object A fitted model, as returned by [fit_linear_regression()],
#'   [fit_logistic_regression()], [fit_elastic_net()], [fit_rf()] or
#'   [fit_svm()].
#' @param newdata Rows to predict, a data.frame or matrix carrying the predictor
#'   columns, or `NULL` for the rows the model was fitted on.
#' @param type `"raw"` for the prediction itself, a fitted value for a regression
#'   and a class label for a classification; `"response"` for the prediction on
#'   the scale of the outcome, which for a classification is the probability of
#'   `design$outcome_lv[2]`, the class the coefficients describe; or `"prob"` for
#'   one column per class.
#' @param ... Passed on to [caret::predict.train()].
#'
#' @return One prediction per row of `newdata`: an unnamed numeric vector for a
#'   regression and for `type = "response"`, a factor at the levels of
#'   `design$outcome_lv` for `type = "raw"` on a classification, and a data.frame
#'   of one column per class for `type = "prob"`. Rows that are not complete
#'   across the predictors are `NA` throughout.
#'
#' @seealso [fit_elastic_net()], [fit_rf()], [fit_svm()], [split_data()], which
#'   draws the rows to predict, and [predict.sa_fit()] for the engine object in
#'   `$fit`, which is `caret`'s own method and takes what `caret` prepared.
#'
#' @examples
#' sp <- split_data(mtcars, seed = 1)
#' train <- sp$datasets[[1]]$train_data
#' test <- sp$datasets[[1]]$test_data
#'
#' fit <- fit_linear_regression(train, outcome = "mpg",
#'                              predictors = c("wt", "hp"), cv = FALSE)
#' sqrt(mean((test$mpg - predict(fit, newdata = test))^2))
#'
#' ## A factor predictor is what `predict(fit$fit, )` cannot be given on a
#' ## penalized fit: the model has one term per level, and the frame has the
#' ## column the levels came from.
#' train$cyl <- factor(train$cyl)
#' test$cyl <- factor(test$cyl, levels = levels(train$cyl))
#' pen <- fit_elastic_net(train, outcome = "mpg",
#'                        predictors = c("wt", "hp", "cyl"),
#'                        penalty = "lasso", lambda = 0.5, cv = FALSE)
#' predict(pen, newdata = test)
#'
#' ## On a classification, `type = "response"` is the probability of the second
#' ## level of `outcome_lv`, the one the coefficients describe.
#' iris2 <- iris[iris$Species != "setosa", ]
#' clf <- fit_logistic_regression(iris2, outcome = "Species",
#'                                predictors = "Petal.Length",
#'                                outcome_lv = c("versicolor", "virginica"),
#'                                cv = FALSE)
#' head(predict(clf, newdata = iris2, type = "response"))
#'
#' @export
predict.sa_model <- function(object,
                             newdata = NULL,
                             type = c("raw", "response", "prob"),
                             ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    # No rows to code: `caret` predicts the ones it kept, which are the
    # `design$n_used` rows the coefficients were estimated from.
    return(sa_predict_shape(stats::predict(object$fit, type = type, ...)))
  }

  x <- sa_predict_frame(newdata, object$design)
  usable <- stats::complete.cases(x)
  if (!any(usable)) {
    stop("no row of `newdata` is complete across the predictor(s) the model ",
         "was fitted on, so there is nothing to predict from.", call. = FALSE)
  }

  # The incomplete rows are held back rather than passed in and patched up
  # afterwards, so that the engine is asked only what it can answer and the
  # answer is put back where the row was.
  # `engine$x_names` is there exactly when the engine was handed a design matrix
  # rather than the predictor frame, and holds the columns it was given in the
  # order it was given them. A model fitted from the frame needs none of this:
  # `caret` reads its columns by name and did the coding itself.
  ready <- x[usable, , drop = FALSE]
  if (!is.null(object$engine$x_names)) {
    ready <- sa_design_matrix(ready,
                              xlev = object$design$predictor_lv,
                              want = object$engine$x_names)
  }

  out <- sa_predict_shape(stats::predict(object$fit, newdata = ready,
                                         type = type, ...))
  sa_predict_scatter(out, usable)
}


#' Give a prediction the shape the model contract promises
#'
#' The engines do not agree on it. `lm()` names its fitted values after the row
#' names of `newdata` and `glmnet` does not name them at all, and a name that
#' comes and goes with the engine is a name a caller cannot rely on. The row a
#' prediction belongs to is its position, as it is everywhere else in this
#' package, so the numbers come back unnamed and a `type = "prob"` table comes
#' back numbered from one.
#'
#' @param value Whatever the engine's `predict()` returned.
#'
#' @return `value`, unnamed.
#'
#' @keywords internal
#' @noRd
sa_predict_shape <- function(value) {
  if (is.data.frame(value)) {
    rownames(value) <- NULL
    return(value)
  }
  if (is.matrix(value) && ncol(value) == 1L) {
    value <- as.numeric(value)
  }
  unname(value)
}


#' Put predictions back beside the rows they were asked about
#'
#' The engine was given the complete rows only, so what comes back is shorter than
#' `newdata` whenever a row had a missing predictor. It is scattered back to full
#' length here, `NA` where the row could not be read.
#'
#' @param value Predictions for the usable rows, a vector, factor or data.frame.
#' @param usable Logical vector over the rows of `newdata`.
#'
#' @return `value` at the length of `usable`.
#'
#' @keywords internal
#' @noRd
sa_predict_scatter <- function(value, usable) {
  if (all(usable)) {
    return(value)
  }
  # Indexing by `NA` is what makes the empty rows: it produces a missing value of
  # the right type, which for a factor keeps the levels and for a table keeps the
  # columns.
  blank <- rep(NA_integer_, length(usable))
  if (is.data.frame(value)) {
    full <- value[blank, , drop = FALSE]
    full[usable, ] <- value
    rownames(full) <- NULL
    return(full)
  }
  full <- value[blank]
  full[usable] <- value
  full
}


# The engine object in `$fit` is a `caret::train` object wearing an `sa_fit`
# class in front, and these are the three methods that class exists for. `caret`
# defines `predict()`, `fitted()`, `residuals()` and `print()` for its own class
# and those keep working by inheritance; what it does not define is `coef()`, so
# `coef()` on a fitted model returns NULL, which is the one answer a fitted model
# should never give. `coef()` and `summary()` read through to `finalModel`, the
# `lm`, `glm` or `glmnet` underneath.
#
# A penalized fit reads one step further. `caret` fits the whole `lambda` path and
# records the chosen value as `lambdaOpt`, so `coef()` on the `glmnet` object
# returns a sparse matrix of every fit on that path. Only one of them is the model
# the result describes, and it is the one `caret` predicts from, so that is the
# one `coef()` answers with, as the same named vector the other two models give.
#
# `predict()` is the odd one out, because `caret` does define it and the method
# below only widens the vocabulary it accepts. `predict.train()` rejects any
# `type` other than "raw" and "prob", so `type = "response"` — the word anyone
# who has used `glm()` reaches for — is an error on a logistic regression. The
# probability it asks for is already what `type = "prob"` computes, so the
# method translates the name and lets `caret` do the predicting.

#' Coefficients, summary and predictions from the model inside a fit
#'
#' The `fit` element of an `sa_model` is the engine object, kept so that the
#' model can [predict()]. These methods make the rest of what anyone does with a
#' fitted model work on it as well: `coef()` and `summary()` read through to the
#' [stats::lm()] or [stats::glm()] fit that [caret::train()] built, so what comes
#' back is what the same call on that model would give.
#'
#' The numbers are also in `$coefficients`, as a data.frame with the confidence
#' interval beside them and, for a logistic regression, the odds ratio, and that
#' table is what [coef.sa_model()] gives for the result as a whole. These methods
#' are for when an `lm`-shaped answer is what is wanted instead: a named vector to
#' index, or the summary layout that is already familiar.
#'
#' On a [fit_elastic_net()] fit, `coef()` is the coefficients at the penalty that
#' was chosen, which is the model the result describes and the one `predict()`
#' predicts from, rather than the whole `lambda` path `glmnet` keeps. `summary()`
#' has nothing of the same shape to return there, since a penalized fit has no
#' standard errors and no residual degrees of freedom to test against; what it
#' does have is in `$coefficients` and `$fit_stats`.
#'
#' A [fit_rf()] or [fit_svm()] fit has neither. `coef()` on one is an error naming
#' the importance table to read instead, rather than the `NULL` `coef()` on a
#' `randomForest` object gives or the S4 indexing error it gives on a `ksvm` one,
#' and `summary()` reaches a model with no summary method of its own.
#' `print(fit$fit)` is the readable answer for both, and everything the result
#' reports is in `$coefficients` and `$fit_stats`.
#'
#' @details
#' `predict()` here takes `newdata` as `caret` prepared it rather than as the fit
#' received it, which for [fit_elastic_net()] and [fit_svm()] is a design matrix
#' and not a data frame of predictors. [predict.sa_model()] is the method that
#' takes the rows
#' themselves, since the result is what knows which columns they were and what
#' levels its factors had.
#'
#' `predict()` is otherwise `caret`'s, with one word added. [caret::train()] knows
#' `type = "raw"` for the prediction itself and `type = "prob"` for the class
#' probabilities, and refuses anything else; `"response"`, which is what
#' [stats::glm()] and [stats::lm()] call the prediction on the scale of the
#' outcome, is therefore an error on an object that is plainly a logistic
#' regression. It is accepted here and means the same thing in both models: the
#' predicted value for a linear regression, and for a logistic one the
#' probability of `outcome_lv[2]`, the class every coefficient and odds ratio in
#' the result already describes. The probability is `caret`'s own, the second
#' column of `type = "prob"`, so nothing is computed twice.
#'
#' @param object The `fit` element of a [fit_linear_regression()],
#'   [fit_logistic_regression()], [fit_elastic_net()], [fit_rf()] or [fit_svm()]
#'   result.
#' @param newdata Rows to predict as `caret` takes them, which for a penalized fit
#'   or a machine is the design matrix rather than the predictor frame, or `NULL`
#'   for the rows the model was fitted on. [predict.sa_model()] takes the frame.
#' @param type `"raw"` for the prediction as `caret` gives it, a fitted value or
#'   a class label; `"response"` for the same thing on the scale of the outcome,
#'   which for a classification is the probability of `outcome_lv[2]`; or
#'   `"prob"` for one column per class.
#' @param ... Passed on to the method of the underlying fit.
#'
#' @return From `coef()`, a named numeric vector of coefficients, `NA` for a term
#'   the fit could not estimate and exactly 0 for one a penalty dropped. From
#'   `summary()`, the `summary.lm` or
#'   `summary.glm` object of the underlying fit. The `Call` it prints belongs to
#'   `caret` rather than to you: `.outcome ~ .` on the predictor frame for a
#'   linear model, and empty for a logistic one, which `caret` fits without a
#'   formula. What was fitted to what is what `print()` on the result says.
#'
#'   From `predict()`, whatever [caret::predict.train()] returns for that `type`,
#'   except that `type = "response"` on a classification gives an unnamed numeric
#'   vector of probabilities, one per row of `newdata`.
#'
#' @examples
#' fit <- fit_linear_regression(mtcars, outcome = "mpg",
#'                              predictors = c("wt", "hp"), cv = FALSE)
#' coef(fit$fit)
#' summary(fit$fit)
#'
#' ## The same numbers the result reports in its own table.
#' fit$coefficients
#'
#' ## On a classification, `type = "response"` is the probability of the second
#' ## level of `outcome_lv`, the one the coefficients describe.
#' iris2 <- iris[iris$Species != "setosa", ]
#' clf <- fit_logistic_regression(iris2, outcome = "Species",
#'                                predictors = "Petal.Length",
#'                                outcome_lv = c("versicolor", "virginica"),
#'                                cv = FALSE)
#' head(predict(clf$fit, newdata = iris2, type = "response"))
#'
#' ## A penalized fit answers with the coefficients at the chosen penalty, the
#' ## ones `$coefficients` reports, rather than with the whole lambda path.
#' pen <- fit_elastic_net(mtcars, outcome = "mpg", penalty = "lasso",
#'                        lambda = 0.5, cv = FALSE)
#' coef(pen$fit)
#'
#' @rdname sa_fit-methods
#' @export
coef.sa_fit <- function(object, ...) {
  model <- object$finalModel
  # `lambdaOpt` is `caret`'s record of the penalty it chose, so only a penalized
  # fit has one. Reading it has to be safe on a fit that is an S4 object, which is
  # what a machine's is: `$` on one is an error rather than the `NULL` it answers
  # for a list that has no such element.
  chosen <- if (isS4(model)) NULL else model$lambdaOpt
  if (is.null(chosen)) {
    # A model that estimates no coefficient says so in more than one way, and the
    # two mean the same thing: `coef.default()` answers a forest with NULL, which
    # is the one answer this method exists to prevent, and cannot index the S4
    # object a machine is at all. Both become the message below, which names what
    # the model holds instead. Nothing else here can raise: `lm` and `glm` have
    # the coefficients this reads.
    estimate <- tryCatch(stats::coef(model, ...), error = function(e) NULL)
    if (is.null(estimate)) {
      stop("a `", class(model)[1L], "` fit estimates no coefficients, so there ",
           "is no vector to return. What it answers with instead is the ",
           "`coefficients` table of the result: coef() on the `sa_model` rather ",
           "than on its `$fit`.", call. = FALSE)
    }
    return(estimate)
  }
  # One column of the path, and as a named vector rather than the sparse matrix
  # `glmnet` returns, so that a caller can index it the way the other two models
  # let them.
  at_chosen <- as.matrix(stats::coef(model, s = chosen, ...))
  stats::setNames(as.numeric(at_chosen), rownames(at_chosen))
}


#' @rdname sa_fit-methods
#' @export
summary.sa_fit <- function(object, ...) {
  summary(object$finalModel, ...)
}


#' @rdname sa_fit-methods
#' @export
predict.sa_fit <- function(object, newdata = NULL, type = "raw", ...) {
  if (identical(type, "response")) {
    # A regression predicts on the scale of the outcome already, so the two
    # words name the same thing and only one of them has to reach `caret`.
    if (!identical(object$modelType, "Classification")) {
      type <- "raw"
      return(NextMethod())
    }
    # Reassigning `type` is what `NextMethod()` passes on, so the prediction is
    # still `caret`'s; all that happens here is that a column is picked out of
    # it. Which column is the direction rule the whole result follows.
    if (length(object$levels) != 2L) {
      stop("`type = \"response\"` names the probability of one class, which ",
           "only a two-class outcome has. Use `type = \"prob\"` for one column ",
           "per class.", call. = FALSE)
    }
    type <- "prob"
    return(NextMethod()[[object$levels[2L]]])
  }
  NextMethod()
}
