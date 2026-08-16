# What `evaluate_regression_models()` and `evaluate_classification_models()`
# share, which is everything up to the point where the scoring begins. Both take
# a baseline and a set of models to hold against it, both read the rows through
# `predict.sa_model()` rather than through any engine object, and both are only
# meaningful if every model was scored on the same rows.
#
# That last one is the reason this file exists rather than the two functions
# each resolving their own input. `predict.sa_model()` returns `NA` for a row
# that is incomplete across *that model's* predictors, so models fitted on
# different predictor sets come back with different rows filled in. Scoring each
# model on whatever it happened to manage would put two AUCs on two samples and
# call their difference an improvement, and DeLong's test, IDI and NRI are all
# paired statistics that have no meaning at all across different rows. The
# intersection is taken once, here, for both.


#' Resolve the baseline and the models held against it into one named list
#'
#' The names are not decoration: they are what every table is keyed on and what
#' the legend of a plot reads, so an unnamed list is refused rather than given
#' positions for names. The baseline comes first, which is the order `models`
#' fixes for the whole result.
#'
#' @param baseline_model The reference model.
#' @param new_models Named list of models to hold against it, or `NULL`.
#' @param baseline_label What to call the baseline in the tables.
#'
#' @return Named list of `sa_model`, the baseline first.
#'
#' @keywords internal
#' @noRd
sa_resolve_models <- function(baseline_model, new_models, baseline_label) {
  if (!inherits(baseline_model, "sa_model")) {
    stop("`baseline_model` must be a fitted model, as returned by ",
         "fit_linear_regression(), fit_logistic_regression(), ",
         "fit_elastic_net(), fit_rf() or fit_svm().", call. = FALSE)
  }
  if (!is.character(baseline_label) || length(baseline_label) != 1L ||
        is.na(baseline_label) || !nzchar(baseline_label)) {
    stop("`baseline_label` must be a single non-empty name.", call. = FALSE)
  }

  # An empty list is a call that named no comparison, which is the same thing
  # `NULL` says and reads more naturally out of a `lapply()` that found nothing.
  if (is.null(new_models) || length(new_models) == 0L) {
    out <- list(baseline_model)
    names(out) <- baseline_label
    return(out)
  }

  if (!is.list(new_models) || inherits(new_models, "sa_model")) {
    stop("`new_models` must be a named list of fitted models, such as ",
         "list(selected = fit_1, penalized = fit_2), or NULL to score the ",
         "baseline on its own.", call. = FALSE)
  }
  labels <- names(new_models)
  if (is.null(labels) || anyNA(labels) || !all(nzchar(labels))) {
    stop("every element of `new_models` must be named: the names are what the ",
         "result tables and the plot legend call the models.", call. = FALSE)
  }
  dup <- unique(labels[duplicated(labels)])
  if (length(dup) > 0L) {
    stop("`new_models` contains duplicated names: ",
         paste(dup, collapse = ", "), ".", call. = FALSE)
  }
  if (baseline_label %in% labels) {
    stop("`new_models` holds a model called `", baseline_label,
         "`, which is what the baseline is called. Rename it, or pass a ",
         "different `baseline_label`.", call. = FALSE)
  }
  not_models <- labels[!vapply(new_models, inherits, logical(1), "sa_model")]
  if (length(not_models) > 0L) {
    stop("every element of `new_models` must be a fitted model. Not a model: ",
         paste(not_models, collapse = ", "), ".", call. = FALSE)
  }

  out <- c(stats::setNames(list(baseline_model), baseline_label), new_models)
  out
}


#' Refuse a model that answers a different kind of question
#'
#' A classification handed to the regression function would be scored by
#' correlating a probability against a class label, which produces a number
#' rather than an error and is the reason this is checked rather than left to
#' fail downstream.
#'
#' @param models Named list of models.
#' @param want `"continuous"` or `"two classes"`.
#' @param other The function that does take the other kind.
#'
#' @keywords internal
#' @noRd
sa_check_model_family <- function(models, want, other) {
  types <- vapply(models, function(m) m$design$outcome_type, character(1))
  wrong <- names(models)[types != want]
  if (length(wrong) > 0L) {
    stop("every model must have been fitted to ", want, " outcome. Not ",
         want, ": ", paste0(wrong, " (", types[wrong], ")", collapse = ", "),
         ". Use ", other, " for those.", call. = FALSE)
  }
  invisible(models)
}


#' Refuse a set of models that are not describing the same question
#'
#' Two models of different outcomes can both be scored, and the scores can be
#' put in one table, and the table means nothing. The same goes for a pair of
#' classifications whose `outcome_lv` point at different classes: both return a
#' probability from `type = "response"` and one of them is the probability of
#' the other class, so every comparison between them is reversed. Neither is
#' quietly repaired, because a re-pointed level order is a different model from
#' the one the caller fitted and printed.
#'
#' @keywords internal
#' @noRd
sa_check_model_agreement <- function(models) {
  outcomes <- vapply(models, function(m) m$design$outcome, character(1))
  if (length(unique(outcomes)) > 1L) {
    stop("every model must have been fitted to the same outcome, since the ",
         "scores are put side by side. Got ",
         paste0(names(models), " = ", outcomes, collapse = ", "), ".",
         call. = FALSE)
  }

  levels_of <- lapply(models, function(m) m$design$outcome_lv)
  named <- !vapply(levels_of, is.null, logical(1))
  if (any(named)) {
    first <- levels_of[named][[1]]
    same <- vapply(levels_of[named], identical, logical(1), first)
    if (!all(same)) {
      disagree <- names(levels_of[named])[!same]
      stop("every model must hold the same `outcome_lv`, in the same order: ",
           "the second level is the class `type = \"response\"` reports the ",
           "probability of, so a model that names them the other way round ",
           "predicts the other class. Expected ",
           paste(first, collapse = ", "), ", but ",
           paste(disagree, collapse = ", "), " disagree(s). Refit with a ",
           "matching `outcome_lv`.", call. = FALSE)
    }
  }

  invisible(models)
}


#' Resolve the observed outcome of the rows being scored
#'
#' `answer = NULL` reads the column the models were fitted to, which is the
#' usual case: the held-out half of a [split_data()] result carries the outcome
#' under the same name as the half the models were fitted on. A model fitted
#' from a vector remembers no name to look up, and says so.
#'
#' @keywords internal
#' @noRd
sa_resolve_answer <- function(answer, newdata, baseline_model) {
  if (!is.null(answer)) {
    return(sa_resolve_row_vector(answer, "answer", newdata, allow_na = TRUE))
  }

  label <- baseline_model$design$outcome
  if (is.na(label) || identical(label, "<vector>")) {
    stop("`answer` is required: `baseline_model` was fitted to an outcome ",
         "passed as a vector, so it remembers no column name to read from ",
         "`newdata`.", call. = FALSE)
  }
  if (!label %in% names(newdata)) {
    stop("`answer` is NULL, so the outcome is read from the `",  label,
         "` column the models were fitted to, and `newdata` has no such ",
         "column. Name the observed values with `answer`.", call. = FALSE)
  }
  list(value = newdata[[label]], label = label)
}


#' Predict every model on the same rows, and say which rows those are
#'
#' The intersection rather than the union, and one message rather than one per
#' model. Which rows a model can predict is a property of its predictors, so a
#' baseline fitted on nine columns and a reduced model fitted on four disagree
#' on any row that is missing one of the extra five. Scoring each on what it
#' managed would compare two models on two samples.
#'
#' @param models Named list of models, the baseline first.
#' @param newdata The rows to score.
#' @param observed The observed outcome, one entry per row of `newdata`.
#'
#' @return A list of `predicted` (a numeric matrix, one column per model in the
#'   order `models` is in, with only the kept rows), `keep` (the row positions
#'   in `newdata` that survived), `n_obs` and `n_dropped`.
#'
#' @keywords internal
#' @noRd
sa_collect_predictions <- function(models, newdata, observed) {
  n_obs <- nrow(newdata)

  columns <- lapply(names(models), function(nm) {
    out <- tryCatch(
      stats::predict(models[[nm]], newdata = newdata, type = "response"),
      # The engine's message names `newdata` and the predictor at fault, which
      # is the useful half. Which of several models asked for it is the half
      # only this loop knows.
      error = function(e) {
        stop("model `", nm, "` cannot be scored on `newdata`: ",
             conditionMessage(e), call. = FALSE)
      }
    )
    as.numeric(out)
  })
  predicted <- matrix(unlist(columns, use.names = FALSE), nrow = n_obs)
  colnames(predicted) <- names(models)

  answered <- !is.na(observed)
  predictable <- stats::complete.cases(predicted)
  keep <- which(answered & predictable)
  n_dropped <- n_obs - length(keep)

  if (n_dropped > 0L) {
    # Named per reason, since the two are fixed by different things: a missing
    # answer is a row that cannot be scored at all, while a row no model could
    # be given is a row whose predictors are incomplete.
    no_answer <- sum(!answered)
    lost <- names(models)[apply(
      is.na(predicted[answered, , drop = FALSE]), 2, any
    )]
    message(
      n_dropped, " of ", n_obs, " row(s) of `newdata` were left out: ",
      paste(c(
        if (no_answer > 0L) paste0(no_answer, " with no observed outcome"),
        if (length(lost) > 0L) {
          paste0("some incomplete across the predictors of ",
                 paste(lost, collapse = ", "))
        }
      ), collapse = ", "),
      ". Every model is scored on the same rows, so a row one model cannot ",
      "predict is left out of all of them."
    )
  }

  if (length(keep) < 2L) {
    stop("only ", length(keep), " row(s) of `newdata` have an observed ",
         "outcome and a prediction from every model; at least 2 are needed.",
         call. = FALSE)
  }

  list(
    predicted = predicted[keep, , drop = FALSE],
    keep      = keep,
    n_obs     = n_obs,
    n_dropped = n_dropped
  )
}


#' Fold the per-model predictions into the long table the contract holds
#'
#' @keywords internal
#' @noRd
sa_prediction_table <- function(predicted, keep, observed) {
  models <- colnames(predicted)
  data.frame(
    model     = rep(models, each = length(keep)),
    row       = rep(keep, times = length(models)),
    observed  = rep(observed, times = length(models)),
    predicted = as.numeric(predicted),
    stringsAsFactors = FALSE
  )
}


#' Read `newdata` into the frame the rest of the evaluation works on
#'
#' @keywords internal
#' @noRd
sa_evaluate_newdata <- function(newdata) {
  if (is.matrix(newdata)) {
    newdata <- as.data.frame(newdata)
  }
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data.frame or a matrix.", call. = FALSE)
  }
  if (nrow(newdata) == 0L) {
    stop("`newdata` has zero rows, so there is nothing to score.",
         call. = FALSE)
  }
  newdata
}
