# Scoring a set of fitted regressions on rows none of them was fitted on. What
# `fit_*()` already reports in `$performance` is a resampled score, measured
# inside the folds of the data the model was fitted to; this is the other kind,
# measured once on a held-out set the caller drew. The two are different
# questions and the names are shared on purpose, so `rmse` here and `RMSE`
# there read in the same unit.
#
# There is no test in this file. A difference of held-out errors has no null
# this function is in a position to state — the rows are not a sample from
# anything the caller described — so what it reports is the two numbers and
# their difference, and the reader supplies the judgement. The classification
# side does carry tests, because the quantities there are functions of a class
# label and a probability and have sampling distributions that do not depend on
# where the rows came from.


#' Per-model scores against the observed outcome
#'
#' The calibration line is `lm(predicted ~ observed)` written out in closed
#' form, which is the same two numbers and leaves no engine object to store.
#'
#' @keywords internal
#' @noRd
sa_regression_scores <- function(observed, predicted, var_observed) {
  n <- length(observed)
  residual <- predicted - observed
  sse <- sum(residual^2)
  sst <- var_observed * (n - 1)

  # An outcome that takes one value over the scored rows has nothing for a
  # proportion of variance or a calibration slope to be measured against. It is
  # a property of the rows rather than of the model, so it is the caller who is
  # told, once, rather than each model reporting its own NaN.
  degenerate <- !is.finite(var_observed) || var_observed <= 0
  if (degenerate) {
    correlation <- NA_real_
    r_squared <- NA_real_
    slope <- NA_real_
    intercept <- NA_real_
  } else {
    var_predicted <- stats::var(predicted)
    # A model that answers the same value for every row ranks nothing, so it
    # has no correlation. Its calibration line is still defined and is flat,
    # which is the honest description of what it did.
    correlation <- if (var_predicted > 0) {
      stats::cor(observed, predicted)
    } else {
      NA_real_
    }
    r_squared <- 1 - sse / sst
    slope <- stats::cov(observed, predicted) / var_observed
    intercept <- mean(predicted) - slope * mean(observed)
  }

  c(n_used          = n,
    cor             = correlation,
    r_squared       = r_squared,
    rmse            = sqrt(sse / n),
    mae             = mean(abs(residual)),
    bias            = mean(residual),
    calib_slope     = slope,
    calib_intercept = intercept)
}


#' Score fitted regressions on held-out rows
#'
#' Predicts one or more fitted regressions on the same rows and reports how each
#' one did, with the differences against a baseline where more than one model
#' was passed. Every model is read through [predict.sa_model()], so the five
#' fitting functions are interchangeable here and a linear model, a penalized
#' one, a forest and a machine are scored by the same arithmetic.
#'
#' @details
#' The rows are the intersection rather than the union. `predict()` on a model
#' answers `NA` for a row that is incomplete across *that model's* predictors,
#' so a baseline fitted on nine columns and a reduced model fitted on four
#' disagree about any row missing one of the extra five. Scoring each model on
#' whatever it managed would put two numbers from two samples in one table and
#' call their difference an improvement, so a row that any model cannot predict
#' is left out of all of them, and a single message reports how many went and
#' why. `design$n_used` and `design$n_dropped` record the outcome of that.
#'
#' Every model must have been fitted to the same outcome, and to a continuous
#' one. A classification is refused by name and pointed at
#' [evaluate_classification_models()] rather than being scored by correlating a
#' probability against a class label, which would produce a number.
#'
#' `$metrics` reports `cor` and `r_squared` side by side because they answer
#' different questions and agree only for predictions that need no calibration.
#' `r_squared` is `1 - SSE/SST` on these rows, the fraction of the variance the
#' predictions actually removed, and it is negative for a model that does worse
#' than the mean of the outcome. `cor^2` is what that would be if the
#' predictions were first rescaled by a line fitted to these same rows, so the
#' gap between the two is what `calib_slope` and `calib_intercept` describe: the
#' line is `lm(predicted ~ observed)`, the same orientation
#' [draw_prediction_plot()] draws, so a slope under one is a model whose
#' predictions are compressed towards their own mean.
#'
#' `rmse` and `mae` carry the names [fit_linear_regression()] and the rest use
#' in `$performance`, so a resampled score and a held-out score read in the same
#' unit under the same name. They are not the same number and are not meant to
#' be: one is measured inside the folds of the training data and the other on
#' rows drawn away from it.
#'
#' `$comparisons` is every column of `$metrics` as `new - baseline`, so a
#' positive `delta_cor` and a **negative** `delta_rmse` both say the new model
#' did better. There is no p-value beside them, since a difference of held-out
#' errors has no null this function is in a position to state.
#'
#' @param baseline_model The reference model, as returned by
#'   [fit_linear_regression()], [fit_elastic_net()], [fit_rf()] or [fit_svm()].
#'   It is the first row of every table and what `$comparisons` subtracts.
#' @param new_models Named list of further models to hold against it, such as
#'   `list(selected = fit_2, penalized = fit_3)`, or `NULL` to score the
#'   baseline on its own. The names are what the tables and the plot legend call
#'   the models, so an unnamed list is refused rather than numbered.
#' @param newdata The rows to score, typically the test half of a
#'   [split_data()] result. Columns no model was fitted on are ignored.
#' @param answer The observed outcome, either the name of a column of `newdata`
#'   or a vector with one entry per row. `NULL` reads the column the models were
#'   fitted to, which is the usual case; a model fitted from a vector remembers
#'   no name to look up and requires this.
#' @param baseline_label What to call the baseline in the tables and the legend.
#'
#' @return An object of class `sa_performance`, a plain list of eight elements.
#'
#'   \describe{
#'     \item{`analysis`}{`"regression_performance"`.}
#'     \item{`models`}{Model names, the baseline first, in the row order every
#'       table follows.}
#'     \item{`design`}{What was scored: the `outcome` label, its type, the
#'       `baseline` name, and the row counts `n_obs`, `n_used` and `n_dropped`.}
#'     \item{`parameters`}{Empty. Scoring a regression takes no choices, which
#'       is the difference between this slot here and in a classification
#'       evaluation.}
#'     \item{`predictions`}{One row per model and scored row: `model`, `row`
#'       (the position in `newdata`), `observed` and `predicted`.}
#'     \item{`metrics`}{One row per model: `n_used`, `cor`, `r_squared`, `rmse`,
#'       `mae`, `bias`, `calib_slope` and `calib_intercept`.}
#'     \item{`comparisons`}{One row per model other than the baseline, holding
#'       `delta_cor`, `delta_r_squared`, `delta_rmse` and `delta_mae` as
#'       `new - baseline`. Absent when nothing was compared.}
#'     \item{`metadata`}{Package version, R version, platform and timestamp.}
#'   }
#'
#' @seealso [evaluate_classification_models()] for the two-class counterpart,
#'   [draw_prediction_plot()] for the picture of this result, and
#'   [predict.sa_model()] for the call every model is read through.
#'
#' @examples
#' ## Two models of the same outcome, scored on rows neither was fitted on.
#' train <- mtcars[1:24, ]
#' test <- mtcars[25:32, ]
#'
#' full <- fit_linear_regression(train, outcome = "mpg",
#'                               predictors = c("wt", "hp", "disp"),
#'                               cv = FALSE)
#' small <- fit_linear_regression(train, outcome = "mpg",
#'                                predictors = "wt", cv = FALSE)
#'
#' res <- evaluate_regression_models(full, list(weight_only = small),
#'                                   newdata = test)
#' res
#' res$metrics
#'
#' @export
evaluate_regression_models <- function(baseline_model,
                                       new_models = NULL,
                                       newdata,
                                       answer = NULL,
                                       baseline_label = "baseline") {
  newdata <- sa_evaluate_newdata(newdata)
  models <- sa_resolve_models(baseline_model, new_models, baseline_label)
  sa_check_model_family(models, "continuous",
                        "evaluate_classification_models()")
  sa_check_model_agreement(models)

  resolved <- sa_resolve_answer(answer, newdata, baseline_model)
  observed <- resolved$value
  if (!is.numeric(observed)) {
    stop("`answer` must be numeric: a regression is scored against the value ",
         "it predicted, not against a label. Got ", class(observed)[1],
         ". Use evaluate_classification_models() for a two-class outcome.",
         call. = FALSE)
  }

  collected <- sa_collect_predictions(models, newdata, observed)
  observed <- observed[collected$keep]
  predicted <- collected$predicted

  var_observed <- stats::var(observed)
  if (!is.finite(var_observed) || var_observed <= 0) {
    warning("the observed outcome takes a single value over the ",
            length(observed), " scored row(s), so `cor`, `r_squared` and the ",
            "calibration line are NA. `rmse`, `mae` and `bias` are still ",
            "reported.", call. = FALSE)
  }

  scores <- vapply(
    seq_along(models),
    function(i) sa_regression_scores(observed, predicted[, i], var_observed),
    numeric(8)
  )
  metrics <- data.frame(model = names(models), t(scores),
                        stringsAsFactors = FALSE)
  # A count is an integer wherever else the package reports one, and it arrives
  # here as a double only because it travelled beside the scores.
  metrics$n_used <- as.integer(metrics$n_used)
  rownames(metrics) <- NULL

  # One warning for the whole run rather than one per model, the way a feature
  # that could not be tested is reported by the comparison functions.
  flat <- names(models)[is.na(metrics$cor) & !is.na(metrics$calib_slope)]
  if (length(flat) > 0L) {
    warning("`cor` is NA for ", length(flat), " model(s) that answered a ",
            "single value over the scored rows: ",
            paste(flat, collapse = ", "),
            ". A prediction that does not vary ranks nothing.", call. = FALSE)
  }

  comparisons <- NULL
  if (length(models) > 1L) {
    against <- seq_along(models)[-1]
    comparisons <- data.frame(
      model           = names(models)[against],
      delta_cor       = metrics$cor[against] - metrics$cor[1],
      delta_r_squared = metrics$r_squared[against] - metrics$r_squared[1],
      delta_rmse      = metrics$rmse[against] - metrics$rmse[1],
      delta_mae       = metrics$mae[against] - metrics$mae[1],
      stringsAsFactors = FALSE
    )
  }

  sa_new_performance(
    analysis    = "regression_performance",
    models      = names(models),
    design      = list(
      outcome      = resolved$label,
      outcome_type = "continuous",
      baseline     = baseline_label,
      n_obs        = collected$n_obs,
      n_used       = length(collected$keep),
      n_dropped    = collected$n_dropped
    ),
    parameters  = list(),
    predictions = sa_prediction_table(predicted, collected$keep, observed),
    metrics     = metrics,
    comparisons = comparisons
  )
}
