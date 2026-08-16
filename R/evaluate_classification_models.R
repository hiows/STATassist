# Scoring a set of fitted classifications on rows none of them was fitted on.
# The regression counterpart reports differences without tests, because a
# difference of held-out errors has no null this package is in a position to
# state. This side does carry tests, and the reason is that its quantities are
# functions of a class label and a probability: DeLong's variance, the IDI and
# the NRI are all built from per-row terms whose sampling distribution follows
# from the two classes being what they are, not from where the rows came from.
#
# All three are paired, and all three are against the baseline rather than
# across every pair of models. A model is proposed as an improvement on
# something, which is one comparison per model and the one the caller named by
# passing a baseline.


#' Read the observed classes as the event indicator the kernels take
#'
#' The direction is the fit's, not this call's: `outcome_lv[2]` is the class
#' `predict(type = "response")` reports the probability of, so it is the class
#' a 1 has to mean here for the two to line up.
#'
#' @keywords internal
#' @noRd
sa_response_code <- function(answer, outcome_lv) {
  labels <- as.character(answer)
  present <- unique(labels[!is.na(labels)])
  extra <- setdiff(present, outcome_lv)
  if (length(extra) > 0L) {
    stop("`answer` holds class(es) the models were not fitted on: ",
         paste(sort(extra), collapse = ", "), ". The models classify ",
         paste(outcome_lv, collapse = " and "),
         ". Reduce `newdata` to those two classes first.", call. = FALSE)
  }
  as.numeric(labels == outcome_lv[2])
}


#' Settle the class order, which the fits have already fixed
#'
#' `outcome_lv` and `control_label` are read here as a statement to be checked
#' rather than as an instruction. A fitted classification predicts the
#' probability of one particular class and cannot be re-pointed after the fact,
#' so naming the other one is a disagreement to report rather than a request to
#' carry out.
#'
#' @keywords internal
#' @noRd
sa_evaluate_levels <- function(model_lv, outcome_lv, control_label) {
  if (!is.null(outcome_lv)) {
    outcome_lv <- as.character(outcome_lv)
    if (length(outcome_lv) != 2L || anyNA(outcome_lv) ||
          anyDuplicated(outcome_lv) > 0L) {
      stop("`outcome_lv` must be two distinct level names, the reference ",
           "first.", call. = FALSE)
    }
    if (!identical(outcome_lv, model_lv)) {
      stop("`outcome_lv` is ", paste(outcome_lv, collapse = ", "),
           " but the models were fitted with ",
           paste(model_lv, collapse = ", "),
           ". A fitted classification predicts the probability of its own ",
           "second level and cannot be re-pointed here; refit to change it.",
           call. = FALSE)
    }
  }
  if (!is.null(control_label)) {
    if (length(control_label) != 1L || is.na(control_label)) {
      stop("`control_label` must be a single level name, the one the models ",
           "hold as the reference.", call. = FALSE)
    }
    if (!identical(as.character(control_label), model_lv[1])) {
      stop("`control_label` is ", control_label,
           " but the models were fitted with ", model_lv[1],
           " as the reference. A fitted classification cannot be re-pointed ",
           "here; refit to change it.", call. = FALSE)
    }
  }
  model_lv
}


#' Score fitted classifications on held-out rows
#'
#' Predicts one or more fitted two-class models on the same rows and reports how
#' each one discriminated, with three tests of each model against a baseline
#' where more than one was passed. Every model is read through
#' [predict.sa_model()] with `type = "response"`, so the five fitting functions
#' are interchangeable here and a logistic regression, a penalized one, a forest
#' and a machine are scored by the same arithmetic.
#'
#' @details
#' The rows are the intersection rather than the union. `predict()` on a model
#' answers `NA` for a row that is incomplete across *that model's* predictors,
#' so models fitted on different predictor sets come back with different rows
#' filled in. Scoring each on whatever it managed would put two AUCs from two
#' samples in one table, and all three comparisons below are **paired**
#' statistics that have no meaning at all across different rows. A row any model
#' cannot predict is therefore left out of all of them, with one message saying
#' how many went and why.
#'
#' The direction is the fits'. `outcome_lv[2]` is the class
#' `predict(type = "response")` reports the probability of, so it is the class
#' every number here is about: `sensitivity` is measured against it, and an AUC
#' above 0.5 means the models rank it above the reference. `outcome_lv` and
#' `control_label` are read as a statement to be checked rather than as an
#' instruction, since a fitted model cannot be re-pointed after the fact, and
#' naming the other class is an error rather than a silent reversal. All the
#' models must agree on it too, which [evaluate_regression_models()] has no
#' counterpart of.
#'
#' `$metrics` reports what needs no threshold first. `auc` comes with the
#' interval its own DeLong standard error gives, which is a Wald interval on a
#' bounded quantity and can therefore run past 1 for a strong classifier;
#' `brier` is the mean squared distance between the probability and the outcome,
#' which an AUC is blind to, since a model that ranks perfectly and predicts
#' every event at 0.6 has an AUC of 1. `accuracy`, `sensitivity` and
#' `specificity` do need one, and it is `threshold`, recorded in `$parameters`
#' so that the table says what it was measured at.
#'
#' `$comparisons` asks three different questions of the same pair of models,
#' which is why all three are reported rather than one being chosen.
#'
#' \describe{
#'   \item{`delta_auc`}{Whether the ranking improved, tested by DeLong's paired
#'     test. Blind to any change that does not reorder rows.}
#'   \item{`idi`}{How much further apart the two classes' predicted
#'     probabilities moved, on the probability scale. Sees exactly the change an
#'     AUC does not.}
#'   \item{`nri`}{How often a probability moved the right way, counting
#'     direction only, so a model that helps many rows slightly and hurts a few
#'     badly scores well here and can score badly on the IDI. Category-free: no
#'     risk strata are named, since their cut points are a clinical convention
#'     rather than a property of the data.}
#' }
#'
#' Every one of them is `new - baseline`, and every one is positive for a new
#' model that did better.
#'
#' @param baseline_model The reference model, as returned by
#'   [fit_logistic_regression()], [fit_elastic_net()], [fit_rf()] or [fit_svm()]
#'   with a two-class outcome. It is the first row of every table and what
#'   `$comparisons` is measured against.
#' @param new_models Named list of further models to hold against it, such as
#'   `list(selected = fit_2, penalized = fit_3)`, or `NULL` to score the
#'   baseline on its own. The names are what the tables and the plot legend call
#'   the models, so an unnamed list is refused rather than numbered.
#' @param newdata The rows to score, typically the test half of a
#'   [split_data()] result. Columns no model was fitted on are ignored.
#' @param answer The observed classes, either the name of a column of `newdata`
#'   or a vector with one entry per row. `NULL` reads the column the models were
#'   fitted to, which is the usual case.
#' @param outcome_lv The two classes, reference first, checked against the order
#'   the models were fitted with rather than used to change it. `NULL` takes
#'   theirs.
#' @param control_label The reference class on its own, checked the same way.
#' @param threshold Where to cut the predicted probability for `accuracy`,
#'   `sensitivity` and `specificity`. A row is called an event when its
#'   probability is greater than or equal to this.
#' @param conf_level Confidence level of every interval in the result.
#' @param baseline_label What to call the baseline in the tables and the legend.
#'
#' @return An object of class `sa_performance`, a plain list of nine elements.
#'
#'   \describe{
#'     \item{`analysis`}{`"classification_performance"`.}
#'     \item{`models`}{Model names, the baseline first, in the row order every
#'       table follows.}
#'     \item{`design`}{What was scored: the `outcome` label, its type, its
#'       `outcome_lv`, the `baseline` name, the row counts `n_obs`, `n_used` and
#'       `n_dropped`, and `n_events`, how many scored rows were
#'       `outcome_lv[2]`.}
#'     \item{`parameters`}{`threshold` and `conf_level`.}
#'     \item{`predictions`}{One row per model and scored row: `model`, `row`
#'       (the position in `newdata`), `observed` (1 for `outcome_lv[2]`, 0 for
#'       the reference) and `predicted` (the probability of `outcome_lv[2]`).}
#'     \item{`metrics`}{One row per model: `n_used`, `n_events`, `auc` with its
#'       interval, `brier`, and `accuracy`, `sensitivity` and `specificity` at
#'       `threshold`.}
#'     \item{`comparisons`}{One row per model other than the baseline, holding
#'       `delta_auc`, `idi` and `nri`, each with its interval and p-value, and
#'       the two class-wise components `nri_event` and `nri_nonevent`. Absent
#'       when nothing was compared.}
#'     \item{`curves`}{The ROC operating points, one row per model and distinct
#'       predicted value: `model`, `threshold`, `sensitivity`, `specificity`.
#'       Each curve opens at `threshold = Inf`, the cut above every prediction,
#'       so that it starts at the corner where nothing is called an event.}
#'     \item{`metadata`}{Package version, R version, platform and timestamp.}
#'   }
#'
#' @seealso [evaluate_regression_models()] for the continuous counterpart,
#'   [draw_roc_curve()] for the picture of this result, and
#'   [predict.sa_model()] for the call every model is read through.
#'
#' @references
#' DeLong, E. R., DeLong, D. M. and Clarke-Pearson, D. L. (1988). Comparing the
#' areas under two or more correlated receiver operating characteristic curves:
#' a nonparametric approach. *Biometrics*, 44(3), 837-845.
#'
#' Pencina, M. J., D'Agostino, R. B., D'Agostino, R. B. and Vasan, R. S. (2008).
#' Evaluating the added predictive ability of a new marker: from area under the
#' ROC curve to reclassification and beyond. *Statistics in Medicine*, 27(2),
#' 157-172.
#'
#' @examples
#' ## Two models of the same two-class outcome, scored on held-out rows.
#' iris2 <- iris[iris$Species != "setosa", ]
#' iris2$Species <- factor(iris2$Species)
#' train <- iris2[c(1:35, 51:85), ]
#' test <- iris2[c(36:50, 86:100), ]
#'
#' full <- fit_logistic_regression(train, outcome = "Species", cv = FALSE)
#' petal <- fit_logistic_regression(train, outcome = "Species",
#'                                  predictors = "Petal.Width", cv = FALSE)
#'
#' res <- evaluate_classification_models(full, list(petal_only = petal),
#'                                       newdata = test)
#' res
#' res$metrics
#'
#' @export
evaluate_classification_models <- function(baseline_model,
                                           new_models = NULL,
                                           newdata,
                                           answer = NULL,
                                           outcome_lv = NULL,
                                           control_label = NULL,
                                           threshold = 0.5,
                                           conf_level = 0.95,
                                           baseline_label = "baseline") {
  sa_check_scalar_num(threshold, "threshold", 0, 1)
  sa_check_scalar_num(conf_level, "conf_level", 0, 1,
                      lower_open = TRUE, upper_open = TRUE)

  newdata <- sa_evaluate_newdata(newdata)
  models <- sa_resolve_models(baseline_model, new_models, baseline_label)
  sa_check_model_family(models, "two classes", "evaluate_regression_models()")
  sa_check_model_agreement(models)

  model_lv <- sa_evaluate_levels(baseline_model$design$outcome_lv,
                                 outcome_lv, control_label)

  resolved <- sa_resolve_answer(answer, newdata, baseline_model)
  response <- sa_response_code(resolved$value, model_lv)

  collected <- sa_collect_predictions(models, newdata, response)
  response <- response[collected$keep]
  predicted <- collected$predicted

  n_events <- sum(response == 1)
  if (n_events == 0L || n_events == length(response)) {
    stop("the ", length(response), " scored row(s) hold a single class, ",
         model_lv[if (n_events == 0L) 1L else 2L],
         ", so there is nothing to discriminate between. Both classes have to ",
         "be present among the rows every model could predict.", call. = FALSE)
  }
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  scores <- vapply(seq_along(models), function(i) {
    p <- predicted[, i]
    area <- sa_auc_delong(response, p)
    c(n_used         = length(response),
      n_events       = n_events,
      auc            = area[["auc"]],
      auc_lower_conf = area[["auc"]] - z * area[["se"]],
      auc_upper_conf = area[["auc"]] + z * area[["se"]],
      brier          = sa_brier(response, p),
      sa_threshold_scores(response, p, threshold))
  }, numeric(9))
  metrics <- data.frame(model = names(models), t(scores),
                        stringsAsFactors = FALSE)
  # Counts are integers wherever else the package reports one, and they arrive
  # here as doubles only because they travelled beside the scores.
  metrics$n_used <- as.integer(metrics$n_used)
  metrics$n_events <- as.integer(metrics$n_events)
  rownames(metrics) <- NULL

  curves <- do.call(rbind, lapply(seq_along(models), function(i) {
    points <- sa_roc_points(response, predicted[, i])
    cbind(model = names(models)[i], points, stringsAsFactors = FALSE)
  }))
  rownames(curves) <- NULL

  comparisons <- NULL
  if (length(models) > 1L) {
    against <- seq_along(models)[-1]
    rows <- lapply(against, function(i) {
      new <- predicted[, i]
      old <- predicted[, 1]
      area <- sa_delong_test(response, new, old)
      discrimination <- sa_idi(response, old, new)
      reclassification <- sa_nri(response, old, new)
      data.frame(
        model                = names(models)[i],
        delta_auc            = area[["delta"]],
        delta_auc_lower_conf = area[["delta"]] - z * area[["se"]],
        delta_auc_upper_conf = area[["delta"]] + z * area[["se"]],
        delta_auc_pval       = area[["pval"]],
        idi                  = discrimination[["idi"]],
        idi_lower_conf       = discrimination[["idi"]] -
          z * discrimination[["se"]],
        idi_upper_conf       = discrimination[["idi"]] +
          z * discrimination[["se"]],
        idi_pval             = discrimination[["pval"]],
        nri                  = reclassification[["nri"]],
        nri_event            = reclassification[["nri_event"]],
        nri_nonevent         = reclassification[["nri_other"]],
        nri_lower_conf       = reclassification[["nri"]] -
          z * reclassification[["se"]],
        nri_upper_conf       = reclassification[["nri"]] +
          z * reclassification[["se"]],
        nri_pval             = reclassification[["pval"]],
        stringsAsFactors = FALSE
      )
    })
    comparisons <- do.call(rbind, rows)
    rownames(comparisons) <- NULL
  }

  sa_new_performance(
    analysis    = "classification_performance",
    models      = names(models),
    design      = list(
      outcome      = resolved$label,
      outcome_type = "two classes",
      outcome_lv   = model_lv,
      baseline     = baseline_label,
      n_obs        = collected$n_obs,
      n_used       = length(collected$keep),
      n_dropped    = collected$n_dropped,
      n_events     = n_events
    ),
    parameters  = list(threshold = threshold, conf_level = conf_level),
    predictions = sa_prediction_table(predicted, collected$keep, response),
    metrics     = metrics,
    comparisons = comparisons,
    curves      = curves
  )
}
