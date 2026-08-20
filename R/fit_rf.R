# The first model here with no coefficients, and most of what is different below
# follows from that. A forest is a few hundred trees grown on resampled rows and
# resampled predictors, and what it holds is splits rather than one number per
# predictor, so there is no effect per unit to put in the table the contract asks
# for. What a forest does answer is how much each predictor was worth to it, and
# that is what `estimate` carries here.
#
# The other two departures are the same fact seen from elsewhere. A tree splits a
# factor directly rather than on dummy columns, so the terms are the predictors
# themselves and not `k - 1` of them per factor. And every tree is grown on a
# bootstrap sample, so about a third of the rows are out of bag for each tree and
# the fit has a held-out score of its own before any fold is drawn. `fit_stats` is
# that score rather than the in-sample counterpart the other models report, which
# on a forest would be a model predicting rows it had memorised.

#' Fit a random forest
#'
#' Grows a forest of regression or classification trees on a set of predictors,
#' and scores it by cross-validation on the data it was fitted to. Which outcome
#' it is decides which forest is grown: a numeric outcome is a regression forest
#' and a two-class one a classification forest, as in [fit_elastic_net()], so one
#' function covers what [fit_linear_regression()] and [fit_logistic_regression()]
#' cover between them.
#'
#' What comes back in place of a coefficient table is an importance table: one row
#' per predictor, in descending order of how much the forest lost when that
#' predictor was shuffled. It says which predictors carried the fit and not which
#' way they pushed it, since a forest is free to use the same predictor in opposite
#' directions in different regions of the data.
#'
#' @details
#' \describe{
#'   \item{`estimate` is permutation importance, and there is no standard error
#'     beside it}{Each tree is scored on its out-of-bag rows, then scored again
#'     with one predictor's values shuffled among them, and `estimate` is the mean
#'     loss over the trees: `%IncMSE` for a regression and
#'     `MeanDecreaseAccuracy` for a classification. `impurity` is the other
#'     measure the same fit reports, the total drop in node impurity that
#'     predictor was responsible for — `IncNodePurity` or `MeanDecreaseGini`.
#'     Neither is an estimate of anything a distribution is defined over, so the
#'     inference columns are absent from the table rather than present and `NA`,
#'     exactly as they are for [fit_elastic_net()], and
#'     `is.null(fit$coefficients$pval)` tells the two kinds of table apart.
#'     The values are the unscaled ones. `randomForest::importance()` divides them
#'     by their standard deviation across trees by default, which makes a t-shaped
#'     ratio that is not referred to any distribution, so what is reported is the
#'     mean loss itself, on the scale of the metric it was measured in.}
#'   \item{A negative importance is a value, not a gap}{Shuffling a predictor that
#'     carried nothing can leave the forest very slightly better than it was, and
#'     the mean loss then comes out below zero. It reads as it should: this
#'     predictor did no better than its own permutation.}
#'   \item{The terms are the predictors}{A tree splits a factor on its levels
#'     directly, so nothing is dummy coded and a `k`-level factor is one term
#'     rather than `k - 1`. This is the one model here where `terms` and
#'     `design$predictors` hold the same names — in a different order, since the
#'     table is sorted by importance and the design records the order the columns
#'     were read in.}
#'   \item{`fit_stats` is measured out of bag}{Every tree is grown on a bootstrap
#'     sample of the rows, which leaves the rest of them out of bag for that tree,
#'     and `randomForest` predicts each row from the trees that did not see it.
#'     That is already a held-out score, and it is the honest one: a forest asked
#'     to predict the rows it was fitted to reports something close to perfect
#'     whatever the data held. The `oob_` prefix says these numbers are not the
#'     in-sample kind the other models report.}
#'   \item{Cross-validation chooses `mtry` when there is more than one to choose
#'     from}{`mtry` is the only argument [caret::train()] tunes for a forest;
#'     `ntree` and `nodesize` are the same for every candidate. A single value is
#'     scored rather than chosen, as in [fit_linear_regression()], and a vector is
#'     compared the way [fit_elastic_net()] compares penalties, `parameters$mtry`
#'     then holding the value that won and `performance` one row per candidate.
#'     Unlike a penalty, `mtry` has a default worth fitting, so `NULL` resolves to
#'     the rule of thumb — the square root of the predictor count for a
#'     classification, a third of it for a regression — rather than to a grid, and
#'     `cv = FALSE` needs nothing added to work.}
#' }
#'
#' Rows with a missing value in the outcome or in any predictor are dropped before
#' the folds are drawn rather than inside each fit, and `design$n_dropped` reports
#' how many went. A predictor that takes a single value is left out with a message.
#'
#' One predictor is enough. There is no budget to divide as there is under a
#' penalty, so `predictors = "wt"` grows a forest of stumps on that one column
#' rather than raising the error [fit_elastic_net()] raises.
#'
#' For a two-class outcome the direction is `outcome_lv`, read as
#' [fit_logistic_regression()] reads it: the first level is the reference, so
#' `oob_sensitivity` is the share of `outcome_lv[2]` the forest found,
#' `oob_specificity` the share of `outcome_lv[1]` it left alone, and
#' `predict(model, newdata, type = "response")` is the probability of
#' `outcome_lv[2]`. The importance table itself has no direction to report, which
#' is why there is no `odds_ratio` column beside it. `control_label` says the
#' same thing with one name instead of two and says it alone, so it too is enough
#' to make a column of zeroes and ones a classification. Naming both and pointing
#' them at different classes is an error rather than a re-pointing, since an
#' `outcome_lv` holds the two classes and nothing else; the comparison functions
#' read the argument the other way, as [compare_two_groups()] describes.
#'
#' A forest is random beyond the folds — the rows of each tree and the predictors
#' of each split are both draws — so `seed` fixes more here than it does in the
#' other models, where it fixes only the fold assignment. Two calls without one
#' give slightly different importance values on the same data.
#'
#' @param data A data.frame (or matrix) in wide format, one row per observation.
#'   Typically the training half of a [split_data()] result.
#' @param outcome The outcome, either the name of a column of `data` or a vector
#'   with one entry per row. A numeric outcome is fitted as a regression, and a
#'   factor, character or logical one as a two-class classification.
#' @param predictors Column names to fit on, or `NULL` for every column of `data`
#'   except the outcome. Numeric, logical, factor and character columns are
#'   accepted; a column that takes a single value is left out with a message,
#'   since it cannot contribute.
#' @param outcome_lv The two classes, reference first, so that the reported
#'   sensitivity and the predicted probability are about the second one. Supplying
#'   it also states that the outcome is to be classified, which is the one thing a
#'   numeric column holding two values cannot say on its own. `NULL` sorts the
#'   classes of an outcome that is not numeric, and leaves a numeric one to the
#'   regression.
#' @param control_label The reference class on its own, for when the other one
#'   needs no saying. Defaults to `outcome_lv[1]`, so a call that names one of the
#'   two names the reference either way; naming both and disagreeing is an error.
#' @param mtry Predictors offered at each split, one value to score or several to
#'   choose between. `NULL` is the rule of thumb: `floor(sqrt(p))` for a
#'   classification and `floor(p / 3)` for a regression, at least 1. A value above
#'   the predictor count is refused rather than quietly reset by the engine.
#' @param ntree Trees to grow. More is steadier rather than more prone to
#'   overfitting, and costs time.
#' @param nodesize Smallest number of rows a leaf may hold, so larger grows
#'   shallower trees. `NULL` is `randomForest`'s own default: 1 for a
#'   classification and 5 for a regression.
#' @param cv Whether to cross-validate. `FALSE` grows the single candidate the
#'   grid names, once, and reports no resampled performance. The out-of-bag
#'   `fit_stats` are there either way.
#' @param cv_method Resampling scheme: `"repeated_kfold"` for `n_repeat` runs of
#'   `n_fold`-fold cross-validation, `"kfold"` for one, or `"loocv"` for
#'   leave-one-out.
#' @param n_fold Folds per run, used by `"repeated_kfold"` and `"kfold"`.
#' @param n_repeat Number of runs, used by `"repeated_kfold"`.
#' @param seed Seed for the forest and the fold assignment, or `NULL` to use the
#'   stream as it stands. Supplying one does not disturb the caller: the previous
#'   random number state is put back when the function returns.
#'
#' @return An object of class `sa_model`, the same eleven elements
#'   [fit_linear_regression()] returns, with these differences:
#'
#'   \describe{
#'     \item{`analysis`}{`"random_forest"`.}
#'     \item{`terms`}{The predictors, in descending order of importance, since a
#'       forest expands no factor into dummy terms.}
#'     \item{`design`}{Holds `outcome_lv`, `n_events` and `event_rate` for a
#'       two-class outcome, as [fit_logistic_regression()] does, and neither for a
#'       continuous one.}
#'     \item{`parameters`}{Holds the `mtry` that was fitted rather than the grid
#'       that was searched, `ntree`, `nodesize`, and `n_candidates`, how many
#'       values of `mtry` were scored. The grid itself is the rows of
#'       `performance`.}
#'     \item{`coefficients`}{`estimate`, the permutation importance, and
#'       `impurity`, the impurity-based one. There is no intercept row, no
#'       `odds_ratio` and no inference column.}
#'     \item{`fit_stats`}{Measured on the out-of-bag predictions: `oob_r_squared`,
#'       `oob_rmse` and `oob_mae` for a regression, and `oob_accuracy`,
#'       `oob_error`, `oob_kappa`, `oob_sensitivity` and `oob_specificity` for a
#'       classification, each with `n_oob`, the rows that were out of bag of at
#'       least one tree and so could be predicted.}
#'     \item{`performance`}{One row per value of `mtry`, the chosen one being the
#'       row that matches `parameters$mtry`, or `NULL` when `cv = FALSE`.}
#'   }
#'
#' @seealso [fit_linear_regression()], [fit_logistic_regression()] and
#'   [fit_elastic_net()] for the models that do report an effect per predictor,
#'   [split_data()], which defines the rows this is fitted on, [coef.sa_model()]
#'   for the importance table, and [predict.sa_model()] for predicting the rows it
#'   was not fitted on.
#'
#' @examples
#' ## A single forest without resampling (fast enough for examples).
#' fit <- fit_rf(mtcars, outcome = "mpg", ntree = 50, cv = FALSE, seed = 1)
#' fit$coefficients
#' fit$fit_stats
#'
#' \donttest{
#' ## Cross-validated fit and an `mtry` search.
#' fit_rf(mtcars, outcome = "mpg", ntree = 200,
#'        cv_method = "kfold", n_fold = 5, seed = 1)
#' tuned <- fit_rf(mtcars, outcome = "mpg", mtry = c(2, 5, 10), ntree = 200,
#'                 cv_method = "kfold", n_fold = 5, seed = 1)
#' tuned$parameters[c("mtry", "n_candidates")]
#'
#' ## Two-class outcome and a simulated regression scored against known truth.
#' iris2 <- iris[iris$Species != "setosa", ]
#' clf <- fit_rf(iris2, outcome = "Species",
#'               outcome_lv = c("versicolor", "virginica"),
#'               ntree = 200, cv = FALSE)
#' clf$fit_stats[c("oob_accuracy", "oob_sensitivity", "oob_specificity")]
#' sim <- simulate_regression(seed = 1)
#' rf <- do.call(fit_rf, c(sim$args, list(ntree = 200, cv = FALSE)))
#' scored <- merge(rf$coefficients, sim$truth_term, by = "terms")
#' tapply(scored$estimate, scored$beta != 0, mean)
#' }
#'
#' @export
fit_rf <- function(data,
                   outcome,
                   predictors = NULL,
                   outcome_lv = NULL,
                   control_label = outcome_lv[1],
                   mtry = NULL,
                   ntree = 500,
                   nodesize = NULL,
                   cv = TRUE,
                   cv_method = c("repeated_kfold", "kfold", "loocv"),
                   n_fold = 5,
                   n_repeat = 5,
                   seed = NULL) {

  cv_method <- match.arg(cv_method)
  ntree <- sa_check_count(ntree, "ntree", 1)

  input <- sa_resolve_model_input(data, outcome, predictors)

  # One argument, two forests, and the outcome decides which, the same way it
  # does in `fit_elastic_net()`. A numeric column is a regression, since that is
  # what a number usually is, and anything else is a set of class labels.
  # `outcome_lv` overrules the guess, and so does `control_label` on its own,
  # which is the only way to say that a column of zeroes and ones is two classes
  # rather than two numbers.
  classify <- !is.null(outcome_lv) || !is.null(control_label) ||
    !is.numeric(input$y)
  if (!classify && length(unique(input$y)) == 2L) {
    message("`outcome` is numeric and takes two values, so it was fitted as a ",
            "regression. Pass `outcome_lv` or `control_label`, or a factor ",
            "column, to model it as a classification.")
  }

  if (classify) {
    y <- sa_outcome_levels(input$y, outcome_lv, control_label,
                           model = "a random forest")
    outcome_lv <- levels(y)
  } else {
    if (!all(is.finite(input$y))) {
      stop("`outcome` holds non-finite value(s), which a leaf of a regression ",
           "tree cannot average.", call. = FALSE)
    }
    y <- input$y
  }

  if (is.null(nodesize)) {
    nodesize <- if (classify) 1L else 5L
  }
  nodesize <- sa_check_count(nodesize, "nodesize", 1)
  grid <- sa_rf_grid(mtry, length(input$predictors), classify, cv)
  ctrl <- sa_train_control(cv, cv_method, n_fold, n_repeat, input$n_used)

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  # The predictor frame goes to the engine as it is, which is where this parts
  # company with `fit_elastic_net()`. `glmnet` reads a numeric matrix by position
  # and had to be handed one built here; a tree splits a factor on its levels, so
  # the frame is what it wants. `importance = TRUE` reaches `randomForest()`
  # through `caret`'s `...`, and is what makes the permutation measure available
  # at all: without it the fit reports impurity alone.
  fit <- sa_fit_engine(
    caret::train(x = input$x, y = y, method = "rf", trControl = ctrl$control,
                 tuneGrid = grid, ntree = ntree, nodesize = nodesize,
                 importance = TRUE),
    "Random forest"
  )
  model <- fit$finalModel

  coefs <- sa_rf_importance(model, classify)
  fit_stats <- if (classify) {
    sa_rf_class_stats(model$predicted, y, outcome_lv)
  } else {
    sa_rf_reg_stats(model$predicted, y)
  }
  fit_stats <- lapply(fit_stats, function(v) unname(as.numeric(v)))

  design <- list(
    outcome      = input$outcome,
    outcome_type = if (classify) "two classes" else "continuous"
  )
  if (classify) {
    n_events <- sum(y == outcome_lv[2])
    design <- c(design, list(outcome_lv = outcome_lv,
                             n_events   = n_events,
                             event_rate = n_events / input$n_used))
  }

  sa_new_model(
    analysis = "random_forest",
    terms    = coefs$terms,
    design   = c(design, list(
      n_obs              = input$n_obs,
      n_used             = input$n_used,
      n_dropped          = input$n_dropped,
      predictors         = input$predictors,
      dropped_predictors = input$dropped_predictors
    ), sa_design_lv(input$predictor_lv)),
    # The `mtry` that ran, not the grid that was asked for: `performance` holds
    # every candidate, so recording the grid here as well would say the same thing
    # twice and leave two places for it to be wrong.
    parameters = list(
      mtry         = as.integer(fit$bestTune$mtry),
      ntree        = ntree,
      nodesize     = nodesize,
      n_candidates = nrow(grid),
      cv           = cv,
      cv_method    = ctrl$cv_method,
      n_fold       = ctrl$n_fold,
      n_repeat     = ctrl$n_repeat,
      seed         = seed
    ),
    coefficients = coefs,
    fit_stats    = fit_stats,
    performance  = sa_model_frame(fit$results),
    resampling   = sa_model_frame(fit$resample),
    engine       = list(
      package = "caret",
      method  = "rf",
      label   = paste("Random forest",
                      if (classify) "classification" else "regression"),
      metrics = fit$perfNames
    ),
    fit = fit
  )
}


#' The importance table a forest answers with instead of coefficients
#'
#' Two measures come out of the same fit and both are kept, since they disagree in
#' a way that is worth seeing: permutation importance is measured out of bag and
#' impurity importance is measured on the splits themselves, so a predictor with
#' many distinct values can earn impurity it does not earn on held-out rows.
#' `estimate` is the permutation one, being the measure that was scored on rows the
#' tree had not seen.
#'
#' Descending rather than in the order the columns arrived. `terms` is the row
#' order of every table in the result, and for this model the order worth reading
#' first is the most important predictor; the order the columns were read in is
#' `design$predictors`, which keeps it.
#'
#' @param model The `randomForest` object inside the fit.
#' @param classify Whether the outcome was classified, which is what decides
#'   which two of the four importance columns to read.
#'
#' @return data.frame with `terms`, `estimate` and `impurity`.
#'
#' @keywords internal
#' @noRd
sa_rf_importance <- function(model, classify) {
  importance <- model$importance
  permutation <- if (classify) "MeanDecreaseAccuracy" else "%IncMSE"
  purity <- if (classify) "MeanDecreaseGini" else "IncNodePurity"

  absent <- setdiff(c(permutation, purity), colnames(importance))
  if (length(absent) > 0L) {
    stop("internal error: the forest reports no ",
         paste(absent, collapse = " and no "), " column, so ",
         "`importance = TRUE` did not reach the engine.", call. = FALSE)
  }

  out <- data.frame(
    terms    = rownames(importance),
    estimate = unname(importance[, permutation]),
    impurity = unname(importance[, purity]),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$estimate, decreasing = TRUE), ]
  rownames(out) <- NULL
  out
}


#' Goodness of fit of a forest, measured on the rows each tree did not see
#'
#' The other models measure themselves on the rows they were fitted to, which for
#' a forest is a model predicting what it memorised. What a forest has instead is
#' the out-of-bag prediction: each row predicted by the trees whose bootstrap
#' sample left it out, which is already held out and comes free of the resampling.
#'
#' A row can in principle be in the bag of every tree, in which case it has no
#' out-of-bag prediction and cannot be scored. Those rows are left out and `n_oob`
#' says how many were scored, rather than the whole statistic coming back `NA`
#' because of one of them.
#'
#' @param oob_prediction `model$predicted`, one entry per fitted row.
#' @param observed The outcome as the fit received it.
#'
#' @return Named list of scalars.
#'
#' @keywords internal
#' @noRd
sa_rf_reg_stats <- function(oob_prediction, observed) {
  at <- !is.na(oob_prediction)
  fitted_value <- oob_prediction[at]
  observed <- observed[at]

  residual <- observed - fitted_value
  sse <- sum(residual^2)
  sst <- sum((observed - mean(observed))^2)

  list(
    oob_r_squared = if (sst > 0) 1 - sse / sst else NA_real_,
    oob_rmse      = sqrt(mean(residual^2)),
    oob_mae       = mean(abs(residual)),
    n_oob         = sum(at)
  )
}


#' @rdname sa_rf_reg_stats
#' @param outcome_lv The two classes, reference first, so that the sensitivity is
#'   the one of `outcome_lv[2]` and the whole result reads in one direction.
#' @keywords internal
#' @noRd
sa_rf_class_stats <- function(oob_prediction, observed, outcome_lv) {
  at <- !is.na(oob_prediction)
  called <- as.character(oob_prediction[at]) == outcome_lv[2]
  positive <- as.character(observed[at]) == outcome_lv[2]

  accuracy <- mean(called == positive)
  # Cohen's kappa on the two-by-two table: the agreement above what the observed
  # and the predicted labels would reach from their marginals alone.
  expected <- mean(called) * mean(positive) + mean(!called) * mean(!positive)

  list(
    oob_accuracy    = accuracy,
    oob_error       = 1 - accuracy,
    oob_kappa       = if (expected < 1) {
      (accuracy - expected) / (1 - expected)
    } else {
      NA_real_
    },
    oob_sensitivity = if (any(positive)) mean(called[positive]) else NA_real_,
    oob_specificity = if (any(!positive)) mean(!called[!positive]) else NA_real_,
    n_oob           = sum(at)
  )
}
