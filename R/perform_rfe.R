# The first function here that searches rather than fits. Every `fit_*()` is
# handed the predictors and answers about them; this one is handed candidates and
# answers which of them to keep, so what it returns is `sa_selection` and its row
# axis is `candidates` rather than `terms`.
#
# The search is `caret`'s recursive feature elimination, wrapped so that the
# ranking it eliminates by is one this package is willing to stand behind. Two of
# `caret`'s own ranking functions are replaced for reasons that are the substance
# of this file:
#
# `lmFuncs$rank` ranks by `abs(coef(object))`. A coefficient is an effect per unit
# of its predictor, so ranking by its size ranks by the units the predictors
# happen to be measured in: the same model with a predictor in millimetres rather
# than metres eliminates in a different order. What replaces it is the absolute t
# statistic, which is the coefficient divided by its own standard error and so has
# no units left, and which is what `caret::lrFuncs` already ranks a logistic
# regression by. The two models therefore rank on the same scale.
#
# Neither `lm` nor `glm` sees the columns that were passed in. A factor with `k`
# levels is `k - 1` coefficients, so a ranking read off the coefficients is a
# ranking of dummy columns and `caret` cannot match it back to `x`. The rank
# function here folds the dummies back into the column they came from and keeps
# the largest statistic among them, so what is eliminated is always a column of
# the input. That is what makes `$selected` something that can be handed straight
# back to `fit_rf(predictors = )`, which is the whole point of running a selection.
#
# `$fit` is the `rfe` object and carries no `sa_fit` class, unlike the `$fit` of a
# model. That class exists to route `coef()` and `summary()` to a `$finalModel`,
# and an elimination has none: what it has is `predict()`, `print()` and
# `varImp()`, all of which `caret` already answers on the object as it is.

#' Select the predictors worth keeping
#'
#' Runs a recursive feature elimination: the predictors are ranked, the weakest is
#' dropped, and the model is scored again, over and over, so that every subset size
#' gets a resampled score. What comes back is the size that scored best, the
#' predictors of that size, and the whole profile the search walked, so that a
#' choice of two predictors over eight can be read against what the other six were
#' worth.
#'
#' The input is the wide format the model functions take, **one row per
#' observation with one column as the outcome**, and is normally the training half
#' of a [split_data()] result. Selecting on data a model is later scored on is how
#' a selection flatters itself; see the details.
#'
#' @details
#' # What is resampled
#'
#' The elimination is inside the resampling, not before it. Each fold ranks the
#' predictors on its own training rows, peels them down to each size in turn, and
#' scores every size on rows it did not rank on, so a predictor that looks useful
#' only on the rows that chose it is caught. Ranking once on all the rows and
#' cross-validating afterwards is the mistake this ordering exists to avoid, and it
#' is why there is no `cv` argument: an elimination with nothing held out has no
#' score to choose a size by, so it is not a shorter version of this function but a
#' different and wrong one.
#'
#' What the resampling cannot do is make the reported score an honest estimate of
#' the selected model. The size was chosen because it scored best, so `profile`
#' reads high at the size it picked, for the same reason a maximum of noisy numbers
#' is above their mean. Score the selection on the test half of [split_data()],
#' which the search never saw, and read `profile` as the shape of the search rather
#' than as a performance claim.
#'
#' The ranking is computed once per fold, at the full set of predictors, and the
#' peeling follows it. Re-ranking after every drop is `caret`'s `rerank` and is a
#' far slower procedure — one refit per remaining predictor per fold — so it is not
#' what this function does.
#'
#' # Which model does the ranking
#'
#' `model` names what is fitted inside the search, and the outcome has to agree
#' with it. `"linear"` is a continuous outcome, `"logistic"` a two-class one, and
#' `"rf"` is either, following the outcome the way [fit_rf()] does. A disagreement
#' is an error naming the model that would have fitted, rather than a silently
#' different analysis.
#'
#' What each ranks by is reported in `engine$importance`:
#'
#' \describe{
#'   \item{`"linear"`, `"logistic"`}{The absolute t or Wald z statistic of each
#'     coefficient, largest first. Not the coefficient itself, which is an effect
#'     per unit of its predictor and so would rank by the units the predictors were
#'     measured in; the statistic divides that by its own standard error and has no
#'     units left. A factor is ranked as one column, by the largest statistic among
#'     its levels, and a term the fit could not estimate ranks at zero, since a
#'     predictor the others already span costs nothing to drop.}
#'   \item{`"rf"`}{Permutation importance, the loss when a predictor's values are
#'     shuffled among the rows, measured out of bag. This is the same measure
#'     [fit_rf()] reports as `estimate`, so a forest's ranking here and its
#'     importance table there can be read together.}
#' }
#'
#' A forest inside the search grows at `randomForest()`'s own `mtry` for each
#' subset — the square root of the predictor count for a classification and a third
#' of it for a regression, which is the rule [fit_rf()] uses — rather than at one
#' value throughout. A fixed `mtry` would exceed the predictor count at the small
#' end of the profile, where the whole question is what a handful of predictors can
#' do.
#'
#' # Which class the direction is fixed on
#'
#' `outcome_lv` is read as it is everywhere else in this package: the first level
#' is the reference, so `outcome_lv = c("control", "case")` searches for the
#' predictors of `case`. `control_label` names that same first level on its own, for
#' the usual case where the sort has it backwards and the other level needs no
#' saying. The two are the same statement, so passing both and disagreeing is an
#' error rather than a precedence rule.
#'
#' For this function the direction changes nothing but the reading. A ranking is a
#' statement about how much a predictor is worth, which is the same number whichever
#' class is called the reference; what the levels decide is which class `design`
#' counts as the events, and which way a later [fit_logistic_regression()] on the
#' selected predictors will read.
#'
#' # What is dropped before anything runs
#'
#' The same listwise deletion the model functions use: rows that are missing the
#' outcome or any candidate go before the folds are drawn, so every fold sees the
#' same rows, and `design$n_dropped` says how many went. A candidate that takes a
#' single value is left out with a message, since a column with nothing in it
#' cannot be eliminated for a reason.
#'
#' # Portability
#'
#' Everything but `$fit` is a scalar, a character vector, a named list or a
#' data.frame, so dropping that one slot leaves an object that writes out as JSON.
#' In a Python transcription this is `sklearn.feature_selection.RFECV`.
#'
#' @param data A data.frame (or matrix) in wide format, one row per observation.
#'   Typically the training half of a [split_data()] result.
#' @param outcome The outcome, either the name of a column of `data` or a vector
#'   with one entry per row.
#' @param predictors Candidate column names, or `NULL` for every column of `data`
#'   except the outcome. Numeric, logical, factor and character columns are all
#'   accepted.
#' @param outcome_lv For a two-class outcome, the two classes with the reference
#'   first. `NULL` sorts them, which puts `"control"` before `"treated"` and `0`
#'   before `1`. Naming it is also what tells this function that a numeric column
#'   of zeroes and ones is two classes rather than two numbers.
#' @param control_label The reference class on its own, for when the other one
#'   needs no saying. Defaults to `outcome_lv[1]`, so a call that names one of the
#'   two names the reference either way; naming both and disagreeing is an error.
#' @param model What is fitted inside the search: `"linear"` for a continuous
#'   outcome, `"logistic"` for a two-class one, or `"rf"` for a random forest over
#'   either.
#' @param subset_sizes Subset sizes to score, or `NULL` for a ladder that is dense
#'   where the answer usually is — every size up to ten, then 15, 20, 30, 50 and
#'   100 — capped at the number of candidates. The full set is always scored, since
#'   keeping everything is the option a selection is being compared against.
#' @param metric Which resampled number chooses the size: `"RMSE"`, `"Rsquared"` or
#'   `"MAE"` for a regression, `"Accuracy"` or `"Kappa"` for a classification, or
#'   `NULL` for the first of each. Whether it is maximised or minimised follows
#'   from the metric and is reported in `parameters`.
#' @param ntree,nodesize Trees to grow and the smallest leaf to split, used by
#'   `model = "rf"` and ignored otherwise. `NULL` leaves `nodesize` at 1 for a
#'   classification and 5 for a regression, which is `randomForest()`'s own rule.
#' @param cv_method Resampling scheme: `"repeated_kfold"` for `n_repeat` runs of
#'   `n_fold`-fold cross-validation, `"kfold"` for one, or `"loocv"` for
#'   leave-one-out.
#' @param n_fold Folds per run, used by `"repeated_kfold"` and `"kfold"`.
#' @param n_repeat Number of runs, used by `"repeated_kfold"`.
#' @param seed Seed for the fold assignment, or `NULL` to use the stream as it
#'   stands. Supplying one does not disturb the caller: the previous random number
#'   state is put back when the function returns.
#'
#' @return An object of class `sa_selection`, a plain list.
#'
#'   \describe{
#'     \item{`analysis`}{`"rfe"`.}
#'     \item{`candidates`}{The predictors that were offered, most important first.
#'       This is the row order `ranking` follows, and what `terms` is to a model.}
#'     \item{`design`}{What the search saw: the `outcome` and its
#'       `outcome_type`, `outcome_lv` with `n_events` and `event_rate` for a
#'       classification, the row counts `n_obs`, `n_used` and `n_dropped`, the
#'       `predictors` in the order they arrived, any `dropped_predictors`, and
#'       `predictor_lv` for those that are factors.}
#'     \item{`parameters`}{The choices as they were used: `model`, `metric` and
#'       `maximize`, the forest's `ntree` and `nodesize` when there was one, the
#'       resampling scheme with `NA` where it used none of an argument, and
#'       `seed`.}
#'     \item{`selected`}{The predictors of the winning size, most important first.
#'       These are the names to hand to `predictors =` in a `fit_*()` call.}
#'     \item{`ranking`}{One row per candidate: `candidates`, the `estimate` it was
#'       ranked by averaged over the resamples, its `rank`, and whether it was
#'       `selected`.}
#'     \item{`profile`}{One row per subset size that was scored: `n_vars`, one
#'       column per metric with its standard deviation over the resamples, and
#'       `chosen`, which is `TRUE` on exactly one row.}
#'     \item{`resampling`}{One row per resample at the chosen size.}
#'     \item{`engine`}{What ran the search, including `importance`, the name of
#'       what `ranking$estimate` measures.}
#'     \item{`fit`}{The [caret::rfe()] object. This is the slot that is not
#'       portable; dropping it leaves an object that writes out as JSON.}
#'     \item{`metadata`}{Package version, R version, platform and timestamp.}
#'   }
#'
#' @seealso [split_data()], which defines the rows a selection should be run on,
#'   [fit_rf()] and [fit_linear_regression()] for fitting the predictors it kept,
#'   and [fit_elastic_net()], which answers the same question from the other end by
#'   shrinking a coefficient to exactly zero rather than by dropping a column.
#'
#' @examples
#' ## Linear RFE on a small candidate set (fast enough for examples).
#' res <- perform_rfe(mtcars, outcome = "mpg",
#'                    predictors = c("wt", "hp", "disp", "qsec"),
#'                    cv_method = "kfold", n_fold = 3, seed = 1)
#' res$selected
#'
#' \donttest{
#' res_full <- perform_rfe(mtcars, outcome = "mpg",
#'                         predictors = c("wt", "hp", "disp", "qsec", "drat",
#'                                        "carb"),
#'                         cv_method = "kfold", n_fold = 3, seed = 1)
#' fit_linear_regression(mtcars, outcome = "mpg", predictors = res_full$selected,
#'                       cv = FALSE)$coefficients
#'
#' ## Forest ranking on a two-class outcome.
#' iris2 <- iris[iris$Species != "setosa", ]
#' perform_rfe(iris2, outcome = "Species", control_label = "versicolor",
#'             model = "rf", ntree = 100, cv_method = "kfold", n_fold = 3,
#'             seed = 1)
#' }
#'
#' @export
perform_rfe <- function(data,
                        outcome,
                        predictors = NULL,
                        outcome_lv = NULL,
                        control_label = outcome_lv[1],
                        model = c("linear", "logistic", "rf"),
                        subset_sizes = NULL,
                        metric = NULL,
                        ntree = 500,
                        nodesize = NULL,
                        cv_method = c("repeated_kfold", "kfold", "loocv"),
                        n_fold = 5,
                        n_repeat = 5,
                        seed = NULL) {

  model <- match.arg(model)
  cv_method <- match.arg(cv_method)
  ntree <- sa_check_count(ntree, "ntree", 1)

  input <- sa_resolve_model_input(data, outcome, predictors)

  # `model` is the first say and the outcome is the second, which is where this
  # parts company with `fit_rf()`: there one function fits either kind and the
  # outcome decides alone, while here the name of a model is already an answer to
  # the question the outcome would have been asked.
  classify <- model == "logistic" || !is.null(outcome_lv) ||
    !is.null(control_label) || !is.numeric(input$y)

  if (model == "linear" && classify) {
    stop("`model = \"linear\"` ranks by the coefficients of a straight line ",
         "through a number, and `outcome` is a set of class labels. Use ",
         "`model = \"logistic\"` or `model = \"rf\"` for a two-class outcome.",
         call. = FALSE)
  }
  if (model == "logistic" && is.numeric(input$y) &&
        length(unique(input$y)) > 2L) {
    stop("`model = \"logistic\"` classifies two classes, and `outcome` is a ",
         "numeric column taking ", length(unique(input$y)),
         " values. Use `model = \"linear\"` or `model = \"rf\"` for a ",
         "continuous outcome.", call. = FALSE)
  }
  if (!classify && length(unique(input$y)) == 2L) {
    message("`outcome` is numeric and takes two values, so it was searched as ",
            "a regression. Pass `control_label`, or a factor column, with ",
            "`model = \"logistic\"` or `model = \"rf\"` to treat it as a ",
            "classification.")
  }

  if (classify) {
    y <- sa_outcome_levels(input$y, outcome_lv, control_label,
                           model = "a recursive feature elimination")
    outcome_lv <- levels(y)
  } else {
    if (!all(is.finite(input$y))) {
      stop("`outcome` holds non-finite value(s), which a model fitted inside ",
           "the search cannot be scored against.", call. = FALSE)
    }
    y <- input$y
  }

  # Validated whether or not the chosen model reads them, for the reason
  # `sa_train_control()` validates the folds of a scheme with no folds.
  if (is.null(nodesize)) {
    nodesize <- if (classify) 1L else 5L
  }
  nodesize <- sa_check_count(nodesize, "nodesize", 1)

  p <- length(input$predictors)
  sizes <- sa_rfe_sizes(subset_sizes, p)
  scoring <- sa_rfe_metric(metric, classify)
  funcs <- sa_rfe_funcs(model, classify, outcome_lv, ntree, nodesize)
  ctrl <- sa_rfe_control(funcs, cv_method, n_fold, n_repeat, input$n_used)

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  label <- sa_search_label(model, classify)
  fit <- sa_quiet_engine(
    caret::rfe(x = input$x, y = y, sizes = sizes, metric = scoring$metric,
               maximize = scoring$maximize, rfeControl = ctrl$control),
    label
  )

  ranking <- sa_rfe_ranking(fit, input$predictors)
  profile <- sa_rfe_profile(fit)

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

  sa_new_selection(
    analysis   = "rfe",
    candidates = ranking$candidates,
    design     = c(design, list(
      n_obs              = input$n_obs,
      n_used             = input$n_used,
      n_dropped          = input$n_dropped,
      predictors         = input$predictors,
      dropped_predictors = input$dropped_predictors
    ), sa_design_lv(input$predictor_lv)),
    # No `subset_sizes` here: `profile` holds one row per size that was scored, so
    # recording the ladder as well would say the same thing twice and leave two
    # places for it to be wrong. This is the rule `fit_rf()` follows with `mtry`.
    parameters = c(
      list(model = model, metric = scoring$metric, maximize = scoring$maximize),
      if (model == "rf") list(ntree = ntree, nodesize = nodesize),
      list(cv_method = ctrl$cv_method,
           n_fold    = ctrl$n_fold,
           n_repeat  = ctrl$n_repeat,
           seed      = seed)
    ),
    selected   = fit$optVariables,
    ranking    = ranking,
    profile    = profile,
    resampling = sa_model_frame(fit$resample),
    engine     = list(
      package    = "caret",
      method     = "rfe",
      label      = label,
      metrics    = fit$perfNames,
      importance = sa_rfe_importance_label(model)
    ),
    fit = fit
  )
}


#' The subset sizes to score, and the ladder that is scored by default
#'
#' Dense at the small end because that is where the answer usually is: the
#' difference between three predictors and four is a different model, while the
#' difference between sixty and seventy is the same model with noise in it. The
#' ladder is capped at the number of candidates rather than reaching past it, and
#' `caret` scores the full set whichever sizes are asked for, since keeping
#' everything is what a selection is being compared against.
#'
#' @param subset_sizes Sizes as requested, or `NULL` for the ladder.
#' @param p Number of candidates.
#'
#' @return Increasing integer vector.
#'
#' @keywords internal
#' @noRd
sa_rfe_sizes <- function(subset_sizes, p) {
  if (is.null(subset_sizes)) {
    return(sort(unique(pmin(p, c(1:10, 15L, 20L, 30L, 50L, 100L)))))
  }

  sa_check_num_vector(subset_sizes, "subset_sizes", 1, p)
  fractional <- unique(subset_sizes[subset_sizes != trunc(subset_sizes)])
  if (length(fractional) > 0L) {
    stop("`subset_sizes` counts predictors, so it must hold whole numbers, ",
         "but holds ", paste(fractional, collapse = ", "), ".", call. = FALSE)
  }
  sort(unique(as.integer(subset_sizes)))
}


#' The metric the size is chosen by, and which way it is read
#'
#' A metric belongs to a kind of outcome: there is no accuracy of a continuous
#' prediction and no root mean squared error of a class label. Naming one from the
#' other list is refused here rather than by `caret`, which scores the folds first
#' and then reports that the metric it was asked for is not among the columns it
#' produced.
#'
#' `maximize` follows from the metric and is not an argument. An error is better
#' when it is smaller and a rate is better when it is larger, and letting a caller
#' say otherwise would let them ask for the worst subset by accident.
#'
#' @param metric Metric name as requested, or `NULL` for the first of the list.
#' @param classify Whether the outcome is being classified.
#'
#' @return List with `metric` and `maximize`.
#'
#' @keywords internal
#' @noRd
sa_rfe_metric <- function(metric, classify) {
  available <- if (classify) {
    c("Accuracy", "Kappa")
  } else {
    c("RMSE", "Rsquared", "MAE")
  }

  if (is.null(metric)) {
    metric <- available[1L]
  } else if (!is.character(metric) || length(metric) != 1L ||
               !metric %in% available) {
    stop("`metric` must be one of ", paste(available, collapse = ", "),
         " for a ", if (classify) "classification" else "regression", ".",
         call. = FALSE)
  }

  list(metric = metric, maximize = !metric %in% c("RMSE", "MAE"))
}


#' What the search fits and ranks with at each step
#'
#' `caret`'s own `functions` lists are not used as they are, for the two reasons
#' at the head of this file: a ranking by raw coefficient size ranks by the units
#' the predictors were measured in, and a ranking of dummy columns cannot be
#' matched back to the columns that were passed in. What is kept from `caret` is
#' the pair that decides the answer once the ranking exists — the best size and the
#' variables of that size — since those are choices about the search rather than
#' about the model, and they are the same for all three.
#'
#' @param model Which model runs inside the search, already resolved.
#' @param classify Whether the outcome is being classified.
#' @param outcome_lv The two classes, reference first, or `NULL` for a regression.
#'   Held so that a logistic prediction comes back as a label rather than as a
#'   probability that would have to guess which class it belongs to.
#' @param ntree,nodesize Forest arguments, ignored by the other two models.
#'
#' @return A `functions` list for [caret::rfeControl()].
#'
#' @keywords internal
#' @noRd
sa_rfe_funcs <- function(model, classify, outcome_lv, ntree, nodesize) {
  force(outcome_lv)
  force(ntree)
  force(nodesize)

  picks <- list(
    summary    = caret::defaultSummary,
    selectSize = caret::lmFuncs$selectSize,
    selectVar  = caret::lmFuncs$selectVar
  )

  specific <- switch(
    model,
    linear = list(
      fit = function(x, y, first, last, ...) {
        stats::lm(sa_rfe_formula(), data = sa_search_frame(x, y))
      },
      pred = function(object, x) {
        unname(stats::predict(object, newdata = sa_search_frame(x)))
      },
      rank = function(object, x, y) sa_rfe_wald_rank(object)
    ),
    logistic = list(
      fit = function(x, y, first, last, ...) {
        stats::glm(sa_rfe_formula(), data = sa_search_frame(x, y),
                   family = stats::binomial())
      },
      # `glm()` models the probability of the last level, which is the direction
      # rule the whole package follows, so the cut is at the second level and
      # nothing is reversed afterwards.
      pred = function(object, x) {
        prob <- stats::predict(object, newdata = sa_search_frame(x),
                               type = "response")
        factor(ifelse(prob > 0.5, outcome_lv[2L], outcome_lv[1L]),
               levels = outcome_lv)
      },
      rank = function(object, x, y) sa_rfe_wald_rank(object)
    ),
    rf = list(
      # `importance = first` is `caret`'s own arrangement and is what makes the
      # permutation measure available on the one fit that is ranked. The forests
      # grown at the smaller sizes are scored rather than ranked, so computing it
      # for them would be paid for and thrown away.
      fit = function(x, y, first, last, ...) {
        randomForest::randomForest(x, y, importance = first, ntree = ntree,
                                   nodesize = nodesize)
      },
      pred = function(object, x) stats::predict(object, x),
      rank = function(object, x, y) sa_rfe_forest_rank(object, classify)
    ),
    stop("internal error: unhandled `model` ", model, ".", call. = FALSE)
  )

  c(specific, picks)
}


#' The formula the search fits by, built rather than written
#'
#' `.outcome ~ .` written out is a symbol that exists only inside a model frame,
#' which is the one thing `R CMD check` cannot tell from a typo. Building it from
#' a string says the same thing to `lm()` and nothing at all to the checker.
#'
#' @keywords internal
#' @noRd
sa_rfe_formula <- function() {
  stats::as.formula(".outcome ~ .")
}


#' Rank the columns of a linear or logistic fit by their own statistic
#'
#' The statistic rather than the estimate, since the estimate carries the units of
#' its predictor and the statistic has divided them out. The column rather than the
#' coefficient, since a factor is several coefficients and only the column can be
#' eliminated; the largest statistic among a factor's levels stands for it, so a
#' factor is kept as long as one of its levels is worth keeping.
#'
#' A term the fit could not estimate ranks at zero rather than at `NA`. It is a
#' column the others already span, so dropping it costs nothing, which is exactly
#' what a rank of zero says.
#'
#' @param object The `lm` or `glm` fitted inside one step of the search.
#'
#' @return data.frame with `var` and `Overall`, most important first, in the shape
#'   [caret::rfeControl()] expects of a ranking.
#'
#' @keywords internal
#' @noRd
sa_rfe_wald_rank <- function(object) {
  estimate <- stats::coef(object)
  summ <- summary(object)$coefficients
  statistic <- abs(as.numeric(summ[match(names(estimate), rownames(summ)), 3L]))
  statistic[is.na(statistic)] <- 0

  # `assign` maps each column of the model matrix to the term it came from, 0
  # being the intercept, which is the only way back from a dummy column to the
  # factor that produced it.
  assign <- attr(stats::model.matrix(object), "assign")
  labels <- attr(stats::terms(object), "term.labels")
  overall <- vapply(seq_along(labels), function(i) {
    at <- statistic[assign == i]
    if (length(at) == 0L) 0 else max(at)
  }, numeric(1))

  out <- data.frame(var = labels, Overall = overall, stringsAsFactors = FALSE)
  out <- out[order(out$Overall, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- out$var
  out
}


#' Rank the columns of a forest by what shuffling them costs
#'
#' The same measure `fit_rf()` reports as `estimate`, read off the same matrix by
#' the same helper, so that a forest's ranking here and its importance table there
#' are the same number and can be read together.
#'
#' @param object The `randomForest` fitted inside one step of the search.
#' @param classify Whether the outcome was classified, which decides which of the
#'   four importance columns is the permutation one.
#'
#' @return data.frame with `var` and `Overall`, most important first.
#'
#' @keywords internal
#' @noRd
sa_rfe_forest_rank <- function(object, classify) {
  importance <- sa_rf_importance(object, classify)
  out <- data.frame(var = importance$terms, Overall = importance$estimate,
                    stringsAsFactors = FALSE)
  rownames(out) <- out$var
  out
}


#' The ranking table, one row per candidate
#'
#' Averaged over the resamples the way `caret` averages them when it picks the
#' variables of the winning size, so that the top of this table and `$selected` are
#' the same predictors in the same order rather than two answers that usually
#' agree. Ties are broken by name, which is what makes a search of a data frame of
#' constant columns reproducible.
#'
#' @param fit The `rfe` object.
#' @param predictors The candidates, in the order they arrived.
#'
#' @return data.frame with `candidates`, `estimate`, `rank` and `selected`.
#'
#' @keywords internal
#' @noRd
sa_rfe_ranking <- function(fit, predictors) {
  averaged <- stats::aggregate(fit$variables["Overall"],
                               by = list(var = as.character(fit$variables$var)),
                               FUN = mean, na.rm = TRUE)

  # Alphabetical first, so that the stable sort below leaves ties in a defined
  # order; a candidate no fold ever ranked has no importance rather than a zero
  # one, and sorts last.
  candidates <- sort(predictors)
  estimate <- averaged$Overall[match(candidates, averaged$var)]
  at <- order(estimate, decreasing = TRUE, na.last = TRUE)

  out <- data.frame(
    candidates = candidates[at],
    estimate   = estimate[at],
    rank       = seq_along(candidates),
    selected   = candidates[at] %in% fit$optVariables,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}


#' The profile table, one row per subset size that was scored
#'
#' `Variables` is renamed because the count is of predictors rather than of
#' variables in the sense the rest of the package uses the word, and because
#' `n_vars` reads as the count it is beside the other `n_` fields of a result.
#'
#' @param fit The `rfe` object.
#'
#' @return data.frame with `n_vars`, the metrics with their standard deviations,
#'   and `chosen`.
#'
#' @keywords internal
#' @noRd
sa_rfe_profile <- function(fit) {
  out <- fit$results
  names(out)[names(out) == "Variables"] <- "n_vars"
  out$chosen <- out$n_vars == fit$bestSubset
  rownames(out) <- NULL
  out
}


#' What `ranking$estimate` measures, for `engine$importance`
#'
#' @param model Which model did the ranking, already resolved.
#'
#' @keywords internal
#' @noRd
sa_rfe_importance_label <- function(model) {
  switch(
    model,
    linear = "absolute t statistic",
    logistic = "absolute Wald z",
    rf = "permutation importance"
  )
}
