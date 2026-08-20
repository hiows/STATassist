# The penalized counterpart of `fit_linear_regression()` and
# `fit_logistic_regression()`, and the first model here that is tuned rather than
# merely scored. Two things follow from the penalty and account for most of what
# is different below. The estimates are shrunk, so there is no standard error and
# no p-value to report beside them, and the amount of shrinkage is a choice, so
# the resampling that used to score the model now also picks it.
#
# One outcome argument stands for two engines. A continuous outcome is a Gaussian
# elastic net and a two-class one is a binomial elastic net; `glmnet` fits both
# and the direction rule of the classification is the one `glm()` already follows,
# so the two share every line except the family and the fit statistics.

#' Fit an elastic net, lasso or ridge regression
#'
#' Fits a penalized linear model of one continuous or two-class outcome on a set
#' of predictors, choosing how much to penalize by cross-validation, and scores
#' the chosen model on the same folds. Which outcome it is decides which model is
#' fitted: a numeric outcome is a linear regression and a two-class one is a
#' logistic regression, both with the elastic net penalty on their coefficients.
#'
#' The penalty is what makes this different from [fit_linear_regression()]. It
#' pulls the coefficients towards zero and, when `alpha` is above zero, sets some
#' of them exactly to zero, so the fit selects predictors as well as estimating
#' them. `selected` in the coefficient table is that answer.
#'
#' @details
#' `penalty` names a corner of one model rather than three different models.
#' `alpha` is the weight on the L1 penalty against the L2 one, so `"lasso"` is
#' `alpha = 1`, `"ridge"` is `alpha = 0`, and `"elastic_net"` is everything
#' between, `alpha` then being tuned like `lambda`. The first two ignore the
#' `alpha` argument, and `parameters$alpha` reports what was used either way.
#'
#' \describe{
#'   \item{Cross-validation chooses the model here, and does not merely score
#'     it}{This is the one place where this family departs from
#'     [fit_linear_regression()], where `cv` decides nothing about the fit.
#'     `lambda` and `alpha` are chosen by the resampled metric, which is `RMSE`
#'     for a regression and `Accuracy` for a classification, and the final model
#'     is then fitted on all usable rows at that pair. `performance` therefore
#'     holds one row per candidate, and `parameters$alpha` and
#'     `parameters$lambda` are the pair that won. With `cv = FALSE` there is
#'     nothing to choose between, so the grid must name exactly one candidate.}
#'   \item{The estimates have no standard error, and the table says so by not
#'     having the column}{A penalized estimate is deliberately biased, and the
#'     usual standard error assumes an unbiased one, so there is no honest number
#'     to put under `stderr`, `statistic`, `df`, `pval` or the confidence limits.
#'     The table carries only what a penalized fit answers — `estimate` and
#'     `selected` — rather than those six columns filled with `NA` down their
#'     whole length, which would read as a table that lost its values instead of a
#'     model that never had them. `is.null(fit$coefficients$pval)` is how a
#'     consumer tells the two kinds of table apart. Use
#'     [fit_linear_regression()] or [fit_logistic_regression()] when the p-value
#'     is the point.}
#'   \item{The predictors are penalized on a standardised scale}{`glmnet`
#'     standardises each column before penalizing it, since otherwise a predictor
#'     measured in millimetres would be shrunk less than the same predictor in
#'     metres, and returns the coefficients on their original scale. They are
#'     therefore comparable with each other in what they mean per unit, but not
#'     with an unpenalized coefficient of the same predictor, which was not
#'     shrunk at all.}
#'   \item{Terms are not predictors}{As in the unpenalized models, a factor or
#'     character predictor with `k` levels becomes `k - 1` terms named after the
#'     levels, and the coding is the one [stats::lm()] would have used, so the
#'     two coefficient tables can be read side by side. Each of those terms is
#'     penalized on its own, so a factor can have one level selected and another
#'     dropped.}
#' }
#'
#' Rows with a missing value in the outcome or in any predictor are dropped
#' before the folds are drawn rather than inside each fit, and
#' `design$n_dropped` reports how many went. A predictor that takes a single
#' value is left out with a message, which is one of the two ways a call can end
#' up with fewer terms than it named.
#'
#' Two terms are the fewest that can be fitted. A penalty divides a budget
#' between coefficients, and `glmnet` refuses a design matrix of one column
#' outright, so a single predictor is an error here rather than a model with
#' nothing to trade off. [fit_linear_regression()] and
#' [fit_logistic_regression()] fit one predictor as readily as ten.
#'
#' For a two-class outcome the direction is `outcome_lv`, read exactly as
#' [fit_logistic_regression()] reads it: the first level is the reference, every
#' coefficient is the change in the log odds of `outcome_lv[2]`, `odds_ratio` is
#' above 1 for a predictor that raises the chance of it, and
#' `predict(model, newdata, type = "response")` is its probability. That is
#' also what the engine does unaided, `glmnet` modelling the last level of a
#' factor as `glm()` does. `control_label` names the reference on its own and is
#' enough by itself to make a column of zeroes and ones a classification; naming
#' it alongside an `outcome_lv` that puts the other class first is an error, as
#' it is in [fit_logistic_regression()].
#'
#' New rows are predicted through the result rather than through `$fit`, and on
#' this model that is the only route: `glmnet` was handed a design matrix and
#' reads one by position, so the frame the fit was given cannot be handed to it
#' again. [predict.sa_model()] rebuilds the terms from the levels the fit recorded
#' and puts them in the model's order by name.
#'
#' A note on `$fit`: `caret` fits the whole `lambda` path and records the chosen
#' value on it, so `coef()` and `predict()` on `$fit` interpolate the path at
#' that value rather than reading a model fitted at it alone. The two agree to
#' within the resolution of the path, and everything the result reports comes
#' from the same interpolation, so the object is consistent with itself.
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
#' @param outcome_lv The two classes, reference first, so that the coefficients
#'   describe the odds of the second one. Supplying it also states that the
#'   outcome is to be classified, which is the one thing a numeric column holding
#'   two values cannot say on its own. `NULL` sorts the classes of an outcome
#'   that is not numeric, and leaves a numeric one to the regression.
#' @param control_label The reference class on its own, for when the other one
#'   needs no saying. Defaults to `outcome_lv[1]`, so a call that names one of the
#'   two names the reference either way; naming both and disagreeing is an error.
#' @param penalty Which penalty to apply: `"elastic_net"` for a tuned mixture of
#'   the L1 and L2 penalties, `"lasso"` for the L1 penalty alone (`alpha = 1`),
#'   or `"ridge"` for the L2 penalty alone (`alpha = 0`). A lasso and an elastic
#'   net can set a coefficient to zero; a ridge only shrinks.
#' @param alpha Mixing weights to search, between 0 and 1, read only when
#'   `penalty = "elastic_net"`.
#' @param lambda Penalty sizes to search. Larger shrinks more, and 0 is no
#'   penalty at all.
#' @param cv Whether to cross-validate. `FALSE` fits the single candidate the
#'   grid names, once, and reports no resampled performance.
#' @param cv_method Resampling scheme: `"repeated_kfold"` for `n_repeat` runs of
#'   `n_fold`-fold cross-validation, `"kfold"` for one, or `"loocv"` for
#'   leave-one-out.
#' @param n_fold Folds per run, used by `"repeated_kfold"` and `"kfold"`.
#' @param n_repeat Number of runs, used by `"repeated_kfold"`.
#' @param seed Seed for the fold assignment, or `NULL` to use the stream as it
#'   stands. Supplying one does not disturb the caller: the previous random number
#'   state is put back when the function returns.
#'
#' @return An object of class `sa_model`, the same eleven elements
#'   [fit_linear_regression()] returns, with these differences:
#'
#'   \describe{
#'     \item{`analysis`}{`"elastic_net"`, whichever corner `penalty` named.}
#'     \item{`design`}{Holds `outcome_lv`, `n_events` and `event_rate` for a
#'       two-class outcome, as [fit_logistic_regression()] does, and neither for
#'       a continuous one.}
#'     \item{`parameters`}{Holds `penalty`, the `alpha` and `lambda` that were
#'       chosen rather than the grids that were searched, and `n_candidates`, how
#'       many pairs were scored. The grids themselves are the rows of
#'       `performance`.}
#'     \item{`coefficients`}{`estimate` at the chosen penalty and `selected`,
#'       whether the term survived it. The intercept is never penalized, so it is
#'       always selected. The inference columns are absent rather than `NA`; a
#'       two-class outcome also gets `odds_ratio`, the exponentiated estimate,
#'       but no interval to go with it.}
#'     \item{`fit_stats`}{Measured on the rows the model was fitted to:
#'       `r_squared`, `rmse` and `mae` for a regression, and `null_deviance`,
#'       `residual_deviance` and `mcfadden_r2` for a classification, each with
#'       `n_selected` and `n_zero`, the terms the penalty kept and dropped. These
#'       describe the fit rather than what it would do next; `performance` is the
#'       one that was measured on rows the model had not seen.}
#'     \item{`performance`}{One row per candidate, the chosen one being the row
#'       that matches `parameters$alpha` and `parameters$lambda`, or `NULL` when
#'       `cv = FALSE`.}
#'   }
#'
#' @seealso [fit_linear_regression()] and [fit_logistic_regression()] for the
#'   unpenalized fits and the inference they can report, [split_data()], which
#'   defines the rows this is fitted on, [coef.sa_model()] for the coefficient
#'   table, [predict.sa_model()] for predicting the rows it was not fitted on, and
#'   [coef.sa_fit()] for the methods `$fit` answers to, among them the estimates
#'   as a named vector.
#'
#' @examples
#' ## A single lasso without resampling (fast enough for examples).
#' fit <- fit_elastic_net(mtcars, outcome = "mpg", penalty = "lasso",
#'                        lambda = 0.5, cv = FALSE)
#' fit$parameters[c("penalty", "alpha", "lambda")]
#' fit$coefficients[c("terms", "estimate", "selected")]
#'
#' \donttest{
#' ## Cross-validated grid search, ridge vs hard lasso, and known-truth scoring.
#' fit_elastic_net(mtcars, outcome = "mpg", penalty = "lasso",
#'                 lambda = c(0.01, 0.1, 0.5, 1, 2),
#'                 cv_method = "kfold", n_fold = 5, seed = 1)
#' hard <- fit_elastic_net(mtcars, outcome = "mpg", penalty = "lasso",
#'                         lambda = 2, cv = FALSE)
#' ridge <- fit_elastic_net(mtcars, outcome = "mpg", penalty = "ridge",
#'                          lambda = 2, cv = FALSE)
#' c(lasso = hard$fit_stats$n_selected, ridge = ridge$fit_stats$n_selected)
#' sim <- simulate_regression(seed = 1)
#' lasso <- do.call(fit_elastic_net,
#'                  c(sim$args, list(penalty = "lasso", lambda = 0.5,
#'                                   cv = FALSE)))
#' scored <- merge(lasso$coefficients, sim$truth_term, by = "terms")
#' table(planted = scored$beta != 0, selected = scored$selected)
#'
#' iris2 <- iris[iris$Species != "setosa", ]
#' clf <- fit_elastic_net(iris2, outcome = "Species",
#'                        outcome_lv = c("versicolor", "virginica"),
#'                        penalty = "lasso", lambda = c(0.001, 0.01, 0.1),
#'                        cv_method = "kfold", n_fold = 5, seed = 1)
#' clf$coefficients[c("terms", "estimate", "odds_ratio", "selected")]
#' }
#'
#' @export
fit_elastic_net <- function(data,
                            outcome,
                            predictors = NULL,
                            outcome_lv = NULL,
                            control_label = outcome_lv[1],
                            penalty = c("elastic_net", "lasso", "ridge"),
                            alpha = seq(0, 1, by = 0.1),
                            lambda = 10^seq(-4, 1, length.out = 50),
                            cv = TRUE,
                            cv_method = c("repeated_kfold", "kfold", "loocv"),
                            n_fold = 5,
                            n_repeat = 5,
                            seed = NULL) {

  penalty <- match.arg(penalty)
  cv_method <- match.arg(cv_method)

  input <- sa_resolve_model_input(data, outcome, predictors)
  grid <- sa_enet_grid(penalty, alpha, lambda, cv)
  ctrl <- sa_train_control(cv, cv_method, n_fold, n_repeat, input$n_used)

  # The predictor frame is dummy coded here rather than left to the engine: what
  # `caret` hands `glmnet` is `as.matrix()` of the frame, which turns a factor
  # into its integer codes and fits without complaint. See `sa_design_matrix()`.
  x <- sa_design_matrix(input$x)
  if (ncol(x) < 2L) {
    stop("a penalty divides its budget between terms, and `glmnet` refuses a ",
         "single column, but the model has ", ncol(x), " term(s): ",
         paste(colnames(x), collapse = ", "), ". Add a predictor, or use ",
         "fit_linear_regression() or fit_logistic_regression(), which fit one ",
         "predictor as readily as ten.", call. = FALSE)
  }

  # One argument, two models, and the outcome decides which. A numeric column is
  # a regression, since that is what a number usually is, and anything else is a
  # set of class labels. `outcome_lv` overrules the guess, and so does
  # `control_label` on its own, which is the only way to say that a column of
  # zeroes and ones is two classes rather than two numbers; the guess is
  # announced when it is the ambiguous one so that nobody has to infer it from
  # the output.
  classify <- !is.null(outcome_lv) || !is.null(control_label) ||
    !is.numeric(input$y)
  if (!classify && length(unique(input$y)) == 2L) {
    message("`outcome` is numeric and takes two values, so it was fitted as a ",
            "regression. Pass `outcome_lv` or `control_label`, or a factor ",
            "column, to model it as a classification.")
  }

  if (classify) {
    y <- sa_outcome_levels(input$y, outcome_lv, control_label,
                           model = "an elastic net")
    outcome_lv <- levels(y)
    family <- "binomial"
  } else {
    if (!all(is.finite(input$y))) {
      stop("`outcome` holds non-finite value(s), which least squares has no ",
           "residual for.", call. = FALSE)
    }
    y <- input$y
    family <- "gaussian"
  }

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  fit <- sa_fit_engine(
    caret::train(x = x, y = y, method = "glmnet", family = family,
                 trControl = ctrl$control, tuneGrid = grid),
    switch(penalty, lasso = "Lasso", ridge = "Ridge", "Elastic net")
  )
  model <- fit$finalModel
  chosen <- model$lambdaOpt

  estimate <- as.numeric(as.matrix(stats::coef(model, s = chosen)))
  terms <- rownames(stats::coef(model, s = chosen))
  intercept <- terms == "(Intercept)"

  coefs <- data.frame(
    terms      = terms,
    estimate   = estimate,
    # The intercept is not penalized, so it is not the penalty that would have
    # set it to zero and it counts as kept whatever its value.
    selected   = estimate != 0 | intercept,
    stringsAsFactors = FALSE
  )
  if (classify) {
    coefs$odds_ratio <- exp(coefs$estimate)
  }

  penalised <- coefs[!intercept, ]
  n_selected <- sum(penalised$selected)
  fitted_value <- as.numeric(stats::predict(model, newx = x, s = chosen,
                                            type = "response"))
  fit_stats <- c(
    if (classify) {
      sa_enet_binomial_stats(fitted_value, y == outcome_lv[2])
    } else {
      sa_enet_gaussian_stats(fitted_value, y)
    },
    list(n_selected = n_selected, n_zero = nrow(penalised) - n_selected)
  )
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
    analysis = "elastic_net",
    terms    = coefs$terms,
    design   = c(design, list(
      n_obs              = input$n_obs,
      n_used             = input$n_used,
      n_dropped          = input$n_dropped,
      predictors         = input$predictors,
      dropped_predictors = input$dropped_predictors
    ), sa_design_lv(input$predictor_lv)),
    # The penalty that ran, not the grid that was asked for: `performance` holds
    # every candidate, so recording the grid here as well would say the same
    # thing twice and leave two places for it to be wrong.
    parameters = list(
      penalty      = penalty,
      alpha        = unname(fit$bestTune$alpha),
      lambda       = unname(fit$bestTune$lambda),
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
      method  = "glmnet",
      family  = family,
      label   = paste(
        switch(penalty,
               lasso       = "Lasso (L1 penalty)",
               ridge       = "Ridge (L2 penalty)",
               elastic_net = "Elastic net (L1 and L2 penalties)"),
        if (classify) "binomial classification" else "linear regression"
      ),
      metrics = fit$perfNames,
      # `glmnet` was handed a design matrix and reads one by position, so the
      # columns it was given have to be rebuilt in that order to predict new rows.
      # These are `finalModel$xNames`, kept here so that `predict.sa_model()` can
      # ask the result rather than the engine object, which is the slot that does
      # not survive being written out.
      x_names = colnames(x)
    ),
    fit = fit
  )
}


#' Goodness of fit of a penalized model on the rows it was fitted to
#'
#' The unpenalized models take these from `summary()`, which a `glmnet` object
#' has none of: it holds a whole path of fits and no inference about any of them.
#' They are computed from the fitted values at the chosen penalty instead, which
#' is the same arithmetic `summary.lm()` does and the same deviance
#' `summary.glm()` reports, minus the degrees of freedom. A penalized fit has no
#' honest count of those, which is why no adjusted R-squared, `AIC` or `BIC`
#' appears beside them.
#'
#' @param fitted_value Fitted outcome values, or fitted probabilities of the
#'   modelled class.
#' @param observed The outcome, or a logical vector that is `TRUE` on the
#'   modelled class.
#'
#' @return Named list of scalars.
#'
#' @keywords internal
#' @noRd
sa_enet_gaussian_stats <- function(fitted_value, observed) {
  residual <- observed - fitted_value
  sse <- sum(residual^2)
  sst <- sum((observed - mean(observed))^2)

  list(
    r_squared = if (sst > 0) 1 - sse / sst else NA_real_,
    rmse      = sqrt(mean(residual^2)),
    mae       = mean(abs(residual))
  )
}


#' @rdname sa_enet_gaussian_stats
#' @keywords internal
#' @noRd
sa_enet_binomial_stats <- function(fitted_value, observed) {
  # A probability of exactly 0 or 1 makes the deviance infinite, which reads as a
  # failure rather than as the perfect in-sample separation it is. It is held off
  # the boundary instead, so the deviance stays finite and large.
  eps <- .Machine$double.eps
  p <- pmin(pmax(fitted_value, eps), 1 - eps)
  event <- as.numeric(observed)

  deviance <- -2 * sum(event * log(p) + (1 - event) * log(1 - p))
  rate <- mean(event)
  null_deviance <- if (rate > 0 && rate < 1) {
    -2 * length(event) * (rate * log(rate) + (1 - rate) * log(1 - rate))
  } else {
    NA_real_
  }

  list(
    null_deviance     = null_deviance,
    residual_deviance = deviance,
    mcfadden_r2       = if (is.na(null_deviance) || null_deviance == 0) {
      NA_real_
    } else {
      1 - deviance / null_deviance
    }
  )
}
