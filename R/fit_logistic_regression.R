# Binomial logistic regression, the classification counterpart of
# `fit_linear_regression()`. The two are deliberately near-identical below the
# roxygen: the same input resolution, the same resampling control, the same
# coefficient table. What is specific to this one is the direction rule. A
# classification has a class of interest, and which of the two levels that is
# decides the sign of every coefficient in the table.

#' Fit a logistic regression
#'
#' Fits a binomial logistic regression of one two-class outcome on a set of
#' predictors, and scores it by cross-validation on the data it was fitted to. The
#' coefficient table describes the fit on every usable row, while `performance`
#' describes how the same procedure did on rows it had not seen inside each fold.
#'
#' @details
#' `outcome_lv` fixes the direction, and it does so by the rule the rest of the
#' package follows: the first level is the reference. Every coefficient is the
#' change in the log odds of `outcome_lv[2]` per unit of its predictor, and
#' `odds_ratio` is above 1 for a predictor that raises the chance of it. With
#' `outcome_lv = c("control", "case")` the table therefore reads as a statement
#' about `case`, and the same vector handed to [compare_two_groups()] as
#' `group_lv` would put `control` in the denominator of its fold change, so the
#' two point the same way.
#'
#' `control_label` states the same direction with one name instead of two, which
#' is the shorter thing to say when the other class needs no saying. Naming both
#' and pointing them at different classes is an error rather than a re-pointing:
#' `outcome_lv` holds the two classes and nothing else, so there is no reading
#' under which a different reference leaves any of it standing. The comparison
#' functions take the argument the other way round, since their `group_lv`
#' carries the display order too — see [compare_two_groups()].
#'
#' This is also what the engine does unaided. [stats::glm()] models the
#' probability of the last level of a factor, so ordering the levels
#' reference-first is the whole of the implementation and nothing is reversed
#' afterwards.
#'
#' The rule reaches the predictions too:
#' `predict(model, newdata, type = "response")` is the probability of
#' `outcome_lv[2]`, the same class the coefficients describe. See
#' [predict.sa_model()], and [coef.sa_fit()] for the same word on the engine
#' object, since `"response"` is one this package adds to the ones
#' [caret::train()] accepts.
#'
#' A third class is an error rather than a silently dropped set of rows: two
#' classes are what this model is, and quietly fitting a subset of the data that
#' was passed in would answer a question nobody asked. Reduce `data` to the two
#' classes first; naming two of three with `outcome_lv` is refused for the same
#' reason.
#'
#' The rest matches [fit_linear_regression()]: the model is fitted by
#' [caret::train()] with `method = "glm"`, `cv` decides whether the model is
#' scored rather than how it is fitted, a factor predictor with `k` levels becomes
#' `k - 1` terms, and rows missing anything the model needs are dropped before
#' the folds are drawn. The interval is the Wald interval, the one matching the z
#' statistic and standard error reported beside it, rather than the profile
#' likelihood interval [stats::confint()] would give.
#'
#' Perfect separation is reported rather than hidden. A predictor that splits the
#' two classes exactly gives an estimate that grows until the fit stops, with a
#' standard error to match, and the engine says so once per fold; those notes come
#' back as one message with a count.
#'
#' @param data A data.frame (or matrix) in wide format, one row per observation.
#'   Typically the training half of a [split_data()] result.
#' @param outcome The two-class outcome, either the name of a column of `data` or
#'   a vector with one entry per row. Factor, character, logical and numeric
#'   columns are all read as class labels.
#' @param predictors Column names to fit on, or `NULL` for every column of `data`
#'   except the outcome. Numeric, logical, factor and character columns are
#'   accepted; a column that takes a single value is left out with a message,
#'   since it cannot contribute.
#' @param outcome_lv The two classes, reference first, so that the coefficients
#'   describe the odds of the second one. `NULL` sorts the classes, which puts
#'   `"control"` before `"treated"` and `0` before `1`.
#' @param control_label The reference class on its own, for when the other one
#'   needs no saying. Defaults to `outcome_lv[1]`, so a call that names one of the
#'   two names the reference either way; naming both and disagreeing is an error.
#' @param cv Whether to cross-validate. `FALSE` fits the model once and reports no
#'   resampled performance.
#' @param cv_method Resampling scheme: `"repeated_kfold"` for `n_repeat` runs of
#'   `n_fold`-fold cross-validation, `"kfold"` for one, or `"loocv"` for
#'   leave-one-out.
#' @param n_fold Folds per run, used by `"repeated_kfold"` and `"kfold"`.
#' @param n_repeat Number of runs, used by `"repeated_kfold"`.
#' @param conf_level Confidence level of the coefficient intervals.
#' @param seed Seed for the fold assignment, or `NULL` to use the stream as it
#'   stands. Supplying one does not disturb the caller: the previous random number
#'   state is put back when the function returns.
#'
#' @return An object of class `sa_model`, the same eleven elements
#'   [fit_linear_regression()] returns, with these differences:
#'
#'   \describe{
#'     \item{`analysis`}{`"logistic_regression"`.}
#'     \item{`design`}{Also holds `outcome_lv`, and `n_events` with `event_rate`,
#'       the number and proportion of rows in `outcome_lv[2]`.}
#'     \item{`coefficients`}{`statistic` is a Wald z rather than a t and `df` is
#'       `NA`, since the z is not referred to any. `odds_ratio`,
#'       `or_lower_conf` and `or_upper_conf` are added, being the exponentiated
#'       estimate and its limits.}
#'     \item{`fit_stats`}{`null_deviance` and `residual_deviance` with their
#'       degrees of freedom, `mcfadden_r2`, the likelihood ratio test of the model
#'       against the intercept alone, `aic` and `bic`.}
#'     \item{`performance`}{Resampled `Accuracy` and `Kappa` rather than the
#'       regression metrics.}
#'   }
#'
#' @seealso [split_data()], which defines the rows this is fitted on,
#'   [fit_linear_regression()] for a continuous outcome, [coef.sa_model()] for
#'   the coefficient table, [predict.sa_model()] for predicting the rows it was
#'   not fitted on, and [coef.sa_fit()] for the methods `$fit` answers to, among
#'   them the estimates as a named vector.
#'
#' @examples
#' ## Two of the three iris species, so that the outcome has two classes.
#' ## `versicolor` first makes it the reference, so every odds ratio is the odds
#' ## of `virginica` rather than of the other way round.
#' iris2 <- iris[iris$Species != "setosa", ]
#' fit <- fit_logistic_regression(iris2, outcome = "Species",
#'                                predictors = c("Petal.Length", "Sepal.Width"),
#'                                outcome_lv = c("versicolor", "virginica"),
#'                                cv_method = "kfold", seed = 1)
#' fit
#' fit$coefficients[c("terms", "estimate", "odds_ratio", "pval")]
#' fit$fit_stats$mcfadden_r2
#'
#' ## Swapping the two levels turns every coefficient around and inverts every
#' ## odds ratio, which is the same rule `group_lv` follows in a comparison.
#' other_way <- fit_logistic_regression(
#'   iris2, outcome = "Species", predictors = c("Petal.Length", "Sepal.Width"),
#'   outcome_lv = c("virginica", "versicolor"), cv = FALSE
#' )
#' cbind(versicolor_ref = fit$coefficients$odds_ratio,
#'       virginica_ref  = other_way$coefficients$odds_ratio)
#'
#' ## The training half of a split is what a model is normally fitted on, and
#' ## stratifying on the outcome keeps both halves able to see both classes.
#' sp <- split_data(iris2, stratified = "Species", seed = 1)
#' train <- sp$datasets[[1]]$train_data
#' fit_logistic_regression(train, outcome = "Species",
#'                         predictors = "Petal.Width", cv = FALSE)$design$n_events
#'
#' @export
fit_logistic_regression <- function(data,
                                    outcome,
                                    predictors = NULL,
                                    outcome_lv = NULL,
                                    control_label = outcome_lv[1],
                                    cv = TRUE,
                                    cv_method = c("repeated_kfold", "kfold",
                                                  "loocv"),
                                    n_fold = 5,
                                    n_repeat = 5,
                                    conf_level = 0.95,
                                    seed = NULL) {

  cv_method <- match.arg(cv_method)
  sa_check_scalar_num(conf_level, "conf_level", 0, 1,
                      lower_open = TRUE, upper_open = TRUE)

  input <- sa_resolve_model_input(data, outcome, predictors)
  y <- sa_outcome_levels(input$y, outcome_lv, control_label)
  outcome_lv <- levels(y)
  ctrl <- sa_train_control(cv, cv_method, n_fold, n_repeat, input$n_used)

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  fit <- sa_fit_engine(
    caret::train(x = input$x, y = y, method = "glm",
                 family = stats::binomial(), trControl = ctrl$control),
    "Logistic regression"
  )
  model <- fit$finalModel
  summ <- summary(model)

  # No df column: a Wald z is referred to the normal distribution, so reporting
  # the residual degrees of freedom next to it would suggest a t test.
  coefs <- sa_coef_table(
    model,
    sa_wald_interval(summ$coefficients, conf_level, df = NULL),
    df = NA_real_
  )
  coefs$odds_ratio <- exp(coefs$estimate)
  coefs$or_lower_conf <- exp(coefs$lower_conf)
  coefs$or_upper_conf <- exp(coefs$upper_conf)

  aliased <- names(which(summ$aliased))
  if (length(aliased) > 0L) {
    warning("term(s) could not be estimated because the other predictors ",
            "already span them; their rows are NA: ",
            paste(aliased, collapse = ", "), ".", call. = FALSE)
  }

  lr_stat <- model$null.deviance - model$deviance
  lr_df <- model$df.null - model$df.residual
  fit_stats <- list(
    null_deviance     = model$null.deviance,
    residual_deviance = model$deviance,
    df_null           = model$df.null,
    df_residual       = model$df.residual,
    mcfadden_r2       = if (model$null.deviance > 0) {
      1 - model$deviance / model$null.deviance
    } else {
      NA_real_
    },
    lr_stat           = lr_stat,
    lr_df             = lr_df,
    lr_pval           = if (lr_df > 0) {
      stats::pchisq(lr_stat, lr_df, lower.tail = FALSE)
    } else {
      NA_real_
    },
    aic               = stats::AIC(model),
    bic               = stats::BIC(model)
  )
  fit_stats <- lapply(fit_stats, function(v) unname(as.numeric(v)))

  n_events <- sum(y == outcome_lv[2])
  sa_new_model(
    analysis = "logistic_regression",
    terms    = coefs$terms,
    design   = c(list(
      outcome            = input$outcome,
      outcome_type       = "two classes",
      outcome_lv         = outcome_lv,
      n_events           = n_events,
      event_rate         = n_events / input$n_used,
      n_obs              = input$n_obs,
      n_used             = input$n_used,
      n_dropped          = input$n_dropped,
      predictors         = input$predictors,
      dropped_predictors = input$dropped_predictors
    ), sa_design_lv(input$predictor_lv)),
    parameters = list(
      cv         = cv,
      cv_method  = ctrl$cv_method,
      n_fold     = ctrl$n_fold,
      n_repeat   = ctrl$n_repeat,
      conf_level = conf_level,
      seed       = seed
    ),
    coefficients = coefs,
    fit_stats    = fit_stats,
    performance  = sa_model_frame(fit$results),
    resampling   = sa_model_frame(fit$resample),
    engine       = list(
      package = "caret",
      method  = "glm",
      family  = "binomial",
      label   = "Binomial logistic regression",
      metrics = fit$perfNames
    ),
    fit = fit
  )
}
