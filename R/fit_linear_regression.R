# Ordinary least squares, fitted through `caret` rather than through `lm()`
# directly. On its own that looks like a detour, since the coefficients are the
# ones `lm()` would give either way. What `caret` adds is the resampling: the same
# `trControl` argument that scores this model scores the elastic net and the
# random forest, so the numbers in `performance` are comparable across models
# that have nothing else in common.

#' Fit a linear regression
#'
#' Fits an ordinary least squares model of one continuous outcome on a set of
#' predictors, and scores it by cross-validation on the data it was fitted to.
#' Both halves of that sentence matter: the coefficient table describes the fit on
#' every usable row, while `performance` describes how the same procedure did on
#' rows it had not seen inside each fold.
#'
#' Nothing is selected on the outcome here. `predictors` is the set the caller
#' names, and a predictor that turns out not to matter stays in the table with the
#' p-value that says so.
#'
#' @details
#' The model is fitted by [caret::train()] with `method = "lm"`, so the fit itself
#' is [stats::lm()]'s and the resampling is `caret`'s. Two consequences are worth
#' knowing:
#'
#' \describe{
#'   \item{Cross-validation does not change the model}{`cv` decides whether the
#'     model is scored, not how it is fitted. The final model is fitted on all
#'     usable rows either way, so the coefficients of `cv = TRUE` and
#'     `cv = FALSE` are identical and only `performance` and `resampling`
#'     differ.}
#'   \item{Terms are not predictors}{A factor or character predictor with `k`
#'     levels becomes `k - 1` terms, named after the level each one stands for.
#'     `terms` therefore holds the row order of `coefficients`, while
#'     `design$predictors` holds the columns that were read.}
#' }
#'
#' Rows with a missing value in the outcome or in any predictor are dropped before
#' the folds are drawn rather than inside each fit. Left to the engine, deletion
#' would happen once per fold on whatever that fold held, and the folds would then
#' be scored on different subsets of the data; `design$n_dropped` reports how many
#' rows went.
#'
#' The confidence interval is the t interval on the residual degrees of freedom,
#' the one that matches the t statistic and standard error reported beside it. A
#' term the fit could not estimate, because another predictor already spans it,
#' keeps its row with its estimate and inference `NA` and is named in a warning.
#' Dropping the row instead would make the table quietly shorter than the model
#' it describes.
#'
#' @param data A data.frame (or matrix) in wide format, one row per observation.
#'   Typically the training half of a [split_data()] result.
#' @param outcome The continuous outcome, either the name of a column of `data` or
#'   a vector with one entry per row.
#' @param predictors Column names to fit on, or `NULL` for every column of `data`
#'   except the outcome. Numeric, logical, factor and character columns are
#'   accepted; a column that takes a single value is left out with a message,
#'   since it cannot contribute.
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
#' @return An object of class `sa_model`, a plain list of eleven elements.
#'
#'   \describe{
#'     \item{`analysis`}{`"linear_regression"`.}
#'     \item{`terms`}{Coefficient term names, in the row order `coefficients`
#'       follows.}
#'     \item{`design`}{What the model saw: the `outcome` label, its type, the
#'       `predictors` kept and any `dropped_predictors`, and the row counts
#'       `n_obs`, `n_used` and `n_dropped`.}
#'     \item{`parameters`}{The fitting choices, with the resampling arguments as
#'       they were used rather than as they were passed: `n_repeat` is `NA` for a
#'       scheme that does not repeat, and both fold arguments are `NA` for
#'       leave-one-out.}
#'     \item{`coefficients`}{One row per term: `estimate`, `stderr`,
#'       `statistic` (t), `df`, `pval`, `lower_conf` and `upper_conf`.}
#'     \item{`fit_stats`}{The model as a whole: `r_squared`, `adj_r_squared`,
#'       `sigma`, the overall F test, `aic` and `bic`.}
#'     \item{`performance`}{Resampled `RMSE`, `Rsquared` and `MAE` with their
#'       standard deviations across resamples, or `NULL` when `cv = FALSE`.}
#'     \item{`resampling`}{One row per resample, or `NULL`.}
#'     \item{`engine`}{What fitted the model, and which metrics it scored.}
#'     \item{`fit`}{The [caret::train()] object, so that
#'       [predict.sa_model()] has something to predict with — with
#'       `type = "response"` as well as `caret`'s own `"raw"`, the two meaning the
#'       same thing for a regression — and so that `coef(model$fit)` and
#'       `summary(model$fit)` answer as they would for the [stats::lm()] inside
#'       it. This is the one element that is not portable; dropping it leaves an
#'       object that writes out as JSON.}
#'     \item{`metadata`}{Package version, R version, platform and timestamp.}
#'   }
#'
#' @seealso [split_data()], which defines the rows this is fitted on,
#'   [fit_logistic_regression()] for a two-class outcome, [coef.sa_model()] for
#'   the coefficient table, [predict.sa_model()] for predicting the rows it was
#'   not fitted on, and [coef.sa_fit()] for the methods `$fit` answers to, among
#'   them the estimates as a named vector.
#'
#' @examples
#' ## Fitted without resampling (fast enough for examples).
#' fit <- fit_linear_regression(mtcars, outcome = "mpg",
#'                              predictors = c("wt", "hp", "disp"),
#'                              cv = FALSE)
#' fit$coefficients
#'
#' \donttest{
#' ## Cross-validation scores the model; it does not change the coefficients.
#' scored <- fit_linear_regression(mtcars, outcome = "mpg",
#'                                 predictors = c("wt", "hp", "disp"),
#'                                 cv_method = "kfold", seed = 1)
#' all.equal(fit$coefficients, scored$coefficients)
#' }
#'
#' ## A factor predictor becomes one term per level beyond the first.
#' cars <- mtcars
#' cars$cyl <- factor(cars$cyl)
#' fit_linear_regression(cars, outcome = "mpg", predictors = c("wt", "cyl"),
#'                       cv = FALSE)$terms
#'
#' ## The training half of a split is what a model is normally fitted on.
#' sp <- split_data(mtcars, stratified = "mpg", seed = 1)
#' train <- sp$datasets[[1]]$train_data
#' fit_linear_regression(train, outcome = "mpg", predictors = c("wt", "hp"),
#'                       cv = FALSE)$fit_stats$r_squared
#'
#' @export
fit_linear_regression <- function(data,
                                  outcome,
                                  predictors = NULL,
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
  if (!is.numeric(input$y)) {
    stop("`outcome` must be a numeric column for a linear regression, but is ",
         class(input$y)[1], ". Use fit_logistic_regression() for an outcome ",
         "with two classes.", call. = FALSE)
  }
  if (!all(is.finite(input$y))) {
    stop("`outcome` holds non-finite value(s), which least squares has no ",
         "residual for.", call. = FALSE)
  }
  ctrl <- sa_train_control(cv, cv_method, n_fold, n_repeat, input$n_used)

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  fit <- sa_fit_engine(
    caret::train(x = input$x, y = input$y, method = "lm",
                 trControl = ctrl$control),
    "Linear regression"
  )
  model <- fit$finalModel
  summ <- summary(model)

  coefs <- sa_coef_table(
    model,
    sa_wald_interval(summ$coefficients, conf_level, df = model$df.residual),
    df = model$df.residual
  )
  aliased <- names(which(summ$aliased))
  if (length(aliased) > 0L) {
    warning("term(s) could not be estimated because the other predictors ",
            "already span them; their rows are NA: ",
            paste(aliased, collapse = ", "), ".", call. = FALSE)
  }

  # NULL for an intercept-only model, which has nothing to test against itself.
  f <- summ$fstatistic
  fit_stats <- list(
    r_squared     = summ$r.squared,
    adj_r_squared = summ$adj.r.squared,
    sigma         = summ$sigma,
    f_stat        = if (is.null(f)) NA_real_ else unname(f[1]),
    df1           = if (is.null(f)) NA_real_ else unname(f[2]),
    df2           = if (is.null(f)) NA_real_ else unname(f[3]),
    pval          = if (is.null(f)) {
      NA_real_
    } else {
      stats::pf(f[1], f[2], f[3], lower.tail = FALSE)
    },
    aic           = stats::AIC(model),
    bic           = stats::BIC(model)
  )
  fit_stats <- lapply(fit_stats, function(v) unname(as.numeric(v)))

  sa_new_model(
    analysis = "linear_regression",
    terms    = coefs$terms,
    design   = c(list(
      outcome            = input$outcome,
      outcome_type       = "continuous",
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
      method  = "lm",
      label   = "Ordinary least squares linear regression",
      metrics = fit$perfNames
    ),
    fit = fit
  )
}
