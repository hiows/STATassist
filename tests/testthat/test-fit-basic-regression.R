# These two functions are wrappers, so what is worth pinning is not arithmetic
# they perform but arithmetic they must not disturb: the coefficients have to be
# the ones `lm()` and `glm()` give on the same rows, and cross-validating must
# not move them. Everything else here is about the boundary — which rows and
# columns entered the model, which direction the classification reads in, and
# what the object says happened rather than what was asked for.
#
# No external package is used for the reference values. `stats` fits both models
# already, and the Wald interval is one line of arithmetic on the standard error
# the summary reports.

sa_fit_cars <- function() mtcars[c("mpg", "wt", "hp", "disp")]

# Two species rather than three, and two predictors that are not jointly
# separable, so the fit converges and the engine has nothing to warn about.
sa_fit_iris2 <- function() {
  out <- iris[iris$Species != "setosa",
              c("Species", "Petal.Length", "Sepal.Width")]
  rownames(out) <- NULL
  out
}


test_that("the result has the shape the contract promises", {
  fit <- fit_linear_regression(sa_fit_cars(), outcome = "mpg", cv = FALSE)

  expect_s3_class(fit, "sa_model")
  expect_s3_class(fit, "sa_result")
  expect_named(fit, c("analysis", "terms", "design", "parameters",
                      "coefficients", "fit_stats", "performance", "resampling",
                      "engine", "fit", "metadata"))
  expect_identical(fit$analysis, "linear_regression")

  # The slot names do not depend on whether anything was resampled: the two
  # resampling slots are present and NULL rather than absent, since the question
  # they answer was asked and came back empty.
  cv_fit <- fit_linear_regression(sa_fit_cars(), outcome = "mpg",
                                  cv_method = "kfold", n_fold = 3, seed = 1)
  expect_named(cv_fit, names(fit))

  # The inference columns are optional in the contract, because a penalized fit
  # has nothing to put in them, but a model that can report them must: these two
  # are the reason the columns exist.
  expect_true(all(c(sa_model_coef_columns(), sa_model_inference_columns()) %in%
                    names(fit$coefficients)))
  expect_false(anyNA(fit$coefficients$pval))
  expect_identical(fit$coefficients$terms, fit$terms)
  expect_identical(fit$terms, c("(Intercept)", "wt", "hp", "disp"))
  expect_identical(fit$design$predictors, c("wt", "hp", "disp"))

  # The engine object is kept so that the model can predict, which is the one
  # thing a plain list cannot stand in for.
  expect_s3_class(fit$fit, "train")
  expect_length(stats::predict(fit$fit, newdata = mtcars), nrow(mtcars))
})


test_that("`$fit` answers coef() and summary() as the model inside it does", {
  fit <- fit_linear_regression(sa_fit_cars(), outcome = "mpg", cv = FALSE)
  ref <- stats::lm(mpg ~ wt + hp + disp, data = mtcars)

  # `caret` defines no `coef()` for its own class, so without `sa_fit` in front
  # this returns NULL, which is the one answer a fitted model must not give.
  expect_equal(coef(fit$fit), coef(ref))
  expect_s3_class(summary(fit$fit), "summary.lm")
  expect_equal(summary(fit$fit)$coefficients, summary(ref)$coefficients)

  # The result object answers the same generic with its own table, which is the
  # richer answer and the one it was assembled to hold.
  expect_identical(coef(fit), fit$coefficients)

  # The class is prepended rather than substituted, so what `caret` already
  # answered to is still found by inheritance.
  expect_s3_class(fit$fit, "sa_fit")
  expect_s3_class(fit$fit, "train")
  expect_length(stats::predict(fit$fit, newdata = mtcars), nrow(mtcars))

  d <- sa_fit_iris2()
  clf <- fit_logistic_regression(d, outcome = "Species",
                                 outcome_lv = c("versicolor", "virginica"),
                                 cv = FALSE)
  ref_glm <- stats::glm(
    factor(Species, levels = c("versicolor", "virginica")) ~ Petal.Length +
      Sepal.Width,
    data = d, family = stats::binomial()
  )

  expect_equal(coef(clf$fit), coef(ref_glm))
  expect_identical(coef(clf), clf$coefficients)
  expect_s3_class(summary(clf$fit), "summary.glm")
  expect_equal(summary(clf$fit)$coefficients, summary(ref_glm)$coefficients)

  # A classification still predicts class labels and their probabilities, which
  # is what `predict.train` adds over the `glm` it wraps.
  expect_s3_class(stats::predict(clf$fit, newdata = d), "factor")
  expect_named(stats::predict(clf$fit, newdata = d, type = "prob"),
               c("versicolor", "virginica"))
})


test_that("`type = \"response\"` is the probability of the second level", {
  d <- sa_fit_iris2()
  clf <- fit_logistic_regression(d, outcome = "Species",
                                 outcome_lv = c("versicolor", "virginica"),
                                 cv = FALSE)
  ref_glm <- stats::glm(
    factor(Species, levels = c("versicolor", "virginica")) ~ Petal.Length +
      Sepal.Width,
    data = d, family = stats::binomial()
  )

  # `caret` refuses any type but "raw" and "prob", so this is the word this
  # package adds. The number behind it is `caret`'s own, not a second
  # calculation: the probability column the engine already produces.
  resp <- stats::predict(clf$fit, newdata = d, type = "response")
  expect_equal(resp, stats::predict(clf$fit, newdata = d,
                                    type = "prob")[["virginica"]])
  expect_equal(resp, unname(stats::predict(ref_glm, newdata = d,
                                           type = "response")))
  expect_null(names(resp))

  # Which class it is the probability of is the direction rule the rest of the
  # result follows, so swapping the levels turns the probability over.
  other <- fit_logistic_regression(d, outcome = "Species",
                                   outcome_lv = c("virginica", "versicolor"),
                                   cv = FALSE)
  expect_equal(stats::predict(other$fit, newdata = d, type = "response"),
               1 - resp)

  # A regression predicts on the scale of the outcome either way, so the two
  # words name the same numbers there.
  fit <- fit_linear_regression(sa_fit_cars(), outcome = "mpg", cv = FALSE)
  expect_equal(stats::predict(fit$fit, newdata = mtcars, type = "response"),
               stats::predict(fit$fit, newdata = mtcars, type = "raw"))

  # What `caret` already accepted is untouched, and what it never accepted is
  # still its own error rather than a silent guess.
  expect_s3_class(stats::predict(clf$fit, newdata = d, type = "raw"), "factor")
  expect_named(stats::predict(clf$fit, newdata = d, type = "prob"),
               c("versicolor", "virginica"))
  expect_error(stats::predict(clf$fit, newdata = d, type = "link"),
               "raw")
})


test_that("the coefficients are the ones stats::lm gives", {
  fit <- fit_linear_regression(mtcars, outcome = "mpg",
                               predictors = c("wt", "hp", "disp"), cv = FALSE)
  ref <- stats::lm(mpg ~ wt + hp + disp, data = mtcars)
  ref_coefs <- summary(ref)$coefficients

  expect_equal(fit$coefficients$estimate, unname(ref_coefs[, 1]))
  expect_equal(fit$coefficients$stderr, unname(ref_coefs[, 2]))
  expect_equal(fit$coefficients$statistic, unname(ref_coefs[, 3]))
  expect_equal(fit$coefficients$pval, unname(ref_coefs[, 4]))
  expect_equal(unique(fit$coefficients$df), as.numeric(ref$df.residual))

  # The t interval, which is the one that agrees with the t statistic and the
  # standard error in the same row.
  expect_equal(fit$coefficients$lower_conf, unname(stats::confint(ref)[, 1]))
  expect_equal(fit$coefficients$upper_conf, unname(stats::confint(ref)[, 2]))

  expect_equal(fit$fit_stats$r_squared, summary(ref)$r.squared)
  expect_equal(fit$fit_stats$adj_r_squared, summary(ref)$adj.r.squared)
  expect_equal(fit$fit_stats$sigma, summary(ref)$sigma)
  expect_equal(fit$fit_stats$f_stat, unname(summary(ref)$fstatistic[1]))
  expect_equal(fit$fit_stats$df1, unname(summary(ref)$fstatistic[2]))
  expect_equal(fit$fit_stats$df2, unname(summary(ref)$fstatistic[3]))
  expect_equal(fit$fit_stats$aic, stats::AIC(ref))
  expect_equal(fit$fit_stats$bic, stats::BIC(ref))
})


test_that("`conf_level` is the only thing that moves the interval", {
  wide <- fit_linear_regression(mtcars, outcome = "mpg", predictors = "wt",
                                conf_level = 0.99, cv = FALSE)
  narrow <- fit_linear_regression(mtcars, outcome = "mpg", predictors = "wt",
                                  conf_level = 0.8, cv = FALSE)

  expect_equal(wide$coefficients$estimate, narrow$coefficients$estimate)
  expect_true(all(wide$coefficients$lower_conf <
                    narrow$coefficients$lower_conf))
  expect_equal(wide$coefficients$lower_conf,
               unname(stats::confint(stats::lm(mpg ~ wt, data = mtcars),
                                     level = 0.99)[, 1]))
  expect_error(
    fit_linear_regression(mtcars, outcome = "mpg", conf_level = 1),
    "`conf_level` must be in"
  )
})


test_that("the coefficients are the ones stats::glm gives", {
  d <- sa_fit_iris2()
  fit <- fit_logistic_regression(d, outcome = "Species",
                                 outcome_lv = c("versicolor", "virginica"),
                                 cv = FALSE)
  ref <- stats::glm(
    factor(Species, levels = c("versicolor", "virginica")) ~ Petal.Length +
      Sepal.Width,
    data = d, family = stats::binomial()
  )
  ref_coefs <- summary(ref)$coefficients

  expect_identical(fit$analysis, "logistic_regression")
  expect_equal(fit$coefficients$estimate, unname(ref_coefs[, 1]))
  expect_equal(fit$coefficients$stderr, unname(ref_coefs[, 2]))
  expect_equal(fit$coefficients$statistic, unname(ref_coefs[, 3]))
  expect_equal(fit$coefficients$pval, unname(ref_coefs[, 4]))

  # A Wald z is referred to the normal distribution, so there are no degrees of
  # freedom to report and the interval is the normal one rather than a profile.
  expect_true(all(is.na(fit$coefficients$df)))
  crit <- stats::qnorm(0.975)
  expect_equal(fit$coefficients$upper_conf,
               unname(ref_coefs[, 1] + crit * ref_coefs[, 2]))
  expect_equal(fit$coefficients$odds_ratio, exp(fit$coefficients$estimate))
  expect_equal(fit$coefficients$or_lower_conf, exp(fit$coefficients$lower_conf))

  expect_equal(fit$fit_stats$null_deviance, ref$null.deviance)
  expect_equal(fit$fit_stats$residual_deviance, ref$deviance)
  expect_equal(fit$fit_stats$mcfadden_r2, 1 - ref$deviance / ref$null.deviance)
  expect_equal(fit$fit_stats$lr_stat, ref$null.deviance - ref$deviance)
  expect_equal(fit$fit_stats$lr_df, as.numeric(ref$df.null - ref$df.residual))
  expect_equal(fit$fit_stats$aic, stats::AIC(ref))
})


test_that("cross-validation scores the model, it does not change it", {
  args <- list(data = mtcars, outcome = "mpg",
               predictors = c("wt", "hp", "disp"))
  plain <- do.call(fit_linear_regression, c(args, list(cv = FALSE)))
  kfold <- do.call(fit_linear_regression,
                   c(args, list(cv_method = "kfold", n_fold = 4, seed = 1)))
  loocv <- do.call(fit_linear_regression,
                   c(args, list(cv_method = "loocv")))

  expect_identical(kfold$coefficients, plain$coefficients)
  expect_identical(loocv$coefficients, plain$coefficients)
  expect_identical(kfold$fit_stats, plain$fit_stats)
})


test_that("the resampling record counts the fits that actually happened", {
  args <- list(data = mtcars, outcome = "mpg", predictors = "wt")

  repeated <- do.call(fit_linear_regression,
                      c(args, list(n_fold = 4, n_repeat = 3, seed = 1)))
  expect_identical(nrow(repeated$resampling), 12L)
  expect_identical(nrow(repeated$performance), 1L)

  kfold <- do.call(fit_linear_regression,
                   c(args, list(cv_method = "kfold", n_fold = 4, seed = 1)))
  expect_identical(nrow(kfold$resampling), 4L)

  # Leave-one-out reports its score but no per-fold table, and no resampling at
  # all leaves both slots empty.
  loocv <- do.call(fit_linear_regression, c(args, list(cv_method = "loocv")))
  expect_null(loocv$resampling)
  expect_identical(nrow(loocv$performance), 1L)

  plain <- do.call(fit_linear_regression, c(args, list(cv = FALSE)))
  expect_null(plain$resampling)
  expect_null(plain$performance)

  # Which metrics exist is a property of the model type rather than of the
  # resampling, so a model that was not scored still names them.
  expect_identical(repeated$engine$metrics, c("RMSE", "Rsquared", "MAE"))
  expect_identical(plain$engine$metrics, c("RMSE", "Rsquared", "MAE"))
  expect_identical(
    fit_logistic_regression(sa_fit_iris2(), outcome = "Species",
                            cv = FALSE)$engine$metrics,
    c("Accuracy", "Kappa")
  )
})


test_that("engine notes are reported once with a count, not per fit", {
  # A predictor that splits the classes exactly has no finite maximum likelihood
  # estimate, and the engine says so. The note is grouped and re-emitted rather
  # than muffled: a model that did not converge has to be able to say it.
  separable <- data.frame(y = rep(c("a", "b"), each = 10), x = c(1:10, 21:30))

  expect_message(
    fit <- fit_logistic_regression(separable, outcome = "y", cv = FALSE),
    "engine note\\(s\\) while fitting"
  )
  expect_message(
    fit_logistic_regression(separable, outcome = "y", cv = FALSE),
    "\\[1 time\\(s\\)\\]"
  )
  # The fit still comes back, with the enormous standard error that says why.
  expect_gt(fit$coefficients$stderr[2], 100)
})


test_that("`parameters` records the scheme that ran, not the one asked for", {
  args <- list(data = mtcars, outcome = "mpg", predictors = "wt",
               n_fold = 4, n_repeat = 3)

  repeated <- do.call(fit_linear_regression, c(args, list(seed = 1)))
  expect_identical(repeated$parameters$cv_method, "repeated_kfold")
  expect_identical(repeated$parameters$n_fold, 4L)
  expect_identical(repeated$parameters$n_repeat, 3L)

  # Plain k-fold does not repeat and leave-one-out has no fold count, so the
  # arguments they ignore are reported as NA rather than as the number that was
  # passed in and never used.
  kfold <- do.call(fit_linear_regression,
                   c(args, list(cv_method = "kfold", seed = 1)))
  expect_identical(kfold$parameters$n_fold, 4L)
  expect_identical(kfold$parameters$n_repeat, NA_integer_)

  loocv <- do.call(fit_linear_regression, c(args, list(cv_method = "loocv")))
  expect_identical(loocv$parameters$cv_method, "loocv")
  expect_identical(loocv$parameters$n_fold, NA_integer_)
  expect_identical(loocv$parameters$n_repeat, NA_integer_)

  plain <- do.call(fit_linear_regression, c(args, list(cv = FALSE)))
  expect_identical(plain$parameters$cv_method, NA_character_)
  expect_identical(plain$parameters$n_fold, NA_integer_)

  # Invalid resampling arguments are refused whether or not the chosen scheme
  # would have read them, so what is accepted does not depend on `cv_method`.
  expect_error(
    fit_linear_regression(mtcars, outcome = "mpg", cv = FALSE, n_fold = 1),
    "`n_fold` must be in"
  )
  expect_error(
    fit_linear_regression(mtcars, outcome = "mpg", cv_method = "loocv",
                          n_repeat = 0),
    "`n_repeat` must be in"
  )
  expect_error(
    fit_linear_regression(mtcars, outcome = "mpg", predictors = "wt",
                          cv_method = "kfold", n_fold = 40),
    "exceeds the 32 usable observation"
  )
})


test_that("the first level of `outcome_lv` is the reference", {
  d <- sa_fit_iris2()
  versi <- fit_logistic_regression(d, outcome = "Species",
                                   outcome_lv = c("versicolor", "virginica"),
                                   cv = FALSE)
  virgi <- fit_logistic_regression(d, outcome = "Species",
                                   outcome_lv = c("virginica", "versicolor"),
                                   cv = FALSE)

  expect_identical(versi$design$outcome_lv, c("versicolor", "virginica"))
  expect_equal(versi$coefficients$estimate, -virgi$coefficients$estimate)
  expect_equal(versi$coefficients$odds_ratio,
               1 / virgi$coefficients$odds_ratio)
  expect_equal(versi$coefficients$stderr, virgi$coefficients$stderr)
  expect_equal(versi$coefficients$pval, virgi$coefficients$pval)

  # `n_events` counts the class the coefficients are about, so it moves with
  # the reference rather than staying with one label.
  expect_identical(versi$design$n_events, 50L)
  expect_identical(virgi$design$n_events, 50L)
  expect_equal(versi$design$event_rate, 0.5)

  # Longer petals make `virginica` more likely, so the odds ratio is above 1
  # when `virginica` is the level being modelled and below 1 when it is the
  # reference. This is the direction `compare_two_groups()` reads `group_lv` in.
  petal <- versi$coefficients$terms == "Petal.Length"
  expect_gt(versi$coefficients$odds_ratio[petal], 1)
  expect_lt(virgi$coefficients$odds_ratio[petal], 1)

  # Sorting is what `NULL` does, and it puts `versicolor` first here.
  sorted <- fit_logistic_regression(d, outcome = "Species", cv = FALSE)
  expect_identical(sorted$coefficients, versi$coefficients)
})


test_that("an outcome is read as its classes, whatever type it arrived as", {
  d <- mtcars[c("am", "wt")]
  numeric_am <- fit_logistic_regression(d, outcome = "am", cv = FALSE)
  expect_identical(numeric_am$design$outcome_lv, c("0", "1"))
  expect_identical(numeric_am$design$n_events, 13L)

  d$am <- factor(ifelse(d$am == 1, "manual", "auto"),
                 levels = c("auto", "manual", "unobserved"))
  # An unused factor level is a property of how the column was built, not a
  # third class, so it does not make this a call for a different model.
  labelled <- fit_logistic_regression(d, outcome = "am", cv = FALSE)
  expect_identical(labelled$design$outcome_lv, c("auto", "manual"))
  expect_equal(labelled$coefficients$estimate, numeric_am$coefficients$estimate)
})


test_that("an outcome the model cannot use is refused by name", {
  expect_error(
    fit_linear_regression(iris, outcome = "Species", cv = FALSE),
    "must be a numeric column for a linear regression"
  )
  expect_error(
    fit_logistic_regression(iris, outcome = "Species", cv = FALSE),
    "holds 3 classes, but a logistic regression models two"
  )
  # Naming two of three would fit a subset of the rows that were passed in.
  expect_error(
    fit_logistic_regression(iris, outcome = "Species",
                            outcome_lv = c("setosa", "virginica"), cv = FALSE),
    "would be silently left out"
  )
  expect_error(
    fit_logistic_regression(iris[iris$Species == "setosa", ],
                            outcome = "Species", cv = FALSE),
    "nothing to classify"
  )
  expect_error(
    fit_logistic_regression(sa_fit_iris2(), outcome = "Species",
                            outcome_lv = c("versicolor", "setosa"),
                            cv = FALSE),
    "absent from `outcome`"
  )
  expect_error(
    fit_logistic_regression(sa_fit_iris2(), outcome = "Species",
                            outcome_lv = "versicolor", cv = FALSE),
    "two distinct level names"
  )
})


test_that("a factor predictor becomes one term per level beyond the first", {
  cars <- mtcars
  cars$cyl <- factor(cars$cyl)
  fit <- fit_linear_regression(cars, outcome = "mpg",
                               predictors = c("wt", "cyl"), cv = FALSE)

  expect_identical(fit$terms, c("(Intercept)", "wt", "cyl6", "cyl8"))
  expect_identical(fit$design$predictors, c("wt", "cyl"))
  expect_equal(fit$coefficients$estimate,
               unname(stats::coef(stats::lm(mpg ~ wt + cyl, data = cars))))

  # A character column is the same model as the factor it stands for, since it
  # is turned into one before the folds are drawn rather than inside the engine.
  cars$cyl <- as.character(cars$cyl)
  chr <- fit_linear_regression(cars, outcome = "mpg",
                               predictors = c("wt", "cyl"), cv = FALSE)
  expect_identical(chr$coefficients, fit$coefficients)
})


test_that("rows missing anything the model needs are dropped once", {
  d <- mtcars[c("mpg", "wt", "hp")]
  d$wt[1:3] <- NA
  d$mpg[4] <- NA
  fit <- fit_linear_regression(d, outcome = "mpg", cv = FALSE)

  expect_identical(fit$design$n_obs, 32L)
  expect_identical(fit$design$n_used, 28L)
  expect_identical(fit$design$n_dropped, 4L)

  # Dropped before the folds are drawn, so every fold sees the same rows. The
  # fit is therefore the one the complete rows give on their own.
  ref <- stats::lm(mpg ~ wt + hp, data = d[stats::complete.cases(d), ])
  expect_equal(fit$coefficients$estimate, unname(stats::coef(ref)))

  expect_error(
    fit_linear_regression(data.frame(y = c(1, NA, NA), x = c(1, 2, 3)),
                          outcome = "y", cv = FALSE),
    "at least 2 are needed"
  )
})


test_that("a predictor that takes one value is left out with a message", {
  d <- mtcars[c("mpg", "wt")]
  d$flat <- 1

  expect_message(
    fit <- fit_linear_regression(d, outcome = "mpg", cv = FALSE),
    "single value cannot contribute"
  )
  expect_identical(fit$design$predictors, "wt")
  expect_identical(fit$design$dropped_predictors, "flat")
  expect_identical(fit$terms, c("(Intercept)", "wt"))

  expect_error(
    suppressMessages(
      fit_linear_regression(data.frame(y = c(1, 2, 3), x = c(1, 1, 1)),
                            outcome = "y", cv = FALSE)
    ),
    "nothing to fit"
  )
})


test_that("a term the fit could not estimate keeps its row", {
  d <- mtcars[c("mpg", "wt")]
  d$wt_again <- d$wt * 2

  expect_warning(
    fit <- fit_linear_regression(d, outcome = "mpg", cv = FALSE),
    "could not be estimated"
  )
  # Present and NA, not absent: the term was in the model and has no answer,
  # which is a different fact from a table that never mentioned it.
  expect_identical(fit$terms, c("(Intercept)", "wt", "wt_again"))
  aliased <- fit$coefficients[fit$coefficients$terms == "wt_again", ]
  expect_true(is.na(aliased$estimate))
  expect_true(is.na(aliased$pval))
  expect_true(is.na(aliased$lower_conf))
  # The vector `$fit` gives keeps it too, so it is as long as the model has terms.
  expect_identical(names(coef(fit$fit)), fit$terms)
  expect_true(is.na(coef(fit$fit)[["wt_again"]]))
})


test_that("the outcome cannot also be a predictor", {
  expect_error(
    fit_linear_regression(mtcars, outcome = "mpg",
                          predictors = c("wt", "mpg"), cv = FALSE),
    "predict from the answer"
  )
  # `NULL` takes every other column, which is the same rule stated once.
  fit <- fit_linear_regression(mtcars[c("mpg", "wt", "hp")], outcome = "mpg",
                               cv = FALSE)
  expect_identical(fit$design$predictors, c("wt", "hp"))
})


test_that("an outcome given as a vector is fitted the same as one named", {
  by_name <- fit_linear_regression(mtcars[c("mpg", "wt", "hp")],
                                   outcome = "mpg", cv = FALSE)
  by_vector <- fit_linear_regression(mtcars[c("wt", "hp")],
                                     outcome = mtcars$mpg, cv = FALSE)

  expect_identical(by_vector$coefficients, by_name$coefficients)
  # A resolved vector no longer remembers where it came from, and the design
  # says so rather than inventing a column name.
  expect_identical(by_name$design$outcome, "mpg")
  expect_identical(by_vector$design$outcome, "<vector>")

  expect_error(
    fit_linear_regression(mtcars, outcome = mtcars$mpg[1:5]),
    "one entry per row"
  )
  expect_error(
    fit_linear_regression(mtcars, outcome = "mpg", predictors = "nope"),
    "`predictors` not found in `data`"
  )
  expect_error(fit_linear_regression(1:5, outcome = "mpg"),
               "must be a data.frame or a matrix")
})


test_that("`seed` fixes the folds and leaves the stream as it was", {
  args <- list(data = mtcars, outcome = "mpg", predictors = c("wt", "hp"),
               cv_method = "kfold", n_fold = 4)
  once <- do.call(fit_linear_regression, c(args, list(seed = 11)))
  twice <- do.call(fit_linear_regression, c(args, list(seed = 11)))
  other <- do.call(fit_linear_regression, c(args, list(seed = 12)))

  expect_identical(once$resampling, twice$resampling)
  expect_false(isTRUE(all.equal(once$resampling$RMSE, other$resampling$RMSE)))

  set.seed(99)
  before <- .Random.seed
  invisible(do.call(fit_linear_regression, c(args, list(seed = 11))))
  expect_identical(before, .Random.seed)
})


test_that("print says what was fitted to what and how it did", {
  fit <- fit_linear_regression(mtcars, outcome = "mpg",
                               predictors = c("wt", "hp", "disp"),
                               cv_method = "kfold", n_fold = 4, seed = 1)

  expect_output(print(fit), "<sa_model> linear_regression")
  expect_output(print(fit), "outcome  : mpg")
  expect_output(print(fit), "kfold, 4 fold\\(s\\)")
  expect_output(print(fit), "r_squared")
  expect_output(print(fit), "RMSE")
  expect_output(print(fit), "\\(Intercept\\)")

  # The table is summarised rather than dumped, and what was left out is counted.
  expect_output(print(fit, n = 1), "and 3 more term\\(s\\)")
  expect_output(expect_invisible(print(fit)))

  plain <- fit_linear_regression(mtcars, outcome = "mpg", predictors = "wt",
                                 cv = FALSE)
  expect_output(print(plain), "no resampling")

  logistic <- fit_logistic_regression(sa_fit_iris2(), outcome = "Species",
                                      cv = FALSE)
  expect_output(print(logistic), "two classes")
  expect_output(print(logistic),
                "modelling the odds of virginica against versicolor")
  expect_output(print(logistic), "mcfadden_r2")
})


test_that("a model is fitted on the training half of a split", {
  sp <- split_data(mtcars, stratified = "mpg", seed = 1)
  train <- sp$datasets[[1]]$train_data
  test <- sp$datasets[[1]]$test_data
  fit <- fit_linear_regression(train, outcome = "mpg",
                               predictors = c("wt", "hp"), cv = FALSE)

  expect_identical(fit$design$n_used, nrow(train))
  expect_length(stats::predict(fit$fit, newdata = test), nrow(test))

  # Everything but the engine object survives the trip a JSON export would take,
  # which is the whole reason the rest of the object is a plain list.
  skip_if_not_installed("jsonlite")
  portable <- fit[setdiff(names(fit), "fit")]
  round_trip <- jsonlite::fromJSON(
    jsonlite::toJSON(portable, na = "string", digits = NA)
  )
  expect_identical(round_trip$terms, fit$terms)
  expect_equal(round_trip$coefficients$estimate, fit$coefficients$estimate)
})
