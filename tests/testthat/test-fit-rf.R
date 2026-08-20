# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# A forest is random twice over — the rows of each tree and the predictors of each
# split — so almost nothing here can be pinned against a second fit of the same
# data. Two forests grown from the same seed with the class levels in the other
# order are not the same forest either, which is why the direction rule is checked
# by recomputing each fit's own numbers from its own out-of-bag predictions rather
# than by comparing one fit with another.
#
# What that leaves worth pinning is the boundary. Which of the four importance
# columns `estimate` reads, that the reported fit statistics are the out-of-bag
# ones and not the in-sample ones a forest would answer near perfectly, that the
# terms are the predictors rather than dummy columns, and that the table has no
# inference columns at all — the same promise `fit_elastic_net()` makes, for a
# different reason.
#
# The reference values are the forest's own. It is what grew the trees, so this
# checks that the right column of the right matrix is being read, not that a
# forest is built correctly.

sa_rf_cars <- function() mtcars[c("mpg", "wt", "hp", "disp")]

sa_rf_iris2 <- function() {
  out <- iris[iris$Species != "setosa", ]
  rownames(out) <- NULL
  out
}


test_that("the result has the shape the contract promises", {
  fit <- fit_rf(sa_rf_cars(), outcome = "mpg", ntree = 100, cv = FALSE,
                seed = 1)

  expect_s3_class(fit, "sa_model")
  expect_s3_class(fit, "sa_result")
  expect_named(fit, c("analysis", "terms", "design", "parameters",
                      "coefficients", "fit_stats", "performance", "resampling",
                      "engine", "fit", "metadata"))
  expect_identical(fit$analysis, "random_forest")
  expect_identical(fit$coefficients$terms, fit$terms)
  expect_setequal(fit$terms, c("wt", "hp", "disp"))
  expect_identical(fit$design$predictors, c("wt", "hp", "disp"))
  expect_identical(fit$design$outcome_type, "continuous")
  expect_identical(fit$engine$method, "rf")
  expect_identical(fit$engine$label, "Random forest regression")

  # The table holds what a forest answers and nothing else. There is no effect per
  # unit of a predictor, so no standard error, so no inference columns at all
  # rather than six columns of `NA`; and no intercept, since a tree has none.
  expect_named(fit$coefficients, c("terms", "estimate", "impurity"))
  expect_length(intersect(sa_model_inference_columns(),
                          names(fit$coefficients)), 0L)
  expect_null(fit$coefficients$selected)
  expect_null(fit$coefficients$odds_ratio)
  expect_false("(Intercept)" %in% fit$terms)

  # The slot names do not depend on whether anything was resampled.
  cv_fit <- fit_rf(sa_rf_cars(), outcome = "mpg", ntree = 100,
                   cv_method = "kfold", n_fold = 3, seed = 1)
  expect_named(cv_fit, names(fit))

  expect_s3_class(fit$fit, "train")
  expect_length(stats::predict(fit$fit, newdata = mtcars), nrow(mtcars))
})


test_that("`estimate` is the permutation importance the forest reports", {
  fit <- fit_rf(sa_rf_cars(), outcome = "mpg", ntree = 200, cv = FALSE,
                seed = 1)
  reported <- fit$fit$finalModel$importance
  at <- match(fit$terms, rownames(reported))

  # Which of the four columns is read is the whole of what this fixes. The
  # permutation measure is scored on out-of-bag rows and the impurity one on the
  # splits themselves, so a predictor can rank differently under the two.
  expect_equal(fit$coefficients$estimate, unname(reported[at, "%IncMSE"]))
  expect_equal(fit$coefficients$impurity,
               unname(reported[at, "IncNodePurity"]))

  # Unscaled: `randomForest::importance()` divides by the standard deviation
  # across trees by default, which is a different number.
  expect_equal(fit$coefficients$estimate,
               unname(randomForest::importance(fit$fit$finalModel,
                                               scale = FALSE)[at, "%IncMSE"]))

  # Descending, so the first rows of the table are the predictors the forest
  # leaned on. The order the columns were read in is `design$predictors`.
  expect_identical(fit$terms, fit$terms[order(fit$coefficients$estimate,
                                              decreasing = TRUE)])
  expect_setequal(fit$terms, fit$design$predictors)

  # The result answers `coef()` with its own table, as the other models do.
  expect_identical(coef(fit), fit$coefficients)
  expect_s3_class(coef(fit), "data.frame")

  # `coef()` on the engine object is where a forest differs: `coef.default()`
  # would answer NULL, which is the one answer a fitted model must not give.
  expect_error(coef(fit$fit), "estimates no coefficients")
})


test_that("the classification importance is the other pair of columns", {
  clf <- fit_rf(sa_rf_iris2(), outcome = "Species",
                outcome_lv = c("versicolor", "virginica"), ntree = 200,
                cv = FALSE, seed = 1)
  reported <- clf$fit$finalModel$importance
  at <- match(clf$terms, rownames(reported))

  expect_equal(clf$coefficients$estimate,
               unname(reported[at, "MeanDecreaseAccuracy"]))
  expect_equal(clf$coefficients$impurity,
               unname(reported[at, "MeanDecreaseGini"]))
  expect_identical(clf$engine$label, "Random forest classification")

  # The per-class columns are there too and are not what is reported: the table
  # holds one importance per predictor, not one per predictor and class.
  expect_true(all(c("versicolor", "virginica") %in% colnames(reported)))
  expect_identical(nrow(clf$coefficients), 4L)
})


test_that("the fit statistics are the out-of-bag ones", {
  d <- sa_rf_cars()
  fit <- fit_rf(d, outcome = "mpg", ntree = 200, cv = FALSE, seed = 1)
  oob <- fit$fit$finalModel$predicted
  residual <- d$mpg - oob

  expect_named(fit$fit_stats,
               c("oob_r_squared", "oob_rmse", "oob_mae", "n_oob"))
  expect_equal(fit$fit_stats$oob_rmse, sqrt(mean(residual^2)))
  expect_equal(fit$fit_stats$oob_mae, mean(abs(residual)))
  expect_equal(fit$fit_stats$oob_r_squared,
               1 - sum(residual^2) / sum((d$mpg - mean(d$mpg))^2))
  expect_identical(fit$fit_stats$n_oob, as.numeric(nrow(d)))

  # In-sample would be a forest predicting rows it had memorised, and would come
  # back far better than this. That it does not is what the prefix promises.
  in_sample <- stats::predict(fit$fit, newdata = d)
  expect_lt(sqrt(mean((d$mpg - in_sample)^2)), fit$fit_stats$oob_rmse)
})


test_that("a classification's out-of-bag statistics read `outcome_lv`", {
  d <- sa_rf_iris2()
  lv <- c("versicolor", "virginica")
  clf <- fit_rf(d, outcome = "Species", outcome_lv = lv, ntree = 200,
                cv = FALSE, seed = 1)
  oob <- as.character(clf$fit$finalModel$predicted)

  expect_named(clf$fit_stats,
               c("oob_accuracy", "oob_error", "oob_kappa", "oob_sensitivity",
                 "oob_specificity", "n_oob"))
  expect_equal(clf$fit_stats$oob_accuracy, mean(oob == d$Species))
  expect_equal(clf$fit_stats$oob_error, 1 - clf$fit_stats$oob_accuracy)

  # The forest's own out-of-bag error rate, reached by another route.
  expect_equal(clf$fit_stats$oob_accuracy,
               1 - unname(clf$fit$finalModel$err.rate[clf$parameters$ntree,
                                                      "OOB"]))

  # Sensitivity is the share of `outcome_lv[2]` found and specificity the share of
  # `outcome_lv[1]` left alone, which is the direction the whole family reads in.
  expect_equal(clf$fit_stats$oob_sensitivity,
               mean(oob[d$Species == lv[2]] == lv[2]))
  expect_equal(clf$fit_stats$oob_specificity,
               mean(oob[d$Species == lv[1]] == lv[1]))

  called <- oob == lv[2]
  positive <- d$Species == lv[2]
  expected <- mean(called) * mean(positive) + mean(!called) * mean(!positive)
  expect_equal(clf$fit_stats$oob_kappa,
               (clf$fit_stats$oob_accuracy - expected) / (1 - expected))

  # The same rule read the other way round: with the levels swapped, what was the
  # sensitivity is now measured on the other class.
  other <- fit_rf(d, outcome = "Species", outcome_lv = rev(lv), ntree = 200,
                  cv = FALSE, seed = 1)
  other_oob <- as.character(other$fit$finalModel$predicted)
  expect_equal(other$fit_stats$oob_sensitivity,
               mean(other_oob[d$Species == lv[1]] == lv[1]))
  expect_identical(other$design$outcome_lv, rev(lv))
})


test_that("`mtry` is the rule of thumb unless the caller names one", {
  # A third of the predictors for a regression, their square root for a
  # classification, and at least one however few there are.
  regression <- fit_rf(mtcars, outcome = "mpg", ntree = 100, cv = FALSE,
                       seed = 1)
  expect_identical(regression$parameters$mtry, 3L)
  expect_identical(regression$parameters$n_candidates, 1L)

  classification <- fit_rf(sa_rf_iris2(), outcome = "Species", ntree = 100,
                           cv = FALSE, seed = 1)
  expect_identical(classification$parameters$mtry, 2L)

  one <- fit_rf(mtcars, outcome = "mpg", predictors = "wt", ntree = 100,
                cv = FALSE, seed = 1)
  expect_identical(one$parameters$mtry, 1L)

  # The defaults `randomForest` uses, and both are recorded as they were used.
  expect_identical(regression$parameters$nodesize, 5L)
  expect_identical(classification$parameters$nodesize, 1L)
  expect_identical(regression$parameters$ntree, 100L)
})


test_that("the reported candidate is the one the resampling chose", {
  fit <- fit_rf(mtcars, outcome = "mpg", mtry = c(2, 5, 10), ntree = 100,
                cv_method = "kfold", n_fold = 5, seed = 1)

  expect_identical(nrow(fit$performance), 3L)
  expect_identical(fit$parameters$n_candidates, 3L)
  expect_identical(fit$parameters$mtry, as.integer(fit$fit$bestTune$mtry))
  best <- which(fit$performance$mtry == fit$parameters$mtry)
  expect_equal(fit$performance$RMSE[best], min(fit$performance$RMSE))

  # `mtry` is the one argument the resampling tunes; `ntree` and `nodesize` are
  # the same for every candidate, so they are scalars in `parameters`.
  expect_null(fit$performance$ntree)
  expect_length(fit$parameters$ntree, 1L)

  plain <- fit_rf(mtcars, outcome = "mpg", ntree = 100, cv = FALSE, seed = 1)
  expect_null(plain$performance)
  expect_null(plain$resampling)
  # The out-of-bag score is there whether or not anything was resampled, which is
  # what a forest has that the other models do not.
  expect_false(is.null(plain$fit_stats$oob_r_squared))
})


test_that("a `mtry` the forest could not honour is refused by name", {
  expect_error(
    fit_rf(mtcars, outcome = "mpg", mtry = c(2, 5), cv = FALSE),
    "must hold one candidate"
  )
  # `randomForest()` resets an impossible `mtry` to the valid range and fits, so
  # the forest would be at a different `mtry` than the one reported.
  expect_error(
    fit_rf(mtcars, outcome = "mpg", mtry = 99, cv = FALSE),
    "cannot exceed the 10 predictor\\(s\\)"
  )
  expect_error(
    fit_rf(mtcars, outcome = "mpg", mtry = 2.5, cv = FALSE),
    "must hold whole numbers"
  )
  expect_error(fit_rf(mtcars, outcome = "mpg", mtry = 0, cv = FALSE),
               "`mtry` must be in")
  expect_error(fit_rf(mtcars, outcome = "mpg", ntree = 0, cv = FALSE),
               "`ntree` must be in")
  expect_error(fit_rf(mtcars, outcome = "mpg", nodesize = 0, cv = FALSE),
               "`nodesize` must be in")
  expect_error(
    fit_rf(mtcars, outcome = "mpg", predictors = "wt", cv_method = "kfold",
           n_fold = 40),
    "exceeds the 32 usable observation"
  )
})


test_that("the outcome decides which forest is grown", {
  numeric_outcome <- fit_rf(sa_rf_cars(), outcome = "mpg", ntree = 100,
                            cv = FALSE, seed = 1)
  expect_identical(numeric_outcome$design$outcome_type, "continuous")
  expect_identical(numeric_outcome$engine$metrics,
                   c("RMSE", "Rsquared", "MAE"))
  expect_null(numeric_outcome$design$outcome_lv)

  labels <- fit_rf(sa_rf_iris2(), outcome = "Species", ntree = 100,
                   cv = FALSE, seed = 1)
  expect_identical(labels$design$outcome_type, "two classes")
  expect_identical(labels$engine$metrics, c("Accuracy", "Kappa"))
  expect_identical(labels$design$outcome_lv, c("versicolor", "virginica"))
  expect_identical(labels$design$n_events, 50L)
  expect_equal(labels$design$event_rate, 0.5)

  # A numeric column holding two values is the one case the outcome cannot settle
  # on its own, so the guess is announced and `outcome_lv` overrules it.
  d <- mtcars[c("am", "wt", "hp")]
  expect_message(
    as_numbers <- fit_rf(d, outcome = "am", ntree = 100, cv = FALSE, seed = 1),
    "fitted as a regression"
  )
  expect_identical(as_numbers$design$outcome_type, "continuous")

  as_classes <- fit_rf(d, outcome = "am", outcome_lv = c("0", "1"),
                       ntree = 100, cv = FALSE, seed = 1)
  expect_identical(as_classes$design$outcome_type, "two classes")
  expect_identical(as_classes$design$n_events, 13L)

  # More than two classes is the same refusal the other classifications make, in
  # this model's name.
  expect_error(
    fit_rf(iris, outcome = "Species", ntree = 100, cv = FALSE),
    "holds 3 classes, but a random forest models two"
  )
  expect_error(
    fit_rf(iris, outcome = "Species", outcome_lv = c("setosa", "virginica"),
           ntree = 100, cv = FALSE),
    "would be silently left out"
  )
  expect_error(
    fit_rf(mtcars, outcome = "mpg", predictors = c("wt", "mpg"), ntree = 100,
           cv = FALSE),
    "predict from the answer"
  )
})


test_that("a factor predictor is one term rather than one per level", {
  cars <- mtcars[c("mpg", "wt")]
  cars$cyl <- factor(mtcars$cyl)
  fit <- fit_rf(cars, outcome = "mpg", ntree = 100, cv = FALSE, seed = 1)

  # A tree splits a factor on its levels directly, so nothing is dummy coded.
  # This is what the linear models cannot do: they report `cyl6` and `cyl8`.
  expect_setequal(fit$terms, c("wt", "cyl"))
  expect_identical(fit$design$predictors, c("wt", "cyl"))
  expect_identical(fit$fit$finalModel$xNames, c("wt", "cyl"))
  expect_identical(
    fit_linear_regression(cars, outcome = "mpg", cv = FALSE)$terms,
    c("(Intercept)", "wt", "cyl6", "cyl8")
  )
  expect_identical(fit$design$predictor_lv, list(cyl = c("4", "6", "8")))

  # A character column is the same model as the factor it stands for, since it is
  # turned into one before the forest is grown.
  cars$cyl <- as.character(cars$cyl)
  chr <- fit_rf(cars, outcome = "mpg", ntree = 100, cv = FALSE, seed = 1)
  expect_identical(chr$coefficients, fit$coefficients)
})


test_that("a model of one predictor is grown rather than refused", {
  # The contrast with `fit_elastic_net()`, which refuses a single term because a
  # penalty has no budget to divide. A forest of stumps on one column is a model.
  fit <- fit_rf(mtcars, outcome = "mpg", predictors = "wt", ntree = 100,
                cv = FALSE, seed = 1)
  expect_identical(fit$terms, "wt")
  expect_identical(nrow(fit$coefficients), 1L)
  expect_error(
    fit_elastic_net(mtcars, outcome = "mpg", predictors = "wt",
                    penalty = "ridge", lambda = 1, cv = FALSE),
    "the model has 1 term"
  )
})


test_that("the rows and columns that entered the model are the reported ones", {
  d <- mtcars[c("mpg", "wt", "hp")]
  d$wt[1:3] <- NA
  d$mpg[4] <- NA
  fit <- fit_rf(d, outcome = "mpg", ntree = 100, cv = FALSE, seed = 1)

  expect_identical(fit$design$n_obs, 32L)
  expect_identical(fit$design$n_used, 28L)
  expect_identical(fit$design$n_dropped, 4L)
  # Dropped before the trees are grown, so the out-of-bag score is over the rows
  # that were kept and nothing else.
  expect_identical(fit$fit_stats$n_oob, 28)

  flat <- mtcars[c("mpg", "wt", "hp")]
  flat$flat <- 1
  expect_message(
    dropped <- fit_rf(flat, outcome = "mpg", ntree = 100, cv = FALSE, seed = 1),
    "single value cannot contribute"
  )
  expect_identical(dropped$design$dropped_predictors, "flat")
  expect_setequal(dropped$terms, c("wt", "hp"))
})


test_that("`seed` fixes the forest as well as the folds", {
  # More is fixed here than in the other models, where the seed fixes only the
  # fold assignment: the rows of every tree and the predictors of every split are
  # draws, so two forests of the same data differ without one.
  args <- list(data = mtcars, outcome = "mpg", predictors = c("wt", "hp", "disp"),
               ntree = 100, cv = FALSE)
  once <- do.call(fit_rf, c(args, list(seed = 11)))
  twice <- do.call(fit_rf, c(args, list(seed = 11)))
  other <- do.call(fit_rf, c(args, list(seed = 12)))

  expect_identical(once$coefficients, twice$coefficients)
  expect_identical(once$fit_stats, twice$fit_stats)
  expect_false(isTRUE(all.equal(once$coefficients$estimate,
                                other$coefficients$estimate)))

  set.seed(99)
  before <- .Random.seed
  invisible(do.call(fit_rf, c(args, list(seed = 11))))
  expect_identical(before, .Random.seed)

  # The resampling arguments are recorded as the scheme used them, the way the
  # rest of the family records them.
  cv_args <- list(data = mtcars, outcome = "mpg", predictors = "wt",
                  ntree = 100, n_fold = 4, n_repeat = 3, seed = 1)
  repeated <- do.call(fit_rf, cv_args)
  expect_identical(repeated$parameters$cv_method, "repeated_kfold")
  expect_identical(nrow(repeated$resampling), 12L)

  kfold <- do.call(fit_rf, c(cv_args, list(cv_method = "kfold")))
  expect_identical(kfold$parameters$n_fold, 4L)
  expect_identical(kfold$parameters$n_repeat, NA_integer_)

  plain <- do.call(fit_rf, c(args, list(seed = 1)))
  expect_identical(plain$parameters$cv_method, NA_character_)
  expect_identical(plain$parameters$n_fold, NA_integer_)
})


test_that("print says what was grown, at what settings, and how it did", {
  fit <- fit_rf(mtcars, outcome = "mpg", predictors = c("wt", "hp", "disp"),
                mtry = c(1, 2, 3), ntree = 100, cv_method = "kfold",
                n_fold = 4, seed = 1)

  expect_output(print(fit), "<sa_model> random_forest")
  expect_output(print(fit), "forest   : 100 tree\\(s\\), mtry = ")
  expect_output(print(fit), "nodesize = 5")
  expect_output(print(fit), "mtry chosen from 3 candidate\\(s\\)")
  expect_output(print(fit), "kfold, 4 fold\\(s\\)")
  expect_output(print(fit), "oob_r_squared")
  expect_output(print(fit), "RMSE")
  expect_output(expect_invisible(print(fit)))

  # The block is not called `coefficients`, because these are not coefficients,
  # and there is no interval and no p-value to print beside them.
  # The heading also names which of the two importances the printed column is,
  # since `impurity` is in the table and is not the number beside each term.
  output <- utils::capture.output(print(fit))
  expect_true(any(grepl("importance  (permutation)", output, fixed = TRUE)))
  expect_false(any(grepl("coefficients", output, fixed = TRUE)))
  expect_false(any(grepl("p = NA", output, fixed = TRUE)))
  expect_false(any(grepl("selected", output, fixed = TRUE)))
  expect_false(any(grepl("conf_level", output, fixed = TRUE)))

  # A single `mtry` was scored rather than chosen between, so nothing won.
  single <- utils::capture.output(
    print(fit_rf(mtcars, outcome = "mpg", predictors = "wt", ntree = 100,
                 cv_method = "kfold", n_fold = 4, seed = 1))
  )
  expect_true(any(grepl("forest   : ", single, fixed = TRUE)))
  expect_false(any(grepl("chosen from", single, fixed = TRUE)))

  # A forest models the same class in the same direction as the logistic
  # regression and reports no odds, so the header does not claim any.
  clf <- utils::capture.output(
    print(fit_rf(sa_rf_iris2(), outcome = "Species", ntree = 100, cv = FALSE,
                 seed = 1))
  )
  expect_true(any(grepl("modelling virginica against versicolor", clf,
                        fixed = TRUE)))
  expect_true(any(grepl("oob_accuracy", clf, fixed = TRUE)))
  expect_false(any(grepl("the odds of", clf, fixed = TRUE)))

  # The models that do report odds still say so.
  logistic <- utils::capture.output(
    print(fit_logistic_regression(sa_rf_iris2(), outcome = "Species",
                                  predictors = "Petal.Length", cv = FALSE))
  )
  expect_true(any(grepl("modelling the odds of virginica", logistic,
                        fixed = TRUE)))
})


test_that("what the forest leaned on can be scored against a known answer", {
  # A simulated regression plants some coefficients at exactly zero, so the null
  # predictors are the ones the importance ought to leave at the bottom. Nothing
  # about a forest guarantees a perfect ordering, so what is pinned is that the
  # planted predictors are the more important as a group.
  sim <- simulate_regression(seed = 1)
  rf <- do.call(fit_rf, c(sim$args, list(ntree = 200, cv = FALSE, seed = 1)))
  scored <- merge(rf$coefficients, sim$truth_term, by = "terms")

  planted <- scored$estimate[scored$beta != 0]
  null <- scored$estimate[scored$beta == 0]
  expect_gt(length(planted), 0L)
  expect_gt(min(planted), max(null))

  # A predictor that carried nothing can be left very slightly better by its own
  # permutation, so a negative importance is a value rather than a gap.
  expect_false(anyNA(scored$estimate))

  # Everything but the engine object survives the trip a JSON export would take.
  skip_if_not_installed("jsonlite")
  portable <- rf[setdiff(names(rf), "fit")]
  round_trip <- jsonlite::fromJSON(
    jsonlite::toJSON(portable, na = "string", digits = NA)
  )
  expect_identical(round_trip$terms, rf$terms)
  expect_equal(round_trip$coefficients$estimate, rf$coefficients$estimate)
  expect_equal(round_trip$fit_stats$oob_rmse, rf$fit_stats$oob_rmse)
})
