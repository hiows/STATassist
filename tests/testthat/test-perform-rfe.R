# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# A search is scored by resampling and the fold assignment is random, so what is
# pinned here is the boundary rather than which predictors won. That the ranking
# is the columns that were passed in and not the dummy columns a factor becomes,
# that it does not move when a predictor is measured in other units, that the
# three tables of the result describe the same search, and that `control_label`
# and `outcome_lv` are two ways of saying the same thing and cannot be used to say
# two.
#
# The reference values are the `rfe` object's own. It ran the search, so what
# these check is that the right thing is being read out of it and not that caret
# eliminates correctly.

sa_rfe_cars <- function() {
  mtcars[c("mpg", "wt", "hp", "disp", "qsec", "drat", "carb")]
}

sa_rfe_iris2 <- function() {
  out <- iris[iris$Species != "setosa", ]
  rownames(out) <- NULL
  out
}

sa_rfe_quick <- function(...) {
  perform_rfe(..., cv_method = "kfold", n_fold = 3, seed = 1)
}


test_that("the result has the shape the contract promises", {
  res <- sa_rfe_quick(sa_rfe_cars(), outcome = "mpg")

  expect_s3_class(res, "sa_selection")
  expect_s3_class(res, "sa_result")
  expect_named(res, c("analysis", "candidates", "design", "parameters",
                      "selected", "ranking", "profile", "resampling", "engine",
                      "fit", "metadata"))
  expect_identical(res$analysis, "rfe")
  expect_identical(res$ranking$candidates, res$candidates)
  expect_setequal(res$candidates, c("wt", "hp", "disp", "qsec", "drat", "carb"))
  # The row axis is the ranking and `design` keeps the order the columns arrived
  # in, which is the same division `fit_rf()` makes between `terms` and
  # `design$predictors`.
  expect_identical(res$design$predictors,
                   c("wt", "hp", "disp", "qsec", "drat", "carb"))
  expect_identical(res$design$outcome_type, "continuous")
  expect_null(res$design$outcome_lv)
  expect_identical(res$engine$method, "rfe")
  expect_identical(res$engine$label, "Linear regression")
  expect_identical(res$engine$importance, "absolute t statistic")

  expect_named(res$ranking, c("candidates", "estimate", "rank", "selected"))
  expect_identical(res$ranking$rank, seq_along(res$candidates))
  expect_false(is.unsorted(rev(res$ranking$estimate)))
  expect_identical(res$ranking$selected, res$candidates %in% res$selected)

  # A forest's two arguments belong to a search that grew one, so they are absent
  # rather than present and unused.
  expect_null(res$parameters$ntree)
  expect_null(res$parameters$nodesize)
  expect_identical(res$parameters$model, "linear")
  expect_identical(res$parameters$metric, "RMSE")
  expect_false(res$parameters$maximize)
})


test_that("the three tables describe the same search", {
  res <- sa_rfe_quick(sa_rfe_cars(), outcome = "mpg")

  expect_identical(res$selected, res$fit$optVariables)
  expect_identical(sum(res$profile$chosen), 1L)
  expect_identical(res$profile$n_vars[res$profile$chosen], res$fit$bestSubset)
  expect_identical(res$profile$n_vars[res$profile$chosen],
                   length(res$selected))
  # The top of the ranking is the selection, rather than a second answer that
  # usually agrees with it.
  expect_setequal(utils::head(res$candidates, length(res$selected)),
                  res$selected)

  expect_identical(res$profile$n_vars, res$fit$results$Variables)
  expect_null(res$profile$Variables)
  expect_identical(res$resampling$Variables,
                   rep(res$fit$bestSubset, nrow(res$resampling)))
  expect_identical(nrow(res$resampling), 3L)
})


test_that("the ranking is of the columns that were passed in", {
  cars <- sa_rfe_cars()
  cars$shape <- factor(rep(c("wide", "narrow", "tall", "flat"), length.out = 32))
  res <- sa_rfe_quick(cars, outcome = "mpg")

  # A four-level factor is three dummy columns to `lm()`, and none of the three
  # is something `fit_rf(predictors = )` would accept.
  expect_true("shape" %in% res$candidates)
  expect_false(any(grepl("^shape", setdiff(res$candidates, "shape"))))
  expect_identical(sort(res$candidates), sort(res$design$predictors))
  expect_identical(res$design$predictor_lv$shape,
                   c("flat", "narrow", "tall", "wide"))

  # Which is the point: what the search returns goes straight back into a fit.
  refit <- fit_linear_regression(cars, outcome = "mpg",
                                 predictors = res$selected, cv = FALSE)
  expect_true(all(res$selected %in% refit$design$predictors))
})


test_that("the ranking does not move when a predictor changes units", {
  # The reason the rank is the t statistic and not the coefficient. Measuring
  # weight in grams rather than in thousands of pounds divides its coefficient by
  # a million and leaves its statistic exactly where it was, so a ranking by
  # coefficient size would reorder here and this one does not.
  cars <- sa_rfe_cars()
  rescaled <- cars
  rescaled$wt <- rescaled$wt * 1e6

  as_is <- sa_rfe_quick(cars, outcome = "mpg")
  in_grams <- sa_rfe_quick(rescaled, outcome = "mpg")

  expect_identical(in_grams$candidates, as_is$candidates)
  expect_equal(in_grams$ranking$estimate, as_is$ranking$estimate)
  expect_identical(in_grams$selected, as_is$selected)
})


test_that("a two-class outcome is searched in the direction it was given", {
  res <- sa_rfe_quick(sa_rfe_iris2(), outcome = "Species", model = "rf",
                      ntree = 100)

  expect_identical(res$design$outcome_type, "two classes")
  expect_identical(res$design$outcome_lv, c("versicolor", "virginica"))
  expect_identical(res$design$n_events, 50L)
  expect_equal(res$design$event_rate, 0.5)
  expect_identical(res$parameters$metric, "Accuracy")
  expect_true(res$parameters$maximize)
  expect_identical(res$engine$label, "Random forest classification")
  expect_identical(res$engine$importance, "permutation importance")
  expect_identical(res$parameters$ntree, 100L)
  expect_identical(res$parameters$nodesize, 1L)

  # `control_label` names the reference on its own, and naming the other one
  # turns the result round: the events are the class it is not.
  other_way <- sa_rfe_quick(sa_rfe_iris2(), outcome = "Species",
                            control_label = "virginica", model = "rf",
                            ntree = 100)
  expect_identical(other_way$design$outcome_lv, c("virginica", "versicolor"))
  expect_identical(other_way$design$n_events, 50L)

  # Naming the reference twice and agreeing is the default rather than a
  # conflict, since `control_label` defaults to `outcome_lv[1]`.
  named_both <- sa_rfe_quick(sa_rfe_iris2(), outcome = "Species",
                             outcome_lv = c("versicolor", "virginica"),
                             model = "rf", ntree = 100)
  expect_identical(named_both$design$outcome_lv, res$design$outcome_lv)
  expect_identical(named_both$selected, res$selected)
})


test_that("a disagreement about the reference class is refused", {
  expect_error(
    sa_rfe_quick(sa_rfe_iris2(), outcome = "Species",
                 outcome_lv = c("versicolor", "virginica"),
                 control_label = "virginica", model = "rf", ntree = 50),
    "disagree"
  )
  expect_error(
    sa_rfe_quick(sa_rfe_iris2(), outcome = "Species",
                 control_label = "setosa", model = "rf", ntree = 50),
    "does not hold"
  )
  expect_error(
    sa_rfe_quick(sa_rfe_iris2(), outcome = "Species",
                 control_label = c("versicolor", "virginica"), model = "rf",
                 ntree = 50),
    "single level name"
  )
})


test_that("the model and the outcome have to agree", {
  expect_error(
    sa_rfe_quick(sa_rfe_iris2(), outcome = "Species"),
    "model = \"logistic\""
  )
  expect_error(
    sa_rfe_quick(sa_rfe_cars(), outcome = "mpg", model = "logistic"),
    "model = \"linear\""
  )
  # A numeric two-valued outcome is a regression unless something says otherwise,
  # and the something is named rather than left to be discovered.
  expect_message(
    sa_rfe_quick(mtcars[c("am", "wt", "hp", "qsec")], outcome = "am"),
    "searched as a regression"
  )
  expect_identical(
    suppressMessages(
      sa_rfe_quick(mtcars[c("am", "wt", "hp", "qsec")], outcome = "am",
                   control_label = "0", model = "rf",
                   ntree = 50)
    )$design$outcome_lv,
    c("0", "1")
  )
})


test_that("the sizes and the metric are the ones that were asked for", {
  res <- sa_rfe_quick(sa_rfe_cars(), outcome = "mpg",
                      subset_sizes = c(2, 4), metric = "MAE")

  # caret scores the full set whatever sizes are named, since keeping everything
  # is what the search is being compared against.
  expect_identical(res$profile$n_vars, c(2L, 4L, 6L))
  expect_identical(res$parameters$metric, "MAE")
  expect_false(res$parameters$maximize)
  expect_true(res$parameters$metric %in% res$engine$metrics)

  expect_true(
    sa_rfe_quick(sa_rfe_cars(), outcome = "mpg",
                 metric = "Rsquared")$parameters$maximize
  )

  # The default ladder is dense where the answer usually is and stops at the
  # number of candidates.
  expect_identical(res$profile$n_vars[nrow(res$profile)], 6L)
  expect_error(sa_rfe_quick(sa_rfe_cars(), outcome = "mpg",
                            subset_sizes = c(1, 99)),
               "must be in \\[1, 6\\]")
  expect_error(sa_rfe_quick(sa_rfe_cars(), outcome = "mpg",
                            subset_sizes = 1.5),
               "whole numbers")
  expect_error(sa_rfe_quick(sa_rfe_cars(), outcome = "mpg",
                            metric = "Accuracy"),
               "RMSE, Rsquared, MAE")
  expect_error(perform_rfe(sa_rfe_cars(), outcome = "mpg", cv_method = "kfold",
                           n_fold = 100),
               "exceeds the 32 usable")
})


test_that("a candidate with nothing in it is left out with a message", {
  cars <- sa_rfe_cars()
  cars$flat <- 1

  expect_message(res <- sa_rfe_quick(cars, outcome = "mpg"),
                 "single value")
  expect_identical(res$design$dropped_predictors, "flat")
  expect_false("flat" %in% res$candidates)
})


test_that("the same seed searches the same way and the stream is put back", {
  args <- list(sa_rfe_cars(), outcome = "mpg", cv_method = "kfold", n_fold = 3)

  once <- do.call(perform_rfe, c(args, list(seed = 1)))
  twice <- do.call(perform_rfe, c(args, list(seed = 1)))
  other <- do.call(perform_rfe, c(args, list(seed = 2)))

  expect_identical(once$selected, twice$selected)
  expect_identical(once$ranking, twice$ranking)
  expect_identical(once$profile, twice$profile)
  expect_false(isTRUE(all.equal(once$resampling$RMSE, other$resampling$RMSE)))

  set.seed(99)
  before <- .Random.seed
  invisible(do.call(perform_rfe, c(args, list(seed = 11))))
  expect_identical(before, .Random.seed)

  expect_identical(once$parameters$cv_method, "kfold")
  expect_identical(once$parameters$n_fold, 3L)
  expect_identical(once$parameters$n_repeat, NA_integer_)
})


test_that("everything but the fit writes out as JSON", {
  res <- sa_rfe_quick(sa_rfe_iris2(), outcome = "Species", model = "rf",
                      ntree = 100)

  portable <- res[setdiff(names(res), "fit")]
  unportable <- rapply(portable, function(v) is.function(v) || is.environment(v),
                       how = "unlist")
  expect_false(any(unportable))

  skip_if_not_installed("jsonlite")
  round_trip <- jsonlite::fromJSON(
    jsonlite::toJSON(portable, na = "string", digits = NA)
  )
  expect_identical(round_trip$analysis, "rfe")
  expect_identical(round_trip$candidates, res$candidates)
  expect_identical(round_trip$selected, res$selected)
  expect_identical(round_trip$design$outcome_lv, res$design$outcome_lv)
  expect_equal(round_trip$ranking$estimate, res$ranking$estimate)
  expect_identical(round_trip$profile$chosen, res$profile$chosen)
})


test_that("printing says what was searched and what was kept", {
  res <- sa_rfe_quick(sa_rfe_cars(), outcome = "mpg")

  expect_output(print(res), "<sa_selection> rfe")
  expect_output(print(res), "RMSE minimised")
  expect_output(print(res), "absolute t statistic")
  expect_output(print(res),
                paste0("selected : ", length(res$selected), " of 6"))
  expect_invisible(print(res))

  clf <- sa_rfe_quick(sa_rfe_iris2(), outcome = "Species", model = "rf",
                      ntree = 100)
  expect_output(print(clf), "modelling virginica against versicolor")
  expect_output(print(clf), "permutation importance")

  # One line per candidate up to `n`, and a count of what is left.
  expect_output(print(res, n = 2), "more candidate\\(s\\) in \\$ranking")
})


test_that("the contract refuses a result that disagrees with itself", {
  # These fire on a mistake inside the package rather than on a call, so they are
  # tested through the assembler directly.
  ranking <- data.frame(candidates = c("a", "b"), estimate = c(2, 1),
                        rank = 1:2, selected = c(TRUE, FALSE))
  profile <- data.frame(n_vars = 1:2, RMSE = c(3, 4), chosen = c(TRUE, FALSE))
  engine <- list(package = "caret", method = "rfe", label = "Linear regression",
                 metrics = "RMSE", importance = "absolute t statistic")
  ok <- function(...) {
    args <- list(analysis = "rfe", candidates = c("a", "b"),
                 design = list(), parameters = list(), selected = "a",
                 ranking = ranking, profile = profile, engine = engine,
                 fit = NULL)
    # Replaced rather than merged: modifyList() would fold a shortened `engine`
    # back into the full one and put the missing element back.
    changed <- list(...)
    args[names(changed)] <- changed
    do.call(sa_new_selection, args)
  }

  expect_s3_class(ok(), "sa_selection")
  expect_error(ok(analysis = "boruta"), "must be one of")
  expect_error(ok(candidates = c("b", "a")), "aligned with")
  expect_error(ok(selected = "c"), "not candidates")
  expect_error(ok(selected = c("a", "b")), "disagrees with")
  expect_error(ok(profile = data.frame(n_vars = 1:2, chosen = c(TRUE, TRUE))),
               "exactly one row")
  expect_error(ok(engine = engine[setdiff(names(engine), "importance")]),
               "missing `importance`")
})
