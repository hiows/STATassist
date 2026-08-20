# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# Nothing here is resampled, so unlike the elimination there is no fold assignment
# to pin and the path is a fact about the data. What these check is that it is read
# out of `step()` correctly, that the three tables describe the same search,
# and that the two claims the documentation makes about the criterion hold: that a
# difference of criteria is the same number on either scale, and that both criteria
# are reported whichever one is searching.
#
# The reference values are recomputed with `lm()` rather than written down, since
# what is in question is the arithmetic this package does on top of `step()` and
# not whether `step()` walks the path correctly.

sa_step_cars <- function() {
  mtcars[c("mpg", "wt", "hp", "disp", "qsec", "drat", "carb")]
}

sa_step_iris2 <- function() {
  out <- iris[iris$Species != "setosa", ]
  rownames(out) <- NULL
  out
}

sa_step_orthogonal <- function() {
  # `a` is exactly uncorrelated with `y`: its pattern repeats every four rows and
  # `y` cycles with it, so the fit that keeps it has the residual sum of squares of
  # the fit that does not and the criterion is charging for a parameter that bought
  # nothing at all. Whichever charge is levied, the search drops it.
  data.frame(y = rep(c(1, 2, 3, 4), 5), a = rep(c(0, 1, 1, 0), 5))
}


test_that("the result has the shape the contract promises", {
  res <- perform_stepwise(sa_step_cars(), outcome = "mpg")

  expect_s3_class(res, "sa_selection")
  expect_s3_class(res, "sa_result")
  expect_named(res, c("analysis", "candidates", "design", "parameters",
                      "selected", "ranking", "profile", "resampling", "engine",
                      "fit", "metadata"))
  expect_identical(res$analysis, "stepwise")
  expect_identical(res$ranking$candidates, res$candidates)
  expect_setequal(res$candidates, c("wt", "hp", "disp", "qsec", "drat", "carb"))
  expect_identical(res$design$predictors,
                   c("wt", "hp", "disp", "qsec", "drat", "carb"))
  expect_identical(res$design$outcome_type, "continuous")
  expect_null(res$design$outcome_lv)

  expect_identical(res$engine$package, "stats")
  expect_identical(res$engine$method, "step")
  expect_identical(res$engine$label, "Linear regression")
  expect_identical(res$engine$metrics, c("AIC", "BIC"))
  expect_identical(res$engine$importance,
                   "AIC increase when the predictor is left out")

  expect_named(res$ranking, c("candidates", "estimate", "rank", "selected"))
  expect_identical(res$ranking$rank, seq_along(res$candidates))
  expect_false(is.unsorted(rev(res$ranking$estimate)))
  expect_identical(res$ranking$selected, res$candidates %in% res$selected)

  expect_named(res$profile, c("n_vars", "AIC", "BIC", "step", "chosen"))
  expect_identical(res$parameters$model, "linear")
  expect_identical(res$parameters$criterion, "AIC")
  expect_false(res$parameters$maximize)
  expect_identical(res$parameters$direction, "backward")
})


test_that("the path and the two tables describe the same search", {
  res <- perform_stepwise(sa_step_cars(), outcome = "mpg")

  # `step()` is the authority on where the search stopped; `selected` is the same
  # set in the ranking's order rather than the formula's.
  expect_setequal(res$selected, attr(stats::terms(res$fit), "term.labels"))
  expect_identical(res$selected,
                   res$ranking$candidates[res$ranking$selected])

  expect_identical(sum(res$profile$chosen), 1L)
  expect_true(res$profile$chosen[nrow(res$profile)])
  expect_identical(res$profile$n_vars[res$profile$chosen],
                   length(res$selected))
  # A backward search only drops, so the path is a ladder down from every
  # candidate, and each step lowered the criterion it was searching on.
  expect_identical(res$profile$n_vars[1], length(res$candidates))
  expect_true(all(diff(res$profile$n_vars) == -1L))
  expect_true(all(diff(res$profile$AIC) < 0))
  expect_identical(res$profile$step[1], "")
  expect_true(all(grepl("^- ", res$profile$step[-1])))

  # The top of the ranking is the selection rather than a second answer that
  # usually agrees with it, and the sign says which side of it a candidate is on.
  expect_setequal(utils::head(res$candidates, length(res$selected)),
                  res$selected)
  expect_true(all(res$ranking$estimate[res$ranking$selected] > 0))
  expect_true(all(res$ranking$estimate[!res$ranking$selected] < 0))
})


test_that("the ranking prices leaving one candidate out", {
  res <- perform_stepwise(sa_step_cars(), outcome = "mpg")
  chosen <- stats::lm(
    stats::reformulate(res$selected, response = "mpg"), data = mtcars
  )
  priced <- stats::setNames(res$ranking$estimate, res$ranking$candidates)

  # The documented claim: the estimate is a difference of criteria, so it is the
  # same number whether it is taken on `AIC()`'s scale or on the `extractAIC()`
  # scale `drop1()` reports, the two differing by a constant on these rows.
  for (nm in res$selected) {
    without <- stats::lm(
      stats::reformulate(setdiff(res$selected, nm), response = "mpg"),
      data = mtcars
    )
    expect_equal(priced[[nm]], stats::AIC(without) - stats::AIC(chosen))
  }
  for (nm in setdiff(res$candidates, res$selected)) {
    with_it <- stats::lm(
      stats::reformulate(c(res$selected, nm), response = "mpg"), data = mtcars
    )
    expect_equal(priced[[nm]], stats::AIC(chosen) - stats::AIC(with_it))
  }

  # `profile` reports the fitted model's own AIC and BIC, which is what makes a
  # step of the path and a `fit_*()` result the same number.
  best <- res$profile[res$profile$chosen, ]
  expect_equal(best$AIC, stats::AIC(chosen))
  expect_equal(best$BIC, stats::BIC(chosen))
  expect_equal(
    best$AIC,
    fit_linear_regression(mtcars, outcome = "mpg", predictors = res$selected,
                          cv = FALSE)$fit_stats$aic
  )
})


test_that("the criterion is one search at two charges", {
  aic <- perform_stepwise(sa_step_cars(), outcome = "mpg")
  bic <- perform_stepwise(sa_step_cars(), outcome = "mpg", criterion = "BIC")

  expect_identical(aic$parameters$k, 2)
  expect_equal(bic$parameters$k, log(32))
  expect_identical(bic$parameters$criterion, "BIC")
  expect_identical(bic$engine$importance,
                   "BIC increase when the predictor is left out")

  # Both criteria are reported whichever one is searching, so a path chosen by one
  # can be read against the other. What each path is monotone in is its own.
  expect_true(all(diff(bic$profile$BIC) < 0))
  expect_named(aic$profile, names(bic$profile))
  # A heavier charge per parameter cannot keep more of them at any step it shares.
  expect_true(length(bic$selected) <= length(aic$selected))
  expect_true(all(bic$selected %in% aic$selected))
})


test_that("the direction decides where the search starts and how it moves", {
  cars <- sa_step_cars()
  fwd <- perform_stepwise(cars, outcome = "mpg", direction = "forward")
  back <- perform_stepwise(cars, outcome = "mpg", direction = "backward")
  both <- perform_stepwise(cars, outcome = "mpg", direction = "both")

  expect_identical(fwd$parameters$direction, "forward")
  expect_identical(fwd$profile$n_vars[1], 0L)
  expect_true(all(diff(fwd$profile$n_vars) == 1L))
  expect_true(all(grepl("^\\+ ", fwd$profile$step[-1])))

  expect_identical(back$profile$n_vars[1], length(back$candidates))

  # Every direction is bounded by the same two models, the intercept below and all
  # of `predictors` above.
  for (res in list(fwd, back, both)) {
    expect_true(all(res$profile$n_vars >= 0L))
    expect_true(all(res$profile$n_vars <= length(res$candidates)))
  }
})


test_that("a search that keeps nothing is an error saying so", {
  expect_error(perform_stepwise(sa_step_orthogonal(), outcome = "y"),
               "walked back to the intercept")
  # The charge that was levied is named, and a BIC search is told what AIC charges.
  expect_error(perform_stepwise(sa_step_orthogonal(), outcome = "y",
                            criterion = "BIC"),
               "`criterion = \"AIC\"` charges 2 per parameter")
  expect_error(perform_stepwise(sa_step_orthogonal(), outcome = "y"),
               "1 candidate\\(s\\)")
})


test_that("the model and the outcome have to agree", {
  expect_error(perform_stepwise(sa_step_iris2(), outcome = "Species"),
               "model = \"logistic\"")
  expect_error(perform_stepwise(sa_step_cars(), outcome = "mpg",
                            model = "logistic"),
               "model = \"linear\"")

  clf <- perform_stepwise(sa_step_iris2(), outcome = "Species",
                      control_label = "versicolor", model = "logistic")
  expect_identical(clf$design$outcome_type, "two classes")
  expect_identical(clf$design$outcome_lv, c("versicolor", "virginica"))
  expect_identical(clf$design$n_events, 50L)
  expect_identical(clf$engine$label, "Binomial logistic regression")
  expect_error(
    perform_stepwise(sa_step_iris2(), outcome = "Species", model = "logistic",
                 outcome_lv = c("versicolor", "virginica"),
                 control_label = "virginica"),
    "disagree"
  )

  # A numeric two-valued outcome is a regression unless something says otherwise,
  # and the something is named rather than left to be discovered.
  expect_message(
    perform_stepwise(mtcars[c("am", "wt", "hp", "qsec")], outcome = "am"),
    "searched as a regression"
  )
  expect_identical(
    perform_stepwise(mtcars[c("am", "wt", "hp", "qsec")], outcome = "am",
                 control_label = "0", model = "logistic")$design$outcome_lv,
    c("0", "1")
  )
})


test_that("a factor is one candidate however many terms it becomes", {
  cars <- mtcars[c("mpg", "wt", "qsec", "carb", "cyl")]
  cars$cyl <- factor(cars$cyl)
  res <- perform_stepwise(cars, outcome = "mpg",
                      predictors = c("wt", "cyl", "qsec", "carb"))

  expect_identical(length(res$candidates), 4L)
  expect_true("cyl" %in% res$candidates)
  expect_identical(res$design$predictor_lv$cyl, c("4", "6", "8"))
  # A column of the input rather than a dummy column of the model frame, which is
  # what lets the selection be handed straight back to a fit.
  expect_true(all(res$selected %in% names(cars)))
  expect_identical(
    sort(fit_linear_regression(cars, outcome = "mpg",
                               predictors = res$selected,
                               cv = FALSE)$design$predictors),
    sort(res$selected)
  )
})


test_that("a candidate with nothing in it is left out with a message", {
  cars <- sa_step_cars()
  cars$flat <- 1

  expect_message(res <- perform_stepwise(cars, outcome = "mpg"), "single value")
  expect_identical(res$design$dropped_predictors, "flat")
  expect_false("flat" %in% res$candidates)
})


test_that("nothing is resampled, so nothing is random", {
  once <- perform_stepwise(sa_step_cars(), outcome = "mpg")
  twice <- perform_stepwise(sa_step_cars(), outcome = "mpg")

  expect_identical(once$selected, twice$selected)
  expect_identical(once$ranking, twice$ranking)
  expect_identical(once$profile, twice$profile)

  # `resampling` is the slot that tells the two searches apart, and a search with
  # no random part has no seed to record.
  expect_null(once$resampling)
  expect_null(once$parameters$seed)
  expect_null(once$parameters$cv_method)

  set.seed(99)
  before <- .Random.seed
  invisible(perform_stepwise(sa_step_cars(), outcome = "mpg"))
  expect_identical(before, .Random.seed)
})


test_that("everything but the fit writes out as JSON", {
  res <- perform_stepwise(sa_step_iris2(), outcome = "Species",
                      control_label = "versicolor", model = "logistic")

  portable <- res[setdiff(names(res), "fit")]
  unportable <- rapply(portable, function(v) is.function(v) || is.environment(v),
                       how = "unlist")
  expect_false(any(unportable))

  skip_if_not_installed("jsonlite")
  round_trip <- jsonlite::fromJSON(
    jsonlite::toJSON(portable, na = "string", digits = NA)
  )
  expect_identical(round_trip$analysis, "stepwise")
  expect_identical(round_trip$candidates, res$candidates)
  expect_identical(round_trip$selected, res$selected)
  expect_identical(round_trip$design$outcome_lv, res$design$outcome_lv)
  expect_equal(round_trip$ranking$estimate, res$ranking$estimate)
  expect_identical(round_trip$profile$chosen, res$profile$chosen)
})


test_that("printing says what was searched and what was kept", {
  res <- perform_stepwise(sa_step_cars(), outcome = "mpg")

  expect_output(print(res), "<sa_selection> stepwise")
  expect_output(print(res), "backward search, AIC minimised at 2 per parameter")
  expect_output(print(res), "step\\(s\\)")
  expect_output(print(res), "AIC increase when the predictor is left out")
  expect_output(print(res),
                paste0("selected : ", length(res$selected), " of 6"))
  expect_invisible(print(res))

  # The resampling line of an elimination has nothing to say here, so it says
  # nothing rather than reporting a scheme that was never used.
  expect_output(print(res), "search   : Linear regression over 6 candidate")
  out <- utils::capture.output(print(res))
  expect_false(any(grepl("fold|resample|NA", out)))

  bic <- perform_stepwise(sa_step_cars(), outcome = "mpg", criterion = "BIC")
  expect_output(print(bic), "BIC minimised at 3.47 per parameter")
})
