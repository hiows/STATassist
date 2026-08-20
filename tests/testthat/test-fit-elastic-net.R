# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# Like the two unpenalized fits, this is a wrapper, so what is pinned is what it
# must not disturb: the coefficients have to be the ones `glmnet` gives at the
# same penalty on the same columns. Three things are specific to this one and none
# of them is arithmetic.
#
# The columns, first. `caret` hands `glmnet` `as.matrix()` of the predictor frame,
# which turns a three-level factor into one integer-coded numeric predictor and
# fits it without complaint, so the terms are checked against the terms
# `fit_linear_regression()` produces on the same data rather than merely counted.
#
# The tuning, second. This is the first model here that the resampling picks
# rather than merely scores, so what is reported has to be the candidate that won
# and not the one that happened to be scored first.
#
# The absence, third. A penalized estimate has no standard error, and the
# inference columns being missing from the table rather than present and `NA` is a
# promise about the contract rather than an accident of the fit: it is what tells
# a consumer, in one test, which kind of table it is holding.
#
# The reference values are `glmnet`'s own. It is what fitted the model, so this
# checks that the right column of the right path is being read, not that the
# arithmetic of the elastic net is right.

sa_enet_cars <- function() mtcars[c("mpg", "wt", "hp", "disp")]

sa_enet_iris2 <- function() {
  out <- iris[iris$Species != "setosa",
              c("Species", "Petal.Length", "Sepal.Width")]
  rownames(out) <- NULL
  out
}

# The matrix the fit is given, built the way `sa_design_matrix()` builds it, so a
# reference fit can be asked for on the same columns.
sa_enet_matrix <- function(data, predictors) {
  stats::model.matrix(~ ., data = data[predictors])[, -1L, drop = FALSE]
}


test_that("the result has the shape the contract promises", {
  fit <- fit_elastic_net(sa_enet_cars(), outcome = "mpg", penalty = "lasso",
                         lambda = 0.5, cv = FALSE)

  expect_s3_class(fit, "sa_model")
  expect_s3_class(fit, "sa_result")
  expect_named(fit, c("analysis", "terms", "design", "parameters",
                      "coefficients", "fit_stats", "performance", "resampling",
                      "engine", "fit", "metadata"))
  expect_identical(fit$analysis, "elastic_net")
  expect_identical(fit$terms, c("(Intercept)", "wt", "hp", "disp"))
  expect_identical(fit$coefficients$terms, fit$terms)
  expect_identical(fit$design$predictors, c("wt", "hp", "disp"))
  expect_identical(fit$design$outcome_type, "continuous")
  expect_identical(fit$engine$method, "glmnet")
  expect_identical(fit$engine$family, "gaussian")

  # The table carries what a penalized fit answers and nothing else: no standard
  # error means no inference columns at all, rather than six columns of `NA`.
  expect_named(fit$coefficients, c("terms", "estimate", "selected"))
  expect_length(intersect(sa_model_inference_columns(),
                          names(fit$coefficients)), 0L)

  # What it answers instead. The intercept is not penalized, so it is kept
  # whatever its value, and every other term is kept exactly when it is not zero.
  expect_type(fit$coefficients$selected, "logical")
  penalised <- fit$coefficients[fit$coefficients$terms != "(Intercept)", ]
  expect_identical(penalised$selected, penalised$estimate != 0)
  expect_true(fit$coefficients$selected[fit$coefficients$terms ==
                                          "(Intercept)"])
  expect_identical(fit$fit_stats$n_selected + fit$fit_stats$n_zero,
                   as.numeric(nrow(penalised)))

  expect_s3_class(fit$fit, "train")
  expect_length(stats::predict(fit$fit, newdata = mtcars), nrow(mtcars))
})


test_that("the coefficients are the ones glmnet gives at the same penalty", {
  d <- sa_enet_cars()
  fit <- fit_elastic_net(d, outcome = "mpg", penalty = "lasso", lambda = 0.5,
                         cv = FALSE)
  ref <- glmnet::glmnet(sa_enet_matrix(d, c("wt", "hp", "disp")), d$mpg,
                        family = "gaussian", alpha = 1)
  ref_coef <- as.matrix(stats::coef(ref, s = 0.5))

  expect_equal(fit$coefficients$estimate, unname(ref_coef[, 1L]))
  expect_identical(fit$coefficients$terms, rownames(ref_coef))

  # `coef()` on the fit answers with the chosen penalty as a named vector, not
  # with the sparse matrix of every fit on the lambda path.
  from_fit <- coef(fit$fit)
  expect_type(from_fit, "double")
  expect_null(dim(from_fit))
  expect_equal(from_fit, stats::setNames(fit$coefficients$estimate, fit$terms))

  # The fit statistics are the fitted values, so they are recomputable from them.
  fitted_value <- stats::predict(fit$fit, newdata = d)
  residual <- d$mpg - fitted_value
  expect_equal(fit$fit_stats$rmse, sqrt(mean(residual^2)))
  expect_equal(fit$fit_stats$mae, mean(abs(residual)))
  expect_equal(fit$fit_stats$r_squared,
               1 - sum(residual^2) / sum((d$mpg - mean(d$mpg))^2))
})


test_that("the result answers `coef()` with the table and `$fit` with the vector", {
  # A penalty hard enough to drop two of the three terms, so that what happens to
  # a dropped one is visible rather than assumed.
  fit <- fit_elastic_net(sa_enet_cars(), outcome = "mpg", penalty = "lasso",
                         lambda = 5, cv = FALSE)

  # Two objects, two answers. The result holds the table, so that is what it
  # gives; the engine object gives what the engine gives.
  expect_identical(coef(fit), fit$coefficients)
  expect_s3_class(coef(fit), "data.frame")

  estimate <- coef(fit$fit)
  expect_type(estimate, "double")
  expect_null(dim(estimate))
  expect_identical(names(estimate), fit$terms)

  # A dropped term keeps its row and its place in the vector, as an exact zero.
  # Leaving it out would make either one shorter than the model it came from, and
  # which terms survived is `selected` rather than the length.
  expect_length(estimate, nrow(fit$coefficients))
  expect_true(any(!fit$coefficients$selected))
  expect_identical(unname(estimate == 0), !fit$coefficients$selected)

  # `complete = FALSE` would have `coef.default()` index the table with a logical
  # matrix. The rows are what say which terms the model has, so it is ignored.
  expect_identical(coef(fit, complete = FALSE), fit$coefficients)
})


test_that("`penalty` names the corner of the model that is fitted", {
  args <- list(data = mtcars, outcome = "mpg", lambda = 2, cv = FALSE)
  lasso <- do.call(fit_elastic_net, c(args, list(penalty = "lasso")))
  ridge <- do.call(fit_elastic_net, c(args, list(penalty = "ridge")))

  expect_identical(lasso$parameters$penalty, "lasso")
  expect_identical(lasso$parameters$alpha, 1)
  expect_identical(ridge$parameters$alpha, 0)

  # The L1 penalty sets coefficients to exactly zero and the L2 one only shrinks
  # them, which is the whole difference between selecting and not selecting.
  expect_gt(lasso$fit_stats$n_zero, 0)
  expect_identical(ridge$fit_stats$n_zero, 0)
  expect_identical(ridge$fit_stats$n_selected,
                   as.numeric(length(ridge$terms) - 1L))
  expect_true(all(abs(ridge$coefficients$estimate) > 0))

  # `alpha` is read only by the mixture, but it is validated either way, so what
  # is accepted does not depend on `penalty`.
  expect_error(
    fit_elastic_net(mtcars, outcome = "mpg", penalty = "lasso", alpha = 2),
    "`alpha` must be in \\[0, 1\\]"
  )
  expect_error(
    fit_elastic_net(mtcars, outcome = "mpg", lambda = -1),
    "`lambda` must be in"
  )
  expect_error(
    fit_elastic_net(mtcars, outcome = "mpg", lambda = numeric(0)),
    "non-empty numeric vector"
  )
})


test_that("the reported candidate is the one the resampling chose", {
  fit <- fit_elastic_net(mtcars, outcome = "mpg", penalty = "lasso",
                         lambda = c(0.01, 0.1, 0.5, 1, 2),
                         cv_method = "kfold", n_fold = 5, seed = 1)

  expect_identical(nrow(fit$performance), 5L)
  expect_identical(fit$parameters$n_candidates, 5L)
  expect_identical(fit$parameters$lambda, fit$fit$bestTune$lambda)
  expect_identical(fit$parameters$alpha, fit$fit$bestTune$alpha)

  # The table is in the order the candidates were scored, so the chosen one is
  # not the first row. Nothing may quietly report the first row as the winner.
  best_row <- which(fit$performance$lambda == fit$parameters$lambda)
  expect_gt(best_row, 1L)
  expect_equal(fit$performance$RMSE[best_row], min(fit$performance$RMSE))
  expect_output(print(fit), paste0("lambda = ", fit$parameters$lambda))
  expect_output(print(fit),
                paste0("RMSE = ",
                       sa_fmt_num(fit$performance$RMSE[best_row], 3)))

  # `fit` is the one slot that cannot be serialised, so nothing the object
  # reports about itself may be reachable only through it. Without this, the
  # penalty line named the chosen lambda while the resample line beneath it
  # reported the first candidate's score.
  stripped <- fit
  stripped$fit <- NULL
  expect_output(print(stripped), paste0("lambda = ", fit$parameters$lambda))
  expect_output(print(stripped),
                paste0("RMSE = ",
                       sa_fmt_num(fit$performance$RMSE[best_row], 3)))

  # Both grids are searched together, so the count is their product.
  mixed <- fit_elastic_net(mtcars, outcome = "mpg", alpha = c(0, 0.5, 1),
                           lambda = c(0.1, 1), cv_method = "kfold", n_fold = 3,
                           seed = 1)
  expect_identical(nrow(mixed$performance), 6L)
  expect_true(mixed$parameters$alpha %in% c(0, 0.5, 1))

  # Nothing scores the candidates when nothing is resampled, so there is nothing
  # to choose between and the grid has to name one.
  expect_error(
    fit_elastic_net(mtcars, outcome = "mpg", cv = FALSE),
    "must hold one candidate"
  )
  expect_error(
    fit_elastic_net(mtcars, outcome = "mpg", penalty = "ridge",
                    lambda = c(0.1, 1), cv = FALSE),
    "must hold one candidate"
  )
  plain <- fit_elastic_net(mtcars, outcome = "mpg", penalty = "ridge",
                           lambda = 1, cv = FALSE)
  expect_null(plain$performance)
  expect_null(plain$resampling)
  expect_identical(plain$parameters$n_candidates, 1L)
})


test_that("a factor predictor becomes one term per level beyond the first", {
  cars <- mtcars
  cars$cyl <- factor(cars$cyl)
  args <- list(data = cars, outcome = "mpg", predictors = c("wt", "cyl"),
               cv = FALSE)
  fit <- do.call(fit_elastic_net, c(args, list(penalty = "ridge",
                                               lambda = 0.1)))

  # The terms of the unpenalized fit on the same data. What this guards against
  # is `caret` passing the frame to `glmnet` as `as.matrix()`, which codes a
  # three-level factor as one evenly spaced numeric predictor: that fit succeeds
  # and reports a term called `cyl`.
  expect_identical(fit$terms,
                   fit_linear_regression(cars, outcome = "mpg",
                                         predictors = c("wt", "cyl"),
                                         cv = FALSE)$terms)
  expect_identical(fit$terms, c("(Intercept)", "wt", "cyl6", "cyl8"))
  expect_identical(fit$design$predictors, c("wt", "cyl"))

  ref <- glmnet::glmnet(sa_enet_matrix(cars, c("wt", "cyl")), cars$mpg,
                        family = "gaussian", alpha = 0)
  expect_equal(fit$coefficients$estimate,
               unname(as.matrix(stats::coef(ref, s = 0.1))[, 1L]))

  # A character column is the same model as the factor it stands for, since it is
  # turned into one before the design matrix is built.
  cars$cyl <- as.character(cars$cyl)
  chr <- do.call(fit_elastic_net, c(args[c("outcome", "predictors", "cv")],
                                    list(data = cars, penalty = "ridge",
                                         lambda = 0.1)))
  expect_identical(chr$coefficients, fit$coefficients)

  # Each term is penalized on its own, so one level of a factor can be dropped
  # while another survives.
  cars$cyl <- factor(cars$cyl)
  lasso <- do.call(fit_elastic_net, c(args, list(penalty = "lasso",
                                                 lambda = 1.5)))
  expect_identical(lasso$terms, fit$terms)
  expect_gt(lasso$fit_stats$n_zero, 0)
})


test_that("the outcome decides which of the two models is fitted", {
  # A number is a regression, and a set of labels is a classification. Nothing
  # else is asked.
  numeric_outcome <- fit_elastic_net(sa_enet_cars(), outcome = "mpg",
                                     penalty = "ridge", lambda = 1,
                                     cv = FALSE)
  expect_identical(numeric_outcome$design$outcome_type, "continuous")
  expect_identical(numeric_outcome$engine$family, "gaussian")
  expect_null(numeric_outcome$design$outcome_lv)
  expect_null(numeric_outcome$coefficients$odds_ratio)

  labels <- fit_elastic_net(sa_enet_iris2(), outcome = "Species",
                            penalty = "ridge", lambda = 0.1, cv = FALSE)
  expect_identical(labels$design$outcome_type, "two classes")
  expect_identical(labels$engine$family, "binomial")
  expect_identical(labels$design$outcome_lv, c("versicolor", "virginica"))
  expect_identical(labels$engine$metrics, c("Accuracy", "Kappa"))
  # The odds ratio is the one column the classification adds, and it does not
  # bring the inference columns back with it.
  expect_named(labels$coefficients,
               c("terms", "estimate", "selected", "odds_ratio"))

  # A numeric column holding two values is the one case the outcome cannot
  # settle on its own, so the guess is announced and `outcome_lv` overrules it.
  d <- mtcars[c("am", "wt", "hp")]
  expect_message(
    as_numbers <- fit_elastic_net(d, outcome = "am", penalty = "ridge",
                                  lambda = 0.1, cv = FALSE),
    "fitted as a regression"
  )
  expect_identical(as_numbers$design$outcome_type, "continuous")

  as_classes <- fit_elastic_net(d, outcome = "am", outcome_lv = c("0", "1"),
                                penalty = "ridge", lambda = 0.1, cv = FALSE)
  expect_identical(as_classes$design$outcome_type, "two classes")
  expect_identical(as_classes$design$n_events, 13L)

  # An outcome with more than two classes is the same refusal the logistic
  # regression makes, in this model's name.
  expect_error(
    fit_elastic_net(iris, outcome = "Species", penalty = "ridge", lambda = 1,
                    cv = FALSE),
    "holds 3 classes, but an elastic net models two"
  )
  expect_error(
    fit_elastic_net(iris, outcome = "Species", penalty = "ridge", lambda = 1,
                    outcome_lv = c("setosa", "virginica"), cv = FALSE),
    "would be silently left out"
  )
  expect_error(
    fit_elastic_net(mtcars, outcome = "mpg", predictors = c("wt", "mpg"),
                    lambda = 1, penalty = "ridge", cv = FALSE),
    "predict from the answer"
  )
})


test_that("a two-class fit reads in the direction `outcome_lv` states", {
  d <- sa_enet_iris2()
  args <- list(data = d, outcome = "Species", penalty = "lasso",
               lambda = 0.01, cv = FALSE)
  versi <- do.call(fit_elastic_net,
                   c(args, list(outcome_lv = c("versicolor", "virginica"))))
  virgi <- do.call(fit_elastic_net,
                   c(args, list(outcome_lv = c("virginica", "versicolor"))))

  ref <- glmnet::glmnet(
    sa_enet_matrix(d, c("Petal.Length", "Sepal.Width")),
    factor(d$Species, levels = c("versicolor", "virginica")),
    family = "binomial", alpha = 1
  )
  expect_equal(versi$coefficients$estimate,
               unname(as.matrix(stats::coef(ref, s = 0.01))[, 1L]))

  # Swapping the levels turns every coefficient around and inverts every odds
  # ratio, which is the rule `group_lv` follows in a comparison and `outcome_lv`
  # follows in a logistic regression.
  expect_equal(versi$coefficients$estimate, -virgi$coefficients$estimate)
  expect_equal(versi$coefficients$odds_ratio, exp(versi$coefficients$estimate))
  expect_equal(versi$coefficients$odds_ratio,
               1 / virgi$coefficients$odds_ratio)
  expect_identical(versi$design$n_events, 50L)
  expect_equal(versi$design$event_rate, 0.5)

  # Longer petals make `virginica` more likely, so its odds ratio is above 1
  # when `virginica` is the level being modelled.
  petal <- versi$coefficients$terms == "Petal.Length"
  expect_gt(versi$coefficients$odds_ratio[petal], 1)
  expect_lt(virgi$coefficients$odds_ratio[petal], 1)

  # Sorting is what `NULL` does, and it puts `versicolor` first here.
  sorted <- do.call(fit_elastic_net, args)
  expect_identical(sorted$coefficients, versi$coefficients)

  # `type = "response"` is the probability of the class the coefficients
  # describe, as it is for the unpenalized fit.
  resp <- stats::predict(versi$fit, newdata = d, type = "response")
  expect_equal(resp, stats::predict(versi$fit, newdata = d,
                                    type = "prob")[["virginica"]])
  expect_equal(stats::predict(virgi$fit, newdata = d, type = "response"),
               1 - resp)

  # The deviance is the one those probabilities give, and McFadden's R-squared
  # is it against the intercept-only model.
  event <- as.numeric(d$Species == "virginica")
  expect_equal(versi$fit_stats$residual_deviance,
               -2 * sum(event * log(resp) + (1 - event) * log(1 - resp)))
  expect_equal(versi$fit_stats$mcfadden_r2,
               1 - versi$fit_stats$residual_deviance /
                 versi$fit_stats$null_deviance)
})


test_that("the rows and columns that entered the model are the reported ones", {
  d <- mtcars[c("mpg", "wt", "hp")]
  d$wt[1:3] <- NA
  d$mpg[4] <- NA
  fit <- fit_elastic_net(d, outcome = "mpg", penalty = "ridge", lambda = 0.1,
                         cv = FALSE)

  expect_identical(fit$design$n_obs, 32L)
  expect_identical(fit$design$n_used, 28L)
  expect_identical(fit$design$n_dropped, 4L)

  # Dropped before the folds are drawn, so the fit is the one the complete rows
  # give on their own.
  complete <- d[stats::complete.cases(d), ]
  ref <- glmnet::glmnet(sa_enet_matrix(complete, c("wt", "hp")), complete$mpg,
                        family = "gaussian", alpha = 0)
  expect_equal(fit$coefficients$estimate,
               unname(as.matrix(stats::coef(ref, s = 0.1))[, 1L]))

  flat <- mtcars[c("mpg", "wt", "hp")]
  flat$flat <- 1
  expect_message(
    dropped <- fit_elastic_net(flat, outcome = "mpg", penalty = "ridge",
                               lambda = 0.1, cv = FALSE),
    "single value cannot contribute"
  )
  expect_identical(dropped$design$dropped_predictors, "flat")
  expect_identical(dropped$terms, c("(Intercept)", "wt", "hp"))
})


test_that("a model of one term is refused rather than handed to the engine", {
  # A penalty has nothing to trade off between one coefficient and no other, and
  # `glmnet` refuses a one-column matrix in a message about matrices. The
  # refusal is made here, in terms of the model that was asked for.
  expect_error(
    fit_elastic_net(mtcars, outcome = "mpg", predictors = "wt",
                    penalty = "ridge", lambda = 1, cv = FALSE),
    "the model has 1 term"
  )

  # It counts terms rather than predictors, so a single factor predictor with
  # three levels is two terms and fits.
  cars <- mtcars
  cars$cyl <- factor(cars$cyl)
  fit <- fit_elastic_net(cars, outcome = "mpg", predictors = "cyl",
                         penalty = "ridge", lambda = 1, cv = FALSE)
  expect_identical(fit$terms, c("(Intercept)", "cyl6", "cyl8"))

  # And a predictor dropped for taking a single value can be what leaves one
  # term behind, which is why the message about it comes first.
  flat <- mtcars[c("mpg", "wt")]
  flat$flat <- 1
  expect_error(
    suppressMessages(
      fit_elastic_net(flat, outcome = "mpg", penalty = "ridge", lambda = 1,
                      cv = FALSE)
    ),
    "the model has 1 term"
  )
})


test_that("`seed` fixes the folds and leaves the stream as it was", {
  args <- list(data = mtcars, outcome = "mpg", predictors = c("wt", "hp"),
               penalty = "lasso", lambda = c(0.1, 1), cv_method = "kfold",
               n_fold = 4)
  once <- do.call(fit_elastic_net, c(args, list(seed = 11)))
  twice <- do.call(fit_elastic_net, c(args, list(seed = 11)))
  other <- do.call(fit_elastic_net, c(args, list(seed = 12)))

  expect_identical(once$resampling, twice$resampling)
  expect_identical(once$parameters, twice$parameters)
  expect_false(isTRUE(all.equal(once$resampling$RMSE, other$resampling$RMSE)))

  set.seed(99)
  before <- .Random.seed
  invisible(do.call(fit_elastic_net, c(args, list(seed = 11))))
  expect_identical(before, .Random.seed)

  # The resampling arguments are recorded as the scheme used them, the same way
  # the unpenalized fits record them.
  expect_identical(once$parameters$cv_method, "kfold")
  expect_identical(once$parameters$n_fold, 4L)
  expect_identical(once$parameters$n_repeat, NA_integer_)
  expect_identical(
    fit_elastic_net(mtcars, outcome = "mpg", predictors = c("wt", "hp"),
                    penalty = "ridge", lambda = 1,
                    cv = FALSE)$parameters$cv_method,
    NA_character_
  )
})


test_that("print says what was fitted, at what penalty, and how it did", {
  fit <- fit_elastic_net(mtcars, outcome = "mpg", predictors = c("wt", "hp"),
                         penalty = "lasso", lambda = c(0.1, 0.5, 1),
                         cv_method = "kfold", n_fold = 4, seed = 1)

  expect_output(print(fit), "<sa_model> elastic_net")
  expect_output(print(fit), "penalty  : lasso")
  expect_output(print(fit), "chosen from 3 candidate\\(s\\)")
  expect_output(print(fit), "r_squared")
  expect_output(print(fit), "n_selected")
  expect_output(print(fit), "RMSE")
  expect_output(expect_invisible(print(fit)))

  # There is no interval and no p-value to print, so the line says what the
  # penalty did with the term instead of printing three NAs.
  expect_output(print(fit), "selected")
  expect_output(print(fit, n = 1), "and 2 more term\\(s\\)")
  output <- utils::capture.output(print(fit))
  expect_false(any(grepl("p = NA", output, fixed = TRUE)))
  expect_false(any(grepl("conf_level", output, fixed = TRUE)))

  dropped <- fit_elastic_net(mtcars, outcome = "mpg", penalty = "lasso",
                             lambda = 2, cv = FALSE)
  expect_output(print(dropped), "dropped")
  expect_output(print(dropped), "no resampling")

  # The unpenalized fits still print their interval and their confidence level,
  # which is what the two branches of the coefficient line are for.
  linear <- utils::capture.output(
    print(fit_linear_regression(mtcars, outcome = "mpg", predictors = "wt",
                                cv = FALSE))
  )
  expect_true(any(grepl("p = ", linear, fixed = TRUE)))
  expect_true(any(grepl("conf_level", linear, fixed = TRUE)))
})


test_that("what the penalty selected can be scored against a known answer", {
  # A simulated regression plants some coefficients at exactly zero, so both
  # kinds of mistake are defined: a null predictor the lasso kept and a planted
  # one it dropped.
  sim <- simulate_regression(seed = 1)
  lasso <- do.call(fit_elastic_net,
                   c(sim$args, list(penalty = "lasso", lambda = 0.5,
                                    cv = FALSE)))
  scored <- merge(lasso$coefficients, sim$truth_term, by = "terms")

  expect_identical(nrow(scored), nrow(sim$truth_term))
  expect_true(all(scored$selected[scored$beta != 0]))
  expect_gt(sum(!scored$selected[scored$beta == 0]), 0)

  # Everything but the engine object survives the trip a JSON export would take.
  skip_if_not_installed("jsonlite")
  portable <- lasso[setdiff(names(lasso), "fit")]
  round_trip <- jsonlite::fromJSON(
    jsonlite::toJSON(portable, na = "string", digits = NA)
  )
  expect_identical(round_trip$terms, lasso$terms)
  expect_equal(round_trip$coefficients$estimate, lasso$coefficients$estimate)
  expect_identical(round_trip$coefficients$selected, lasso$coefficients$selected)
})
