# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# A machine is deterministic given its data and its two arguments, so more can be
# pinned here than for a forest. What cannot be pinned by comparing one fit with
# another is the importance, because it is measured by shuffling and the shuffles
# are draws. It is checked against its own definition instead: `sa_svm_importance()`
# is called with the stream at a known point and the same arithmetic is written out
# by hand beside it, so what this fixes is which metric is measured, in which
# direction, and against which baseline.
#
# The rest of the boundary is where this model differs from the two it sits between.
# The terms are dummy columns and not predictors, as in `fit_elastic_net()` and
# unlike `fit_rf()`; `fit_stats` is measured on the rows the machine was fitted to,
# as in `fit_elastic_net()` and unlike `fit_rf()`; and the table carries no
# inference column and no `selected` either, which is what tells an importance
# table from a penalized one.
#
# `C` is small and named throughout. The default grid is eight values, and nothing
# here is about how long a search takes.

sa_svm_cars <- function() mtcars[c("mpg", "wt", "hp", "disp")]

sa_svm_iris2 <- function() {
  out <- iris[iris$Species != "setosa", ]
  rownames(out) <- NULL
  out
}


test_that("the result has the shape the contract promises", {
  fit <- fit_svm(sa_svm_cars(), outcome = "mpg", C = 1, cv = FALSE, seed = 1)

  expect_s3_class(fit, "sa_model")
  expect_s3_class(fit, "sa_result")
  expect_named(fit, c("analysis", "terms", "design", "parameters",
                      "coefficients", "fit_stats", "performance", "resampling",
                      "engine", "fit", "metadata"))
  expect_identical(fit$analysis, "svm")
  expect_identical(fit$coefficients$terms, fit$terms)
  expect_setequal(fit$terms, c("wt", "hp", "disp"))
  expect_identical(fit$design$predictors, c("wt", "hp", "disp"))
  expect_identical(fit$design$outcome_type, "continuous")
  expect_identical(fit$engine$method, "svmRadial")
  expect_identical(fit$engine$kernel, "radial")
  expect_identical(fit$engine$label,
                   "Support vector machine regression (radial basis kernel)")

  # The table holds what a machine answers and nothing else. There is no effect per
  # unit of a term, so no standard error, so no inference columns at all rather
  # than six columns of `NA`; no `selected`, since nothing was penalized; and no
  # intercept, since the surface is not a line through one.
  expect_named(fit$coefficients, c("terms", "estimate"))
  expect_length(intersect(sa_model_inference_columns(),
                          names(fit$coefficients)), 0L)
  expect_null(fit$coefficients$selected)
  expect_null(fit$coefficients$impurity)
  expect_null(fit$coefficients$odds_ratio)
  expect_false("(Intercept)" %in% fit$terms)

  # The slot names do not depend on whether anything was resampled.
  cv_fit <- fit_svm(sa_svm_cars(), outcome = "mpg", C = 1, cv_method = "kfold",
                    n_fold = 3, seed = 1)
  expect_named(cv_fit, names(fit))

  expect_s3_class(fit$fit, "train")
  expect_s4_class(fit$fit$finalModel, "ksvm")
})


test_that("`estimate` is the loss when the term is shuffled", {
  d <- sa_svm_cars()
  fit <- fit_svm(d, outcome = "mpg", C = 1, cv = FALSE, seed = 1)
  x <- sa_design_matrix(d[c("wt", "hp", "disp")])
  y <- d$mpg

  # The same arithmetic written out beside the helper, with the stream at the same
  # point: the baseline is scored before any shuffle, so the draws are the
  # `n_permute` samples of each column in the order the columns arrive.
  set.seed(7)
  got <- sa_svm_importance(fit$fit, x, y, classify = FALSE, n_permute = 3)
  set.seed(7)
  as_fitted <- sqrt(mean((y - stats::predict(fit$fit, newdata = x))^2))
  manual <- vapply(colnames(x), function(nm) {
    rose_to <- vapply(seq_len(3), function(i) {
      permuted <- x
      permuted[, nm] <- sample(permuted[, nm])
      sqrt(mean((y - stats::predict(fit$fit, newdata = permuted))^2))
    }, numeric(1))
    mean(rose_to) - as_fitted
  }, numeric(1))

  expect_equal(got$estimate, unname(manual[got$terms]))
  # Measured in the metric the resampling tuned on, and against the machine's own
  # in-sample score, which is the one the result reports.
  expect_equal(as_fitted, fit$fit_stats$rmse)

  # Descending, so the first rows of the table are the terms the machine leaned on.
  # The order the columns were read in is `engine$x_names`.
  expect_identical(fit$terms, fit$terms[order(fit$coefficients$estimate,
                                              decreasing = TRUE)])
  expect_setequal(fit$terms, fit$engine$x_names)
  expect_identical(fit$parameters$n_permute, 10L)
  expect_identical(
    fit_svm(d, outcome = "mpg", C = 1, n_permute = 3,
            cv = FALSE, seed = 1)$parameters$n_permute,
    3L
  )

  # The result answers `coef()` with its own table, as the other models do.
  expect_identical(coef(fit), fit$coefficients)
  expect_s3_class(coef(fit), "data.frame")

  # `coef()` on the engine object is where a machine differs, and it fails one step
  # earlier than a forest does: a `ksvm` object is S4, so `coef.default()` cannot
  # even look for the element it would have found empty.
  expect_error(coef(fit$fit), "estimates no coefficients")
})


test_that("a classification's importance is the accuracy it lost", {
  d <- sa_svm_iris2()
  lv <- c("versicolor", "virginica")
  clf <- fit_svm(d, outcome = "Species", outcome_lv = lv, C = 1, cv = FALSE,
                 seed = 1)
  x <- sa_design_matrix(d[setdiff(names(d), "Species")])
  y <- factor(as.character(d$Species), levels = lv)

  set.seed(3)
  got <- sa_svm_importance(clf$fit, x, y, classify = TRUE, n_permute = 2)
  set.seed(3)
  as_fitted <- mean(as.character(stats::predict(clf$fit, newdata = x)) ==
                      as.character(y))
  manual <- vapply(colnames(x), function(nm) {
    fell_to <- vapply(seq_len(2), function(i) {
      permuted <- x
      permuted[, nm] <- sample(permuted[, nm])
      mean(as.character(stats::predict(clf$fit, newdata = permuted)) ==
             as.character(y))
    }, numeric(1))
    as_fitted - mean(fell_to)
  }, numeric(1))

  # Accuracy falls where RMSE rose, so the subtraction is the other way round and a
  # term the machine needed is positive in both models.
  expect_equal(got$estimate, unname(manual[got$terms]))
  expect_equal(as_fitted, clf$fit_stats$accuracy)
  expect_identical(clf$engine$label,
                   "Support vector machine classification (radial basis kernel)")
})


test_that("the fit statistics are the in-sample ones and the support vectors", {
  d <- sa_svm_cars()
  fit <- fit_svm(d, outcome = "mpg", C = 1, cv = FALSE, seed = 1)
  fitted_value <- stats::predict(fit$fit,
                                 newdata = sa_design_matrix(d[-1L]))
  residual <- d$mpg - fitted_value

  expect_named(fit$fit_stats, c("r_squared", "rmse", "mae", "n_support_vector",
                               "support_vector_rate"))
  expect_equal(fit$fit_stats$rmse, sqrt(mean(residual^2)))
  expect_equal(fit$fit_stats$mae, mean(abs(residual)))
  expect_equal(fit$fit_stats$r_squared,
               1 - sum(residual^2) / sum((d$mpg - mean(d$mpg))^2))

  # A machine has no out-of-bag rows to score itself on, which is the one place it
  # is weaker than a forest: these describe the surface it found, and `performance`
  # is the number measured on rows it had not seen. The prefix says so by its
  # absence.
  expect_null(fit$fit_stats$oob_rmse)
  expect_equal(fit$fit_stats$n_support_vector,
               as.numeric(kernlab::nSV(fit$fit$finalModel)))
  expect_equal(fit$fit_stats$support_vector_rate,
               fit$fit_stats$n_support_vector / fit$design$n_used)
})


test_that("a classification's statistics read `outcome_lv`", {
  d <- sa_svm_iris2()
  lv <- c("versicolor", "virginica")
  clf <- fit_svm(d, outcome = "Species", outcome_lv = lv, C = 1, cv = FALSE,
                 seed = 1)
  called <- as.character(stats::predict(clf, newdata = d))

  expect_named(clf$fit_stats, c("accuracy", "error", "kappa", "sensitivity",
                               "specificity", "n_support_vector",
                               "support_vector_rate"))
  expect_equal(clf$fit_stats$accuracy, mean(called == d$Species))
  expect_equal(clf$fit_stats$error, 1 - clf$fit_stats$accuracy)

  # Sensitivity is the share of `outcome_lv[2]` found and specificity the share of
  # `outcome_lv[1]` left alone, which is the direction the whole family reads in.
  expect_equal(clf$fit_stats$sensitivity,
               mean(called[d$Species == lv[2]] == lv[2]))
  expect_equal(clf$fit_stats$specificity,
               mean(called[d$Species == lv[1]] == lv[1]))

  positive <- d$Species == lv[2]
  hit <- called == lv[2]
  expected <- mean(hit) * mean(positive) + mean(!hit) * mean(!positive)
  expect_equal(clf$fit_stats$kappa,
               (clf$fit_stats$accuracy - expected) / (1 - expected))

  # The same rule read the other way round: with the levels swapped, what was the
  # sensitivity is now measured on the other class.
  other <- fit_svm(d, outcome = "Species", outcome_lv = rev(lv), C = 1,
                   cv = FALSE, seed = 1)
  other_called <- as.character(stats::predict(other, newdata = d))
  expect_equal(other$fit_stats$sensitivity,
               mean(other_called[d$Species == lv[1]] == lv[1]))
  expect_identical(other$design$outcome_lv, rev(lv))

  # `type = "response"` is the probability of the class the sensitivity is about,
  # which a machine can only answer because `prob.model` reached the engine.
  resp <- stats::predict(clf, newdata = d, type = "response")
  expect_length(resp, nrow(d))
  expect_true(all(resp >= 0 & resp <= 1))
  expect_equal(resp, stats::predict(clf, newdata = d, type = "prob")[[lv[2]]])
})


test_that("`sigma` is read off the data unless the caller names one", {
  d <- sa_svm_cars()
  x <- sa_design_matrix(d[-1L])

  # The middle of the three widths `sigest()` reports, drawn from the same stream
  # position the fit draws it from, which is why the seed is what makes the two
  # calls agree.
  set.seed(1)
  expected <- unname(kernlab::sigest(x, scaled = TRUE)[2L])
  fit <- fit_svm(d, outcome = "mpg", C = 1, cv = FALSE, seed = 1)
  expect_equal(fit$parameters$sigma, expected)
  expect_identical(fit$parameters$kernel, "radial")
  expect_identical(fit$parameters$n_candidates, 1L)

  named <- fit_svm(d, outcome = "mpg", C = 1, sigma = 0.25, cv = FALSE, seed = 1)
  expect_equal(named$parameters$sigma, 0.25)
  expect_equal(named$fit$bestTune$sigma, 0.25)
})


test_that("the reported candidate is the one the resampling chose", {
  fit <- fit_svm(sa_svm_cars(), outcome = "mpg", C = c(0.5, 2, 8),
                 cv_method = "kfold", n_fold = 5, seed = 1)

  expect_identical(nrow(fit$performance), 3L)
  expect_identical(fit$parameters$n_candidates, 3L)
  expect_equal(fit$parameters$C, fit$fit$bestTune$C)
  expect_equal(fit$parameters$sigma, fit$fit$bestTune$sigma)
  best <- which(fit$performance$C == fit$parameters$C)
  expect_equal(fit$performance$RMSE[best], min(fit$performance$RMSE))

  # Both arguments are tuned, so both are columns of the grid rather than scalars
  # passed through to the engine the way a forest's `ntree` is.
  expect_true(all(c("sigma", "C") %in% names(fit$performance)))
  expect_length(fit$parameters$C, 1L)

  plain <- fit_svm(sa_svm_cars(), outcome = "mpg", C = 1, cv = FALSE, seed = 1)
  expect_null(plain$performance)
  expect_null(plain$resampling)
  # Unlike a forest, nothing held out is reported when nothing was resampled: the
  # machine saw every row it was scored on.
  expect_false(is.null(plain$fit_stats$rmse))
})


test_that("a grid the machine could not honour is refused by name", {
  expect_error(
    fit_svm(mtcars, outcome = "mpg", C = c(1, 2), cv = FALSE),
    "must hold one candidate"
  )
  expect_error(
    fit_svm(mtcars, outcome = "mpg", C = 1, sigma = c(0.1, 0.2), cv = FALSE),
    "must hold one candidate"
  )
  # `kernlab` fits a machine at either zero and answers something that is not a
  # machine, so the refusal has to be here.
  expect_error(fit_svm(mtcars, outcome = "mpg", C = 0, cv = FALSE),
               "`C` must be above 0")
  expect_error(fit_svm(mtcars, outcome = "mpg", C = 1, sigma = 0, cv = FALSE),
               "`sigma` must be above 0")
  expect_error(fit_svm(mtcars, outcome = "mpg", C = -1, cv = FALSE),
               "`C` must be in")
  expect_error(fit_svm(mtcars, outcome = "mpg", C = 1, n_permute = 0,
                       cv = FALSE),
               "`n_permute` must be in")
  expect_error(
    fit_svm(mtcars, outcome = "mpg", predictors = "wt", C = 1,
            cv_method = "kfold", n_fold = 40),
    "exceeds the 32 usable observation"
  )
})


test_that("the outcome decides which machine is fitted", {
  numeric_outcome <- fit_svm(sa_svm_cars(), outcome = "mpg", C = 1, cv = FALSE,
                             seed = 1)
  expect_identical(numeric_outcome$design$outcome_type, "continuous")
  expect_identical(numeric_outcome$engine$metrics,
                   c("RMSE", "Rsquared", "MAE"))
  expect_null(numeric_outcome$design$outcome_lv)

  labels <- fit_svm(sa_svm_iris2(), outcome = "Species", C = 1, cv = FALSE,
                    seed = 1)
  expect_identical(labels$design$outcome_type, "two classes")
  expect_identical(labels$engine$metrics, c("Accuracy", "Kappa"))
  expect_identical(labels$design$outcome_lv, c("versicolor", "virginica"))
  expect_identical(labels$design$n_events, 50L)
  expect_equal(labels$design$event_rate, 0.5)

  # A numeric column holding two values is the one case the outcome cannot settle
  # on its own, so the guess is announced and `outcome_lv` overrules it.
  d <- mtcars[c("am", "wt", "hp")]
  expect_message(
    as_numbers <- fit_svm(d, outcome = "am", C = 1, cv = FALSE, seed = 1),
    "fitted as a regression"
  )
  expect_identical(as_numbers$design$outcome_type, "continuous")

  as_classes <- fit_svm(d, outcome = "am", outcome_lv = c("0", "1"), C = 1,
                        cv = FALSE, seed = 1)
  expect_identical(as_classes$design$outcome_type, "two classes")
  expect_identical(as_classes$design$n_events, 13L)

  # More than two classes is the same refusal the other classifications make, in
  # this model's name.
  expect_error(
    fit_svm(iris, outcome = "Species", C = 1, cv = FALSE),
    "holds 3 classes, but a support vector machine models two"
  )
  expect_error(
    fit_svm(iris, outcome = "Species", outcome_lv = c("setosa", "virginica"),
            C = 1, cv = FALSE),
    "would be silently left out"
  )
  expect_error(
    fit_svm(mtcars, outcome = "mpg", predictors = c("wt", "mpg"), C = 1,
            cv = FALSE),
    "predict from the answer"
  )
})


test_that("a factor predictor is one term per level beyond the first", {
  cars <- mtcars[c("mpg", "wt")]
  cars$cyl <- factor(mtcars$cyl)
  fit <- fit_svm(cars, outcome = "mpg", C = 1, cv = FALSE, seed = 1)

  # A kernel measures a distance along each column, so the levels have to be
  # columns. This is where the machine sides with `fit_elastic_net()` against
  # `fit_rf()`, whose trees split the factor itself.
  expect_setequal(fit$terms, c("wt", "cyl6", "cyl8"))
  expect_identical(fit$engine$x_names, c("wt", "cyl6", "cyl8"))
  expect_identical(fit$design$predictors, c("wt", "cyl"))
  expect_identical(fit$design$predictor_lv, list(cyl = c("4", "6", "8")))
  expect_setequal(fit_rf(cars, outcome = "mpg", ntree = 50, cv = FALSE,
                         seed = 1)$terms, c("wt", "cyl"))

  # A character column is the same machine as the factor it stands for, since it is
  # turned into one before the columns are coded.
  cars$cyl <- as.character(cars$cyl)
  chr <- fit_svm(cars, outcome = "mpg", C = 1, cv = FALSE, seed = 1)
  expect_identical(chr$coefficients, fit$coefficients)
})


test_that("a machine of one term is fitted rather than refused", {
  # The contrast with `fit_elastic_net()`, which refuses a single term because a
  # penalty has no budget to divide. A distance along one axis is a distance.
  fit <- fit_svm(mtcars, outcome = "mpg", predictors = "wt", C = 1, cv = FALSE,
                 seed = 1)
  expect_identical(fit$terms, "wt")
  expect_identical(nrow(fit$coefficients), 1L)
  expect_error(
    fit_elastic_net(mtcars, outcome = "mpg", predictors = "wt",
                    penalty = "ridge", lambda = 1, cv = FALSE),
    "the model has 1 term"
  )
})


test_that("the rows and columns that entered the machine are the reported ones", {
  d <- mtcars[c("mpg", "wt", "hp")]
  d$wt[1:3] <- NA
  d$mpg[4] <- NA
  fit <- fit_svm(d, outcome = "mpg", C = 1, cv = FALSE, seed = 1)

  expect_identical(fit$design$n_obs, 32L)
  expect_identical(fit$design$n_used, 28L)
  expect_identical(fit$design$n_dropped, 4L)
  # Dropped before the machine is fitted, so its score is over the rows that were
  # kept and the support vectors are counted against those.
  expect_lte(fit$fit_stats$n_support_vector, 28)
  expect_equal(fit$fit_stats$support_vector_rate,
               fit$fit_stats$n_support_vector / 28)

  flat <- mtcars[c("mpg", "wt", "hp")]
  flat$flat <- 1
  expect_message(
    dropped <- fit_svm(flat, outcome = "mpg", C = 1, cv = FALSE, seed = 1),
    "single value cannot contribute"
  )
  expect_identical(dropped$design$dropped_predictors, "flat")
  expect_setequal(dropped$terms, c("wt", "hp"))
})


test_that("`seed` fixes the width, the shuffles and the folds", {
  # More is fixed here than in the unpenalized models, where the seed fixes only
  # the fold assignment: the kernel width is estimated from sampled pairs of rows
  # and the importance is measured by shuffling, so two machines fitted on the same
  # data report different numbers without one.
  args <- list(data = mtcars, outcome = "mpg",
               predictors = c("wt", "hp", "disp"), C = 1, cv = FALSE)
  once <- do.call(fit_svm, c(args, list(seed = 11)))
  twice <- do.call(fit_svm, c(args, list(seed = 11)))
  other <- do.call(fit_svm, c(args, list(seed = 12)))

  expect_identical(once$coefficients, twice$coefficients)
  expect_identical(once$fit_stats, twice$fit_stats)
  expect_identical(once$parameters, twice$parameters)
  expect_false(isTRUE(all.equal(once$coefficients$estimate,
                                other$coefficients$estimate)))

  set.seed(99)
  before <- .Random.seed
  invisible(do.call(fit_svm, c(args, list(seed = 11))))
  expect_identical(before, .Random.seed)

  # The resampling arguments are recorded as the scheme used them, the way the rest
  # of the family records them.
  cv_args <- list(data = mtcars, outcome = "mpg", predictors = "wt", C = 1,
                  n_fold = 4, n_repeat = 3, seed = 1)
  repeated <- do.call(fit_svm, cv_args)
  expect_identical(repeated$parameters$cv_method, "repeated_kfold")
  expect_identical(nrow(repeated$resampling), 12L)

  kfold <- do.call(fit_svm, c(cv_args, list(cv_method = "kfold")))
  expect_identical(kfold$parameters$n_fold, 4L)
  expect_identical(kfold$parameters$n_repeat, NA_integer_)

  plain <- do.call(fit_svm, c(args, list(seed = 1)))
  expect_identical(plain$parameters$cv_method, NA_character_)
  expect_identical(plain$parameters$n_fold, NA_integer_)
})


test_that("print says what was fitted, at what settings, and how it did", {
  fit <- fit_svm(mtcars, outcome = "mpg", predictors = c("wt", "hp", "disp"),
                 C = c(0.5, 2, 8), cv_method = "kfold", n_fold = 4, seed = 1)

  expect_output(print(fit), "<sa_model> svm")
  expect_output(print(fit), "kernel   : radial, C = ")
  expect_output(print(fit), "sigma = ")
  expect_output(print(fit), "chosen from 3 candidate\\(s\\)")
  expect_output(print(fit), "kfold, 4 fold\\(s\\)")
  expect_output(print(fit), "n_support_vector")
  expect_output(print(fit), "RMSE")
  expect_output(expect_invisible(print(fit)))

  # The block is not called `coefficients`, because these are not coefficients, and
  # there is no interval, no p-value and no `selected` to print beside them. The
  # heading is the forest's, since the two answer the same kind of question.
  output <- utils::capture.output(print(fit))
  expect_true(any(grepl("importance  (permutation)", output, fixed = TRUE)))
  expect_false(any(grepl("coefficients", output, fixed = TRUE)))
  expect_false(any(grepl("p = NA", output, fixed = TRUE)))
  expect_false(any(grepl("selected", output, fixed = TRUE)))
  expect_false(any(grepl("conf_level", output, fixed = TRUE)))

  # A single pair was scored rather than chosen between, so nothing won.
  single <- utils::capture.output(
    print(fit_svm(mtcars, outcome = "mpg", predictors = "wt", C = 1,
                  cv_method = "kfold", n_fold = 4, seed = 1))
  )
  expect_true(any(grepl("kernel   : radial", single, fixed = TRUE)))
  expect_false(any(grepl("chosen from", single, fixed = TRUE)))

  # A machine models the same class in the same direction as the logistic
  # regression and reports no odds, so the header does not claim any.
  clf <- utils::capture.output(
    print(fit_svm(sa_svm_iris2(), outcome = "Species", C = 1, cv = FALSE,
                  seed = 1))
  )
  expect_true(any(grepl("modelling virginica against versicolor", clf,
                        fixed = TRUE)))
  expect_true(any(grepl("accuracy", clf, fixed = TRUE)))
  expect_false(any(grepl("the odds of", clf, fixed = TRUE)))

  # The penalized fit still prints its own table under its own heading, which is
  # what the new one had to be told apart from.
  pen <- utils::capture.output(
    print(fit_elastic_net(mtcars, outcome = "mpg", penalty = "lasso",
                          lambda = 0.5, cv = FALSE))
  )
  expect_true(any(grepl("coefficients", pen, fixed = TRUE)))
  expect_false(any(grepl("importance", pen, fixed = TRUE)))
})


test_that("what the machine leaned on can be scored against a known answer", {
  # A simulated regression plants some coefficients at exactly zero, so the null
  # terms are the ones the importance ought to leave at the bottom. Nothing about a
  # kernel machine guarantees a perfect ordering, so what is pinned is that the
  # planted terms are the more important as a group.
  sim <- simulate_regression(seed = 1)
  svm <- do.call(fit_svm, c(sim$args, list(C = 1, cv = FALSE, seed = 1)))
  scored <- merge(svm$coefficients, sim$truth_term, by = "terms")

  planted <- scored$estimate[scored$beta != 0]
  null <- scored$estimate[scored$beta == 0]
  expect_gt(length(planted), 0L)
  expect_gt(mean(planted), mean(null))

  # A term that carried nothing can be left very slightly better by its own
  # permutation, so a negative importance is a value rather than a gap.
  expect_false(anyNA(scored$estimate))

  # Everything but the engine object survives the trip a JSON export would take.
  skip_if_not_installed("jsonlite")
  portable <- svm[setdiff(names(svm), "fit")]
  round_trip <- jsonlite::fromJSON(
    jsonlite::toJSON(portable, na = "string", digits = NA)
  )
  expect_identical(round_trip$terms, svm$terms)
  expect_equal(round_trip$coefficients$estimate, svm$coefficients$estimate)
  expect_equal(round_trip$fit_stats$rmse, svm$fit_stats$rmse)
  expect_identical(round_trip$engine$x_names, svm$engine$x_names)
})
