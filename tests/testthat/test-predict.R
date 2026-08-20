# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# What is pinned here is that the rows to predict can be handed over as they came.
# The method exists because the engine object cannot take them: `caret` prepares
# `newdata` for `glmnet` by keeping the columns whose names it recognises, in the
# order they happen to be in, so a factor predictor is dropped whole and a set of
# numeric ones can be matched to the wrong coefficients without any error being
# raised. The two tests that would have failed before this method are therefore
# the ones about a factor predictor and about column order, and the second of them
# is the reason the order is fixed by name rather than trusted.
#
# The arithmetic is checked against a design matrix built here by hand: the point
# of the method is that it builds the same one, not that `glmnet` multiplies
# correctly.
#
# One prediction per row of `newdata` is the third promise, which is why the
# missing values are tested through a predictor the penalty dropped as well as
# through one it kept. The sparse product never reads a cell it multiplies by
# zero, so without a rule of its own the answer would depend on which column the
# hole fell in.

sa_pred_cars <- function() {
  d <- mtcars[c("mpg", "wt", "hp", "disp", "cyl")]
  d$cyl <- factor(d$cyl, labels = c("four", "six", "eight"))
  rownames(d) <- NULL
  d
}

sa_pred_iris2 <- function() {
  d <- iris[iris$Species != "setosa",
            c("Species", "Petal.Length", "Sepal.Width", "Sepal.Length")]
  d$size <- cut(d$Sepal.Length, 3, labels = c("small", "mid", "large"))
  d$Sepal.Length <- NULL
  rownames(d) <- NULL
  d
}

sa_pred_train <- function(d) d[seq(1L, nrow(d), by = 2L), , drop = FALSE]
sa_pred_test <- function(d) d[seq(2L, nrow(d), by = 2L), , drop = FALSE]

sa_pred_enet <- function(d = sa_pred_cars(), lambda = 0.5) {
  fit_elastic_net(sa_pred_train(d), outcome = "mpg",
                  predictors = c("wt", "hp", "disp", "cyl"),
                  penalty = "lasso", lambda = lambda, cv = FALSE)
}


test_that("a penalized fit predicts the design matrix it was fitted on", {
  d <- sa_pred_cars()
  fit <- sa_pred_enet(d)
  test <- sa_pred_test(d)

  # The reference is the matrix `caret` was given, coded here rather than by the
  # method under test, and handed to the engine object directly.
  x <- stats::model.matrix(~ ., data = test[fit$design$predictors])
  x <- x[, fit$fit$finalModel$xNames, drop = FALSE]
  reference <- as.numeric(stats::predict(fit$fit, newdata = x))

  got <- stats::predict(fit, newdata = test)
  expect_equal(got, reference)
  expect_length(got, nrow(test))
  expect_type(got, "double")
  expect_null(names(got))

  # The terms are the factor's levels, which is what the engine object holds and
  # the frame does not.
  expect_identical(fit$fit$finalModel$xNames,
                   c("wt", "hp", "disp", "cylsix", "cyleight"))
})


test_that("the predictors are read by name, not by position", {
  d <- sa_pred_cars()
  fit <- sa_pred_enet(d)
  test <- sa_pred_test(d)
  got <- stats::predict(fit, newdata = test)

  # Reversing the columns is what `caret` cannot survive on its own: it keeps the
  # ones it recognises in the order it was given them, so the coefficients would
  # meet the wrong columns and no error would be raised.
  expect_identical(stats::predict(fit, newdata = test[rev(names(test))]), got)

  # A column the model never saw is ignored rather than refused, which is what
  # lets the test half be handed over with its outcome still in it.
  extra <- test
  extra$note <- "kept"
  extra$gear <- mtcars$gear[seq(2L, nrow(mtcars), by = 2L)]
  expect_identical(stats::predict(fit, newdata = extra), got)

  # A matrix is accepted the way `data` is, as long as it can carry the columns.
  num <- fit_elastic_net(sa_pred_train(d), outcome = "mpg",
                         predictors = c("wt", "hp", "disp"),
                         penalty = "lasso", lambda = 0.5, cv = FALSE)
  as_matrix <- as.matrix(sa_pred_test(d)[c("wt", "hp", "disp")])
  expect_equal(stats::predict(num, newdata = as_matrix),
               stats::predict(num, newdata = sa_pred_test(d)))
})


test_that("what cannot be predicted from is an error naming it", {
  d <- sa_pred_cars()
  fit <- sa_pred_enet(d)
  test <- sa_pred_test(d)

  expect_error(stats::predict(fit, newdata = test[c("mpg", "wt", "cyl")]),
               "missing predictor column\\(s\\).*hp, disp")
  expect_error(stats::predict(fit, newdata = test[0L, ]), "zero rows")
  expect_error(stats::predict(fit, newdata = "test"),
               "must be a data.frame or a matrix")

  # A level the fit never saw has no coefficient to be predicted with, which is a
  # different fact from a level `newdata` happens not to hold.
  unseen <- test
  levels(unseen$cyl) <- c("four", "six", "twelve")
  expect_error(stats::predict(fit, newdata = unseen),
               "level\\(s\\) of `cyl` the model was not fitted on.*twelve")

  # A column whose coding changed between fitting and predicting would be coded
  # into different terms, so it is refused rather than guessed at.
  swapped <- test
  swapped$wt <- as.character(swapped$wt)
  expect_error(stats::predict(fit, newdata = swapped),
               "`wt` was a numeric predictor")

  numeric_cyl <- test
  numeric_cyl$cyl <- as.numeric(numeric_cyl$cyl)
  expect_error(stats::predict(fit, newdata = numeric_cyl),
               "level\\(s\\) of `cyl` the model was not fitted on")
})


test_that("a level `newdata` does not hold keeps its term", {
  d <- sa_pred_cars()
  fit <- sa_pred_enet(d)
  test <- sa_pred_test(d)
  some <- test[test$cyl != "six", , drop = FALSE]

  expect_false("six" %in% as.character(some$cyl))
  got <- stats::predict(fit, newdata = some)
  expect_length(got, nrow(some))

  # The rows are the same rows, so the predictions are the same numbers: the
  # coding did not shift when a level went missing.
  expect_equal(got, stats::predict(fit, newdata = test)[test$cyl != "six"])
})


test_that("a row the model cannot read is predicted as NA", {
  d <- sa_pred_cars()
  fit <- sa_pred_enet(d, lambda = 2)
  test <- sa_pred_test(d)
  full <- stats::predict(fit, newdata = test)

  holed <- test
  holed$wt[1L] <- NA
  holed$cyl[2L] <- NA
  got <- stats::predict(fit, newdata = holed)

  expect_length(got, nrow(test))
  expect_true(all(is.na(got[1:2])))
  # The rows that could be read are answered as if the others were not there.
  expect_equal(got[-(1:2)], full[-(1:2)])

  # Including a hole in a predictor whose every term the penalty dropped. Its
  # coefficients are exactly zero, so the sparse product would never have read the
  # missing cell and the row would have been predicted as if it were complete.
  expect_setequal(fit$coefficients$terms[!fit$coefficients$selected],
                  c("cylsix", "cyleight"))
  zeroed <- test
  zeroed$cyl[3L] <- NA
  expect_true(is.na(stats::predict(fit, newdata = zeroed)[3L]))

  none <- test
  none$wt <- NA_real_
  expect_error(stats::predict(fit, newdata = none), "no row of `newdata`")
})


test_that("a classification predicts in the direction its coefficients read", {
  d <- sa_pred_iris2()
  lv <- c("versicolor", "virginica")
  fit <- fit_elastic_net(sa_pred_train(d), outcome = "Species",
                         outcome_lv = lv, penalty = "lasso", lambda = 0.01,
                         cv = FALSE)
  test <- sa_pred_test(d)

  raw <- stats::predict(fit, newdata = test)
  expect_s3_class(raw, "factor")
  expect_identical(levels(raw), lv)
  expect_length(raw, nrow(test))

  prob <- stats::predict(fit, newdata = test, type = "prob")
  expect_named(prob, lv)
  expect_identical(nrow(prob), nrow(test))
  expect_equal(unname(rowSums(prob)), rep(1, nrow(test)))

  # `"response"` is the probability of the second level, the class every
  # coefficient and odds ratio in the result describes.
  resp <- stats::predict(fit, newdata = test, type = "response")
  expect_equal(resp, prob[[lv[2L]]])
  expect_null(names(resp))

  # The label is the more probable class, so the two answers agree.
  expect_identical(as.character(raw), ifelse(resp > 0.5, lv[2L], lv[1L]))

  # A row that cannot be read is missing in every shape the answer takes.
  holed <- test
  holed$size[1L] <- NA
  expect_true(is.na(stats::predict(fit, newdata = holed)[1L]))
  expect_true(all(is.na(stats::predict(fit, newdata = holed,
                                       type = "prob")[1L, ])))
  expect_identical(levels(stats::predict(fit, newdata = holed)), lv)
})


test_that("the unpenalized fits predict what their engines predict", {
  d <- sa_pred_cars()
  train <- sa_pred_train(d)
  test <- sa_pred_test(d)

  lin <- fit_linear_regression(train, outcome = "mpg",
                               predictors = c("wt", "hp", "disp", "cyl"),
                               cv = FALSE)
  reference <- stats::lm(mpg ~ wt + hp + disp + cyl, data = train)
  expect_equal(stats::predict(lin, newdata = test),
               unname(stats::predict(reference, newdata = test)))

  d2 <- sa_pred_iris2()
  lv <- c("versicolor", "virginica")
  ref_train <- sa_pred_train(d2)
  clf <- fit_logistic_regression(ref_train, outcome = "Species",
                                 predictors = c("Petal.Length", "size"),
                                 outcome_lv = lv, cv = FALSE)
  ref_train$Species <- factor(ref_train$Species, levels = lv)
  ref_glm <- stats::glm(Species ~ Petal.Length + size, data = ref_train,
                        family = stats::binomial())
  expect_equal(
    stats::predict(clf, newdata = sa_pred_test(d2), type = "response"),
    unname(stats::predict(ref_glm, newdata = sa_pred_test(d2),
                          type = "response"))
  )
  expect_identical(levels(stats::predict(clf, newdata = sa_pred_test(d2))), lv)
})


test_that("a forest predicts the frame it was handed", {
  # The other branch of the method: a tree splits a factor on its levels, so the
  # engine took the predictor frame and not a design matrix. What still has to
  # hold is everything the method promises regardless of engine — the columns read
  # by name, the levels put back as the fit coded them, and one prediction per row.
  d <- sa_pred_cars()
  train <- sa_pred_train(d)
  test <- sa_pred_test(d)
  fit <- fit_rf(train, outcome = "mpg",
                predictors = c("wt", "hp", "disp", "cyl"), ntree = 100,
                cv = FALSE, seed = 1)

  expect_identical(fit$fit$finalModel$xNames, c("wt", "hp", "disp", "cyl"))
  got <- stats::predict(fit, newdata = test)
  expect_length(got, nrow(test))
  expect_type(got, "double")
  expect_null(names(got))
  expect_identical(stats::predict(fit, newdata = test[rev(names(test))]), got)

  holed <- test
  holed$cyl[1L] <- NA
  scattered <- stats::predict(fit, newdata = holed)
  expect_true(is.na(scattered[1L]))
  expect_equal(scattered[-1L], got[-1L])

  expect_error(stats::predict(fit, newdata = test[c("mpg", "wt")]),
               "missing predictor column\\(s\\).*hp, disp, cyl")
  unseen <- test
  levels(unseen$cyl) <- c("four", "six", "twelve")
  expect_error(stats::predict(fit, newdata = unseen),
               "level\\(s\\) of `cyl` the model was not fitted on.*twelve")

  # A two-class forest answers in all three shapes, and `"response"` is the
  # probability of `outcome_lv[2]` as it is for the models with coefficients.
  d2 <- sa_pred_iris2()
  lv <- c("versicolor", "virginica")
  clf <- fit_rf(sa_pred_train(d2), outcome = "Species", outcome_lv = lv,
                ntree = 100, cv = FALSE, seed = 1)
  test2 <- sa_pred_test(d2)

  raw <- stats::predict(clf, newdata = test2)
  expect_s3_class(raw, "factor")
  expect_identical(levels(raw), lv)
  prob <- stats::predict(clf, newdata = test2, type = "prob")
  expect_named(prob, lv)
  resp <- stats::predict(clf, newdata = test2, type = "response")
  expect_equal(resp, prob[[lv[2L]]])
  expect_null(names(resp))
})


test_that("a machine predicts the design matrix it was fitted on", {
  # The same branch the penalized fit takes, reached by a different engine: a
  # kernel measures a distance along each column, so `kernlab` was handed a matrix
  # too. What names the branch is `engine$x_names` rather than the engine, since it
  # is the result that has to remember the columns; `caret` records none for a fit
  # it was given a matrix.
  d <- sa_pred_cars()
  train <- sa_pred_train(d)
  test <- sa_pred_test(d)
  fit <- fit_svm(train, outcome = "mpg",
                 predictors = c("wt", "hp", "disp", "cyl"), C = 1, cv = FALSE,
                 seed = 1)

  expect_identical(fit$engine$x_names,
                   c("wt", "hp", "disp", "cylsix", "cyleight"))
  expect_null(fit$fit$coefnames)

  x <- stats::model.matrix(~ ., data = test[fit$design$predictors])
  x <- x[, fit$engine$x_names, drop = FALSE]
  reference <- as.numeric(stats::predict(fit$fit, newdata = x))

  got <- stats::predict(fit, newdata = test)
  expect_equal(got, reference)
  expect_length(got, nrow(test))
  expect_type(got, "double")
  expect_null(names(got))
  expect_identical(stats::predict(fit, newdata = test[rev(names(test))]), got)

  holed <- test
  holed$cyl[1L] <- NA
  scattered <- stats::predict(fit, newdata = holed)
  expect_true(is.na(scattered[1L]))
  expect_equal(scattered[-1L], got[-1L])

  expect_error(stats::predict(fit, newdata = test[c("mpg", "wt")]),
               "missing predictor column\\(s\\).*hp, disp, cyl")
  unseen <- test
  levels(unseen$cyl) <- c("four", "six", "twelve")
  expect_error(stats::predict(fit, newdata = unseen),
               "level\\(s\\) of `cyl` the model was not fitted on.*twelve")

  # A two-class machine answers in all three shapes, and `"response"` is the
  # probability of `outcome_lv[2]` as it is for the models with coefficients. The
  # probability exists at all because `prob.model` reached the engine.
  d2 <- sa_pred_iris2()
  lv <- c("versicolor", "virginica")
  clf <- fit_svm(sa_pred_train(d2), outcome = "Species", outcome_lv = lv, C = 1,
                 cv = FALSE, seed = 1)
  test2 <- sa_pred_test(d2)

  raw <- stats::predict(clf, newdata = test2)
  expect_s3_class(raw, "factor")
  expect_identical(levels(raw), lv)
  prob <- stats::predict(clf, newdata = test2, type = "prob")
  expect_named(prob, lv)
  expect_equal(unname(rowSums(prob)), rep(1, nrow(test2)))
  resp <- stats::predict(clf, newdata = test2, type = "response")
  expect_equal(resp, prob[[lv[2L]]])
  expect_null(names(resp))
})


test_that("no newdata predicts the rows the model was fitted on", {
  d <- sa_pred_cars()
  d$wt[c(1L, 3L)] <- NA
  fit <- sa_pred_enet(d)
  train <- sa_pred_train(d)

  expect_identical(fit$design$n_dropped, nrow(train) - fit$design$n_used)
  got <- stats::predict(fit)
  expect_length(got, fit$design$n_used)
  expect_false(anyNA(got))
})


test_that("`design` records the levels the fit coded against", {
  d <- sa_pred_cars()
  fit <- sa_pred_enet(d)

  expect_identical(fit$design$predictor_lv,
                   list(cyl = c("four", "six", "eight")))

  # A model with nothing to record does not carry the slot, the way a regression
  # does not carry `outcome_lv`.
  numeric_only <- fit_elastic_net(sa_pred_train(d), outcome = "mpg",
                                  predictors = c("wt", "hp", "disp"),
                                  penalty = "lasso", lambda = 0.5, cv = FALSE)
  expect_null(numeric_only$design$predictor_lv)
  expect_false("predictor_lv" %in% names(numeric_only$design))

  # A character predictor is recorded as the factor it was turned into, and an
  # unused level is not recorded, since the fit dropped it before coding.
  chr <- sa_pred_train(d)
  chr$cyl <- as.character(chr$cyl)
  chr <- chr[chr$cyl != "six", , drop = FALSE]
  from_chr <- fit_elastic_net(chr, outcome = "mpg",
                              predictors = c("wt", "hp", "cyl"),
                              penalty = "lasso", lambda = 0.5, cv = FALSE)
  expect_identical(from_chr$design$predictor_lv, list(cyl = c("eight", "four")))
  expect_error(stats::predict(from_chr, newdata = sa_pred_test(d)),
               "level\\(s\\) of `cyl`.*six")
})
