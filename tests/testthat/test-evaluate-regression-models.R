# Scoring fitted regressions on held-out rows. The arithmetic is checked against
# the same arithmetic written out a second way rather than against another
# package, since every quantity here is two lines long and the point of writing
# it out was that the Python port would have the formula rather than a call.

test_that("the result carries the contract slots and nothing else", {
  res <- sa_perf_reg_fixture()
  expect_s3_class(res, "sa_performance")
  expect_s3_class(res, "sa_result")
  expect_identical(names(res),
                   c("analysis", "models", "design", "parameters",
                     "predictions", "metrics", "comparisons", "metadata"))
  expect_identical(res$analysis, "regression_performance")
  # A regression has no ROC curve, so the slot is absent rather than empty.
  expect_null(res$curves)
})

test_that("every table is keyed on the models in one order", {
  res <- sa_perf_reg_fixture()
  expect_identical(res$models, c("baseline", "wt_only", "qsec_only"))
  expect_identical(res$metrics$model, res$models)
  expect_identical(unique(res$predictions$model), res$models)
  # The baseline is what the comparisons are against, so it has no row of its
  # own among them.
  expect_identical(res$comparisons$model, res$models[-1])
})

test_that("the metric table holds every contract column", {
  res <- sa_perf_reg_fixture()
  expect_true(all(sa_regression_metric_columns() %in% names(res$metrics)))
  expect_true(all(sa_regression_comparison_columns() %in%
                    names(res$comparisons)))
  expect_true(all(sa_prediction_table_columns() %in% names(res$predictions)))
})

test_that("the metrics are the arithmetic they claim to be", {
  res <- sa_perf_reg_fixture()
  drawn <- res$predictions[res$predictions$model == "baseline", ]
  observed <- drawn$observed
  predicted <- drawn$predicted
  residual <- predicted - observed
  row <- res$metrics[1, ]

  expect_identical(row$n_used, length(observed))
  expect_equal(row$cor, stats::cor(observed, predicted))
  expect_equal(row$rmse, sqrt(mean(residual^2)))
  expect_equal(row$mae, mean(abs(residual)))
  expect_equal(row$bias, mean(residual))
  expect_equal(row$r_squared,
               1 - sum(residual^2) / sum((observed - mean(observed))^2))
})

test_that("r_squared is the held-out one rather than the squared correlation", {
  res <- sa_perf_reg_fixture()
  # The two answer different questions and agree only for predictions that need
  # no calibration, which is exactly what the calibration columns report. A test
  # that let them be equal would not notice a `cor^2` slipping into the column.
  expect_false(isTRUE(all.equal(res$metrics$r_squared, res$metrics$cor^2)))
  expect_true(all(res$metrics$r_squared <= res$metrics$cor^2 + 1e-8))
})

test_that("the calibration line is lm(predicted ~ observed)", {
  res <- sa_perf_reg_fixture()
  for (i in seq_along(res$models)) {
    drawn <- res$predictions[res$predictions$model == res$models[i], ]
    fitted <- stats::lm(drawn$predicted ~ drawn$observed)
    expect_equal(unname(stats::coef(fitted)),
                 c(res$metrics$calib_intercept[i], res$metrics$calib_slope[i]),
                 info = res$models[i])
  }
})

test_that("a comparison is the new model less the baseline", {
  res <- sa_perf_reg_fixture()
  at <- seq_along(res$models)[-1]
  expect_equal(res$comparisons$delta_cor,
               res$metrics$cor[at] - res$metrics$cor[1])
  expect_equal(res$comparisons$delta_r_squared,
               res$metrics$r_squared[at] - res$metrics$r_squared[1])
  expect_equal(res$comparisons$delta_rmse,
               res$metrics$rmse[at] - res$metrics$rmse[1])
  expect_equal(res$comparisons$delta_mae,
               res$metrics$mae[at] - res$metrics$mae[1])
  # A worse model raises the error and lowers the correlation, which is what
  # fixes the sign convention the print method and the roxygen describe.
  expect_true(all(res$comparisons$delta_rmse > 0))
  expect_true(all(res$comparisons$delta_cor < 0))
})

test_that("no test is reported beside a difference of held-out errors", {
  res <- sa_perf_reg_fixture()
  expect_false(any(grepl("pval", names(res$comparisons))))
  expect_false(any(grepl("conf", names(res$comparisons))))
  # And the evaluation itself records no confidence level, since it reports no
  # interval to have one.
  expect_identical(res$parameters, list())
})

test_that("one model on its own is scored without a comparison", {
  parts <- sa_perf_reg_parts()
  res <- evaluate_regression_models(parts$models$full, newdata = parts$test)
  expect_identical(res$models, "baseline")
  expect_identical(nrow(res$metrics), 1L)
  # An empty table would read as a result that lost its values rather than as
  # one for which the question does not arise.
  expect_null(res$comparisons)
  expect_false("comparisons" %in% names(res))
  # An empty list is the same statement as NULL and reads better out of a
  # lapply() that found nothing.
  expect_identical(evaluate_regression_models(parts$models$full,
                                              list(), newdata = parts$test),
                   res, ignore_attr = "metadata")
})

test_that("the baseline is labelled as the caller asks", {
  parts <- sa_perf_reg_parts()
  res <- evaluate_regression_models(parts$models$full,
                                    list(reduced = parts$models$wt),
                                    newdata = parts$test,
                                    baseline_label = "All features")
  expect_identical(res$models, c("All features", "reduced"))
  expect_identical(res$design$baseline, "All features")
  expect_identical(res$metrics$model[1], "All features")
})

test_that("the observed outcome is found three ways and read the same", {
  parts <- sa_perf_reg_parts()
  by_default <- evaluate_regression_models(parts$models$full,
                                           newdata = parts$test)
  by_name <- evaluate_regression_models(parts$models$full,
                                        newdata = parts$test, answer = "mpg")
  by_vector <- evaluate_regression_models(parts$models$full,
                                          newdata = parts$test,
                                          answer = parts$test$mpg)
  expect_equal(by_default$metrics, by_name$metrics)
  expect_equal(by_default$metrics, by_vector$metrics)
  # Only the label differs, since a vector no longer remembers where it came
  # from and the result records what the analysis was made on.
  expect_identical(by_name$design$outcome, "mpg")
  expect_identical(by_vector$design$outcome, "<vector>")
})

test_that("every model is scored on the rows all of them could predict", {
  parts <- sa_perf_reg_parts()
  holed <- parts$test
  # `disp` is read by the baseline alone and `wt` by two of the three, so the
  # union and the intersection differ and a per-model score would put two
  # numbers from two samples in one table.
  holed$disp[1:2] <- NA
  holed$wt[3] <- NA

  expect_message(
    res <- evaluate_regression_models(parts$models$full,
                                      list(wt_only = parts$models$wt),
                                      newdata = holed),
    "left out"
  )
  expect_identical(res$design$n_obs, nrow(holed))
  expect_identical(res$design$n_used, nrow(holed) - 3L)
  expect_identical(res$design$n_dropped, 3L)
  expect_identical(res$metrics$n_used, rep(nrow(holed) - 3L, 2L))
  # The same rows, named by their position in `newdata` rather than renumbered.
  rows <- split(res$predictions$row, res$predictions$model)
  expect_identical(rows[[1]], rows[[2]])
  expect_false(any(c(1L, 2L, 3L) %in% rows[[1]]))
})

test_that("a row with no observed outcome is dropped and counted too", {
  parts <- sa_perf_reg_parts()
  holed <- parts$test
  holed$mpg[1:2] <- NA
  expect_message(
    res <- evaluate_regression_models(parts$models$full, newdata = holed),
    "no observed outcome"
  )
  expect_identical(res$design$n_dropped, 2L)
})

test_that("a classification is refused by name and redirected", {
  cls <- sa_perf_cls_parts()
  reg <- sa_perf_reg_parts()
  expect_error(
    evaluate_regression_models(cls$models$both, newdata = cls$test),
    "evaluate_classification_models"
  )
  # And when it is one of several, the message says which one.
  expect_error(
    evaluate_regression_models(reg$models$full,
                               list(wrong = cls$models$both),
                               newdata = reg$test),
    "wrong \\(two classes\\)"
  )
})

test_that("models of different outcomes are refused rather than tabulated", {
  parts <- sa_perf_reg_parts()
  other <- fit_linear_regression(parts$train, outcome = "qsec",
                                 predictors = c("wt", "hp"), cv = FALSE)
  expect_error(
    evaluate_regression_models(parts$models$full, list(other = other),
                               newdata = parts$test),
    "same outcome"
  )
})

test_that("the models have to arrive named, once each, and as models", {
  parts <- sa_perf_reg_parts()
  expect_error(
    evaluate_regression_models(parts$models$full, list(parts$models$wt),
                               newdata = parts$test),
    "must be named"
  )
  expect_error(
    evaluate_regression_models(parts$models$full,
                               list(a = parts$models$wt, a = parts$models$qsec),
                               newdata = parts$test),
    "duplicated names"
  )
  expect_error(
    evaluate_regression_models(parts$models$full,
                               list(baseline = parts$models$wt),
                               newdata = parts$test),
    "what the baseline is called"
  )
  expect_error(
    evaluate_regression_models(parts$models$full, list(a = 1),
                               newdata = parts$test),
    "Not a model: a"
  )
  expect_error(
    evaluate_regression_models(parts$models$full, parts$models$wt,
                               newdata = parts$test),
    "named list"
  )
  expect_error(evaluate_regression_models(mtcars, newdata = parts$test),
               "must be a fitted model")
})

test_that("the rows to score have to be rows", {
  parts <- sa_perf_reg_parts()
  expect_error(evaluate_regression_models(parts$models$full,
                                          newdata = parts$test[0, ]),
               "zero rows")
  expect_error(evaluate_regression_models(parts$models$full, newdata = 1:5),
               "data.frame or a matrix")
  expect_error(
    evaluate_regression_models(parts$models$full,
                               newdata = parts$test[c("wt", "hp", "disp")]),
    "no such column"
  )
})

test_that("a label where a number belongs is refused with the way out", {
  parts <- sa_perf_reg_parts()
  expect_error(
    evaluate_regression_models(parts$models$full, newdata = parts$test,
                               answer = rep("a", nrow(parts$test))),
    "evaluate_classification_models"
  )
})

test_that("an outcome that does not vary is reported once, not per model", {
  parts <- sa_perf_reg_parts()
  flat <- parts$test
  flat$mpg <- 20
  expect_warning(
    res <- evaluate_regression_models(parts$models$full,
                                      list(wt_only = parts$models$wt),
                                      newdata = flat),
    "single value"
  )
  # Everything measured against the spread of the outcome is unavailable, and
  # everything measured against the outcome itself is still reported.
  expect_true(all(is.na(res$metrics$cor)))
  expect_true(all(is.na(res$metrics$r_squared)))
  expect_true(all(is.na(res$metrics$calib_slope)))
  expect_false(any(is.na(res$metrics$rmse)))
  expect_false(any(is.na(res$metrics$bias)))
})

test_that("the result survives a JSON round trip", {
  skip_if_not_installed("jsonlite")
  res <- sa_perf_reg_fixture()
  # No engine object anywhere: the calibration line is two numbers rather than
  # the `lm` that produced them, which is what leaves every slot a scalar, a
  # character vector, a named list or a data.frame.
  expect_false(any(vapply(res, function(s) is.object(s) && !is.data.frame(s),
                          logical(1))))

  back <- jsonlite::fromJSON(
    jsonlite::toJSON(unclass(res), digits = NA, null = "null")
  )
  expect_identical(back$analysis, res$analysis)
  expect_identical(back$models, res$models)
  expect_identical(back$design$outcome, res$design$outcome)
  expect_identical(back$metrics$model, res$metrics$model)
  expect_equal(back$metrics$rmse, res$metrics$rmse)
  expect_equal(back$metrics$calib_slope, res$metrics$calib_slope)
  expect_equal(back$comparisons$delta_rmse, res$comparisons$delta_rmse)
  expect_equal(back$predictions$predicted, res$predictions$predicted)
})

test_that("printing summarises without dumping the tables", {
  res <- sa_perf_reg_fixture()
  expect_output(print(res), "<sa_performance> regression_performance")
  expect_output(print(res), "outcome  : mpg")
  expect_output(print(res), "baseline = baseline")
  expect_output(print(res), "r_squared")
  expect_output(print(res), "comparisons")
  expect_output(print(res), "delta_rmse")
  # A regression has no threshold to have measured anything at, and no
  # confidence level either, so neither line is printed against nothing.
  printed <- utils::capture.output(print(res))
  expect_false(any(grepl("threshold", printed)))
  expect_false(any(grepl("conf_level", printed)))
  expect_identical(withr::with_output_sink(tempfile(), print(res)), res)
})

test_that("printing counts the models it did not show", {
  res <- sa_perf_reg_fixture()
  expect_output(print(res, n = 1), "more model\\(s\\) in \\$metrics")
  expect_output(print(res, n = 1), "more model\\(s\\) in \\$comparisons")
})
