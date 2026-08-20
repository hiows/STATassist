# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# Scoring fitted classifications on held-out rows. The kernels are checked in
# test-kernel-performance.R; what is checked here is that the right kernel is
# handed the right vectors in the right direction, that the three comparisons
# are paired against the baseline, and that the class order the fits fixed is
# not quietly reversed.

test_that("the result carries the contract slots, curves included", {
  res <- sa_perf_cls_fixture()
  expect_s3_class(res, "sa_performance")
  expect_identical(names(res),
                   c("analysis", "models", "design", "parameters",
                     "predictions", "metrics", "comparisons", "curves",
                     "metadata"))
  expect_identical(res$analysis, "classification_performance")
  expect_true(all(sa_classification_metric_columns() %in% names(res$metrics)))
  expect_true(all(sa_classification_comparison_columns() %in%
                    names(res$comparisons)))
  expect_true(all(sa_roc_curve_columns() %in% names(res$curves)))
})

test_that("the design records the class the numbers are about", {
  res <- sa_perf_cls_fixture()
  expect_identical(res$design$outcome, "Species")
  expect_identical(res$design$outcome_type, "two classes")
  expect_identical(res$design$outcome_lv, c("versicolor", "virginica"))
  expect_identical(res$design$n_events,
                   sum(sa_perf_cls_parts()$test$Species == "virginica"))
  expect_identical(res$parameters, list(threshold = 0.5, conf_level = 0.95))
})

test_that("the stored outcome is the event indicator, not the label", {
  res <- sa_perf_cls_fixture()
  drawn <- res$predictions[res$predictions$model == "baseline", ]
  expect_true(all(drawn$observed %in% c(0, 1)))
  # 1 is `outcome_lv[2]`, which is the class `type = "response"` reports the
  # probability of, so the two line up row for row.
  parts <- sa_perf_cls_parts()
  expect_equal(drawn$observed,
               as.numeric(parts$test$Species == "virginica"))
  expect_equal(drawn$predicted,
               as.numeric(stats::predict(parts$models$both,
                                         newdata = parts$test,
                                         type = "response")))
})

test_that("the metrics are the kernels applied to those two vectors", {
  res <- sa_perf_cls_fixture()
  for (i in seq_along(res$models)) {
    drawn <- res$predictions[res$predictions$model == res$models[i], ]
    area <- sa_auc_delong(drawn$observed, drawn$predicted)
    expect_equal(res$metrics$auc[i], area[["auc"]], info = res$models[i])
    expect_equal(res$metrics$auc[i], sa_auc(drawn$observed, drawn$predicted),
                 info = res$models[i])
    expect_equal(res$metrics$brier[i],
                 sa_brier(drawn$observed, drawn$predicted),
                 info = res$models[i])
    # The interval is the Wald one on the DeLong standard error, which is the
    # one that matches the statistic beside it.
    expect_equal(res$metrics$auc_lower_conf[i],
                 area[["auc"]] - stats::qnorm(0.975) * area[["se"]])
    expect_equal(res$metrics$auc_upper_conf[i],
                 area[["auc"]] + stats::qnorm(0.975) * area[["se"]])
  }
})

test_that("the last three metrics are measured at the stated threshold", {
  parts <- sa_perf_cls_parts()
  low <- evaluate_classification_models(parts$models$both,
                                        newdata = parts$test, threshold = 0.2)
  high <- evaluate_classification_models(parts$models$both,
                                         newdata = parts$test, threshold = 0.8)
  expect_identical(low$parameters$threshold, 0.2)
  # A lower cut calls more rows events, so it can only raise sensitivity and
  # lower specificity.
  expect_gte(low$metrics$sensitivity, high$metrics$sensitivity)
  expect_lte(low$metrics$specificity, high$metrics$specificity)
  # And nothing that does not need a threshold may move with it.
  expect_equal(low$metrics$auc, high$metrics$auc)
  expect_equal(low$metrics$brier, high$metrics$brier)
})

test_that("the confidence level widens every interval it is asked to", {
  parts <- sa_perf_cls_parts()
  narrow <- evaluate_classification_models(parts$models$both,
                                           list(w = parts$models$width),
                                           newdata = parts$test,
                                           conf_level = 0.8)
  wide <- evaluate_classification_models(parts$models$both,
                                         list(w = parts$models$width),
                                         newdata = parts$test,
                                         conf_level = 0.99)
  expect_true(all(wide$metrics$auc_lower_conf < narrow$metrics$auc_lower_conf))
  expect_true(all(wide$metrics$auc_upper_conf > narrow$metrics$auc_upper_conf))
  expect_lt(wide$comparisons$idi_lower_conf, narrow$comparisons$idi_lower_conf)
  expect_gt(wide$comparisons$nri_upper_conf, narrow$comparisons$nri_upper_conf)
  # The estimates and the p-values are not the caller's to move.
  expect_equal(wide$metrics$auc, narrow$metrics$auc)
  expect_equal(wide$comparisons$idi_pval, narrow$comparisons$idi_pval)
})

test_that("the three comparisons are paired against the baseline", {
  res <- sa_perf_cls_fixture()
  drawn <- split(res$predictions$predicted, res$predictions$model)
  observed <- res$predictions$observed[res$predictions$model == "baseline"]

  for (i in seq_len(nrow(res$comparisons))) {
    nm <- res$comparisons$model[i]
    new <- drawn[[nm]]
    old <- drawn[["baseline"]]
    expect_equal(res$comparisons$delta_auc[i],
                 sa_delong_test(observed, new, old)[["delta"]], info = nm)
    expect_equal(res$comparisons$delta_auc[i],
                 res$metrics$auc[res$metrics$model == nm] -
                   res$metrics$auc[1], info = nm)
    expect_equal(res$comparisons$idi[i],
                 sa_idi(observed, old, new)[["idi"]], info = nm)
    expect_equal(res$comparisons$nri[i],
                 sa_nri(observed, old, new)[["nri"]], info = nm)
  }
  # The NRI is the sum of the two components reported beside it.
  expect_equal(res$comparisons$nri,
               res$comparisons$nri_event + res$comparisons$nri_nonevent)
})

test_that("a model that is genuinely worse is called worse three ways", {
  res <- sa_perf_cls_fixture()
  worse <- res$comparisons[res$comparisons$model == "width_only", ]
  expect_lt(worse$delta_auc, 0)
  expect_lt(worse$idi, 0)
  expect_lt(worse$nri, 0)
  expect_lt(worse$delta_auc_pval, 0.05)
  expect_lt(worse$idi_pval, 0.05)
})

test_that("the curves run corner to corner and are monotone", {
  res <- sa_perf_cls_fixture()
  for (nm in res$models) {
    points <- res$curves[res$curves$model == nm, ]
    expect_identical(points$threshold[1], Inf)
    expect_equal(points$sensitivity[1], 0)
    expect_equal(points$specificity[1], 1)
    expect_equal(points$sensitivity[nrow(points)], 1)
    expect_equal(points$specificity[nrow(points)], 0)
    expect_false(is.unsorted(points$sensitivity))
    expect_false(is.unsorted(rev(points$specificity)))
    # The area under the drawn curve has to be the AUC reported beside it, or
    # the picture and the table describe different curves.
    fpr <- 1 - points$specificity
    tpr <- points$sensitivity
    trapezoid <- sum(diff(fpr) *
                       (utils::head(tpr, -1) + utils::tail(tpr, -1)) / 2)
    expect_equal(trapezoid, res$metrics$auc[res$metrics$model == nm],
                 info = nm)
  }
})

test_that("one model on its own is scored with a curve and no comparison", {
  parts <- sa_perf_cls_parts()
  res <- evaluate_classification_models(parts$models$both,
                                        newdata = parts$test)
  expect_null(res$comparisons)
  expect_false("comparisons" %in% names(res))
  # The curve is a property of the model rather than of a comparison, so it
  # stays.
  expect_s3_class(res$curves, "data.frame")
  expect_identical(unique(res$curves$model), "baseline")
})

test_that("the class order is the fits', and disagreeing with it is an error", {
  parts <- sa_perf_cls_parts()
  expect_error(
    evaluate_classification_models(parts$models$both, newdata = parts$test,
                                   outcome_lv = c("virginica", "versicolor")),
    "cannot be re-pointed"
  )
  expect_error(
    evaluate_classification_models(parts$models$both, newdata = parts$test,
                                   control_label = "virginica"),
    "cannot be re-pointed"
  )
  # Agreeing with it changes nothing, which is what makes naming it a way of
  # saying in the call which class the numbers are about.
  stated <- evaluate_classification_models(
    parts$models$both, newdata = parts$test,
    outcome_lv = c("versicolor", "virginica"), control_label = "versicolor"
  )
  expect_equal(stated$metrics, sa_perf_cls_fixture()$metrics[1, ])
})

test_that("models pointed at different classes are refused", {
  parts <- sa_perf_cls_parts()
  # Both return a probability from `type = "response"` and one of them is the
  # probability of the other class, so every comparison between them would be
  # reversed. Silently turning one round would be a different model from the
  # one the caller fitted and printed.
  expect_error(
    evaluate_classification_models(parts$models$both,
                                   list(flipped = parts$models$flipped),
                                   newdata = parts$test),
    "same `outcome_lv`"
  )
})

test_that("a regression is refused by name and redirected", {
  reg <- sa_perf_reg_parts()
  cls <- sa_perf_cls_parts()
  expect_error(
    evaluate_classification_models(reg$models$full, newdata = reg$test),
    "evaluate_regression_models"
  )
  expect_error(
    evaluate_classification_models(cls$models$both,
                                   list(wrong = reg$models$full),
                                   newdata = cls$test),
    "wrong \\(continuous\\)"
  )
})

test_that("scored rows holding one class are refused, not scored", {
  parts <- sa_perf_cls_parts()
  one_class <- parts$test[parts$test$Species == "virginica", ]
  expect_error(
    evaluate_classification_models(parts$models$both, newdata = one_class),
    "single class, virginica"
  )
})

test_that("a class the models were not fitted on is named and refused", {
  parts <- sa_perf_cls_parts()
  mixed <- parts$test
  mixed$Species <- as.character(mixed$Species)
  mixed$Species[1:3] <- "setosa"
  expect_error(
    evaluate_classification_models(parts$models$both, newdata = mixed),
    "setosa"
  )
})

test_that("the threshold and the confidence level are checked", {
  parts <- sa_perf_cls_parts()
  expect_error(evaluate_classification_models(parts$models$both,
                                              newdata = parts$test,
                                              threshold = 1.5),
               "`threshold` must be in \\[0, 1\\]")
  expect_error(evaluate_classification_models(parts$models$both,
                                              newdata = parts$test,
                                              conf_level = 1),
               "`conf_level` must be in \\(0, 1\\)")
})

test_that("every model is scored on the rows all of them could predict", {
  parts <- sa_perf_cls_parts()
  holed <- parts$test
  holed$Sepal.Width[1:2] <- NA
  expect_message(
    res <- evaluate_classification_models(parts$models$both,
                                          list(len = parts$models$length),
                                          newdata = holed),
    "left out"
  )
  expect_identical(res$design$n_dropped, 2L)
  # `length_only` could have predicted those two rows and is not allowed to,
  # because DeLong, IDI and NRI are paired and have no meaning across
  # different rows.
  expect_identical(res$metrics$n_used, rep(nrow(holed) - 2L, 2L))
  rows <- split(res$predictions$row, res$predictions$model)
  expect_identical(rows[[1]], rows[[2]])
})

test_that("the result survives a JSON round trip", {
  skip_if_not_installed("jsonlite")
  res <- sa_perf_cls_fixture()
  # The curves are a rectangular table rather than a curve object, which is
  # what lets the whole result be written out and read back.
  expect_false(any(vapply(res, function(s) is.object(s) && !is.data.frame(s),
                          logical(1))))

  back <- jsonlite::fromJSON(
    jsonlite::toJSON(unclass(res), digits = NA, null = "null")
  )
  expect_identical(back$analysis, res$analysis)
  expect_identical(back$design$outcome_lv, res$design$outcome_lv)
  expect_equal(back$metrics$auc, res$metrics$auc)
  expect_equal(back$comparisons$idi_pval, res$comparisons$idi_pval)
  expect_identical(back$curves$model, res$curves$model)
  expect_equal(back$curves$sensitivity, res$curves$sensitivity)
  expect_equal(back$curves$specificity, res$curves$specificity)

  # The cut above every prediction is the one cell of the object that JSON has
  # no number for, and no finite value can stand in for it: a curve reaching
  # both corners needs one cut outside the range of the predictions. What it
  # labels survives, so a reader that draws the curve loses nothing.
  opens <- res$curves$threshold == Inf
  expect_identical(sum(opens), length(res$models))
  expect_true(all(is.na(back$curves$threshold[opens])))
  expect_equal(back$curves$threshold[!opens], res$curves$threshold[!opens])
})

test_that("printing summarises without dumping the tables", {
  res <- sa_perf_cls_fixture()
  expect_output(print(res), "<sa_performance> classification_performance")
  expect_output(print(res), "probability of virginica against versicolor")
  expect_output(print(res), "threshold: 0.5")
  expect_output(print(res), "auc = ")
  expect_output(print(res), "delta_auc")
  expect_output(print(res), "IDI")
  expect_output(print(res), "NRI")
  expect_identical(withr::with_output_sink(tempfile(), print(res)), res)
})
