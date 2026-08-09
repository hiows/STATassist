# The result contract is the one thing every other part of the package and the
# planned Python port agree on, so it is tested across all four scenarios at
# once rather than scenario by scenario.

test_that("every scenario records the scale its effect table was built on", {
  for (res in sa_all_scenarios()) {
    expect_identical(res$parameters$input_scale, "raw", info = res$analysis)
    expect_true(res$parameters$fc_mean %in% c("arith", "geom"),
                info = res$analysis)
  }
})

test_that("every scenario carries the full set of contract slots", {
  slots <- c("analysis", "features", "design", "parameters", "effect", "tests",
             "test_info", "diagnostics", "metadata")
  for (res in sa_all_scenarios()) {
    expect_true(all(slots %in% names(res)), info = res$analysis)
    expect_s3_class(res, "sa_comparison")
    expect_s3_class(res, "sa_result")
  }
})

test_that("no scenario reports a schema version", {
  # The field described the shape of the object to a consumer in another
  # language. Nothing in R ever read it, so it was noise in a printed result.
  for (res in sa_all_scenarios()) {
    expect_false("schema_version" %in% names(res), info = res$analysis)
  }
  expect_false("schema_version" %in% names(diagnose_distribution(
    iris, "Sepal.Length", iris$Species)))
})

test_that("every test table carries the guaranteed columns", {
  for (res in sa_all_scenarios()) {
    for (nm in names(res$tests)) {
      expect_true(all(sa_test_table_columns() %in% names(res$tests[[nm]])),
                  info = paste(res$analysis, nm))
    }
  }
})

test_that("every table has one row per feature, in the same order", {
  for (res in sa_all_scenarios()) {
    expect_identical(res$effect$features, res$features)
    for (nm in names(res$tests)) {
      expect_identical(res$tests[[nm]]$features, res$features,
                       info = paste(res$analysis, nm))
    }
  }
})

test_that("tests and test_info name exactly the same tests", {
  for (res in sa_all_scenarios()) {
    expect_setequal(names(res$tests), names(res$test_info))
  }
})

test_that("post-hoc tables carry the guaranteed columns and known features", {
  for (res in sa_all_scenarios()) {
    for (nm in names(res$posthoc)) {
      tbl <- res$posthoc[[nm]]
      expect_true(all(sa_posthoc_table_columns() %in% names(tbl)),
                  info = paste(res$analysis, nm))
      expect_true(all(tbl$features %in% res$features))
      expect_true(all(names(res$posthoc) %in% names(res$tests)))
    }
  }
})

test_that("scenarios without a pairwise stage carry no posthoc slot at all", {
  # Absent rather than empty. A two-group comparison is already the only
  # contrast there is, and a one-sample comparison has no pair of levels, so
  # there is no question for the slot to answer emptily.
  for (res in list(sa_two_group_fixture(), sa_one_sample_fixture())) {
    expect_false("posthoc" %in% names(res), info = res$analysis)
    expect_false("pairwise" %in% names(res), info = res$analysis)
    # Reading the slot anyway still has to be safe, since that is what
    # draw_forest_plot() and print() do before knowing which scenario ran.
    expect_null(res$posthoc[["t_test"]])
  }
})

test_that("pairwise and posthoc are populated together, test for test", {
  for (res in sa_all_scenarios()) {
    # identical() rather than setequal(): both are NULL where the scenario has
    # no pairwise stage, and setequal() has nothing to compare there.
    expect_identical(names(res$pairwise), names(res$posthoc))
  }
})

test_that("a post-hoc estimate reads as group1 minus group2", {
  res <- sa_multi_group_fixture()
  tbl <- res$posthoc$anova_test
  # Setosa is the reference, so it is the level being subtracted rather than the
  # one doing the subtracting.
  row <- tbl[tbl$features == "Petal.Length" &
               tbl$group1 == "virginica" & tbl$group2 == "setosa", ]
  expect_equal(row$estimate,
               mean(iris$Petal.Length[iris$Species == "virginica"]) -
                 mean(iris$Petal.Length[iris$Species == "setosa"]))
  expect_identical(row$contrast, "virginica - setosa")
})

test_that("an interval always contains its own estimate", {
  res <- sa_multi_group_fixture()
  for (tbl in res$posthoc) {
    expect_true(all(tbl$lower_conf <= tbl$estimate))
    expect_true(all(tbl$estimate <= tbl$upper_conf))
  }
})

test_that("the result survives a JSON round trip", {
  skip_if_not_installed("jsonlite")
  res <- sa_multi_group_fixture()
  # unclass() because jsonlite dispatches on class and the S3 class carries no
  # data. Everything below it is a plain list, vector or data.frame, which is
  # the property the port depends on.
  back <- jsonlite::fromJSON(
    jsonlite::toJSON(unclass(res), digits = NA, null = "null")
  )

  expect_identical(back$analysis, res$analysis)
  expect_identical(back$features, res$features)
  expect_equal(back$effect$log2fc, res$effect$log2fc)
  expect_equal(back$tests$anova_test$pval, res$tests$anova_test$pval)
  expect_equal(back$posthoc$anova_test$estimate, res$posthoc$anova_test$estimate)
  expect_identical(back$posthoc$anova_test$contrast,
                   res$posthoc$anova_test$contrast)

  # The pairwise slot is two levels of map deep rather than one, so it is the
  # slot most likely to come back as an array and quietly lose its keys.
  ct <- names(res$pairwise$anova_test)[1]
  expect_identical(names(back$pairwise$anova_test), names(res$pairwise$anova_test))
  expect_equal(back$pairwise$anova_test[[ct]]$log2fc,
               res$pairwise$anova_test[[ct]]$log2fc)
})

test_that("no scenario stores a fitted model or other R-only object", {
  # The whole result has to be writable as JSON, so a htest or an lm hiding in
  # a slot would break the port without breaking any other test.
  allowed <- c("list", "character", "numeric", "integer", "logical",
               "data.frame", "NULL")
  walk <- function(x, path) {
    cls <- class(x)[1]
    expect_true(cls %in% allowed, info = paste(path, "is a", cls))
    if (is.list(x) && !is.data.frame(x)) {
      for (nm in names(x)) walk(x[[nm]], paste0(path, "$", nm))
    }
  }
  for (res in sa_all_scenarios()) {
    walk(unclass(res), res$analysis)
  }
})

test_that("print works for every scenario", {
  for (res in sa_all_scenarios()) {
    expect_output(print(res), "features")
    expect_invisible(print(res))
  }
  # The one-sample header has no group levels to show, which is what broke the
  # 0.1.0 print method.
  expect_output(print(sa_one_sample_fixture()), "mu")
})

test_that("estimate_significance still reads a v0.1.0 two-group result", {
  res <- sa_two_group_fixture()
  sig <- estimate_significance(res, test = "t_test")
  expect_identical(names(sig), c("analysis_type", "significance"))
  expect_identical(sig$analysis_type, res$analysis)
  expect_identical(names(sig$significance),
                   c("features", "log2fc", "pvalue", "adj_pvalue", "is_signif"))
  expect_identical(sig$significance$features, res$features)
})

test_that("estimate_significance reads a multi-group result unchanged", {
  # The point of keeping the effect table at one row per feature with a signed
  # log2fc: the volcano path needs no knowledge of how many groups there were.
  res <- sa_multi_group_fixture()
  sig <- estimate_significance(res, test = "kruskal_test",
                               log2fc_cutoff = 1)$significance
  expect_identical(sig$features, res$features)
  expect_true(all(sig$is_signif[sig$features %in%
                                  c("Petal.Length", "Petal.Width")]))
})

test_that("every scenario's verdict names the analysis it came from", {
  for (res in sa_all_scenarios()) {
    sig <- estimate_significance(res)
    expect_s3_class(sig, c("sa_significance", "sa_result"), exact = TRUE)
    expect_identical(names(sig), c("analysis_type", "significance"))
    expect_identical(sig$analysis_type, res$analysis)
    expect_s3_class(sig$significance, "data.frame")
  }
})

test_that("the cutoffs stay on the table the wrapper holds", {
  # draw_volcano_plot() reads them from there, so promoting the analysis name out
  # of the attributes must not have taken the rest with it.
  sig <- estimate_significance(sa_two_group_fixture(), test = "wilcox_test",
                               log2fc_cutoff = 0.25, pval_cutoff = 0.01)
  expect_identical(attr(sig$significance, "log2fc_cutoff"), 0.25)
  expect_identical(attr(sig$significance, "pval_cutoff"), 0.01)
  expect_identical(attr(sig$significance, "test"), "wilcox_test")
  expect_identical(attr(sig$significance, "analysis"), sig$analysis_type)
})

test_that("printing a verdict summarises the rule rather than the table", {
  sig <- estimate_significance(sa_two_group_fixture(), log2fc_cutoff = 0.1)
  expect_output(print(sig), "two_group_comparison")
  expect_output(print(sig), "abs\\(log2fc\\) >= 0.1")
  expect_output(print(sig), "significant")

  by_pair <- estimate_significance(sa_multi_group_fixture(), by = "contrast")
  expect_output(print(by_pair), "one table per contrast")
  expect_output(print(by_pair), "virginica - setosa")
})
