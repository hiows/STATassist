# The result contract is the one thing every other part of the package and the
# planned Python port agree on, so it is tested across all four scenarios at
# once rather than scenario by scenario.

test_that("every scenario reports schema version 0.2.0", {
  for (res in sa_all_scenarios()) {
    expect_identical(res$schema_version, "0.2.0")
  }
})

test_that("every scenario carries the full set of contract slots", {
  slots <- c("schema_version", "analysis", "features", "design", "parameters",
             "effect", "tests", "posthoc", "test_info", "diagnostics",
             "metadata")
  for (res in sa_all_scenarios()) {
    expect_true(all(slots %in% names(res)), info = res$analysis)
    expect_s3_class(res, "sa_comparison")
    expect_s3_class(res, "sa_result")
  }
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

test_that("scenarios without a pairwise stage report an empty posthoc slot", {
  # Not NULL: a consumer should be able to ask any result how many contrasts it
  # holds without first asking which scenario produced it. Named rather than
  # bare, so that it stays a JSON object once serialised.
  empty <- structure(list(), names = character(0))
  expect_identical(sa_two_group_fixture()$posthoc, empty)
  expect_identical(sa_one_sample_fixture()$posthoc, empty)
})

test_that("a post-hoc estimate reads as group1 minus group2", {
  res <- sa_multi_group_fixture()
  tbl <- res$posthoc$anova_test
  row <- tbl[tbl$features == "Petal.Length" &
               tbl$group1 == "setosa" & tbl$group2 == "virginica", ]
  expect_equal(row$estimate,
               mean(iris$Petal.Length[iris$Species == "setosa"]) -
                 mean(iris$Petal.Length[iris$Species == "virginica"]))
  expect_identical(row$contrast, "setosa - virginica")
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

  expect_identical(back$schema_version, res$schema_version)
  expect_identical(back$features, res$features)
  expect_equal(back$effect$log2fc, res$effect$log2fc)
  expect_equal(back$tests$anova_test$pval, res$tests$anova_test$pval)
  expect_equal(back$posthoc$anova_test$estimate, res$posthoc$anova_test$estimate)
  expect_identical(back$posthoc$anova_test$contrast,
                   res$posthoc$anova_test$contrast)
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
  expect_identical(names(sig),
                   c("features", "log2fc", "pvalue", "adj_pvalue", "is_signif"))
  expect_identical(sig$features, res$features)
})

test_that("estimate_significance reads a multi-group result unchanged", {
  # The point of keeping the effect table at one row per feature with a signed
  # log2fc: the volcano path needs no knowledge of how many groups there were.
  res <- sa_multi_group_fixture()
  sig <- estimate_significance(res, test = "kruskal_test", log2fc_cutoff = 1)
  expect_identical(sig$features, res$features)
  expect_true(all(sig$is_signif[sig$features %in%
                                  c("Petal.Length", "Petal.Width")]))
})
