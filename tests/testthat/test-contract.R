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
    # A factorial comparison is the exception, and deliberately: its contrast
    # axis is factor by stratum, which a list keyed by contrast label cannot
    # hold, so it has contrasts without a second view of them.
    if (!is.null(res$terms)) {
      expect_false("pairwise" %in% names(res), info = res$analysis)
      next
    }
    # identical() rather than setequal(): both are NULL where the scenario has
    # no pairwise stage, and setequal() has nothing to compare there.
    expect_identical(names(res$pairwise), names(res$posthoc))
  }
})

test_that("the term axis exists exactly where a model has several terms", {
  for (res in sa_all_scenarios()) {
    if (identical(res$analysis, "factorial_comparison")) {
      expect_true(all(sa_term_table_columns() %in% names(res$terms)))
      expect_true(all(res$terms$features %in% res$features))
    } else {
      expect_false("terms" %in% names(res), info = res$analysis)
    }
  }
})

test_that("the cell axis exists exactly where a design has cells", {
  for (res in sa_all_scenarios()) {
    if (identical(res$analysis, "factorial_comparison")) {
      expect_true(all(sa_cell_table_columns() %in% names(res$cells)))
      expect_identical(unique(res$cells$features), res$features)
      expect_identical(nrow(res$cells),
                       length(res$features) * length(res$design$group_lv))
      # The same cells in the same order the design lists them in, so a plot of
      # this table and a post-hoc stratum cannot disagree about what a cell is.
      expect_identical(unique(res$cells$cell), res$design$group_lv)
      # And one column per factor, holding the level name rather than an index,
      # which is what lets a reader subset on a factor without parsing a label.
      for (nm in names(res$design$factor_lv)) {
        expect_setequal(unique(res$cells[[nm]]), res$design$factor_lv[[nm]])
      }
    } else {
      expect_false("cells" %in% names(res), info = res$analysis)
    }
  }
})

test_that("a cell mean is the mean the crossed model was fitted on", {
  res <- sa_factorial_fixture()
  cells <- res$cells
  by_hand <- tapply(warpbreaks$breaks,
                    list(warpbreaks$wool, warpbreaks$tension), mean)
  for (i in seq_len(nrow(cells))) {
    expect_equal(cells$mean[i],
                 by_hand[cells$wool[i], cells$tension[i]],
                 info = cells$cell[i])
  }
  # With no missing value the per-feature count is the row count of the cell.
  expect_identical(cells$n, as.integer(res$design$cell_n))
  expect_equal(cells$sd, tapply(warpbreaks$breaks,
                               list(warpbreaks$wool, warpbreaks$tension),
                               stats::sd)[cbind(cells$wool, cells$tension)],
               ignore_attr = TRUE)
})

test_that("a marginal mean of the cell table is the one Tukey contrasts", {
  # `cells` stores its `se` pooled over the whole model rather than within the
  # cell so that the mean and the variance of any marginal mean are both
  # recoverable from it. That is what an interaction plot's points and bars rest
  # on, so it is checked against the contrasts the post-hoc stage reports rather
  # than against a formula written out a second time.
  res <- sa_factorial_fixture()
  cells <- res$cells[res$cells$features == res$features[1], ]
  ph <- res$posthoc$anova_test
  ph <- ph[is.na(ph$stratum), ]
  expect_gt(nrow(ph), 0L)

  for (i in seq_len(nrow(ph))) {
    fct <- ph$factor[i]
    mean_of <- tapply(cells$mean, cells[[fct]], mean)
    var_of <- function(level) {
      at <- cells[[fct]] == level
      sum(cells$se[at]^2) / sum(at)^2
    }
    expect_equal(unname(mean_of[ph$group1[i]] - mean_of[ph$group2[i]]),
                 ph$estimate[i], info = ph$contrast[i])
    expect_equal(sqrt((var_of(ph$group1[i]) + var_of(ph$group2[i])) / 2),
                 ph$stderr[i], info = ph$contrast[i])
  }
})

test_that("the cell table survives a JSON round trip", {
  skip_if_not_installed("jsonlite")
  res <- sa_factorial_fixture()
  back <- jsonlite::fromJSON(
    jsonlite::toJSON(unclass(res), digits = NA, null = "null")
  )
  # The factor columns are named after the factors rather than fixed, so they
  # are the ones a column-oriented writer is most likely to drop.
  expect_identical(names(back$cells), names(res$cells))
  expect_identical(back$cells$cell, res$cells$cell)
  expect_equal(back$cells$mean, res$cells$mean)
  expect_equal(back$cells$se, res$cells$se)
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
  expect_identical(sig$extreme_level, res$effect$extreme_level)
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
    expect_identical(sig$significance$features, res$features)
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

test_that("a categorical result is outside the comparison contract", {
  res <- sa_categorical_fixture()
  expect_s3_class(res, "sa_categorical")
  expect_false(inherits(res, "sa_comparison"))
  expect_false("features" %in% names(res))
  expect_false("effect" %in% names(res))
  expect_false(res$analysis %in% vapply(sa_all_scenarios(), `[[`, "", "analysis"))
})

test_that("every categorical design names a null the contract defines", {
  designs <- list(
    independent = sa_categorical_fixture(),
    matched     = sa_categorical_matched_fixture(),
    repeated    = sa_categorical_repeated_fixture()
  )
  for (nm in names(designs)) {
    null <- designs[[nm]]$design$null
    expect_true(is.character(null) && length(null) == 1L, info = nm)
    expect_true(null %in% sa_categorical_nulls(), info = nm)
  }
  # Three designs, three nulls, so `expected` never means the same thing twice
  # by accident.
  expect_length(unique(vapply(designs, function(r) r$design$null, "")), 3L)
})

test_that("estimate_significance refuses a categorical result", {
  expect_error(
    estimate_significance(sa_categorical_fixture()),
    "no feature axis"
  )
  # And names the function that does read one, so the refusal is a redirection
  # rather than a dead end.
  expect_error(
    estimate_significance(sa_categorical_fixture()),
    "estimate_categorical_significance"
  )
})

test_that("a categorical result survives a JSON round trip whole", {
  skip_if_not_installed("jsonlite")
  res <- sa_categorical_fixture()
  # Nothing has to be dropped first. The table is reached through `as.table()`
  # rather than stored, which is what leaves every slot a scalar, a character
  # vector, a named list or a data.frame.
  expect_false(any(vapply(res, inherits, logical(1), "table")))
  expect_false(any(vapply(res, is.matrix, logical(1))))

  back <- jsonlite::fromJSON(
    jsonlite::toJSON(unclass(res), digits = NA, null = "null")
  )
  expect_identical(back$analysis, res$analysis)
  expect_identical(back$variables, res$variables)
  expect_identical(back$design$null, res$design$null)
  expect_identical(back$design$dim, res$design$dim)
  expect_equal(back$cells$observed, res$cells$observed)
  expect_equal(back$cells$expected, res$cells$expected)
  expect_equal(back$tests$chisq_test$pval, res$tests$chisq_test$pval)
  expect_identical(back$association$measure, res$association$measure)

  # And with the approximation check attached, which is the one slot whose shape
  # depends on the design.
  diagnosed <- compare_categorical_groups(sa_categorical_frame())
  round_tripped <- jsonlite::fromJSON(
    jsonlite::toJSON(unclass(diagnosed), digits = NA, null = "null")
  )
  expect_identical(round_tripped$diagnostics$rule, "expected_count_min")
  expect_equal(round_tripped$diagnostics$min_expected,
               diagnosed$diagnostics$min_expected)
})

test_that("the table of a categorical result is rebuilt from the cells alone", {
  res <- sa_categorical_fixture()
  # Which is what makes the round trip above lossless: the cells are the table,
  # so a reader in another language loses nothing by not having `as.table()`.
  expect_equal(sa_categorical_table(res$cells, "smoker", "grade"),
               as.table(res))
})
