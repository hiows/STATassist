# The v0.0.1 surface, pinned so that the 0.2.0 contract extension and the new
# `diagnose` argument cannot change any of it without a test saying so.

test_that("the two-group tests still match their base R equivalents", {
  res <- sa_two_group_fixture()
  d <- iris[iris$Species != "setosa", ]
  x <- d$Petal.Length[d$Species == "virginica"]
  y <- d$Petal.Length[d$Species == "versicolor"]
  row <- res$features == "Petal.Length"

  reference <- stats::t.test(x, y)
  expect_equal(res$tests$t_test$t_stat[row], unname(reference$statistic))
  expect_equal(res$tests$t_test$pval[row], reference$p.value)

  ranked <- suppressWarnings(stats::wilcox.test(x, y, conf.int = TRUE))
  expect_equal(res$tests$wilcox_test$pval[row], ranked$p.value)
})

test_that("the fold change still puts group_lv[1] over group_lv[2]", {
  res <- sa_two_group_fixture()
  d <- iris[iris$Species != "setosa", ]
  row <- res$effect[res$effect$features == "Petal.Length", ]
  expect_equal(row$x_center, mean(d$Petal.Length[d$Species == "virginica"]))
  expect_equal(row$y_center, mean(d$Petal.Length[d$Species == "versicolor"]))
  expect_equal(row$fold_change, row$x_center / row$y_center)
  expect_gt(row$log2fc, 0)
})

test_that("the paired two-group path is unchanged", {
  res <- compare_two_groups(
    data.frame(v = sleep$extra), "v", paste0("g", sleep$group),
    c("g2", "g1"), id = sleep$ID, paired = TRUE, diagnose = FALSE
  )
  wide <- stats::reshape(sleep, idvar = "ID", timevar = "group",
                         direction = "wide")
  reference <- stats::t.test(wide$extra.2, wide$extra.1, paired = TRUE)
  expect_equal(res$tests$t_test$t_stat, unname(reference$statistic))
  expect_equal(res$tests$t_test$pval, reference$p.value)
  expect_identical(res$test_info$robust_test$id, "yuen_paired")
})

test_that("adding the diagnose argument changed no other output", {
  d <- iris[iris$Species != "setosa", ]
  without <- compare_two_groups(d, sa_feats(), d$Species,
                                c("virginica", "versicolor"), diagnose = FALSE)
  with <- compare_two_groups(d, sa_feats(), d$Species,
                             c("virginica", "versicolor"), diagnose = TRUE)

  expect_equal(without$effect, with$effect)
  expect_equal(without$tests, with$tests)
  expect_identical(without$test_info, with$test_info)
  expect_null(without$diagnostics)
  expect_false(is.null(with$diagnostics))
})

test_that("the descriptive summary is unchanged", {
  out <- summarize_descriptive_stats(iris, sa_feats(), iris$Species)
  expect_true(all(c("features", "group", "n", "mean", "sd", "median",
                    "skewness", "excess_kurtosis") %in% names(out)))
  setosa <- iris$Sepal.Length[iris$Species == "setosa"]
  row <- out[out$features == "Sepal.Length" & out$group == "setosa", ]
  expect_equal(row$mean, mean(setosa))
  expect_equal(row$sd, stats::sd(setosa))
  expect_equal(row$median, stats::median(setosa))
})

test_that("the significance table and the volcano input are unchanged", {
  res <- sa_two_group_fixture()
  sig <- estimate_significance(res, test = "t_test", log2fc_cutoff = 0.5,
                               pval_cutoff = 0.01)
  expect_identical(names(sig),
                   c("features", "log2fc", "pvalue", "adj_pvalue", "is_signif"))
  expect_identical(sig$log2fc, res$effect$log2fc)
  expect_identical(sig$pvalue, res$tests$t_test$pval)
  expect_identical(sig$is_signif,
                   abs(sig$log2fc) >= 0.5 & sig$adj_pvalue <= 0.01)
})

test_that("the drawing functions still return their summary statistics", {
  skip_if_not_installed("withr")
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 700, height = 500)
  withr::defer(grDevices::dev.off())

  stats_out <- draw_grouped_boxplot(iris, sa_feats(), iris$Species,
                                    levels(iris$Species))
  expect_named(stats_out, c("box_summary_stats", "median_confidence_stats"))
  expect_identical(names(stats_out$box_summary_stats), sa_feats())
})

test_that("the registry still ships and still parses", {
  skip_if_not_installed("jsonlite")
  path <- system.file("extdata", "registry.json", package = "STATassist")
  skip_if(path == "")
  reg <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  expect_true(length(reg$methods) > 0L)
  ids <- vapply(reg$methods, function(m) m$id, character(1))
  expect_false(anyDuplicated(ids) > 0L)
})
