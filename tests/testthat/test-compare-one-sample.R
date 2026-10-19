test_that("the t-test matches stats::t.test", {
  res <- sa_one_sample_fixture()
  reference <- stats::t.test(iris$Sepal.Length, mu = 3)
  row <- res$tests$t_test[res$features == "Sepal.Length", ]

  expect_equal(row$t_stat, unname(reference$statistic))
  expect_equal(row$df, unname(reference$parameter))
  expect_equal(row$pval, reference$p.value)
  expect_equal(row$lower_conf, reference$conf.int[1])
  expect_equal(row$upper_conf, reference$conf.int[2])
})

test_that("the Wilcoxon test matches stats::wilcox.test", {
  res <- sa_one_sample_fixture()
  reference <- suppressWarnings(
    stats::wilcox.test(iris$Petal.Width, mu = 3, conf.int = TRUE)
  )
  row <- res$tests$wilcox_test[res$features == "Petal.Width", ]

  expect_equal(row$v_stat, unname(reference$statistic))
  expect_equal(row$pval, reference$p.value)
  expect_equal(row$hl_shift, unname(reference$estimate))
})

test_that("the proportion test matches stats::prop.test on a binary feature", {
  cars <- mtcars
  cars$am <- as.numeric(cars$am)
  res <- suppressWarnings(
    compare_one_sample(cars, c("am", "mpg"), mu = 0.5, p = 0.4,
                       diagnose = FALSE)
  )
  reference <- stats::prop.test(sum(mtcars$am == 1), nrow(mtcars), p = 0.4)
  row <- res$tests$prop_test[res$features == "am", ]

  expect_equal(row$n_success, sum(mtcars$am == 1))
  expect_equal(row$chi_sq, unname(reference$statistic))
  expect_equal(row$pval, reference$p.value)
  expect_equal(row$lower_conf, reference$conf.int[1])
  # The Wilson interval cannot leave the unit interval, which is the reason for
  # preferring it over a Wald interval.
  expect_gte(row$lower_conf, 0)
  expect_lte(row$upper_conf, 1)
})

test_that("a non-binary feature comes back NA rather than coerced", {
  expect_warning(
    res <- compare_one_sample(iris, "Sepal.Length", mu = 3, diagnose = FALSE),
    "binary"
  )
  expect_true(is.na(res$tests$prop_test$pval))
  # The other two tests are unaffected by the one that could not run.
  expect_false(is.na(res$tests$t_test$pval))
  expect_false(is.na(res$tests$wilcox_test$pval))
})

test_that("the design records mu instead of group levels", {
  res <- sa_one_sample_fixture()
  expect_identical(res$design$mu, 3)
  expect_null(res$design$group_lv)
  expect_false(res$design$paired)
})

test_that("the effect table reports a distance from mu and no fold change", {
  res <- sa_one_sample_fixture()
  expect_identical(names(res$effect),
                   c("features", "n_used", "center", "mu", "diff"))
  row <- res$effect[res$effect$features == "Sepal.Length", ]
  expect_equal(row$center, mean(iris$Sepal.Length))
  expect_equal(row$mu, 3)
  expect_equal(row$diff, mean(iris$Sepal.Length) - 3)
})

test_that("mu = 0 is an ordinary reference for a distance", {
  expect_no_message(
    res <- suppressWarnings(
      compare_one_sample(iris, "Sepal.Length", mu = 0, diagnose = FALSE)
    )
  )
  expect_equal(res$effect$diff, mean(iris$Sepal.Length))
  expect_false(is.na(res$tests$t_test$pval))
})

test_that("a one-sample verdict carries each test's own distance", {
  res <- sa_one_sample_fixture()

  t_sig <- estimate_significance(res, test = "t_test")$significance
  expect_identical(names(t_sig),
                   c("features", "delta_mu", "pvalue", "adj_pvalue",
                     "is_signif"))
  expect_equal(t_sig$delta_mu, res$tests$t_test$diff)
  expect_identical(attr(t_sig, "reference"), 3)
  expect_identical(attr(t_sig, "estimate_label"), "mean")
  expect_null(attr(t_sig, "log2fc_cutoff"))

  w_sig <- estimate_significance(res, test = "wilcox_test")$significance
  expect_equal(w_sig$delta_mu, res$tests$wilcox_test$hl_shift - 3)
  expect_identical(attr(w_sig, "estimate_label"), "pseudo-median")

  flags <- data.frame(am = mtcars$am, vs = mtcars$vs)
  prop <- suppressWarnings(
    compare_one_sample(flags, c("am", "vs"), mu = 0.5, p = 0.3,
                       diagnose = FALSE)
  )
  p_sig <- estimate_significance(prop, test = "prop_test")$significance
  expect_true("delta_p" %in% names(p_sig))
  expect_false("delta_mu" %in% names(p_sig))
  expect_equal(p_sig$delta_p, colMeans(flags) - 0.3, ignore_attr = TRUE)
  expect_identical(attr(p_sig, "reference"), 0.3)
})

test_that("a one-sample verdict is decided by the p-value alone", {
  res <- sa_one_sample_fixture()
  loose <- estimate_significance(res, log2fc_cutoff = 0)$significance
  strict <- estimate_significance(res, log2fc_cutoff = 100)$significance
  expect_identical(loose$is_signif, strict$is_signif)
  expect_identical(strict$is_signif, strict$adj_pvalue <= 0.05)

  # A test that could not run leaves the verdict undecided, not negative.
  na_sig <- estimate_significance(res, test = "prop_test")$significance
  expect_true(all(is.na(na_sig$is_signif)))
})

test_that("printing a one-sample verdict shows the p-value cutoff only", {
  sig <- estimate_significance(sa_one_sample_fixture())
  expect_output(print(sig), "cutoffs  : adj_pvalue <= 0.05  \\(BH\\)")
  out <- capture.output(print(sig))
  expect_false(any(grepl("log2fc|delta", out)))
})

test_that("a one-sided alternative leaves the untested side open", {
  res <- suppressWarnings(
    compare_one_sample(mtcars, "mpg", mu = 20, alternative = "greater",
                       diagnose = FALSE)
  )
  expect_equal(res$tests$t_test$upper_conf, Inf)
  expect_equal(res$tests$t_test$pval,
               stats::t.test(mtcars$mpg, mu = 20,
                             alternative = "greater")$p.value)
})

test_that("missing values are dropped per feature", {
  d <- iris
  d$Sepal.Length[1:10] <- NA
  res <- suppressWarnings(compare_one_sample(d, sa_feats(), mu = 3,
                                             diagnose = FALSE))
  expect_equal(res$tests$t_test$n_used, c(140, 150, 150, 150))
})
