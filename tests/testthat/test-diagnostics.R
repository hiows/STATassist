test_that("the normality tests match their base R equivalents", {
  d <- diagnose_distribution(iris, "Sepal.Length", iris$Species)
  setosa <- iris$Sepal.Length[iris$Species == "setosa"]
  row <- d$normality[d$normality$group == "setosa", ]

  reference <- stats::shapiro.test(setosa)
  expect_equal(row$shapiro_stat, unname(reference$statistic))
  expect_equal(row$shapiro_pval, reference$p.value)

  ks <- suppressWarnings(
    stats::ks.test(setosa, "pnorm", mean(setosa), stats::sd(setosa))
  )
  expect_equal(row$ks_stat, unname(ks$statistic))
})

test_that("Bartlett's test matches stats::bartlett.test", {
  d <- diagnose_distribution(iris, "Sepal.Length", iris$Species)
  reference <- stats::bartlett.test(Sepal.Length ~ Species, data = iris)
  expect_equal(d$variance$bartlett_stat, unname(reference$statistic))
  expect_equal(d$variance$bartlett_pval, reference$p.value)
})

test_that("Levene's test is an ANOVA on distances from the group median", {
  d <- diagnose_distribution(iris, "Sepal.Length", iris$Species)
  deviation <- abs(iris$Sepal.Length -
                     stats::ave(iris$Sepal.Length, iris$Species,
                                FUN = stats::median))
  reference <- summary(stats::aov(deviation ~ iris$Species))[[1]]

  expect_equal(d$variance$levene_stat, reference[["F value"]][1])
  expect_equal(d$variance$levene_pval, reference[["Pr(>F)"]][1])
  expect_equal(d$variance$levene_df1, reference[["Df"]][1])
})

test_that("centring on the mean gives a different Levene result", {
  # Not a formality: the median-centred Brown-Forsythe variant exists because
  # the mean-centred one is itself sensitive to the skew it is meant to survive.
  median_centred <- diagnose_distribution(iris, "Petal.Length", iris$Species)
  mean_centred <- diagnose_distribution(iris, "Petal.Length", iris$Species,
                                        center = "mean")
  expect_false(isTRUE(all.equal(median_centred$variance$levene_stat,
                                mean_centred$variance$levene_stat)))
})

test_that("normality is reported per feature and level", {
  d <- diagnose_distribution(iris, sa_feats(), iris$Species)
  expect_equal(nrow(d$normality), 4L * 3L)
  expect_identical(unique(d$normality$group),
                   c("setosa", "versicolor", "virginica"))
  expect_equal(nrow(d$variance), 4L)
  expect_equal(nrow(d$summary), 4L)
})

test_that("without a grouping there is nothing to compare variances across", {
  d <- diagnose_distribution(iris, "Petal.Length")
  expect_equal(nrow(d$normality), 1L)
  expect_equal(nrow(d$variance), 0L)
  expect_true(is.na(d$summary$variance_ok))
  # Pooled over species, petal length is strongly bimodal and every normality
  # test says so.
  expect_lt(d$normality$shapiro_pval, 0.001)
})

test_that("the summary flags take the worst level of each feature", {
  d <- diagnose_distribution(iris, sa_feats(), iris$Species, alpha = 0.05)
  for (f in sa_feats()) {
    block <- d$normality[d$normality$features == f, ]
    row <- d$summary[d$summary$features == f, ]
    expect_equal(row$min_shapiro_pval, min(block$shapiro_pval))
    expect_identical(row$normal_ok, min(block$shapiro_pval) > 0.05)
  }
  # Setosa petal width is the one clear departure in the data set.
  expect_false(d$summary$normal_ok[d$summary$features == "Petal.Width"])
})

test_that("a level too small to test yields NA rather than disappearing", {
  d <- iris[c(1:2, 51:100, 101:150), ]
  diag <- diagnose_distribution(d, "Sepal.Length", d$Species)
  row <- diag$normality[diag$normality$group == "setosa", ]
  expect_equal(row$n_used, 2)
  expect_true(is.na(row$shapiro_stat))
  expect_equal(nrow(diag$normality), 3L)
})

test_that("the IQR screen flags exactly what the Tukey fences do", {
  flagged <- screen_outliers(iris, "Sepal.Width")
  q <- unname(stats::quantile(iris$Sepal.Width, c(0.25, 0.75)))
  fence <- 1.5 * (q[2] - q[1])
  expected <- which(iris$Sepal.Width < q[1] - fence |
                      iris$Sepal.Width > q[2] + fence)
  expect_identical(flagged$row, expected)
  expect_identical(attr(flagged, "criterion"), "iqr")
})

test_that("the IQR score is measured in IQR units and clears the multiplier", {
  flagged <- screen_outliers(iris, "Sepal.Width")
  expect_true(all(flagged$score > 1.5))
  wider <- screen_outliers(iris, "Sepal.Width", iqr_multiplier = 3)
  # The score does not depend on the multiplier, only which rows survive it.
  expect_true(all(wider$score > 3))
  expect_lt(nrow(wider), nrow(flagged))
})

test_that("screening within groups is not the same as screening pooled", {
  pooled <- screen_outliers(iris, sa_feats())
  grouped <- screen_outliers(iris, sa_feats(), iris$Species)
  # The between-species spread widens the pooled fences enough to hide flags
  # that are obvious inside a species, which is why the grouped call exists.
  expect_gt(nrow(grouped), nrow(pooled))
  expect_true(all(is.na(pooled$group)))
  expect_false(any(is.na(grouped$group)))
})

test_that("row numbers refer to the data that was passed in", {
  flagged <- screen_outliers(iris, "Sepal.Width", iris$Species,
                             group_lv = c("versicolor", "virginica"))
  # Setosa rows were dropped, so a naive index would be off by 50.
  expect_true(all(flagged$row > 50))
  expect_true(all(iris$Sepal.Width[flagged$row] == flagged$value))
})

test_that("nothing is ever removed, only flagged", {
  flagged <- screen_outliers(iris, sa_feats(), iris$Species)
  expect_s3_class(flagged, "data.frame")
  expect_identical(names(flagged),
                   c("features", "group", "row", "value", "score"))
})

test_that("Grubbs flags at most one observation per group", {
  flagged <- screen_outliers(iris, sa_feats(), iris$Species,
                             criterion = "grubbs", alpha = 0.5)
  counted <- table(flagged$features, flagged$group)
  expect_true(all(counted <= 1L))
})

test_that("an unknown criterion is rejected before anything is computed", {
  expect_error(screen_outliers(iris, "Sepal.Width", criterion = "zscore"),
               "should be one of")
})

test_that("comparisons attach diagnostics computed on what they tested", {
  res <- compare_multiple_groups(iris, sa_feats(), iris$Species,
                                 levels(iris$Species), diagnose = TRUE)
  expect_named(res$diagnostics, c("normality", "variance", "summary"))
  expect_equal(nrow(res$diagnostics$normality), 4L * 3L)

  standalone <- diagnose_distribution(iris, sa_feats(), iris$Species)
  expect_equal(res$diagnostics$normality$shapiro_pval,
               standalone$normality$shapiro_pval)
})

test_that("a within-subject design is not given a variance table", {
  # Homogeneity across independent groups is not the assumption a repeated
  # measures test makes; sphericity is, and it is already in the ANOVA row.
  d <- sa_repeated_frame()
  res <- compare_multiple_groups(d["value"], "value", d$cond,
                                 sort(unique(d$cond)), id = d$subj,
                                 paired = TRUE, diagnose = TRUE)
  expect_equal(nrow(res$diagnostics$variance), 0L)
  expect_equal(nrow(res$diagnostics$normality), 4L)
  expect_false(is.na(res$tests$anova_test$mauchly_pval))
})

test_that("diagnose = FALSE leaves the slot empty", {
  res <- sa_multi_group_fixture()
  expect_null(res$diagnostics)
})

test_that("print works for a diagnosis", {
  d <- diagnose_distribution(iris, sa_feats(), iris$Species)
  expect_output(print(d), "sa_diagnosis")
  expect_output(print(d), "normality")
  expect_invisible(print(d))
})
