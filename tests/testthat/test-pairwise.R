# The pairwise slot is the post-hoc table seen one contrast at a time. It holds
# no number the post-hoc stage did not already produce, apart from the ratio
# columns, so most of what is worth testing is that the two views cannot drift:
# same contrasts, same estimates, same p-values, and a log2fc that points the
# same way as the estimate sitting beside it.

test_that("there is one table per pair of levels, named by contrast", {
  res <- sa_multi_group_fixture()
  k <- length(res$design$group_lv)

  for (nm in names(res$pairwise)) {
    tables <- res$pairwise[[nm]]
    expect_length(tables, choose(k, 2))
    expect_identical(names(tables), unique(res$posthoc[[nm]]$contrast),
                     info = nm)
  }
})

test_that("every contrast table holds every feature, in order", {
  res <- sa_multi_group_fixture()

  for (nm in names(res$pairwise)) {
    for (ct in names(res$pairwise[[nm]])) {
      tbl <- res$pairwise[[nm]][[ct]]
      expect_identical(tbl$features, res$features, info = paste(nm, ct))
      expect_identical(unique(tbl$contrast), ct)
      expect_true(all(sa_pairwise_table_columns() %in% names(tbl)),
                  info = paste(nm, ct))
    }
  }
})

test_that("the numbers are the post-hoc numbers, rearranged", {
  res <- sa_multi_group_fixture()
  long <- res$posthoc$anova_test

  for (ct in names(res$pairwise$anova_test)) {
    wide <- res$pairwise$anova_test[[ct]]
    rows <- long[long$contrast == ct, , drop = FALSE]
    at <- match(rows$features, wide$features)
    expect_equal(wide$estimate[at], rows$estimate, info = ct)
    expect_equal(wide$pval[at], rows$pval, info = ct)
    expect_equal(wide$pval_adj[at], rows$pval_adj, info = ct)
    expect_equal(wide$lower_conf[at], rows$lower_conf, info = ct)
  }
})

test_that("a feature that never entered the pairwise stage is NA, not absent", {
  # posthoc_alpha low enough that only the loudest features qualify, so the
  # difference between "absent" and "present with NA" is actually exercised.
  res <- sa_multi_group_fixture(posthoc_alpha = 1e-30)
  asked <- unique(res$posthoc$anova_test$features)
  expect_lt(length(asked), length(res$features))

  tbl <- res$pairwise$anova_test[[1]]
  never <- !tbl$features %in% asked
  expect_true(any(never))
  expect_true(all(is.na(tbl$pval[never])))
  expect_true(all(is.na(tbl$estimate[never])))
  # The ratio is a property of the data, not of a test having been run.
  expect_true(all(is.finite(tbl$log2fc[never])))
})

test_that("log2fc and estimate read in the same direction", {
  res <- sa_multi_group_fixture()

  for (nm in names(res$pairwise)) {
    for (tbl in res$pairwise[[nm]]) {
      both <- is.finite(tbl$log2fc) & is.finite(tbl$estimate) &
        tbl$estimate != 0
      expect_identical(sign(tbl$log2fc[both]), sign(tbl$estimate[both]),
                       info = nm)
    }
  }
})

test_that("fold_change divides group1 by group2", {
  res <- sa_multi_group_fixture()
  tbl <- res$pairwise$anova_test[["virginica - setosa"]]
  row <- tbl[tbl$features == "Petal.Length", ]

  expect_equal(row$fold_change,
               mean(iris$Petal.Length[iris$Species == "virginica"]) /
                 mean(iris$Petal.Length[iris$Species == "setosa"]))
  expect_equal(row$log2fc, log2(row$fold_change))
})

test_that("a repeated design gets the same treatment", {
  res <- sa_repeated_fixture()
  k <- length(res$design$group_lv)
  expect_length(res$pairwise$anova_test, choose(k, 2))
  expect_identical(res$pairwise$anova_test[[1]]$features, res$features)
})

test_that("turning the pairwise stage off removes the slot", {
  res <- sa_multi_group_fixture(posthoc = FALSE)
  expect_false("pairwise" %in% names(res))
  expect_false("posthoc" %in% names(res))
  # Reading it anyway is still safe, which is what every consumer that does not
  # know how the comparison was called depends on.
  expect_null(res$pairwise[["anova_test"]])
})

test_that("the effect table is unchanged by the shared centre calculation", {
  # sa_group_centers() moved the centres out of sa_multi_fold_change(), so the
  # extreme level and its ratio have to come back exactly as they did.
  res <- sa_multi_group_fixture()
  row <- res$effect[res$effect$features == "Petal.Length", ]
  centers <- vapply(levels(iris$Species), function(lv) {
    mean(iris$Petal.Length[iris$Species == lv])
  }, numeric(1))

  expect_identical(row$extreme_level, "virginica")
  expect_equal(row$ref_center, unname(centers["setosa"]))
  expect_equal(row$fold_change, unname(centers["virginica"] / centers["setosa"]))
})

test_that("a feature whose centres cannot be taken still fails as one NA row", {
  # The geometric mean is undefined at or below zero, and the failure has to
  # arrive as an NA row plus one warning rather than aborting the run.
  d <- iris
  d$Sepal.Width[1] <- -1
  expect_warning(
    res <- compare_multiple_groups(d, c("Sepal.Width", "Petal.Length"),
                                   d$Species, levels(d$Species),
                                   fc_mean = "geom", diagnose = FALSE),
    "fold change could not be computed"
  )
  bad <- res$effect[res$effect$features == "Sepal.Width", ]
  expect_true(is.na(bad$log2fc))
  expect_true(is.na(bad$ref_center))
  # The tests themselves are untouched: only the ratio needs positive centres.
  expect_false(is.na(res$tests$anova_test$pval[
    res$tests$anova_test$features == "Sepal.Width"]))
  expect_true(all(is.na(
    res$pairwise$anova_test[[1]]$log2fc[
      res$pairwise$anova_test[[1]]$features == "Sepal.Width"])))
})


# ---------------------------------------------------------------------------
# estimate_significance(by = "contrast")

test_that("by = 'contrast' returns one verdict table per contrast", {
  res <- sa_multi_group_fixture()
  sig <- estimate_significance(res, test = "anova_test", by = "contrast",
                               log2fc_cutoff = 0.1)

  expect_identical(sig$analysis_type, res$analysis)
  expect_identical(names(sig$significance), names(res$pairwise$anova_test))
  for (one in sig$significance) {
    expect_identical(names(one),
                     c("features", "log2fc", "pvalue", "adj_pvalue",
                       "is_signif"))
    expect_identical(one$features, res$features)
  }
})

test_that("a contrast table describes which contrast it is", {
  res <- sa_multi_group_fixture()
  sig <- estimate_significance(res, test = "anova_test", by = "contrast")
  one <- sig$significance[["virginica - setosa"]]

  expect_identical(attr(one, "contrast"), "virginica - setosa")
  expect_identical(attr(one, "group1"), "virginica")
  expect_identical(attr(one, "group2"), "setosa")
  # The cutoffs travel too, so a single element can be plotted on its own.
  expect_identical(attr(one, "pval_cutoff"), 0.05)
  expect_identical(attr(one, "test"), "anova_test")
})

test_that("the contrast verdict uses that contrast's own two axes", {
  res <- sa_multi_group_fixture()
  sig <- estimate_significance(res, test = "anova_test", by = "contrast",
                               log2fc_cutoff = 0.1)
  ct <- "virginica - setosa"
  tbl <- res$pairwise$anova_test[[ct]]
  one <- sig$significance[[ct]]

  expect_equal(one$log2fc, tbl$log2fc)
  expect_equal(one$pvalue, tbl$pval)
  expect_equal(one$adj_pvalue, tbl$pval_adj)
  expect_identical(attr(one, "adj_type"), res$parameters$posthoc_p_adjust)
})

test_that("a contrast verdict agrees in direction with the omnibus verdict", {
  # Both readings feed the same volcano plot, so a feature drawn as raised under
  # one must not be drawn as lowered under the other.
  res <- sa_multi_group_fixture()
  omnibus <- estimate_significance(res, test = "anova_test")$significance
  by_pair <- estimate_significance(res, test = "anova_test",
                                   by = "contrast")$significance
  reference <- res$design$group_lv[1]

  for (i in seq_along(res$features)) {
    ct <- paste(res$effect$extreme_level[i], reference, sep = " - ")
    one <- by_pair[[ct]]
    expect_false(is.null(one), info = ct)
    expect_identical(sign(one$log2fc[i]), sign(omnibus$log2fc[i]),
                     info = res$features[i])
  }
})

test_that("naming an adjustment re-adjusts along the feature axis", {
  res <- sa_multi_group_fixture()
  sig <- estimate_significance(res, test = "anova_test", by = "contrast",
                               adj_type = "bonferroni")
  ct <- "virginica - setosa"
  raw <- res$pairwise$anova_test[[ct]]$pval
  one <- sig$significance[[ct]]

  expect_equal(one$adj_pvalue, stats::p.adjust(raw, "bonferroni"))
  expect_identical(attr(one, "adj_type"), "bonferroni")
})

test_that("by = 'omnibus' is what the default already did", {
  res <- sa_multi_group_fixture()
  expect_identical(estimate_significance(res, by = "omnibus"),
                   estimate_significance(res))
})

test_that("by = 'contrast' says so when there is no pairwise stage", {
  expect_error(
    estimate_significance(sa_two_group_fixture(), test = "t_test",
                          by = "contrast"),
    "needs a pairwise stage"
  )
  expect_error(
    estimate_significance(sa_multi_group_fixture(posthoc = FALSE),
                          by = "contrast"),
    "posthoc = TRUE"
  )
})

test_that("by must name one of the two readings", {
  expect_error(estimate_significance(sa_multi_group_fixture(), by = "pairwise"))
})
