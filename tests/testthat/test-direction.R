# One rule fixes direction everywhere: `group_lv[1]` is the reference, and every
# difference and ratio is taken against it. The tests below pin that rule on the
# two-group side, where it used to run the other way, and pin that the two-group
# and multi-group readings of the same data agree on which way a feature moved.

sa_two_species <- function() {
  iris[iris$Species != "setosa", ]
}

test_that("differences and ratios are taken against the first level", {
  d <- sa_two_species()
  res <- compare_two_groups(d, "Petal.Length", d$Species,
                            c("versicolor", "virginica"), diagnose = FALSE)
  ref <- d$Petal.Length[d$Species == "versicolor"]
  against <- d$Petal.Length[d$Species == "virginica"]

  expect_equal(res$tests$t_test$mean_diff, mean(against) - mean(ref))
  expect_equal(res$effect$x_center, mean(against))
  expect_equal(res$effect$y_center, mean(ref))
  expect_equal(res$effect$fold_change, mean(against) / mean(ref))
})

test_that("swapping the two levels flips every direction at once", {
  d <- sa_two_species()
  forward <- compare_two_groups(d, sa_feats(), d$Species,
                                c("versicolor", "virginica"),
                                diagnose = FALSE)
  reversed <- compare_two_groups(d, sa_feats(), d$Species,
                                 c("virginica", "versicolor"),
                                 diagnose = FALSE)

  expect_equal(reversed$tests$t_test$mean_diff, -forward$tests$t_test$mean_diff)
  # The Hodges-Lehmann shift is found by root-finding rather than in closed
  # form, so the two directions agree to the solver's tolerance, not to the
  # last bit.
  expect_equal(reversed$tests$wilcox_test$hl_shift,
               -forward$tests$wilcox_test$hl_shift, tolerance = 1e-4)
  expect_equal(reversed$effect$log2fc, -forward$effect$log2fc)
  # A probability rather than a difference, so it mirrors around 0.5.
  expect_equal(reversed$tests$robust_test$relative_effect,
               1 - forward$tests$robust_test$relative_effect)
  # The two-sided p-value is the same question asked from either side.
  expect_equal(reversed$tests$t_test$pval, forward$tests$t_test$pval)
})

test_that("alternative = greater asks about the second level", {
  d <- sa_two_species()
  greater <- compare_two_groups(d, "Petal.Length", d$Species,
                                c("versicolor", "virginica"),
                                alternative = "greater", diagnose = FALSE)
  less <- compare_two_groups(d, "Petal.Length", d$Species,
                             c("versicolor", "virginica"),
                             alternative = "less", diagnose = FALSE)

  # Virginica has the longer petals, so it is the one-sided test naming it as
  # the larger group that finds them.
  expect_lt(greater$tests$t_test$pval, 0.05)
  expect_gt(less$tests$t_test$pval, 0.95)
})

test_that("the two-group and multi-group readings point the same way", {
  d <- sa_two_species()
  two <- compare_two_groups(d, sa_feats(), d$Species,
                            c("versicolor", "virginica"), diagnose = FALSE)
  multi <- sa_multi_group_fixture()

  pair <- multi$posthoc$anova_test
  pair <- pair[pair$contrast == "virginica - versicolor", ]
  pair <- pair[match(two$features, pair$features), ]

  expect_equal(sign(two$effect$log2fc), sign(pair$estimate))
  expect_equal(sign(two$tests$t_test$mean_diff), sign(pair$estimate))
})

test_that("the diagnostics keep the display order rather than the test order", {
  d <- sa_two_species()
  res <- compare_two_groups(d, "Petal.Length", d$Species,
                            c("versicolor", "virginica"), diagnose = TRUE)

  # `x` is virginica here, but the per-level tables are built reference first so
  # that they line up with `group_lv` and with a boxplot of the same input.
  expect_identical(res$diagnostics$normality$group,
                   c("versicolor", "virginica"))
})
