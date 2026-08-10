# One rule fixes direction everywhere: `group_lv[1]` is the reference, and every
# difference and ratio is taken against it. The tests below pin that rule on the
# two-group side, where it used to run the other way, and pin that the two-group
# and multi-group readings of the same data agree on which way a feature moved.
#
# `control_label` is the second way of naming that reference, and the two
# families read it differently on purpose: it re-points a `group_lv`, which
# carries the display order of every level besides, and it must agree with an
# `outcome_lv`, which holds the two classes and nothing else. Both readings are
# pinned here, together with the fact that they end up pointing the same way.

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


test_that("`control_label` re-points a two-group comparison", {
  d <- sa_two_species()
  pointed <- compare_two_groups(d, sa_feats(), d$Species,
                                c("versicolor", "virginica"),
                                control_label = "virginica", diagnose = FALSE)
  written_out <- compare_two_groups(d, sa_feats(), d$Species,
                                    c("virginica", "versicolor"),
                                    diagnose = FALSE)

  # Naming the reference and writing the levels in that order are two ways of
  # saying the one thing, so nothing but the timestamp separates the results.
  expect_identical(pointed$design$group_lv, c("virginica", "versicolor"))
  expect_equal(pointed$effect, written_out$effect)
  expect_equal(pointed$tests, written_out$tests)

  # The default names the level `group_lv` already put first, so it is a
  # no-operation and every existing call keeps its direction.
  plain <- compare_two_groups(d, sa_feats(), d$Species,
                              c("versicolor", "virginica"), diagnose = FALSE)
  defaulted <- compare_two_groups(d, sa_feats(), d$Species,
                                  c("versicolor", "virginica"),
                                  control_label = "versicolor",
                                  diagnose = FALSE)
  expect_equal(defaulted$effect, plain$effect)
  expect_equal(defaulted$tests, plain$tests)

  expect_error(
    compare_two_groups(d, sa_feats(), d$Species,
                       c("versicolor", "virginica"), control_label = "setosa",
                       diagnose = FALSE),
    "names a level `group_lv` does not hold"
  )
})


test_that("`control_label` carries the post-hoc contrasts with it", {
  lv <- levels(iris$Species)
  pointed <- compare_multiple_groups(iris, sa_feats(), iris$Species, lv,
                                     control_label = "versicolor",
                                     diagnose = FALSE)
  written_out <- compare_multiple_groups(iris, sa_feats(), iris$Species,
                                         c("versicolor", "setosa",
                                           "virginica"),
                                         diagnose = FALSE)

  # The move puts the named level first and leaves the rest in the order they
  # arrived, which is what the two calls have in common.
  expect_identical(pointed$design$group_lv,
                   c("versicolor", "setosa", "virginica"))
  expect_equal(pointed$effect, written_out$effect)
  expect_equal(pointed$posthoc, written_out$posthoc)

  # The reference is the denominator of the fold change and the subtracted side
  # of every contrast at once, so the contrast names move with it.
  reference <- sa_multi_group_fixture()
  expect_true("versicolor - setosa" %in% reference$posthoc$anova_test$contrast)
  expect_true("setosa - versicolor" %in% pointed$posthoc$anova_test$contrast)

  moved <- pointed$posthoc$anova_test
  moved <- moved[moved$contrast == "setosa - versicolor", ]
  stayed <- reference$posthoc$anova_test
  stayed <- stayed[stayed$contrast == "versicolor - setosa", ]
  stayed <- stayed[match(moved$features, stayed$features), ]
  expect_equal(moved$estimate, -stayed$estimate)
})


# Fifty of one species against thirty of the other, so that `n_events` is a
# different number depending on which class the model is about.
sa_unbalanced_species <- function() {
  d <- rbind(iris[iris$Species == "versicolor", ][1:50, ],
             iris[iris$Species == "virginica", ][1:30, ])
  d[c("Species", "Petal.Length", "Petal.Width")]
}

sa_direction_models <- function(d) {
  list(
    logistic = function(...) {
      fit_logistic_regression(d, outcome = "Species", cv = FALSE, ...)
    },
    elastic_net = function(...) {
      fit_elastic_net(d, outcome = "Species", penalty = "lasso",
                      lambda = 0.01, cv = FALSE, ...)
    },
    rf = function(...) {
      fit_rf(d, outcome = "Species", ntree = 100, cv = FALSE, seed = 1, ...)
    },
    svm = function(...) {
      fit_svm(d, outcome = "Species", C = 1, cv = FALSE, seed = 1, ...)
    }
  )
}


test_that("`control_label` alone fixes the reference of every model", {
  d <- sa_unbalanced_species()

  for (nm in names(sa_direction_models(d))) {
    fit_it <- sa_direction_models(d)[[nm]]
    versi <- fit_it(control_label = "versicolor")
    virgi <- fit_it(control_label = "virginica")

    expect_identical(versi$design$outcome_lv, c("versicolor", "virginica"),
                     info = nm)
    expect_identical(virgi$design$outcome_lv, c("virginica", "versicolor"),
                     info = nm)
    # `n_events` counts the class the model is about, so it follows the
    # reference to the other side rather than staying with one label.
    expect_identical(versi$design$n_events, 30L, info = nm)
    expect_identical(virgi$design$n_events, 50L, info = nm)
  }
})


test_that("a `control_label` that disagrees with `outcome_lv` is refused", {
  d <- sa_unbalanced_species()

  for (nm in names(sa_direction_models(d))) {
    fit_it <- sa_direction_models(d)[[nm]]
    expect_error(
      fit_it(outcome_lv = c("versicolor", "virginica"),
             control_label = "virginica"),
      "disagree",
      info = nm
    )
    # Naming both and agreeing is the default, so it is no conflict at all.
    both <- fit_it(outcome_lv = c("versicolor", "virginica"),
                   control_label = "versicolor")
    expect_identical(both$design$outcome_lv, c("versicolor", "virginica"),
                     info = nm)
    expect_error(
      fit_it(control_label = "setosa"),
      "names a class `outcome` does not hold",
      info = nm
    )
  }
})


test_that("`control_label` makes a column of zeroes and ones a classification", {
  d <- mtcars[c("am", "wt", "hp")]

  # Without it the column is two numbers and the announcement says what it
  # would take to make it two classes.
  expect_message(fit_rf(d, outcome = "am", ntree = 100, cv = FALSE, seed = 1),
                 "Pass `outcome_lv` or `control_label`")

  clf <- fit_rf(d, outcome = "am", control_label = "1", ntree = 100,
                cv = FALSE, seed = 1)
  expect_identical(clf$design$outcome_lv, c("1", "0"))
  expect_identical(clf$design$n_events, 19L)

  enet <- fit_elastic_net(d, outcome = "am", control_label = "1",
                          penalty = "lasso", lambda = 0.01, cv = FALSE)
  expect_identical(enet$design$outcome_lv, c("1", "0"))

  svm <- fit_svm(d, outcome = "am", control_label = "1", C = 1, cv = FALSE,
                 seed = 1)
  expect_identical(svm$design$outcome_lv, c("1", "0"))
})


test_that("a comparison and a model held at the same control agree", {
  d <- sa_two_species()
  petal <- d[c("Species", "Petal.Length")]

  for (control in c("versicolor", "virginica")) {
    cmp <- compare_two_groups(d, "Petal.Length", d$Species,
                              c("versicolor", "virginica"),
                              control_label = control, diagnose = FALSE)
    fit <- fit_logistic_regression(petal, outcome = "Species",
                                   control_label = control, cv = FALSE)
    at <- fit$coefficients$terms == "Petal.Length"

    # Both are statements about the class that is not the control, so a fold
    # change above 1 and an odds ratio above 1 are the same finding twice.
    expect_equal(sign(cmp$effect$log2fc),
                 sign(log(fit$coefficients$odds_ratio[at])))
  }
})
