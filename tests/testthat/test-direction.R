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


# A crossed design has one reference per factor rather than one in total, so the
# argument holds one name per factor it points and the reference cell is where
# all of them land.
sa_fact_factors <- function() list(wool = "wool", tension = "tension")

sa_fact_pointed <- function(...) {
  compare_factorial_groups(warpbreaks, "breaks", sa_fact_factors(),
                           posthoc = FALSE, diagnose = FALSE, ...)
}

test_that("`control_label` re-points one factor of a crossed design", {
  lv <- list(wool = c("A", "B"), tension = c("L", "M", "H"))
  pointed <- sa_fact_pointed(factor_lv = lv,
                             control_label = list(tension = "M"))
  written_out <- sa_fact_pointed(
    factor_lv = list(wool = c("A", "B"), tension = c("M", "L", "H"))
  )

  # Naming the level and retyping every level in that order are two ways of
  # saying the one thing, as they are for a single factor.
  expect_identical(pointed$design$factor_lv, written_out$design$factor_lv)
  expect_equal(pointed$effect, written_out$effect)
  expect_equal(pointed$terms, written_out$terms)

  # The factor it says nothing about is left as it arrived, which is what makes
  # pointing one factor of several a sentence rather than a rewrite of all.
  expect_identical(pointed$design$factor_lv$wool, c("A", "B"))
  expect_identical(pointed$design$group_lv[1], "A.M")

  # The two shapes of the argument are the same argument.
  expect_equal(sa_fact_pointed(factor_lv = lv,
                               control_label = c(tension = "M"))$effect,
               pointed$effect)

  both <- sa_fact_pointed(factor_lv = lv,
                          control_label = list(wool = "B", tension = "H"))
  expect_identical(both$design$group_lv[1], "B.H")

  # Naming the level that is already first is the no-operation the default is.
  plain <- sa_fact_pointed(factor_lv = lv)
  named <- sa_fact_pointed(factor_lv = lv,
                           control_label = list(wool = "A", tension = "L"))
  expect_equal(plain$effect, named$effect)
  expect_identical(plain$design$factor_lv, named$design$factor_lv)
})

test_that("re-pointing a factorial reference moves the ratios and not the tests", {
  plain <- sa_fact_pointed()
  pointed <- sa_fact_pointed(control_label = list(tension = "M"))

  # A sum-to-zero coding spans the same space whichever level leads a factor, so
  # the design is being read against a different cell rather than refitted.
  expect_equal(plain$tests$anova_test$f_stat, pointed$tests$anova_test$f_stat)
  expect_equal(plain$tests$anova_test$pval, pointed$tests$anova_test$pval)
  expect_equal(plain$terms$f_stat, pointed$terms$f_stat)
  expect_equal(plain$terms$pval, pointed$terms$pval)
  # A term component is a deviation from the rest of the model, so it does not
  # depend on which cell the ratios are taken against either.
  expect_equal(plain$terms$log2_effect, pointed$terms$log2_effect)
  # Term order follows the order the factors were declared in, not their levels.
  expect_identical(plain$terms$terms, pointed$terms$terms)

  # What does move is everything that is read against the reference cell.
  expect_false(isTRUE(all.equal(plain$effect$ref_center,
                                pointed$effect$ref_center)))
  expect_false(isTRUE(all.equal(plain$effect$log2fc, pointed$effect$log2fc)))
  expect_true(all(pointed$effect$extreme_cell %in% pointed$design$group_lv))
})

test_that("`control_label` carries the factorial contrasts with it", {
  plain <- compare_factorial_groups(warpbreaks, "breaks", sa_fact_factors(),
                                    diagnose = FALSE)
  pointed <- compare_factorial_groups(warpbreaks, "breaks", sa_fact_factors(),
                                      control_label = list(tension = "M"),
                                      diagnose = FALSE)
  expect_true("M - H" %in% plain$posthoc$anova_test$contrast)
  expect_true("H - M" %in% pointed$posthoc$anova_test$contrast)

  stayed <- subset(plain$posthoc$anova_test,
                   contrast == "M - H" & is.na(stratum))
  moved <- subset(pointed$posthoc$anova_test,
                  contrast == "H - M" & is.na(stratum))
  expect_equal(stayed$estimate, -moved$estimate)
  # A Tukey contrast has no direction of its own, so only the sign turns.
  expect_equal(stayed$pval, moved$pval)
})

test_that("`control_label` states the reference a sort did not", {
  # Without `factor_lv` the levels arrive in sorted order, which puts H first
  # although the design calls L the baseline. This is the case a single factor
  # does not have, since `group_lv` is required there.
  sorted <- sa_fact_pointed()
  expect_identical(sorted$design$factor_lv$tension, c("H", "L", "M"))
  expect_identical(sorted$design$group_lv[1], "A.H")

  stated <- sa_fact_pointed(control_label = list(tension = "L"))
  expect_identical(stated$design$factor_lv$tension, c("L", "H", "M"))
  expect_identical(stated$design$group_lv[1], "A.L")
})

test_that("a `control_label` that does not fit the design is refused", {
  expect_error(sa_fact_pointed(control_label = list(dose = "low")),
               "names factor\\(s\\) the design does not hold")
  # The message names where the levels came from, which is not `factor_lv` when
  # the caller left it out.
  expect_error(sa_fact_pointed(control_label = list(tension = "X")),
               "names a level `factors\\$tension` does not hold")
  expect_error(
    sa_fact_pointed(factor_lv = list(wool = c("A", "B"),
                                     tension = c("L", "M", "H")),
                    control_label = list(tension = "X")),
    "names a level `factor_lv\\$tension` does not hold"
  )
  expect_error(sa_fact_pointed(control_label = list(tension = c("L", "M"))),
               "one level name per factor")
  expect_error(sa_fact_pointed(control_label = "L"),
               "named list or named character vector")
  expect_error(sa_fact_pointed(control_label = c(tension = "L", tension = "M")),
               "named list or named character vector")
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
