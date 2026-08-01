# Statistics are checked against the base R equivalent wherever one exists, and
# against a hand-written formula where it does not. No external comparison
# package is used here; those live in Test/ and are not a package dependency.

test_that("the omnibus statistics match their base R equivalents", {
  res <- sa_multi_group_fixture()
  y <- iris$Petal.Length
  g <- iris$Species

  anova_row <- res$tests$anova_test[res$features == "Petal.Length", ]
  base_anova <- stats::oneway.test(y ~ g, var.equal = TRUE)
  expect_equal(anova_row$f_stat, unname(base_anova$statistic))
  expect_equal(anova_row$pval, base_anova$p.value)

  welch_row <- res$tests$welch_test[res$features == "Petal.Length", ]
  base_welch <- stats::oneway.test(y ~ g, var.equal = FALSE)
  expect_equal(welch_row$f_stat, unname(base_welch$statistic))
  expect_equal(welch_row$df2, unname(base_welch$parameter[2]))
  expect_equal(welch_row$pval, base_welch$p.value)

  kruskal_row <- res$tests$kruskal_test[res$features == "Petal.Length", ]
  base_kruskal <- stats::kruskal.test(y ~ g)
  expect_equal(kruskal_row$h_stat, unname(base_kruskal$statistic))
  expect_equal(kruskal_row$pval, base_kruskal$p.value)
})

test_that("eta squared is the between-group share of the sum of squares", {
  res <- sa_multi_group_fixture()
  fit <- stats::aov(Petal.Length ~ Species, data = iris)
  ss <- summary(fit)[[1]][["Sum Sq"]]
  row <- res$tests$anova_test[res$features == "Petal.Length", ]
  expect_equal(row$eta_sq, ss[1] / sum(ss))
  # Omega squared corrects for chance-explained variance, so it is the smaller
  # of the two whenever the residual mean square is above zero.
  expect_lt(row$omega_sq, row$eta_sq)
})

test_that("the trimmed mean ANOVA reduces to Welch when nothing is trimmed", {
  # tr = 0 leaves the trimmed mean equal to the mean and the winsorised variance
  # equal to the variance, so the two constructions have to agree exactly.
  res <- sa_multi_group_fixture(tr = 0)
  expect_equal(res$tests$robust_test$f_stat, res$tests$welch_test$f_stat)
  expect_equal(res$tests$robust_test$df2, res$tests$welch_test$df2)
})

test_that("the trimmed mean ANOVA ignores a value already in the tail", {
  d <- iris
  contaminated <- d
  # The largest setosa sepal, blown up. Trimming removes a proportion rather
  # than everything past a threshold, so inflating a value that was already
  # inside the trimmed tail leaves the sort order, the trimmed mean and the
  # winsorised variance all exactly as they were.
  setosa_rows <- which(d$Species == "setosa")
  contaminated$Sepal.Length[setosa_rows[which.max(d$Sepal.Length[setosa_rows])]] <- 500

  clean <- compare_multiple_groups(d, "Sepal.Length", d$Species,
                                   levels(d$Species), diagnose = FALSE,
                                   posthoc = FALSE)
  dirty <- compare_multiple_groups(contaminated, "Sepal.Length",
                                   contaminated$Species, levels(d$Species),
                                   diagnose = FALSE, posthoc = FALSE)

  expect_equal(dirty$tests$robust_test$f_stat, clean$tests$robust_test$f_stat)
  expect_equal(dirty$tests$robust_test$pval, clean$tests$robust_test$pval)
  # The parametric test, judging the same data, is thrown right off.
  expect_gt(abs(dirty$tests$anova_test$f_stat - clean$tests$anova_test$f_stat) /
              clean$tests$anova_test$f_stat, 0.9)
})

test_that("Tukey contrasts match stats::TukeyHSD up to the pair direction", {
  res <- sa_multi_group_fixture()
  reference <- stats::TukeyHSD(stats::aov(Petal.Width ~ Species,
                                          data = iris))$Species
  tbl <- res$posthoc$anova_test
  tbl <- tbl[tbl$features == "Petal.Width", ]

  # TukeyHSD labels a pair "later-earlier"; the package orders it by group_lv,
  # so every estimate and both bounds come out negated and swapped.
  expect_equal(tbl$estimate, unname(-reference[, "diff"]))
  expect_equal(tbl$lower_conf, unname(-reference[, "upr"]))
  expect_equal(tbl$upper_conf, unname(-reference[, "lwr"]))
  expect_equal(tbl$pval, unname(reference[, "p adj"]), tolerance = 1e-7)
})

test_that("Tukey and Games-Howell p-values are left unadjusted", {
  # Both control the error rate over the same set of contrasts through the
  # studentised range, so adjusting again would correct twice.
  res <- sa_multi_group_fixture(posthoc_p_adjust = "bonferroni")
  expect_identical(res$posthoc$anova_test$pval_adj,
                   res$posthoc$anova_test$pval)
  expect_identical(res$posthoc$welch_test$pval_adj,
                   res$posthoc$welch_test$pval)
  # Dunn's test is not family-wise, so it is adjusted.
  dunn <- res$posthoc$kruskal_test
  expect_true(all(dunn$pval_adj >= dunn$pval))
  expect_equal(dunn$pval_adj[1:3],
               stats::p.adjust(dunn$pval[1:3], "bonferroni"))
})

test_that("the post-hoc adjustment is applied within features, not across", {
  res <- sa_multi_group_fixture(posthoc_p_adjust = "bonferroni")
  tbl <- res$posthoc$kruskal_test
  for (f in unique(tbl$features)) {
    block <- tbl[tbl$features == f, ]
    expect_equal(block$pval_adj,
                 stats::p.adjust(block$pval, "bonferroni"), info = f)
  }
})

test_that("only features clearing posthoc_alpha enter the pairwise stage", {
  strict <- sa_multi_group_fixture(posthoc_alpha = 1e-90)
  qualified <- unique(strict$posthoc$anova_test$features)
  expect_true(length(qualified) < length(strict$features))
  expect_identical(strict$parameters$n_posthoc[["anova_test"]],
                   length(qualified))
  # A feature that did not qualify is absent, not present with NA: "never
  # asked" and "asked and unanswerable" are different facts.
  expect_false("Sepal.Width" %in% qualified)
  expect_equal(nrow(strict$posthoc$anova_test), 3L * length(qualified))
})

test_that("posthoc = FALSE runs no pairwise stage at all", {
  res <- sa_multi_group_fixture(posthoc = FALSE)
  expect_identical(res$posthoc, structure(list(), names = character(0)))
  expect_true(all(res$parameters$n_posthoc == 0L))
})

test_that("the effect table points at the level furthest from the reference", {
  res <- sa_multi_group_fixture()
  row <- res$effect[res$effect$features == "Petal.Length", ]
  expect_identical(row$extreme_level, "virginica")
  expect_equal(row$ref_center, mean(iris$Petal.Length[iris$Species == "setosa"]))
  expect_equal(row$extreme_center,
               mean(iris$Petal.Length[iris$Species == "virginica"]))
  expect_equal(row$fold_change, row$extreme_center / row$ref_center)
  expect_equal(row$log2fc, log2(row$fold_change))
  expect_gt(row$log2fc, 0)

  # Sepal width runs the other way: versicolor is furthest below setosa.
  narrow <- res$effect[res$effect$features == "Sepal.Width", ]
  expect_identical(narrow$extreme_level, "versicolor")
  expect_lt(narrow$log2fc, 0)
})

test_that("repeated measures ANOVA matches aov with an Error stratum", {
  d <- sa_repeated_frame()
  res <- sa_repeated_fixture()
  fit <- summary(stats::aov(value ~ cond + Error(subj / cond), data = d))
  inner <- fit[["Error: subj:cond"]][[1]]

  expect_equal(res$tests$anova_test$f_stat, inner[["F value"]][1])
  expect_equal(res$tests$anova_test$df1, inner[["Df"]][1])
  expect_equal(res$tests$anova_test$df2, inner[["Df"]][2])
  expect_equal(res$tests$anova_test$pval, inner[["Pr(>F)"]][1])
})

test_that("Mauchly's W matches stats::mauchly.test", {
  d <- sa_repeated_frame()
  wide <- stats::reshape(d, idvar = "subj", timevar = "cond",
                         direction = "wide")
  mat <- as.matrix(wide[, -1])
  res <- sa_repeated_fixture()

  reference <- stats::mauchly.test(stats::lm(mat ~ 1), X = ~1)
  expect_equal(res$tests$anova_test$mauchly_w, unname(reference$statistic))
  # The p-value is deliberately not compared: stats evaluates one term of the
  # expansion with the number of conditions where the published formula uses
  # the contrast rank, so the two differ in the fourth decimal place.
  expect_equal(res$tests$anova_test$mauchly_pval, reference$p.value,
               tolerance = 1e-3)
})

test_that("the sphericity epsilons stay inside their bounds", {
  res <- sa_repeated_fixture()
  row <- res$tests$anova_test
  lower_bound <- 1 / (row$n_groups - 1)
  for (eps in c(row$gg_eps, row$hf_eps)) {
    expect_gte(eps, lower_bound)
    expect_lte(eps, 1)
  }
  # Huynh-Feldt corrects less than Greenhouse-Geisser, so its epsilon is the
  # larger one and it keeps more degrees of freedom.
  expect_gte(row$hf_eps, row$gg_eps)
  # Any correction shrinks the degrees of freedom, so both corrected p-values
  # sit above the uncorrected one, and the milder correction sits lower.
  expect_gte(row$pval_gg, row$pval)
  expect_gte(row$pval_hf, row$pval)
  expect_gte(row$pval_gg, row$pval_hf)
})

test_that("the Friedman statistic matches stats::friedman.test", {
  d <- sa_repeated_frame()
  wide <- stats::reshape(d, idvar = "subj", timevar = "cond",
                         direction = "wide")
  mat <- as.matrix(wide[, -1])
  res <- sa_repeated_fixture()
  reference <- stats::friedman.test(mat)

  expect_equal(res$tests$kruskal_test$chi_sq, unname(reference$statistic))
  expect_equal(res$tests$kruskal_test$pval, reference$p.value)
  expect_equal(res$tests$kruskal_test$kendalls_w,
               unname(reference$statistic) /
                 (nrow(mat) * (ncol(mat) - 1)))
})

test_that("pairwise paired t contrasts match stats::t.test", {
  d <- sa_repeated_frame()
  wide <- stats::reshape(d, idvar = "subj", timevar = "cond",
                         direction = "wide")
  res <- sa_repeated_fixture()
  tbl <- res$posthoc$anova_test
  row <- tbl[tbl$group1 == "t1" & tbl$group2 == "t3", ]
  reference <- stats::t.test(wide[["value.t1"]], wide[["value.t3"]],
                             paired = TRUE)

  expect_equal(row$statistic, unname(reference$statistic))
  expect_equal(row$pval, reference$p.value)
  expect_equal(row$lower_conf, reference$conf.int[1])
})

test_that("repeated conditions need an id and refuse row order pairing", {
  d <- sa_repeated_frame()
  expect_error(
    compare_multiple_groups(d["value"], "value", d$cond, sort(unique(d$cond)),
                            paired = TRUE),
    "needs `id`"
  )
})

test_that("subjects missing a condition are dropped whole and counted", {
  d <- sa_repeated_frame()
  incomplete <- d[!(d$subj == "s1" & d$cond == "t2"), ]
  expect_message(
    res <- compare_multiple_groups(incomplete["value"], "value",
                                   incomplete$cond, sort(unique(d$cond)),
                                   id = incomplete$subj, paired = TRUE,
                                   diagnose = FALSE),
    "missing at least one condition"
  )
  expect_identical(res$design$unmatched_ids, "s1")
  expect_equal(res$tests$anova_test$n_used, 13)
})

test_that("a repeated id within one condition is an error, not a guess", {
  d <- sa_repeated_frame()
  d$subj[d$cond == "t1"][1] <- d$subj[d$cond == "t1"][2]
  expect_error(
    compare_multiple_groups(d["value"], "value", d$cond, sort(unique(d$cond)),
                            id = d$subj, paired = TRUE),
    "must be unique within each condition"
  )
})

test_that("fewer than three levels is rejected", {
  d <- iris[iris$Species != "setosa", ]
  expect_error(
    compare_multiple_groups(d, "Sepal.Length", d$Species,
                            c("versicolor", "virginica")),
    "at least 3"
  )
})

test_that("an untestable feature yields an NA row and one warning per test", {
  d <- iris
  d$flat <- 1
  # The convention is one warning per test table naming every feature it could
  # not test, not one per feature: four omnibus tests, four warnings.
  warnings <- capture_warnings(
    res <- compare_multiple_groups(d, c("Sepal.Length", "flat"), d$Species,
                                   levels(d$Species), diagnose = FALSE,
                                   posthoc = FALSE)
  )
  expect_length(warnings, length(res$tests))
  expect_true(all(grepl("flat", warnings)))
  # Including the rank-based test: a constant feature leaves the tie correction
  # undefined, and that has to be reported rather than returned as a NaN.
  expect_true(all(vapply(res$tests,
                         function(t) is.na(t$pval[res$features == "flat"]),
                         logical(1))))
  expect_false(is.na(res$tests$anova_test$pval[1]))
})

test_that("an omnibus table reports no interval", {
  res <- sa_multi_group_fixture()
  for (tbl in res$tests) {
    expect_true(all(is.na(tbl$lower_conf)))
    expect_true(all(is.na(tbl$upper_conf)))
  }
})
