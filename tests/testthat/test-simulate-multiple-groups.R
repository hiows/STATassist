# The value of a simulator is that the answer is known, so what is pinned here
# is mostly that the answer it reports is the answer it planted. Beyond what the
# two-group simulator has to promise, this one has a shape per planted feature
# and a second answer table for the pairwise stage, so the shapes and the
# direction of the contrast answers are pinned too.

# Small enough to be quick, with the size arguments merged rather than passed
# after the defaults, so that a test can override any one of them.
small <- function(...) {
  given <- list(...)
  base <- list(n_feats = 12, n_up = 3, n_down = 3, n_control = 6,
               n_treat = rep(6, 3))
  do.call(simulate_multiple_groups,
          c(given, base[setdiff(names(base), names(given))]))
}

deltas_of <- function(sim, feature) {
  sim$truth_group$delta[sim$truth_group$features == feature]
}


test_that("the returned arguments are the ones compare_multiple_groups takes", {
  sim <- small(seed = 1)

  expect_named(sim, c("args", "truth", "truth_group", "truth_contrast"))
  expect_named(sim$args,
               c("data", "feats", "group", "group_lv", "input_scale"))
  expect_true(all(names(sim$args) %in%
                    names(formals(compare_multiple_groups))))
  expect_identical(sim$args$input_scale, "log2")
  expect_identical(names(sim$args$data), sim$args$feats)
  expect_length(sim$args$group, nrow(sim$args$data))
  # The control is the reference denominator, so it has to lead. This is the
  # opposite arrangement to simulate_two_groups(), where the case leads.
  expect_identical(sim$args$group_lv,
                   c("control", "treat_1", "treat_2", "treat_3"))
})

test_that("a repeated design carries its pairing in the arguments", {
  sim <- small(seed = 1, paired = TRUE)

  expect_named(sim$args, c("data", "feats", "group", "group_lv", "id",
                           "paired", "input_scale"))
  expect_true(all(names(sim$args) %in%
                    names(formals(compare_multiple_groups))))
  expect_true(sim$args$paired)
  # Every subject under every condition, so nothing is dropped for being
  # partial and the within-subject tests see the full rectangle.
  expect_equal(unname(table(sim$args$id)), rep(4L, 6L), ignore_attr = TRUE)
  expect_equal(unname(table(sim$args$group)), rep(6L, 4L), ignore_attr = TRUE)
})

test_that("the answer tables are aligned with the features and the pairs", {
  sim <- small(seed = 1)

  expect_identical(sim$truth$features, sim$args$feats)
  expect_named(sim$truth, c("features", "pattern", "direction",
                            "extreme_level", "extreme_tied", "log2fc",
                            "baseline", "sd_subject"))
  expect_named(sim$truth_group, c("features", "group", "is_ref", "delta",
                                  "center", "sd", "n"))
  expect_named(sim$truth_contrast, c("features", "contrast", "group1",
                                     "group2", "delta", "is_diff"))

  expect_equal(nrow(sim$truth_group), 12L * 4L)
  expect_equal(nrow(sim$truth_contrast), 12L * 6L)
  expect_true(all(sim$truth_group$is_ref == (sim$truth_group$group ==
                                               "control")))
  expect_true(all(sim$truth_group$delta[sim$truth_group$is_ref] == 0))
})

test_that("exactly the requested number of features is planted", {
  sim <- small(n_feats = 40, n_up = 7, n_down = 11, seed = 3)

  expect_equal(sum(sim$truth$direction == "up"), 7L)
  expect_equal(sum(sim$truth$direction == "down"), 11L)
  expect_equal(sum(sim$truth$direction == "none"), 22L)
  expect_equal(sum(sim$truth$pattern == "none"), 22L)
  expect_true(all(sim$truth$log2fc[sim$truth$direction == "up"] > 0))
  expect_true(all(sim$truth$log2fc[sim$truth$direction == "down"] < 0))
})

test_that("the shapes are handed out in the proportions asked for", {
  # The split is by largest remainder rather than by lot, so the counts are a
  # function of the arguments and can be pinned exactly. Nine up features over
  # three equal weights is three each; five is 2 / 2 / 1, the remainder going to
  # the earlier shapes.
  even <- small(n_feats = 30, n_up = 9, n_down = 9, seed = 1)
  expect_equal(unname(table(even$truth$pattern)[c("all", "gradient",
                                                  "single")]),
               c(6L, 6L, 6L), ignore_attr = TRUE)

  odd <- small(n_feats = 30, n_up = 5, n_down = 0, seed = 1)
  expect_equal(sum(odd$truth$pattern == "all"), 2L)
  expect_equal(sum(odd$truth$pattern == "gradient"), 2L)
  expect_equal(sum(odd$truth$pattern == "single"), 1L)

  weighted <- small(n_feats = 30, n_up = 8, n_down = 0, seed = 1,
                    pattern_mix = c(all = 3, gradient = 1, single = 0))
  expect_equal(sum(weighted$truth$pattern == "all"), 6L)
  expect_equal(sum(weighted$truth$pattern == "gradient"), 2L)
  expect_equal(sum(weighted$truth$pattern == "single"), 0L)
})

test_that("each shape moves the treatment groups the way it says it does", {
  sim <- small(n_feats = 30, n_up = 9, n_down = 9, n_treat = rep(6, 4),
               seed = 5)

  for (f in sim$truth$features[sim$truth$pattern == "all"]) {
    d <- deltas_of(sim, f)[-1]
    expect_equal(length(unique(d)), 1L)
    expect_true(d[1] != 0)
  }
  for (f in sim$truth$features[sim$truth$pattern == "gradient"]) {
    d <- abs(deltas_of(sim, f)[-1])
    # Strictly increasing, with the last level carrying the whole effect.
    expect_true(all(diff(d) > 0))
    expect_equal(d[4], max(d))
  }
  for (f in sim$truth$features[sim$truth$pattern == "single"]) {
    d <- deltas_of(sim, f)[-1]
    expect_equal(sum(d != 0), 1L)
  }
})

test_that("an unplanted feature is null in the strict sense", {
  sim <- small(seed = 4)
  # Not "small", exactly zero, and in every group rather than on average.
  # Everything downstream that scores a false positive rate depends on there
  # being no true effect here at all.
  null_feats <- sim$truth$features[sim$truth$direction == "none"]

  expect_true(all(sim$truth$log2fc[sim$truth$direction == "none"] == 0))
  expect_true(all(sim$truth_group$delta[
    sim$truth_group$features %in% null_feats] == 0))
  expect_true(all(sim$truth_contrast$delta[
    sim$truth_contrast$features %in% null_feats] == 0))
  expect_true(!any(sim$truth_contrast$is_diff[
    sim$truth_contrast$features %in% null_feats]))
})

test_that("the planted magnitudes stay inside deg_log2fc", {
  sim <- small(n_feats = 30, n_up = 6, n_down = 6, deg_log2fc = c(1.5, 2),
               seed = 5)
  planted <- abs(sim$truth$log2fc[sim$truth$direction != "none"])

  expect_true(all(planted >= 1.5))
  expect_true(all(planted <= 2))
})

test_that("extreme_level names a level only when one stands out", {
  sim <- small(n_feats = 30, n_up = 9, n_down = 9, seed = 5)

  # "all" moves every treatment group alike, so no level is furthest from the
  # control and the tie is reported rather than broken silently.
  tied <- sim$truth$pattern %in% c("all", "none")
  expect_true(all(sim$truth$extreme_tied[tied]))
  expect_true(!any(sim$truth$extreme_tied[!tied]))
  expect_true(all(is.na(sim$truth$extreme_level[sim$truth$pattern == "none"])))
  expect_true(all(!is.na(sim$truth$extreme_level[sim$truth$pattern != "none"])))

  # Where one does stand out, it is the level whose delta the log2fc reports.
  for (f in sim$truth$features[!tied]) {
    row <- sim$truth[sim$truth$features == f, ]
    lv <- sim$truth_group$group[sim$truth_group$features == f]
    d <- deltas_of(sim, f)
    expect_identical(row$extreme_level, lv[which.max(abs(d))])
    expect_equal(row$log2fc, d[which.max(abs(d))])
  }
})

test_that("truth_contrast follows the post-hoc pair order and direction", {
  sim <- small(seed = 1)
  pairs <- sa_level_pairs(sim$args$group_lv)

  expect_identical(sim$truth_contrast$contrast[1:6], pairs$contrast)
  expect_identical(sim$truth_contrast$group1[1:6], pairs$group1)

  # A post-hoc estimate is `group1 - group2`, and the control is the level being
  # subtracted, so the answer for a feature the treatments raised is positive
  # here just as its log2fc is. Getting this backwards would score every call
  # wrong.
  f <- sim$truth$features[sim$truth$direction == "up"][1]
  con <- sim$truth_contrast[sim$truth_contrast$features == f, ]
  d <- deltas_of(sim, f)
  expect_equal(con$delta, d[pairs$i] - d[pairs$j])
  expect_gt(con$delta[con$contrast == "treat_3 - control"], 0)
  expect_gt(sim$truth$log2fc[sim$truth$features == f], 0)

  expect_identical(con$is_diff, con$delta != 0)
})

test_that("every group gets the size it was given", {
  balanced <- small(seed = 1, n_control = 7, n_treat = rep(7, 3))
  expect_equal(nrow(balanced$args$data), 28L)
  expect_equal(unname(table(balanced$args$group)), rep(7L, 4L),
               ignore_attr = TRUE)

  # The control and each treatment group are named separately, so an unbalanced
  # design is written rather than assembled out of one recycled number.
  uneven <- small(seed = 1, n_control = 5, n_treat = c(9, 12, 7))
  expect_equal(nrow(uneven$args$data), 33L)
  expect_equal(unname(table(uneven$args$group)[uneven$args$group_lv]),
               c(5L, 9L, 12L, 7L), ignore_attr = TRUE)
  # The sizes are carried in the answer too, since an unbalanced design is one
  # of the reasons a feature can be missed.
  expect_equal(uneven$truth_group$n[1:4], c(5L, 9L, 12L, 7L))
})

test_that("the length of n_treat is the number of treatment groups", {
  five <- small(seed = 1, n_treat = c(9, 8, 7, 6, 5))
  expect_length(five$args$group_lv, 6L)
  expect_identical(five$args$group_lv[6], "treat_5")
  expect_equal(nrow(five$truth_group), 12L * 6L)
  expect_equal(nrow(five$truth_contrast), 12L * 15L)
})

test_that("group_lv and n_treat settle against each other", {
  # Labels alone say how many groups there are, so the default size is spread
  # over them rather than insisted on as a count.
  named <- simulate_multiple_groups(
    n_feats = 12, n_up = 3, n_down = 3,
    group_lv = c("dmso", "low", "mid", "high", "max"), seed = 1
  )
  expect_identical(named$args$group_lv,
                   c("dmso", "low", "mid", "high", "max"))
  expect_equal(nrow(named$truth_group), 12L * 5L)
  expect_true(all(named$truth_group$is_ref ==
                    (named$truth_group$group == "dmso")))
  expect_equal(unname(table(named$args$group)), rep(50L, 5L),
               ignore_attr = TRUE)

  # One size with labels is the same story: the labels supply the count.
  spread <- small(seed = 1, n_treat = 9,
                  group_lv = c("dmso", "low", "mid", "high"))
  expect_equal(unname(table(spread$args$group)[spread$args$group_lv]),
               c(6L, 9L, 9L, 9L), ignore_attr = TRUE)

  # Sizes and labels that count differently are a mistake, not a guess.
  expect_silent(small(seed = 1, n_treat = c(4, 5), group_lv = c("a", "b", "c")))
  expect_error(small(seed = 1, n_treat = c(4, 5, 6, 7),
                     group_lv = c("a", "b", "c")),
               "names 2 treatment group")
})

test_that("n_up and n_down follow n_feats instead of being fixed counts", {
  # The reported failure: the old absolute defaults of 15 and 15 left no room
  # in a 10 feature simulation and the call could not be made at all.
  expect_silent(simulate_multiple_groups(n_feats = 10, n_control = 4,
                                         n_treat = rep(4, 3), seed = 1))

  # Called without n_up or n_down, so the defaults are the thing under test.
  planted <- function(n_feats) {
    sim <- simulate_multiple_groups(n_feats = n_feats, n_control = 4,
                                    n_treat = rep(4, 3), seed = 1)
    c(sum(sim$truth$direction == "up"), sum(sim$truth$direction == "down"))
  }

  expect_equal(planted(10), c(2L, 2L))
  expect_equal(planted(20), c(3L, 3L))
  # At the default n_feats they are the 15 and 15 the defaults were tuned at.
  expect_equal(planted(100), c(15L, 15L))
})

test_that("feat_prefix names the features", {
  sim <- small(seed = 1, feat_prefix = "gene")
  expect_identical(sim$args$feats[1:2], c("gene_1", "gene_2"))
  expect_identical(sim$truth$features, sim$args$feats)
})

test_that("no feature is planted when none is asked for", {
  # sample() on an empty selection is the classic place for this to return
  # everything instead of nothing.
  sim <- small(seed = 1, n_up = 0, n_down = 0)
  expect_true(all(sim$truth$direction == "none"))
  expect_true(all(sim$truth$pattern == "none"))
  expect_true(all(sim$truth_group$delta == 0))

  up_only <- small(seed = 1, n_up = 2, n_down = 0)
  expect_equal(sum(up_only$truth$direction == "up"), 2L)
  expect_equal(sum(up_only$truth$direction == "down"), 0L)
})

test_that("a seed makes the draw reproducible without stealing the stream", {
  expect_equal(small(seed = 42), small(seed = 42))
  expect_false(isTRUE(all.equal(small(seed = 1)$args$data,
                                small(seed = 2)$args$data)))

  # Seeding inside the function must not reseed the caller, or every draw after
  # a simulation would silently repeat.
  set.seed(99)
  before <- stats::runif(3)
  set.seed(99)
  invisible(small(seed = 7))
  expect_equal(stats::runif(3), before)
})

test_that("without a seed the draw follows the caller's stream", {
  set.seed(11)
  first <- small()
  set.seed(11)
  expect_equal(small(), first)
})

test_that("the simulation feeds compare_multiple_groups directly", {
  sim <- small(n_feats = 8, n_up = 2, n_down = 2, n_control = 10,
               n_treat = rep(10, 3), seed = 1)
  res <- do.call(compare_multiple_groups, c(sim$args, diagnose = FALSE))

  expect_identical(res$features, sim$truth$features)
  expect_identical(res$design$group_lv, sim$args$group_lv)
  expect_false(res$design$paired)
  expect_identical(res$parameters$input_scale, "log2")
  # log2 input takes the geometric mean by default, so the estimate is directly
  # comparable to the planted difference of log2 means.
  expect_identical(res$parameters$fc_mean, "geom")
  expect_named(res$tests, c("anova_test", "welch_test", "robust_test",
                            "kruskal_test"))
})

test_that("a repeated simulation runs the within-subject tests whole", {
  sim <- small(n_feats = 8, n_up = 2, n_down = 2, n_control = 12,
               n_treat = rep(12, 3), paired = TRUE, seed = 1)
  res <- do.call(compare_multiple_groups, c(sim$args, diagnose = FALSE))

  expect_true(res$design$paired)
  # No subject is missing a condition, so none is dropped.
  expect_length(res$design$unmatched_ids, 0L)
  expect_named(res$tests, c("anova_test", "kruskal_test"))
  expect_true(all(res$tests$anova_test$n_used == 12L))
  # The subject offset is real, so the answer records the spread it was drawn
  # with rather than leaving it to be guessed at.
  expect_true(all(!is.na(sim$truth$sd_subject)))
  expect_true(all(is.na(small(n_feats = 4, n_up = 1, n_down = 1,
                              seed = 1)$truth$sd_subject)))
})

test_that("the estimated log2fc tracks what was planted", {
  sim <- small(n_feats = 30, n_up = 5, n_down = 5, n_control = 30,
               n_treat = rep(30, 3), seed = 1)
  res <- do.call(compare_multiple_groups, c(sim$args, diagnose = FALSE))

  # Sampling error keeps this from being exact, and the null features contribute
  # noise against a truth of zero, so the bound is on the two agreeing well
  # rather than on them matching.
  expect_gt(stats::cor(res$effect$log2fc, sim$truth$log2fc), 0.8)
  planted <- sim$truth$direction != "none"
  expect_true(all(sign(res$effect$log2fc[planted]) ==
                    sign(sim$truth$log2fc[planted])))

  # Where the answer singles out a level, the comparison should find the same
  # one. Where it does not, only the magnitude is scoreable.
  untied <- planted & !sim$truth$extreme_tied
  expect_gt(mean(res$effect$extreme_level[untied] ==
                   sim$truth$extreme_level[untied]), 0.7)
})

test_that("the post-hoc stage can be scored against truth_contrast", {
  sim <- small(n_feats = 30, n_up = 5, n_down = 5, n_control = 30,
               n_treat = rep(30, 3), seed = 1)
  res <- do.call(compare_multiple_groups, c(sim$args, diagnose = FALSE))
  ph <- merge(res$posthoc$anova_test, sim$truth_contrast,
              by = c("features", "contrast"))

  expect_gt(nrow(ph), 0L)
  expect_equal(nrow(ph), nrow(res$posthoc$anova_test))
  expect_gt(stats::cor(ph$estimate, ph$delta), 0.8)

  # A contrast the pairwise stage called had better point the way the answer
  # does. The two are on the same `group1 - group2` footing, so this is a
  # like-for-like comparison rather than a sign convention being asserted.
  called <- ph$pval_adj <= 0.05
  expect_gt(mean(sign(ph$estimate[called & ph$is_diff]) ==
                   sign(ph$delta[called & ph$is_diff])), 0.9)
  # A contrast with no planted difference is a false positive if it is called.
  expect_lt(mean(called[!ph$is_diff]), 0.2)
})

test_that("the defaults leave a comparison most but not all of the answer", {
  # Loose bounds on purpose. The point is that the defaults sit between a
  # simulation that is trivially recovered and one that looks broken, not that
  # the recall is any particular number.
  sim <- simulate_multiple_groups(seed = 1)
  res <- do.call(compare_multiple_groups, c(sim$args, posthoc = FALSE,
                                            diagnose = FALSE))
  sig <- estimate_significance(res, test = "anova_test")$significance

  planted <- sim$truth$direction != "none"
  recall <- mean(sig$is_signif[planted] %in% TRUE)
  false_positives <- mean(sig$is_signif[!planted] %in% TRUE)

  expect_gt(recall, 0.5)
  expect_lt(recall, 1)
  expect_lt(false_positives, 0.2)
})

test_that("impossible or malformed arguments are rejected", {
  expect_error(small(n_feats = 10, n_up = 6, n_down = 6), "more features")
  expect_error(small(n_feats = 2.5), "whole number")
  expect_error(small(n_control = 1), "must be in")
  expect_error(small(n_treat = c(6, 1)), "must be in")
  expect_error(small(n_treat = character(0)), "one or more group sizes")
  # One size and no labels says how big a group is, not how many there are.
  expect_error(simulate_multiple_groups(n_treat = 30), "at least two of them")
  expect_error(small(group_lv = c("a", "b")), "distinct")
  expect_error(small(group_lv = c("a", "a", "b")), "distinct")
  expect_error(small(expr_range = c(12, 2)), "increasing")
  expect_error(small(treat_sd = c(-1, 2)), "must not go below")
  expect_error(small(deg_log2fc = 1), "length 2")
  expect_error(small(feat_prefix = ""), "non-empty")
})

test_that("a malformed pattern_mix is rejected", {
  expect_error(small(pattern_mix = c(1, 1, 1)), "named numeric vector")
  expect_error(small(pattern_mix = c(all = 1, nope = 1)), "unknown shape")
  expect_error(small(pattern_mix = c(all = -1, single = 1)), "not be negative")
  expect_error(small(pattern_mix = c(all = 0, gradient = 0, single = 0)),
               "at least one positive weight")
})

test_that("a repeated design cannot have groups of different sizes", {
  # Every condition is measured on the same subjects, so unequal sizes describe
  # a design this cannot build. Levelling them silently would hide the mistake.
  expect_error(small(paired = TRUE, n_control = 6, n_treat = c(6, 5, 6)),
               "same number of them")
  expect_error(small(paired = TRUE, n_control = 8, n_treat = rep(6, 3)),
               "same number of them")
  expect_silent(small(paired = TRUE, n_control = 6, n_treat = rep(6, 3)))
})
