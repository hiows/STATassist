# The value of a simulator is that the answer is known, so what is pinned here
# is mostly that the answer it reports is the answer it planted: the counts, the
# strict nullity of everything else, and the fact that a comparison run on the
# returned arguments lines up with the truth table row for row.

test_that("the returned arguments are the ones compare_two_groups takes", {
  sim <- simulate_two_groups(seed = 1, n_feats = 12, n_up = 3, n_down = 3)

  expect_named(sim, c("args", "truth"))
  expect_named(sim$args,
               c("data", "feats", "group", "group_lv", "input_scale"))
  expect_true(all(names(sim$args) %in% names(formals(compare_two_groups))))
  expect_identical(sim$args$input_scale, "log2")
  expect_identical(names(sim$args$data), sim$args$feats)
  expect_length(sim$args$group, nrow(sim$args$data))
})

test_that("the truth table is aligned with the features", {
  sim <- simulate_two_groups(seed = 1, n_feats = 12, n_up = 3, n_down = 3)

  expect_identical(sim$truth$features, sim$args$feats)
  expect_named(sim$truth, c("features", "direction", "log2fc", "baseline",
                            "sd_case", "sd_control"))
})

test_that("exactly the requested number of features is planted", {
  sim <- simulate_two_groups(seed = 3, n_feats = 40, n_up = 7, n_down = 11)

  expect_equal(sum(sim$truth$direction == "up"), 7L)
  expect_equal(sum(sim$truth$direction == "down"), 11L)
  expect_equal(sum(sim$truth$direction == "none"), 22L)
  # Up and down are disjoint: a feature cannot be both, which a complement-based
  # split would allow if the two index sets overlapped.
  expect_true(all(sim$truth$log2fc[sim$truth$direction == "up"] > 0))
  expect_true(all(sim$truth$log2fc[sim$truth$direction == "down"] < 0))
})

test_that("an unplanted feature is null in the strict sense", {
  sim <- simulate_two_groups(seed = 4)
  # Not "small", exactly zero. Everything downstream that scores a false
  # positive rate depends on there being no true effect here at all.
  expect_true(all(sim$truth$log2fc[sim$truth$direction == "none"] == 0))
})

test_that("the planted magnitudes stay inside deg_log2fc", {
  sim <- simulate_two_groups(seed = 5, deg_log2fc = c(1.5, 2))
  planted <- abs(sim$truth$log2fc[sim$truth$direction != "none"])

  expect_true(all(planted >= 1.5))
  expect_true(all(planted <= 2))
})

test_that("the two group sizes are honoured separately", {
  # The draft this came from generated the control group with n_case rows, so an
  # unequal design silently produced a group vector of the wrong length.
  sim <- simulate_two_groups(seed = 1, n_feats = 6, n_up = 1, n_down = 1,
                             n_case = 9, n_control = 21)

  expect_equal(nrow(sim$args$data), 30L)
  expect_equal(unname(table(sim$args$group)[["case"]]), 9L)
  expect_equal(unname(table(sim$args$group)[["control"]]), 21L)
})

test_that("group_lv sets both the labels and the direction", {
  sim <- simulate_two_groups(seed = 1, n_feats = 6, n_up = 2, n_down = 2,
                             group_lv = c("untreated", "treated"))

  expect_identical(sim$args$group_lv, c("untreated", "treated"))
  expect_identical(sort(unique(sim$args$group)), c("treated", "untreated"))
})

test_that("the control is the first level and the effect goes on the second", {
  sim <- simulate_two_groups(seed = 1, n_feats = 6, n_up = 2, n_down = 2)

  expect_identical(sim$args$group_lv, c("control", "case"))
  # The rows follow the levels, so a comparison and a boxplot of the same
  # arguments both read control first.
  expect_identical(unique(sim$args$group), c("control", "case"))

  res <- do.call(compare_two_groups, c(sim$args, diagnose = FALSE))
  up <- sim$truth$direction == "up"
  expect_true(all(res$effect$log2fc[up] > 0))
})

test_that("no feature is planted when none is asked for", {
  # sample() on an empty selection is the classic place for this to return
  # everything instead of nothing.
  sim <- simulate_two_groups(seed = 1, n_feats = 5, n_up = 0, n_down = 0)
  expect_true(all(sim$truth$direction == "none"))
  expect_true(all(sim$truth$log2fc == 0))

  up_only <- simulate_two_groups(seed = 1, n_feats = 5, n_up = 2, n_down = 0)
  expect_equal(sum(up_only$truth$direction == "up"), 2L)
  expect_equal(sum(up_only$truth$direction == "down"), 0L)
})

test_that("a seed makes the draw reproducible without stealing the stream", {
  small <- function(...) {
    simulate_two_groups(n_feats = 5, n_up = 1, n_down = 1, ...)
  }
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
  first <- simulate_two_groups(n_feats = 5, n_up = 1, n_down = 1)
  set.seed(11)
  expect_equal(simulate_two_groups(n_feats = 5, n_up = 1, n_down = 1), first)
})

test_that("the simulation feeds compare_two_groups directly", {
  sim <- simulate_two_groups(seed = 1, n_feats = 8, n_up = 2, n_down = 2)
  res <- do.call(compare_two_groups, c(sim$args, diagnose = FALSE))

  expect_identical(res$features, sim$truth$features)
  expect_identical(res$design$group_lv, sim$args$group_lv)
  expect_identical(res$parameters$input_scale, "log2")
  # log2 input takes the geometric mean by default, so the estimate is directly
  # comparable to the planted difference of log2 means.
  expect_identical(res$parameters$fc_mean, "geom")
})

test_that("the estimated log2fc tracks what was planted", {
  sim <- simulate_two_groups(seed = 2)
  res <- do.call(compare_two_groups, c(sim$args, diagnose = FALSE))

  # Sampling error keeps this from being exact, and the null features contribute
  # noise against a truth of zero, so the bound is on the two agreeing well
  # rather than on them matching.
  expect_gt(stats::cor(res$effect$log2fc, sim$truth$log2fc), 0.8)
  planted <- sim$truth$direction != "none"
  expect_true(all(sign(res$effect$log2fc[planted]) ==
                    sign(sim$truth$log2fc[planted])))
})

test_that("the defaults leave a comparison most but not all of the answer", {
  # Loose bounds on purpose. The point is that the defaults sit between a
  # simulation that is trivially recovered and one that looks broken, not that
  # the recall is any particular number.
  sim <- simulate_two_groups(seed = 1)
  res <- do.call(compare_two_groups, c(sim$args, diagnose = FALSE))
  sig <- estimate_significance(res)$significance

  planted <- sim$truth$direction != "none"
  recall <- mean(sig$is_signif[planted] %in% TRUE)
  false_positives <- mean(sig$is_signif[!planted] %in% TRUE)

  expect_gt(recall, 0.5)
  expect_lt(recall, 1)
  expect_lt(false_positives, 0.2)
})

test_that("impossible or malformed arguments are rejected", {
  expect_error(simulate_two_groups(n_feats = 10, n_up = 6, n_down = 6),
               "more features")
  expect_error(simulate_two_groups(n_feats = 2.5), "whole number")
  expect_error(simulate_two_groups(n_case = 1), "must be in")
  expect_error(simulate_two_groups(expr_range = c(12, 2)), "increasing")
  expect_error(simulate_two_groups(case_sd = c(-1, 2)), "must not go below")
  expect_error(simulate_two_groups(deg_log2fc = 1), "length 2")
  expect_error(simulate_two_groups(group_lv = c("a", "a")), "distinct")
  expect_error(simulate_two_groups(group_lv = "case"), "distinct")
})
