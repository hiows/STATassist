# `input_scale = "log2"` exists because dividing two means of logged values
# answers a different question than the fold change it looks like. The tests
# below pin the identity that makes the back-transformation correct, the sign it
# gets right and the boundary between what it touches and what it leaves alone.

sa_logged_iris <- function() {
  d <- iris[iris$Species != "setosa", ]
  d[sa_feats()] <- log2(d[sa_feats()])
  d
}

test_that("log2 input reproduces the geometric mean fold change of raw input", {
  d <- iris[iris$Species != "setosa", ]
  raw_geom <- compare_two_groups(d, sa_feats(), d$Species,
                                 c("virginica", "versicolor"),
                                 fc_mean = "geom", diagnose = FALSE)
  logged <- sa_logged_iris()
  on_log2 <- compare_two_groups(logged, sa_feats(), logged$Species,
                                c("virginica", "versicolor"),
                                input_scale = "log2", diagnose = FALSE)

  # exp(mean(log(2^v))) is 2^mean(v), so the two paths are the same quantity
  # reached from different sides.
  expect_equal(on_log2$effect, raw_geom$effect)
})

test_that("log2 input reduces the fold change to a difference of means", {
  logged <- sa_logged_iris()
  res <- compare_two_groups(logged, sa_feats(), logged$Species,
                            c("virginica", "versicolor"),
                            input_scale = "log2", diagnose = FALSE)
  row <- res$effect$features == "Petal.Length"
  x <- logged$Petal.Length[logged$Species == "virginica"]
  y <- logged$Petal.Length[logged$Species == "versicolor"]

  expect_equal(res$effect$log2fc[row], mean(x) - mean(y))
  expect_equal(res$effect$fold_change[row], 2^(mean(x) - mean(y)))
  # The invariant the whole design turns on: the centres and the ratio stay on
  # one scale, so the effect row remains internally consistent.
  expect_equal(res$effect$fold_change, res$effect$x_center / res$effect$y_center)
})

test_that("two negative log2 centres are an increase, not a decrease", {
  # Overlapping spreads around the two centres, so that every test has an
  # answer and the only thing under examination is the effect table.
  d <- data.frame(v = c(-2, -1, 0, -3, -2, -1),
                  g = rep(c("a", "b"), each = 3))
  args <- list(data = d, feats = "v", group = d$g, group_lv = c("a", "b"),
               diagnose = FALSE)

  on_log2 <- do.call(compare_two_groups, c(args, input_scale = "log2"))
  # Group a sits at log2 -1 and group b at log2 -2, so a is twice b.
  expect_equal(on_log2$effect$log2fc, 1)

  # Treating the same numbers as raw divides -1 by -2 and reports a two-fold
  # decrease instead, which is what the argument exists to avoid.
  as_raw <- do.call(compare_two_groups, c(args, fc_mean = "arith"))
  expect_equal(as_raw$effect$log2fc, -1)
})

test_that("input_scale changes the effect table and nothing else", {
  logged <- sa_logged_iris()
  common <- list(data = logged, feats = sa_feats(), group = logged$Species,
                 group_lv = c("virginica", "versicolor"), diagnose = FALSE)
  untouched <- do.call(compare_two_groups, c(common, fc_mean = "arith"))
  converted <- do.call(compare_two_groups, c(common, input_scale = "log2"))

  expect_equal(untouched$tests, converted$tests)
  expect_false(isTRUE(all.equal(untouched$effect, converted$effect)))
  # The t-test keeps the scale it ran on while the effect table does not, so
  # these two columns no longer agree the way they do for raw input.
  expect_equal(untouched$effect$x_center, untouched$tests$t_test$x_mean)
  expect_false(isTRUE(all.equal(converted$effect$x_center,
                                converted$tests$t_test$x_mean)))
})

test_that("fc_mean defaults to geom on the log2 scale and arith otherwise", {
  logged <- sa_logged_iris()
  fc_mean_of <- function(...) {
    compare_two_groups(logged, "Sepal.Length", logged$Species,
                       c("virginica", "versicolor"), diagnose = FALSE,
                       ...)$parameters$fc_mean
  }

  expect_identical(fc_mean_of(input_scale = "log2"), "geom")
  expect_identical(fc_mean_of(), "arith")
  # An explicit choice still wins over the scale-dependent default.
  expect_identical(fc_mean_of(input_scale = "log2", fc_mean = "arith"), "arith")
})

test_that("raw values passed as log2 are reported rather than silently infinite", {
  d <- data.frame(v = c(4e4, 6e4, 8e4, 5e4, 7e4, 9e4),
                  g = rep(c("a", "b"), each = 3))
  expect_warning(
    res <- compare_two_groups(d, "v", d$g, c("a", "b"),
                              input_scale = "log2", diagnose = FALSE),
    "not on the log2 scale"
  )
  expect_true(is.na(res$effect$log2fc))
  # The t-test does not depend on the fold change, so it still has an answer.
  expect_false(is.na(res$tests$t_test$pval))
})

test_that("the multi-group effect table converts the same way", {
  raw_geom <- compare_multiple_groups(iris, sa_feats(), iris$Species,
                                      levels(iris$Species), fc_mean = "geom",
                                      diagnose = FALSE, posthoc = FALSE)
  logged <- iris
  logged[sa_feats()] <- log2(iris[sa_feats()])
  on_log2 <- compare_multiple_groups(logged, sa_feats(), logged$Species,
                                     levels(logged$Species),
                                     input_scale = "log2", diagnose = FALSE,
                                     posthoc = FALSE)

  expect_equal(on_log2$effect, raw_geom$effect)
  expect_identical(on_log2$parameters$input_scale, "log2")
})

test_that("the one-sample ratio takes mu on the scale of the data", {
  raw_geom <- suppressWarnings(
    compare_one_sample(iris, sa_feats(), mu = 3, fc_mean = "geom",
                       diagnose = FALSE)
  )
  logged <- iris
  logged[sa_feats()] <- log2(iris[sa_feats()])
  on_log2 <- suppressWarnings(
    compare_one_sample(logged, sa_feats(), mu = log2(3), input_scale = "log2",
                       diagnose = FALSE)
  )

  expect_equal(on_log2$effect, raw_geom$effect)
  # The tests were given log2(3) and report it unconverted.
  expect_equal(on_log2$tests$t_test$mu, rep(log2(3), length(sa_feats())))
})

test_that("mu = 0 is an ordinary reference on the log2 scale", {
  logged <- iris
  logged[sa_feats()] <- log2(iris[sa_feats()])
  on_log2 <- suppressWarnings(
    compare_one_sample(logged, sa_feats(), mu = 0, input_scale = "log2",
                       diagnose = FALSE)
  )

  # 2^0 is 1, so the ratio is defined and the raw-scale warning must not fire.
  expect_false(any(is.na(on_log2$effect$fold_change)))
  expect_equal(on_log2$effect$mu, rep(1, length(sa_feats())))
  expect_equal(on_log2$effect$log2fc,
               vapply(sa_feats(), function(f) mean(logged[[f]]), numeric(1)),
               ignore_attr = TRUE)
})

test_that("a mu that cannot be a log2 value is rejected outright", {
  expect_error(
    compare_one_sample(iris, "Sepal.Length", mu = 5e4, input_scale = "log2",
                       diagnose = FALSE),
    "not on the log2 scale"
  )
})
