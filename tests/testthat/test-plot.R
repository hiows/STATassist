# The drawing itself is not compared against a reference image; vdiffr is not a
# dependency and a pixel comparison would fail on every graphics device change.
# What is tested is which rows the method decided to draw and in what order,
# which is where the logic lives.

local_null_device <- function(env = parent.frame()) {
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 700, height = 500)
  withr::defer(grDevices::dev.off(), envir = env)
  path
}

test_that("a two-group result draws its estimates by default", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_two_group_fixture()
  drawn <- plot(res)
  expect_identical(drawn$features, res$features)
  expect_identical(attr(drawn, "view"), "estimate")
})

test_that("the view auto resolved to is reported on the returned rows", {
  skip_if_not_installed("withr")
  local_null_device()
  # Otherwise a caller has no way to find out which of the three views it got,
  # short of reading the axis label off the device.
  multi <- sa_multi_group_fixture()
  expect_identical(attr(plot(multi), "view"), "posthoc")
  expect_identical(attr(plot(multi, type = "pvalue"), "view"), "pvalue")
  expect_identical(attr(plot(sa_multi_group_fixture(posthoc = FALSE)), "view"),
                   "pvalue")
})

test_that("a multi-group result falls through to the pairwise contrasts", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_multi_group_fixture()
  drawn <- plot(res)
  # The omnibus table has no interval, so "auto" cannot draw estimates and
  # picks the post-hoc view of the first feature rather than erroring.
  expect_equal(nrow(drawn), 3L)
  expect_identical(unique(drawn$features), res$features[1])
})

test_that("auto falls all the way to p-values when nothing else is available", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_multi_group_fixture(posthoc = FALSE)
  drawn <- plot(res)
  expect_identical(drawn$features, res$features)
})

test_that("naming a feature selects that feature's contrasts", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_multi_group_fixture()
  drawn <- plot(res, test = "kruskal_test", feature = "Petal.Width")
  expect_identical(unique(drawn$features), "Petal.Width")
  expect_equal(nrow(drawn), 3L)
})

test_that("sorting by p-value reorders the drawn rows", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_two_group_fixture()
  drawn <- plot(res, sort_by = "pvalue")
  expect_identical(drawn$features, res$features[order(res$tests$t_test$pval_adj)])
  expect_false(identical(drawn$features, res$features))
})

test_that("the p-value view keeps every feature", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_multi_group_fixture()
  drawn <- plot(res, type = "pvalue")
  expect_identical(drawn$features, res$features)
})

test_that("asking for estimates from an omnibus table is an explicit error", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_multi_group_fixture()
  expect_error(plot(res, type = "estimate"), "no estimate to draw")
})

test_that("asking for contrasts where there are none is an explicit error", {
  skip_if_not_installed("withr")
  local_null_device()
  expect_error(plot(sa_two_group_fixture(), type = "posthoc"),
               "no contrasts to draw")
})

test_that("an unknown feature names the ones that are available", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_multi_group_fixture()
  expect_error(plot(res, feature = "nope"), "Petal.Length")
})

test_that("an unknown test is rejected", {
  skip_if_not_installed("withr")
  local_null_device()
  expect_error(plot(sa_two_group_fixture(), test = "anova_test"), "t_test")
})

test_that("a one-sample result draws like any other", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_one_sample_fixture()
  drawn <- plot(res, test = "wilcox_test")
  expect_identical(drawn$features, res$features)
})

test_that("the device is left as it was found", {
  skip_if_not_installed("withr")
  local_null_device()
  before <- graphics::par(no.readonly = TRUE)
  plot(sa_two_group_fixture(), dark = TRUE)
  after <- graphics::par(no.readonly = TRUE)
  expect_identical(after$bg, before$bg)
  expect_identical(after$mar, before$mar)
})
