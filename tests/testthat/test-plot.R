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

test_that("a volcano plot with nothing above the cutoffs still draws", {
  skip_if_not_installed("withr")
  local_null_device()
  # A fold change no feature can reach leaves every verdict FALSE, which used to
  # send text() a zero-length `labels` and error out instead of drawing.
  sig <- estimate_significance(sa_two_group_fixture(), log2fc_cutoff = 100)
  expect_false(any(sig$is_signif))
  expect_message(drawn <- draw_volcano_plot(sig), "nothing was labelled")
  expect_null(drawn)
})

test_that("labels are still drawn when features do clear the cutoffs", {
  skip_if_not_installed("withr")
  local_null_device()
  sig <- estimate_significance(sa_two_group_fixture(), log2fc_cutoff = 0.1)
  expect_true(any(sig$is_signif))
  expect_no_message(draw_volcano_plot(sig))
})

sa_butterfly_data <- function() iris[iris$Species != "setosa", ]

sa_butterfly_levels <- c("versicolor", "virginica")

test_that("bars alone leave the density estimate out of the result", {
  skip_if_not_installed("withr")
  local_null_device()
  iris2 <- sa_butterfly_data()
  res <- draw_butterfly_hist(iris2, "Petal.Length", iris2$Species,
                             sa_butterfly_levels)
  expect_identical(names(res),
                   c("bin_summary_stats", "group_summary_stats", "group_hists"))
})

test_that("drawing a density returns one density object per group level", {
  skip_if_not_installed("withr")
  local_null_device()
  iris2 <- sa_butterfly_data()
  res <- draw_butterfly_hist(iris2, "Petal.Length", iris2$Species,
                             sa_butterfly_levels, type = "both")
  expect_identical(names(res$group_densities), sa_butterfly_levels)
  expect_true(all(vapply(res$group_densities, inherits, logical(1), "density")))
  # Titled after the feature and the level, not after the internal argument.
  expect_identical(res$group_densities$virginica$data.name,
                   "Petal.Length (virginica)")
})

test_that("a density curve moves the bars to the density scale", {
  skip_if_not_installed("withr")
  local_null_device()
  iris2 <- sa_butterfly_data()
  res <- draw_butterfly_hist(iris2, "Petal.Length", iris2$Species,
                             sa_butterfly_levels, type = "dens")
  expect_equal(res$bin_summary_stats$virginica,
               res$group_hists$virginica$density)
})

test_that("a scale that cannot be read against a curve is refused", {
  skip_if_not_installed("withr")
  local_null_device()
  iris2 <- sa_butterfly_data()
  expect_error(
    draw_butterfly_hist(iris2, "Petal.Length", iris2$Species,
                        sa_butterfly_levels, scale = "count", type = "both"),
    "density scale"
  )
})

test_that("the density fill opacity has to be a fraction", {
  skip_if_not_installed("withr")
  local_null_device()
  iris2 <- sa_butterfly_data()
  expect_error(
    draw_butterfly_hist(iris2, "Petal.Length", iris2$Species,
                        sa_butterfly_levels, type = "dens", dens_alpha = 1.5),
    "dens_alpha"
  )
})

test_that("a group with one distinct value names itself in the error", {
  skip_if_not_installed("withr")
  local_null_device()
  iris2 <- sa_butterfly_data()
  iris2$Petal.Length[iris2$Species == "virginica"] <- 5
  expect_error(
    draw_butterfly_hist(iris2, "Petal.Length", iris2$Species,
                        sa_butterfly_levels, type = "dens"),
    "virginica"
  )
})
