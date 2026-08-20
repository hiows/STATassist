# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
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

test_that("plot() is the draw_forest_plot() it delegates to", {
  skip_if_not_installed("withr")
  local_null_device()
  # The method carries no logic of its own, so the two entry points cannot draw
  # different rows or resolve "auto" to different views.
  two <- sa_two_group_fixture()
  multi <- sa_multi_group_fixture()
  expect_identical(plot(two), draw_forest_plot(two))
  expect_identical(plot(multi), draw_forest_plot(multi))
  expect_identical(attr(plot(multi), "view"),
                   attr(draw_forest_plot(multi), "view"))
  # Arguments named on the method reach the function through `...`.
  expect_identical(plot(multi, test = "kruskal_test", type = "pvalue"),
                   draw_forest_plot(multi, test = "kruskal_test",
                                    type = "pvalue"))
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
  drawn <- plot(res, test = "kruskal_test", feats = "Petal.Width")
  expect_identical(unique(drawn$features), "Petal.Width")
  expect_equal(nrow(drawn), 3L)
})

test_that("naming several features draws all of their contrasts", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_multi_group_fixture()
  drawn <- plot(res, type = "posthoc",
                feats = c("Petal.Width", "Sepal.Length"))
  # In the order they were named, contrasts of one feature staying together.
  expect_identical(unique(drawn$features), c("Petal.Width", "Sepal.Length"))
  expect_equal(nrow(drawn), 6L)
})

test_that("`feats` selects and orders the rows of the estimate view", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_two_group_fixture()
  drawn <- plot(res, feats = c("Petal.Width", "Sepal.Width"))
  expect_identical(drawn$features, c("Petal.Width", "Sepal.Width"))
  expect_identical(attr(drawn, "view"), "estimate")
})

test_that("`feats` selects the rows of the p-value view too", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_multi_group_fixture()
  drawn <- plot(res, type = "pvalue", feats = "Sepal.Width")
  expect_identical(drawn$features, "Sepal.Width")
})

test_that("a feature without contrasts is reported rather than drawn", {
  skip_if_not_installed("withr")
  local_null_device()
  # An absent post-hoc row means the omnibus test never qualified the feature,
  # which is a fact about the result rather than a mistake in the call.
  res <- sa_multi_group_fixture()
  res$posthoc$anova_test <-
    res$posthoc$anova_test[res$posthoc$anova_test$features != "Sepal.Width", ]
  expect_message(
    drawn <- plot(res, type = "posthoc",
                  feats = c("Petal.Width", "Sepal.Width")),
    "Sepal.Width"
  )
  expect_identical(unique(drawn$features), "Petal.Width")
})

test_that("sorting by p-value reorders the drawn rows", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_two_group_fixture()
  drawn <- plot(res, sort_by = "pvalue")
  expect_identical(drawn$features, res$features[order(res$tests$t_test$pval_adj)])
  expect_false(identical(drawn$features, res$features))
})

test_that("the unadjusted p-value can be the one the plot reads", {
  skip_if_not_installed("withr")
  local_null_device()
  # Colouring, sorting and the axis label all follow `use_adjusted`, so the
  # order of the rows is enough to tell which column was read.
  res <- sa_multi_group_fixture()
  drawn <- plot(res, type = "pvalue", use_adjusted = FALSE,
                sort_by = "pvalue")
  expect_identical(drawn$features,
                   res$features[order(res$tests$anova_test$pval)])
})

test_that("the p-value view keeps every feature", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_multi_group_fixture()
  drawn <- plot(res, type = "pvalue")
  expect_identical(drawn$features, res$features)
})

test_that("an axis range of the caller's own is drawn and checked", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_two_group_fixture()
  # A range narrower than the intervals is a legitimate request: the bounds are
  # clamped to it rather than the plot being widened back out.
  drawn <- plot(res, xlim = c(0, 1))
  expect_identical(drawn$features, res$features)
  expect_identical(plot(res, type = "pvalue", xlim = c(0, 5))$features,
                   res$features)
  expect_error(plot(res, xlim = c(0, 1, 2)), "xlim")
  expect_error(plot(res, xlim = c(0, Inf)), "xlim")
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
  expect_error(plot(res, feats = "nope"), "Petal.Length")
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
  expect_identical(after$fg, before$fg)
  expect_identical(after$col.axis, before$col.axis)
  expect_identical(after$mar, before$mar)
  # The legend sits in a panel of its own, so the layout has to come back too.
  expect_identical(after$mfrow, before$mfrow)
  # Nothing that carries an absolute size may be written, or the next plot is
  # pinned to the size this one was drawn at.
  expect_identical(after$fin, before$fin)
  expect_identical(after$pin, before$pin)
  expect_identical(after$mai, before$mai)
})

test_that("a panel grid the caller set up survives the call", {
  skip_if_not_installed("withr")
  local_null_device()
  # `layout()` overwrites the caller's grid while the legend panel is drawn, so
  # the grid has to be put back rather than left at the one panel `layout(1)`
  # resets to.
  withr::local_par(list(mfrow = c(2, 2)))
  plot(sa_two_group_fixture())
  expect_identical(graphics::par("mfrow"), c(2L, 2L))
})

test_that("a volcano plot with nothing above the cutoffs still draws", {
  skip_if_not_installed("withr")
  local_null_device()
  # A fold change no feature can reach leaves every verdict FALSE, which used to
  # send text() a zero-length `labels` and error out instead of drawing.
  sig <- estimate_significance(sa_two_group_fixture(), log2fc_cutoff = 100)
  expect_false(any(sig$significance$is_signif))
  expect_message(drawn <- draw_volcano_plot(sig), "nothing was labelled")
  expect_null(drawn)
})

test_that("labels are still drawn when features do clear the cutoffs", {
  skip_if_not_installed("withr")
  local_null_device()
  sig <- estimate_significance(sa_two_group_fixture(), log2fc_cutoff = 0.1)
  expect_true(any(sig$significance$is_signif))
  expect_no_message(draw_volcano_plot(sig))
})

test_that("the verdict table can be plotted without its wrapper", {
  skip_if_not_installed("withr")
  local_null_device()
  # Selecting rows from the table produces a table rather than the object, so
  # both have to reach the plot.
  sig <- estimate_significance(sa_two_group_fixture(), log2fc_cutoff = 0.1)
  expect_null(draw_volcano_plot(sig$significance))
  expect_null(draw_volcano_plot(sig$significance[1:2, ]))
})

test_that("a contrast reading says which table to name", {
  skip_if_not_installed("withr")
  local_null_device()
  sig <- estimate_significance(sa_multi_group_fixture(), by = "contrast",
                               log2fc_cutoff = 0.1)
  expect_error(draw_volcano_plot(sig), "one verdict table per contrast")
  expect_null(draw_volcano_plot(sig$significance[["virginica - setosa"]]))
})

test_that("a term reading is drawn as one panel per term", {
  skip_if_not_installed("withr")
  local_null_device()
  sig <- estimate_significance(sa_factorial_fixture(posthoc = FALSE),
                              by = "term", log2fc_cutoff = 0.1)

  # The panel grid is undone when the figure is finished, so the caller's own
  # arrangement survives a plot drawn inside it.
  withr::local_par(list(mfrow = c(2, 2)))
  # One feature is one point, and the panel it does not clear the cutoffs in says
  # so, once per panel.
  expect_null(suppressMessages(draw_volcano_plot(sig, main = "warpbreaks")))
  expect_identical(graphics::par("mfrow"), c(2L, 2L))

  # A single term table is still a table, so it draws on its own.
  expect_null(draw_volcano_plot(sig$significance[["wool:tension"]]))
})

test_that("terms says which panels a design of three factors gets", {
  skip_if_not_installed("withr")
  local_null_device()
  dat <- expand.grid(a = c("a1", "a2"), b = c("b1", "b2"), c = c("c1", "c2"),
                     rep = 1:3, stringsAsFactors = FALSE)
  set.seed(20260813)
  dat$y <- stats::rnorm(nrow(dat), mean = 4 + as.integer(dat$a == "a2"))
  res <- compare_factorial_groups(dat["y"], "y",
                                  list(a = dat$a, b = dat$b, c = dat$c),
                                  posthoc = FALSE, diagnose = FALSE)
  sig <- estimate_significance(res, by = "term", log2fc_cutoff = 0.1)
  tables <- sig$significance
  expect_length(tables, 7L)

  # Three panels of the first two factors by default, and the four terms left
  # out are named rather than silently dropped.
  expect_identical(suppressMessages(sa_volcano_terms(tables, NULL)),
                   c("a", "b", "a:b"))
  expect_message(sa_volcano_terms(tables, NULL),
                 "leaving out c, a:c, b:c, a:b:c")
  expect_no_message(sa_volcano_terms(tables[c("a", "b", "a:b")], NULL))

  expect_identical(sa_volcano_terms(tables, c("c", "a:b:c")), c("c", "a:b:c"))
  expect_error(sa_volcano_terms(tables, "a:d"), "does not hold: a:d")
  expect_error(draw_volcano_plot(sig, terms = 3), "character vector")

  expect_null(draw_volcano_plot(sig, terms = c("a", "b", "c", "a:b:c"),
                                panel_nrow = 2, anno_feats = FALSE))
})

test_that("panels drawn together are drawn on one scale", {
  sig <- estimate_significance(sa_factorial_fixture(posthoc = FALSE),
                              by = "term", log2fc_cutoff = 0.1)
  tables <- sig$significance
  shared <- sa_volcano_lims(tables, "adj_pvalue", 0.1, 0.05, NULL, NULL)
  each <- lapply(tables, function(tbl) {
    sa_volcano_lims(list(tbl), "adj_pvalue", 0.1, 0.05, NULL, NULL)
  })

  # The shared range covers every panel's own, which is what makes a point in one
  # panel comparable with a point in another.
  for (one in each) {
    expect_lte(shared$xlim[1], one$xlim[1])
    expect_gte(shared$xlim[2], one$xlim[2])
    expect_gte(shared$ylim[2], one$ylim[2])
  }
  expect_gt(shared$ylim[2], min(vapply(each, function(l) l$ylim[2], numeric(1))))

  # A supplied range is used as given, whatever the data would have asked for.
  fixed <- sa_volcano_lims(tables, "adj_pvalue", 0.1, 0.05, c(-9, 9), c(0, 3))
  expect_identical(fixed$xlim, c(-9, 9))
  expect_identical(fixed$ylim, c(0, 3))
})

test_that("the x axis names what a multi-group log2fc compares", {
  multi <- estimate_significance(sa_multi_group_fixture(),
                                 log2fc_cutoff = 0.1)$significance
  lab <- sa_volcano_xlab(multi)
  expect_true(grepl("most extreme level vs setosa",
                    deparse(lab), fixed = TRUE))
  expect_true("extreme_level" %in% names(multi))
  expect_identical(multi$extreme_level,
                   sa_multi_group_fixture()$effect$extreme_level)

  fact_res <- sa_factorial_fixture(posthoc = FALSE)
  fact <- estimate_significance(fact_res, log2fc_cutoff = 0.1)$significance
  ref <- fact_res$design$group_lv[1]
  lab <- sa_volcano_xlab(fact)
  expect_true(grepl("most extreme cell vs", deparse(lab), fixed = TRUE))
  expect_true(grepl(ref, deparse(lab), fixed = TRUE))
  expect_true("extreme_cell" %in% names(fact))
  expect_identical(fact$extreme_cell, fact_res$effect$extreme_cell)

  # The label names the cell the ratios were actually divided by, so moving the
  # reference moves the label with it.
  moved_res <- sa_factorial_fixture(posthoc = FALSE,
                                    control_label = list(tension = "M"))
  moved <- estimate_significance(moved_res, log2fc_cutoff = 0.1)$significance
  expect_true(grepl("most extreme cell vs A.M",
                    deparse(sa_volcano_xlab(moved)), fixed = TRUE))

  # A two-group verdict and a single contrast both compare two named levels, so
  # neither earns the qualifier.
  two <- estimate_significance(sa_two_group_fixture())$significance
  expect_identical(sa_volcano_xlab(two), expression(log[2] ~ FC))
  ct <- estimate_significance(sa_multi_group_fixture(),
                              by = "contrast")$significance
  expect_identical(sa_volcano_xlab(ct[[1]]), expression(log[2] ~ FC))

  # A term component is not a ratio of two centres at all, so it is not labelled
  # as a fold change.
  term <- estimate_significance(sa_factorial_fixture(posthoc = FALSE),
                                by = "term")$significance
  lab <- deparse(sa_volcano_xlab(term[["wool:tension"]]))
  expect_true(grepl("effect", lab, fixed = TRUE))
  expect_true(grepl("(wool:tension)", lab, fixed = TRUE))
  expect_false(grepl("FC", lab, fixed = TRUE))
})

test_that("xlab overrides the derived label", {
  skip_if_not_installed("withr")
  local_null_device()
  sig <- estimate_significance(sa_multi_group_fixture(), log2fc_cutoff = 0.1)
  expect_null(draw_volcano_plot(sig, xlab = "fold change"))
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

sa_heatmap <- function(...) {
  draw_heatmap(sa_heatmap_matrix(), sa_heatmap_group(), sa_heatmap_levels(),
               ...)
}

test_that("the input is transposed so that features run down the rows", {
  skip_if_not_installed("withr")
  local_null_device()
  # The wide input has one row per sample, and the plot has one row per feature.
  out <- sa_heatmap()
  expect_identical(dim(out$matrix), c(4L, 6L))
  expect_setequal(rownames(out$matrix), colnames(sa_heatmap_matrix()))
  expect_setequal(colnames(out$matrix), rownames(sa_heatmap_matrix()))
})

test_that("the clustering the picture shows is the clustering it reports", {
  skip_if_not_installed("withr")
  local_null_device()
  out <- sa_heatmap()
  # Both blocks are recovered, which is what makes the drawn order meaningful.
  expect_identical(unname(stats::cutree(out$feat_hclust, 2)),
                   c(1L, 2L, 1L, 2L))
  expect_identical(unname(stats::cutree(out$sample_hclust, 2)),
                   c(1L, 1L, 1L, 2L, 2L, 2L))
  # Cells, labels and dendrogram leaves all follow the same permutation, so a
  # row of the returned matrix is the row that was drawn under its own label.
  expect_identical(out$feat_order, out$feat_hclust$order)
  expect_identical(out$sample_order, out$sample_hclust$order)
  expect_identical(rownames(out$matrix),
                   colnames(sa_heatmap_matrix())[out$feat_order])
  expect_identical(colnames(out$matrix),
                   rownames(sa_heatmap_matrix())[out$sample_order])
})

test_that("the order the cells are drawn in is the order that is reported", {
  skip_if_not_installed("withr")
  local_null_device()
  # The picture is drawn by handing stats::heatmap() the dendrograms, and it
  # takes its own permutation from them with order.dendrogram() rather than being
  # told one. The reported order is hclust()'s, so the two have to agree or the
  # returned matrix is not the matrix on the screen.
  out <- sa_heatmap()
  expect_identical(
    stats::order.dendrogram(stats::as.dendrogram(out$feat_hclust)),
    out$feat_order
  )
  expect_identical(
    stats::order.dendrogram(stats::as.dendrogram(out$sample_hclust)),
    out$sample_order
  )
})

test_that("feats picks the features to draw and, unclustered, their order", {
  skip_if_not_installed("withr")
  local_null_device()
  out <- sa_heatmap(feats = c("f4", "f1"), cluster_feats = FALSE,
                    cluster_samples = FALSE)
  expect_identical(rownames(out$matrix), c("f4", "f1"))
  expect_identical(colnames(out$matrix), paste0("s", 1:6))
  expect_identical(out$feat_order, 1:2)
  expect_identical(out$sample_order, 1:6)
  expect_null(out$feat_hclust)
  expect_null(out$sample_hclust)
})

test_that("scale standardises the margin it names and nothing else", {
  skip_if_not_installed("withr")
  local_null_device()
  feature <- sa_heatmap(scale = "feature")
  expect_equal(unname(rowMeans(feature$matrix)), rep(0, 4))
  expect_equal(unname(apply(feature$matrix, 1, stats::sd)), rep(1, 4))

  sample <- sa_heatmap(scale = "sample")
  expect_equal(unname(colMeans(sample$matrix)), rep(0, 6))

  # "none" draws the values as they arrived, so the drawn matrix is the input
  # transposed and nothing more.
  none <- sa_heatmap(scale = "none", cluster_feats = FALSE,
                     cluster_samples = FALSE)
  expect_equal(none$matrix, t(sa_heatmap_matrix()))
})

test_that("the colour range is symmetric only when the values have both signs", {
  skip_if_not_installed("withr")
  local_null_device()
  # z-scores always straddle zero, and zero is the midpoint of the ramp.
  scaled <- sa_heatmap()
  expect_equal(scaled$zlim[1], -scaled$zlim[2])
  expect_equal(scaled$zlim[2], max(abs(scaled$matrix)))

  # Positive values have no such midpoint, so the range is the range of the data
  # and white lands in the middle of it.
  raw <- sa_heatmap(scale = "none")
  expect_equal(raw$zlim, range(sa_heatmap_matrix()))
})

test_that("a supplied range draws the values outside it at the end of the scale", {
  skip_if_not_installed("withr")
  local_null_device()
  expect_message(out <- sa_heatmap(zlim = c(-0.5, 0.5)), "outside `zlim`")
  expect_identical(out$zlim, c(-0.5, 0.5))
  # The clamp is a decision about colour, so the reported values keep theirs.
  expect_true(any(abs(out$matrix) > 0.5))
})

test_that("a feature with no variance is centred rather than divided by zero", {
  skip_if_not_installed("withr")
  local_null_device()
  m <- sa_heatmap_matrix()
  m[, "f3"] <- 7
  expect_message(
    out <- draw_heatmap(m, sa_heatmap_group(), sa_heatmap_levels()),
    "no variance"
  )
  expect_true(all(out$matrix["f3", ] == 0))
  expect_false(anyNA(out$matrix))
})

test_that("undefined distances leave that axis in its input order", {
  skip_if_not_installed("withr")
  local_null_device()
  # Two features observed on disjoint samples have no distance between them at
  # all, which used to be an error from hclust() halfway through the plot.
  m <- sa_heatmap_matrix()
  m[1:3, "f1"] <- NA
  m[4:6, "f2"] <- NA
  expect_message(
    out <- draw_heatmap(m, sa_heatmap_group(), sa_heatmap_levels()),
    "Not clustering the features"
  )
  expect_null(out$feat_hclust)
  expect_identical(out$feat_order, 1:4)
  expect_identical(rownames(out$matrix), paste0("f", 1:4))
  # The samples still share features, so that axis is clustered as usual.
  expect_false(is.null(out$sample_hclust))
})

test_that("rows outside group_lv are dropped along with their labels", {
  skip_if_not_installed("withr")
  local_null_device()
  group <- c("ctrl", "ctrl", "other", "case", "case", "case")
  expect_message(
    out <- draw_heatmap(sa_heatmap_matrix(), group, sa_heatmap_levels(),
                        cluster_samples = FALSE),
    "Dropped 1 row"
  )
  expect_identical(ncol(out$matrix), 5L)
  expect_identical(colnames(out$matrix), c("s1", "s2", "s4", "s5", "s6"))
})

test_that("supplied labels replace the dimnames and are filtered with the rows", {
  skip_if_not_installed("withr")
  local_null_device()
  out <- sa_heatmap(feat_labels = paste0("F", 1:4),
                    sample_labels = paste0("S", 1:6),
                    cluster_feats = FALSE, cluster_samples = FALSE)
  expect_identical(rownames(out$matrix), paste0("F", 1:4))
  expect_identical(colnames(out$matrix), paste0("S", 1:6))

  # A sample label belongs to a row, so dropping the row drops the label with it.
  group <- c("ctrl", "ctrl", "other", "case", "case", "case")
  dropped <- suppressMessages(
    draw_heatmap(sa_heatmap_matrix(), group, sa_heatmap_levels(),
                 sample_labels = paste0("S", 1:6), cluster_samples = FALSE)
  )
  expect_identical(colnames(dropped$matrix), c("S1", "S2", "S4", "S5", "S6"))
})

test_that("a matrix with no dimnames is labelled rather than refused", {
  skip_if_not_installed("withr")
  local_null_device()
  m <- sa_heatmap_matrix()
  dimnames(m) <- NULL
  out <- draw_heatmap(m, sa_heatmap_group(), sa_heatmap_levels(),
                      cluster_feats = FALSE, cluster_samples = FALSE)
  expect_identical(rownames(out$matrix), paste0("V", 1:4))
  expect_identical(colnames(out$matrix), as.character(1:6))

  # Repeated sample names are a naming choice, not an error, so they survive the
  # data.frame the validator builds.
  rownames(m) <- rep("dup", 6)
  colnames(m) <- paste0("f", 1:4)
  out <- draw_heatmap(m, sa_heatmap_group(), sa_heatmap_levels(),
                      cluster_samples = FALSE)
  expect_identical(colnames(out$matrix), rep("dup", 6))
})

test_that("hiding an axis changes the drawing and nothing else", {
  skip_if_not_installed("withr")
  local_null_device()
  shown <- sa_heatmap()
  hidden <- sa_heatmap(show_feat_names = FALSE, show_sample_names = FALSE)
  expect_identical(shown$matrix, hidden$matrix)
  expect_identical(shown$feat_order, hidden$feat_order)
})

test_that("one group level is a heatmap too", {
  skip_if_not_installed("withr")
  local_null_device()
  # Nothing about the picture needs a comparison, so the two-level minimum the
  # comparison functions impose does not apply here.
  m <- sa_heatmap_matrix()[1:3, ]
  out <- draw_heatmap(m, rep("ctrl", 3), "ctrl")
  expect_identical(ncol(out$matrix), 3L)
  expect_identical(names(out$group_colors), "ctrl")
})

test_that("omitting group and group_lv drops the strip and group legend", {
  skip_if_not_installed("withr")
  local_null_device()
  out <- draw_heatmap(sa_heatmap_matrix(), cluster_samples = FALSE,
                    cluster_feats = FALSE)
  expect_identical(dim(out$matrix), c(4L, 6L))
  expect_null(out$group_colors)
  expect_error(
    draw_heatmap(sa_heatmap_matrix(), sa_heatmap_group(), NULL),
    "both be supplied or both be `NULL`"
  )
})

test_that("anno writes cell values without changing the returned matrix", {
  skip_if_not_installed("withr")
  local_null_device()
  plain <- sa_heatmap(cluster_feats = FALSE, cluster_samples = FALSE)
  labelled <- sa_heatmap(cluster_feats = FALSE, cluster_samples = FALSE,
                         anno = TRUE)
  expect_identical(plain$matrix, labelled$matrix)
})

test_that("the layout and the margins are put back", {
  skip_if_not_installed("withr")
  local_null_device()
  withr::local_par(list(mfrow = c(2, 2)))
  before <- graphics::par(no.readonly = TRUE)
  sa_heatmap(main = "restored")
  after <- graphics::par(no.readonly = TRUE)
  # A layout of its own, a narrowed device region and a panel for the key, all
  # of it temporary.
  expect_identical(after$mfrow, c(2L, 2L))
  expect_identical(after$mar, before$mar)
  expect_identical(after$oma, before$oma)
  expect_identical(after$omd, before$omd)
  expect_identical(after$fig, before$fig)
  # The key panel overlays what is already drawn, and `new` is not cleared by
  # the plot it allows, so leaving it set would stop the caller's next plot from
  # starting a page of its own.
  expect_false(after$new)
  # Nothing that carries an absolute size may be written, or the next plot is
  # pinned to the size this one was drawn at.
  expect_identical(after$fin, before$fin)
  expect_identical(after$mai, before$mai)
})

test_that("fewer than two features or samples is refused", {
  skip_if_not_installed("withr")
  local_null_device()
  # stats::heatmap() has nothing to draw a cell grid or a tree from, and a
  # single row or column is more likely a subsetting slip than a request.
  expect_error(sa_heatmap(feats = "f1"), "at least 2 features and 2 samples")
  expect_error(
    draw_heatmap(sa_heatmap_matrix()[1, , drop = FALSE], "ctrl", "ctrl"),
    "at least 2 features and 2 samples"
  )
})

test_that("the colour key is numbered inside zlim", {
  # The widest of these numbers decides where the group legend starts, so sizing
  # the strip and drawing in it have to be asking the same question.
  expect_true(all(sa_key_ticks(c(-2, 2)) >= -2 & sa_key_ticks(c(-2, 2)) <= 2))
  expect_true(all(sa_key_ticks(c(0.3, 0.7)) >= 0.3))
  expect_length(sa_key_ticks(c(-2, 2)), length(pretty(c(-2, 2), 3)))
})

test_that("the key strip is as wide as what goes in it", {
  skip_if_not_installed("withr")
  local_null_device()
  cols <- c(control = "red", case = "blue")
  narrow <- sa_key_width(c(-2, 2), cols, 0.8, 0.9)
  expect_gt(narrow, 0)
  # A longer level name asks for more of the device, which is the whole reason
  # the strip is measured rather than fixed at a fraction of it.
  longer <- stats::setNames(cols, c("control", "a case with a long name"))
  expect_gt(sa_key_width(c(-2, 2), longer, 0.8, 0.9), narrow)
  # So does a wider set of numbers beside the bar.
  expect_gt(sa_key_width(c(-1000, 1000), cols, 0.8, 0.9), narrow)
  no_group <- sa_key_width(c(-2, 2), NULL, 0.8, 0.9)
  expect_lt(no_group, narrow)
})

test_that("an axis label asks for the margin it draws in", {
  skip_if_not_installed("withr")
  local_null_device()
  # In lines of text, which is the unit stats::heatmap(margins = ) takes, and
  # measured rather than counted off the characters: a margin guessed from
  # nchar() is half again too wide for lowercase text, and every line of it that
  # is not needed reads as a gap between the cells and the key.
  expect_gt(sa_text_lines("gene_15", 1), sa_text_lines("g", 1))
  expect_gt(sa_text_lines("gene_15", 1), sa_text_lines("gene_15", 0.5))
  expect_lt(sa_text_lines("gene_15", 1), 0.6 * nchar("gene_15") + 0.6)
})

test_that("the arguments that would break the drawing are refused", {
  skip_if_not_installed("withr")
  local_null_device()
  expect_error(sa_heatmap(zlim = c(1, 1)), "two different ends")
  expect_error(sa_heatmap(zlim = c(1, -1)), "increasing")
  expect_error(sa_heatmap(n_colors = 2), "n_colors")
  expect_error(sa_heatmap(feat_labels = c("only", "two")), "feat_labels")
  expect_error(sa_heatmap(sample_labels = c("only", "two")), "sample_labels")
  expect_error(sa_heatmap(cluster_feats = NA), "cluster_feats")
  expect_error(sa_heatmap(feats = c("f1", "nope")), "nope")
})

# A 3 x 2 design in the wide layout, with the cell written into the values so
# that a box can be checked against the cell it is supposed to hold.
sa_box_frame <- function(n = 5L) {
  cells <- expand.grid(treatment = c("control", "treat_A", "treat_B"),
                       sex = c("male", "female"), rep = seq_len(n),
                       stringsAsFactors = FALSE)
  shift <- match(cells$treatment, c("control", "treat_A", "treat_B")) +
    10 * (cells$sex == "female")
  set.seed(20260813)
  data.frame(
    f1 = shift + stats::rnorm(nrow(cells), 0, 0.1),
    f2 = -shift + stats::rnorm(nrow(cells), 0, 0.1),
    treatment = cells$treatment,
    sex = cells$sex,
    stringsAsFactors = FALSE
  )
}

sa_box_factors <- function() list(treatment = "treatment", sex = "sex")

sa_box_levels <- function() {
  list(treatment = c("control", "treat_A", "treat_B"),
       sex = c("male", "female"))
}

test_that("one factor is one panel of boxes per feature", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  out <- draw_grouped_boxplot(df, c("f1", "f2"), df$treatment,
                              sa_box_levels()$treatment)
  expect_identical(names(out),
                   c("box_summary_stats", "median_confidence_stats"))
  expect_identical(names(out$box_summary_stats), c("f1", "f2"))
  expect_identical(colnames(out$box_summary_stats$f1),
                   sa_box_levels()$treatment)
  # Ten values per box, split between the two tables and nothing left over.
  expect_identical(nrow(out$box_summary_stats$f1), 7L)
  expect_identical(nrow(out$median_confidence_stats$f1), 3L)
  expect_null(draw_grouped_boxplot(df, "f1", df$treatment,
                                   sa_box_levels()$treatment,
                                   out_statistics = FALSE))
})

test_that("a crossed design is keyed on the cells the comparison fits", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  out <- draw_grouped_boxplot(df, c("f1", "f2"), factors = sa_box_factors(),
                              factor_lv = sa_box_levels())
  res <- suppressWarnings(
    compare_factorial_groups(df, "f1", sa_box_factors(), sa_box_levels(),
                             posthoc = FALSE, diagnose = FALSE)
  )
  # The same labels in the same order, so a picture and a result read against
  # each other without either being renamed.
  expect_identical(colnames(out$box_summary_stats$f1), res$design$group_lv)
  expect_identical(as.numeric(out$median_confidence_stats$f1["n", ]),
                   as.numeric(res$design$cell_n))
  # The values were built as treatment index plus ten for the women, so a median
  # says which cell the box came from.
  expect_equal(as.numeric(out$box_summary_stats$f1["median", ]),
               c(1, 2, 3, 11, 12, 13), tolerance = 0.1)
})

test_that("`control_label` draws the reference cell where the analysis puts it", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  pointed <- list(treatment = "treat_A", sex = "female")
  out <- draw_grouped_boxplot(df, c("f1", "f2"), factors = sa_box_factors(),
                              factor_lv = sa_box_levels(),
                              control_label = pointed)
  res <- suppressWarnings(
    compare_factorial_groups(df, "f1", sa_box_factors(), sa_box_levels(),
                             control_label = pointed,
                             posthoc = FALSE, diagnose = FALSE)
  )
  # The same argument means the same thing to both, so the figure and the result
  # cannot disagree about which cell the others are read against.
  expect_identical(colnames(out$box_summary_stats$f1), res$design$group_lv)
  expect_identical(colnames(out$box_summary_stats$f1)[1], "treat_A.female")

  # The boxes are the same boxes, dealt out in another order: a cell keeps the
  # values it had wherever the re-pointing moved it to.
  plain <- draw_grouped_boxplot(df, c("f1", "f2"), factors = sa_box_factors(),
                                factor_lv = sa_box_levels())
  expect_setequal(colnames(out$box_summary_stats$f1),
                  colnames(plain$box_summary_stats$f1))
  expect_identical(out$box_summary_stats$f1[["treat_A.female"]],
                   plain$box_summary_stats$f1[["treat_A.female"]])

  # One list of levels is one place for the draw order to be decided, so the
  # single-factor path says so rather than taking a second.
  expect_error(
    draw_grouped_boxplot(df, "f1", group = df$treatment,
                         group_lv = sa_box_levels()$treatment,
                         control_label = list(treatment = "treat_A")),
    "names a reference level per factor of a crossed design"
  )
})

test_that("the primary factor is the boxes and the rest are the groups", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  input <- sa_box_input(df, "f1", NULL, NULL, sa_box_factors(),
                        sa_box_levels())
  expect_identical(input$lv, sa_box_levels()$treatment)
  expect_identical(input$legend_title, "treatment")
  expect_identical(vapply(input$groups, function(g) g$label, character(1)),
                   sa_box_levels()$sex)
  # A group covers its own cells, in the level order of the primary factor, and
  # no cell is covered twice or left out.
  for (g in input$groups) {
    expect_identical(names(input$samples$f1)[g$cols],
                     paste(input$lv, g$label, sep = "."))
  }
  expect_identical(sort(unlist(lapply(input$groups, function(g) g$cols))),
                   seq_along(input$samples$f1))

  # A third factor multiplies the groups rather than the boxes. One whole
  # replicate of the 3 x 2 design per time point, so every cell is filled.
  df$time <- rep(c("T0", "T1"), each = 6L, length.out = nrow(df))
  three <- sa_box_input(df, "f1", NULL, NULL,
                        c(sa_box_factors(), list(time = "time")),
                        c(sa_box_levels(), list(time = c("T0", "T1"))))
  expect_identical(three$lv, sa_box_levels()$treatment)
  expect_identical(vapply(three$groups, function(g) g$label, character(1)),
                   c("male.T0", "female.T0", "male.T1", "female.T1"))
})

test_that("the two arrangements group the same boxes differently", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  input <- sa_box_input(df, c("f1", "f2"), NULL, NULL, sa_box_factors(),
                        sa_box_levels())

  by_feat <- sa_box_arrange(input, "feature")
  by_fac <- sa_box_arrange(input, "factor")

  # A panel per feature, the sexes the clusters inside it, and back again.
  expect_identical(vapply(by_feat$panels, function(p) p$label, character(1)),
                   c("f1", "f2"))
  expect_identical(by_feat$panels[[1]]$cluster_labels, sa_box_levels()$sex)
  expect_identical(vapply(by_fac$panels, function(p) p$label, character(1)),
                   sa_box_levels()$sex)
  expect_identical(by_fac$panels[[1]]$cluster_labels, c("f1", "f2"))

  # Every panel holds one box per cluster and level either way, which is what
  # lets one `at` layout serve both.
  n_box <- function(a) {
    vapply(a$panels, function(p) length(p$boxes), integer(1))
  }
  expect_identical(n_box(by_feat), rep(6L, 2L))   # 2 sexes x 3 treatments
  expect_identical(n_box(by_fac), rep(6L, 2L))    # 2 features x 3 treatments

  # The same boxes, regrouped: sorting the values away from the arrangement
  # leaves two identical bags.
  bag <- function(a) {
    sort(vapply(unlist(lapply(a$panels, function(p) p$boxes),
                       recursive = FALSE),
                function(v) paste(format(v), collapse = "|"), character(1)))
  }
  expect_identical(bag(by_feat), bag(by_fac))

  # Feature panels are different quantities on different baselines, so they are
  # the case that cannot share a y axis.
  expect_true(by_feat$free_scale)
  expect_false(by_fac$free_scale)

  # A single factor has one group of cells, so there is no panelling for
  # `panel_by` to choose between and the default reads as the old picture.
  one <- sa_box_input(df, c("f1", "f2"), df$treatment,
                      sa_box_levels()$treatment, NULL, NULL)
  flat <- sa_box_arrange(one, "feature")
  expect_length(flat$panels, 1L)
  expect_identical(flat$panels[[1]]$cluster_labels, c("f1", "f2"))
  expect_false(flat$free_scale)
  expect_identical(flat, sa_box_arrange(one, "factor"))
})

test_that("the arrangement does not change the numbers it returns", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  by_feat <- draw_grouped_boxplot(df, c("f1", "f2"), factors = sa_box_factors(),
                                  factor_lv = sa_box_levels())
  by_fac <- draw_grouped_boxplot(df, c("f1", "f2"), factors = sa_box_factors(),
                                 factor_lv = sa_box_levels(),
                                 panel_by = "factor")
  expect_identical(by_feat, by_fac)
  expect_error(
    draw_grouped_boxplot(df, "f1", factors = sa_box_factors(),
                         panel_by = "cell"),
    "should be one of"
  )
})

test_that("a cell holding nothing is reported and left blank", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  full <- draw_grouped_boxplot(df, "f1", factors = sa_box_factors(),
                               factor_lv = sa_box_levels())
  holed <- df[!(df$treatment == "treat_B" & df$sex == "female"), ]
  expect_message(
    out <- draw_grouped_boxplot(holed, "f1", factors = sa_box_factors(),
                                factor_lv = sa_box_levels()),
    "treat_B.female"
  )
  # The column stays, so the cells keep their positions and the missing one is
  # missing rather than absent.
  expect_identical(colnames(out$box_summary_stats$f1),
                   colnames(full$box_summary_stats$f1))
  expect_true(all(is.na(out$box_summary_stats$f1[["treat_B.female"]])))
  expect_identical(out$box_summary_stats$f1[["control.male"]],
                   full$box_summary_stats$f1[["control.male"]])
})

test_that("a level factor_lv leaves out takes its rows with it", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  expect_message(
    out <- draw_grouped_boxplot(
      df, "f1", factors = sa_box_factors(),
      factor_lv = list(treatment = c("control", "treat_A"),
                       sex = c("male", "female"))
    ),
    "Dropped 10 row"
  )
  expect_identical(colnames(out$box_summary_stats$f1),
                   c("control.male", "treat_A.male", "control.female",
                     "treat_A.female"))
})

test_that("what the boxes are has to be said once", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  expect_error(
    draw_grouped_boxplot(df, "f1", df$treatment, sa_box_levels()$treatment,
                         factors = sa_box_factors()),
    "two ways of saying"
  )
  expect_error(draw_grouped_boxplot(df, "f1"), "nothing says what the boxes")
  expect_error(draw_grouped_boxplot(df, "f1", factor_lv = sa_box_levels()),
               "`factors` holds")
  expect_error(
    draw_grouped_boxplot(df, "f1", factors = sa_box_factors()["treatment"]),
    "at least two crossed factors"
  )
  expect_error(
    draw_grouped_boxplot(df, "f1", factors = sa_box_factors(), panel_nrow = 0),
    "panel_nrow"
  )
})

test_that("the panelled plot puts the device back as it found it", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  # `oma` is in the list because a title over several panels is written there,
  # and `mfrow` because layout() overwrites the caller's grid.
  before <- graphics::par(c("mar", "mfrow", "oma", "bg", "fg"))
  draw_grouped_boxplot(df, c("f1", "f2"), factors = sa_box_factors(),
                       factor_lv = sa_box_levels(), panel_nrow = 2,
                       main = "restored", dark = TRUE)
  expect_identical(graphics::par(c("mar", "mfrow", "oma", "bg", "fg")), before)
})

test_that("a simulated crossed design draws from its own args", {
  skip_if_not_installed("withr")
  local_null_device()
  # The simulator's `args` are the comparison's argument names, and now the
  # plot's as well, so the same list draws the design and analyses it.
  sim <- simulate_factorial_groups(n_feats = 3, n_per_cell = 4, seed = 1)
  out <- draw_grouped_boxplot(data = sim$args$data, feats = sim$args$feats,
                              factors = sim$args$factors,
                              factor_lv = sim$args$factor_lv)
  truth <- subset(sim$truth_cell, features == "prot_1")
  expect_identical(as.numeric(out$median_confidence_stats$prot_1["n", ]),
                   as.numeric(truth$n))
  expect_identical(colnames(out$box_summary_stats$prot_1),
                   paste(truth$treatment, truth$sex, sep = "."))
})

test_that("the bar heights are the summary table's own column", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  lv <- sa_box_levels()$treatment
  drawn <- draw_grouped_barplot(df, c("f1", "f2"), df$treatment, lv)
  summ <- summarize_descriptive_stats(df, c("f1", "f2"), df$treatment, lv)

  # Nothing is computed twice for the picture, so a bar is a row of that table
  # and the two cannot disagree.
  expect_identical(drawn$features, summ$features)
  expect_identical(drawn$group, as.character(summ$group))
  expect_identical(drawn$value, summ$mean)
  expect_identical(drawn$n, summ$n)
  expect_identical(attr(drawn, "mainbar"), "mean")
  expect_identical(attr(drawn, "errorbar"), "none")
  # No interval was asked for, so the bound columns say so rather than repeating
  # the height and drawing a bar of zero width.
  expect_true(all(is.na(drawn$lower)))
  expect_true(all(is.na(drawn$upper)))
  expect_null(draw_grouped_barplot(df, "f1", df$treatment, lv,
                                   out_statistics = FALSE))
})

test_that("`se` and `sd` are one of each either side of the mean", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  lv <- sa_box_levels()$treatment
  summ <- summarize_descriptive_stats(df, "f1", df$treatment, lv)

  se <- draw_grouped_barplot(df, "f1", df$treatment, lv, errorbar = "se")
  expect_equal(se$upper - se$value, summ$se)
  expect_equal(se$value - se$lower, summ$se)

  spread <- draw_grouped_barplot(df, "f1", df$treatment, lv, errorbar = "sd")
  expect_equal(spread$upper - spread$value, summ$sd)
  # The two are not the same width, which is the whole reason both are offered.
  expect_true(all(spread$upper > se$upper))
})

test_that("`ci` is Student's interval at the level it was given", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  lv <- sa_box_levels()$treatment
  summ <- summarize_descriptive_stats(df, "f1", df$treatment, lv)

  ci <- draw_grouped_barplot(df, "f1", df$treatment, lv, errorbar = "ci")
  expect_equal(ci$upper - ci$value, stats::qt(0.975, summ$n - 1) * summ$se)
  # Widening the level is the only thing `conf_level` is allowed to change.
  wider <- draw_grouped_barplot(df, "f1", df$treatment, lv, errorbar = "ci",
                                conf_level = 0.99)
  expect_true(all(wider$upper > ci$upper))
  expect_equal(wider$value, ci$value)
})

test_that("a median's interval is the notch the boxplot draws", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  lv <- sa_box_levels()$treatment
  drawn <- draw_grouped_barplot(df, "f1", df$treatment, lv,
                                mainbar = "median", errorbar = "ci")
  # Same interval, so a bar and the notch of the box beside it are the same
  # width on the same data rather than two conventions for one quantity.
  notch <- draw_grouped_boxplot(df, "f1", df$treatment,
                                lv)$median_confidence_stats$f1
  expect_equal(drawn$lower, as.numeric(notch["lower_conf", ]))
  expect_equal(drawn$upper, as.numeric(notch["upper_conf", ]))
})

test_that("an interval the height has no answer for is refused", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  lv <- sa_box_levels()$treatment
  expect_error(draw_grouped_barplot(df, "f1", df$treatment, lv,
                                    mainbar = "median", errorbar = "se"),
               "not a width to draw either side of a median")
  expect_error(draw_grouped_barplot(df, "f1", df$treatment, lv,
                                    mainbar = "n", errorbar = "ci"),
               "takes errorbar = \"none\"")
  # The pair is settled before `data` is read, so a combination that could only
  # mislead fails on its own rather than after a summary has been computed.
  expect_error(draw_grouped_barplot(df, "nope", df$treatment, lv,
                                    mainbar = "sd", errorbar = "se"),
               "itself a spread")
})

test_that("`group` is what says which bars there are", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  expect_error(draw_grouped_barplot(df, "f1"), "`group` says which bars")
  # One level is one bar per feature, which is a summary rather than a
  # comparison of anything.
  expect_error(draw_grouped_barplot(df, "f1", df$treatment, "control"),
               "at least 2 levels")
})

test_that("unnamed levels are taken the way the summary takes them", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  ordered <- c("treat_B", "control", "treat_A")
  df$treatment <- factor(df$treatment, levels = ordered)
  # A factor's own order is the draw order, so a call that names none draws what
  # summarize_descriptive_stats() reports.
  expect_identical(unique(draw_grouped_barplot(df, "f1", df$treatment)$group),
                   ordered)
  expect_identical(
    unique(draw_grouped_barplot(df, "f1", as.character(df$treatment))$group),
    sort(ordered)
  )
})

test_that("`control_label` draws the reference level first", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  lv <- sa_box_levels()$treatment
  pointed <- draw_grouped_barplot(df, "f1", df$treatment, lv,
                                  control_label = "treat_A")
  expect_identical(unique(pointed$group), c("treat_A", "control", "treat_B"))
  # The bars are the same bars, dealt out in another order: a level keeps the
  # value it had wherever the re-pointing moved it to.
  plain <- draw_grouped_barplot(df, "f1", df$treatment, lv)
  expect_setequal(pointed$group, plain$group)
  expect_identical(pointed$value[pointed$group == "treat_A"],
                   plain$value[plain$group == "treat_A"])
})

test_that("a level left out of `group_lv` drops its rows and says so once", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  expect_message(
    drawn <- draw_grouped_barplot(df, "f1", df$treatment,
                                  c("treat_B", "control")),
    "Dropped 10 row"
  )
  expect_identical(unique(drawn$group), c("treat_B", "control"))
})

test_that("a bar the summary could not compute is left blank", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  lv <- sa_box_levels()$treatment
  # A shape estimate needs three observations, so a group of two has no height.
  # The bar goes blank and the others still draw, rather than the call failing.
  thin <- df[sort(c(which(df$treatment == "control")[1:2],
                    which(df$treatment != "control"))), ]
  drawn <- draw_grouped_barplot(thin, "f1", thin$treatment, lv,
                                mainbar = "skewness")
  expect_true(all(is.na(drawn$value[drawn$group == "control"])))
  expect_true(all(is.finite(drawn$value[drawn$group != "control"])))

  # Every group too small is a different case: there is no bar left to give the
  # panel a height, so it is an error rather than an empty picture.
  every <- df[unlist(lapply(lv, function(l) which(df$treatment == l)[1:2])), ]
  expect_error(draw_grouped_barplot(every, "f1", every$treatment, lv,
                                    mainbar = "skewness"),
               "NA for every feature and group")
})

test_that("the barplot puts the device back as it found it", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  lv <- sa_box_levels()$treatment
  # layout() overwrites the caller's grid while the legend panel is drawn, so
  # the grid has to come back rather than be left at the one panel layout(1)
  # resets to.
  withr::local_par(list(mfrow = c(2, 2)))
  before <- graphics::par(no.readonly = TRUE)
  draw_grouped_barplot(df, c("f1", "f2"), df$treatment, lv, dark = TRUE,
                       errorbar = "se", main = "restored")
  after <- graphics::par(no.readonly = TRUE)
  # `fin`, `pin` and `mai` carry absolute sizes: writing them back would pin the
  # next plot to the size this one happened to be drawn at.
  for (nm in c("bg", "fg", "col.axis", "col.lab", "mar", "mfrow", "oma",
               "fin", "pin", "mai")) {
    expect_identical(after[[nm]], before[[nm]])
  }
  expect_identical(graphics::par("mfrow"), c(2L, 2L))
})

test_that("a range and a palette of the caller's own are drawn and checked", {
  skip_if_not_installed("withr")
  local_null_device()
  df <- sa_box_frame()
  lv <- sa_box_levels()$treatment
  # Two colours for three levels is a recycled palette, and a range narrower
  # than the bars clips them rather than being widened back out.
  drawn <- draw_grouped_barplot(df, c("f1", "f2"), df$treatment, lv,
                                errorbar = "sd", col = c("red", "blue"),
                                ylim = c(-5, 5), xlab = "protein", ylab = "",
                                main = "own scale")
  expect_identical(nrow(drawn), 6L)

  expect_error(draw_grouped_barplot(df, "f1", df$treatment, lv,
                                    ylim = c(0, 1, 2)), "ylim")
  expect_error(draw_grouped_barplot(df, "f1", df$treatment, lv, gap = -1),
               "gap")
  expect_error(draw_grouped_barplot(df, "f1", df$treatment, lv, lwd = 0),
               "lwd")
  expect_error(draw_grouped_barplot(df, "f1", df$treatment, lv,
                                    control_label = "missing"),
               "names a level")
  expect_error(draw_grouped_barplot(df, "f1", df$treatment, lv,
                                    conf_level = 1), "conf_level")
  expect_error(draw_grouped_barplot(df, "f1", df$treatment, lv, dark = NA),
               "dark")
})

test_that("a categorical result draws a mosaic of its own cells", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_categorical_fixture()
  drawn <- draw_mosaic_plot(res)
  expect_identical(drawn$cells$row_level, res$cells$row_level)
  expect_identical(drawn$cells$col_level, res$cells$col_level)
  expect_equal(sum(drawn$widths), 1, tolerance = 1e-12)
  expect_identical(names(drawn$widths), unique(res$cells$row_level))
  expect_identical(rownames(drawn$heights), unique(res$cells$row_level))
  expect_identical(colnames(drawn$heights), unique(res$cells$col_level))
  expect_identical(drawn$null, "independence")
  expect_identical(drawn$residual, "pearson")
})

test_that("plot() on a categorical result is draw_mosaic_plot()", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_categorical_fixture()
  expect_identical(plot(res)$cells$observed,
                   draw_mosaic_plot(res)$cells$observed)
})

test_that("the mosaic tiles follow category_lv", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- compare_categorical_groups(
    sa_categorical_frame(),
    category_lv = list(smoker = c("n", "y"), grade = c("low", "mid", "high")),
    diagnose = FALSE
  )
  drawn <- draw_mosaic_plot(res)
  expect_identical(names(drawn$widths), c("n", "y"))
  expect_identical(colnames(drawn$heights), c("low", "mid", "high"))
})

test_that("a strip is as wide as its level's share of the table", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_categorical_fixture()
  drawn <- draw_mosaic_plot(res, gap = 0)
  by_strip <- unique(drawn$cells[c("row_level", "x1", "x2")])
  share <- table(sa_categorical_frame()$smoker) / 120

  expect_equal(by_strip$x2 - by_strip$x1,
               as.numeric(share[by_strip$row_level]))
  # Which is the whole claim of the picture: area is the cell's share of the
  # table, so with no gaps taken out the tiles fill the unit square.
  area <- with(drawn$cells, (x2 - x1) * (y2 - y1))
  expect_equal(sum(area), 1)
  expect_equal(area, res$cells$prop_total)
})

test_that("gap is the width of one gap rather than of all of them", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_categorical_fixture()
  drawn <- draw_mosaic_plot(res, gap = 0.05)
  by_strip <- unique(drawn$cells[c("row_level", "x1", "x2")])

  # Two levels on x, three on y, so one gap across and two up.
  expect_equal(by_strip$x1[2] - by_strip$x2[1], 0.05)
  expect_equal(sum(by_strip$x2 - by_strip$x1), 1 - 0.05)
  area <- with(drawn$cells, (x2 - x1) * (y2 - y1))
  expect_equal(sum(area), (1 - 0.05) * (1 - 2 * 0.05))
})

test_that("many levels cannot be gapped out of existence", {
  skip_if_not_installed("withr")
  local_null_device()
  grid <- expand.grid(g = letters[1:6], h = c("w", "x", "y", "z"),
                      stringsAsFactors = FALSE)
  wide <- compare_categorical_groups(grid[rep(seq_len(nrow(grid)), each = 4), ],
                                     diagnose = FALSE)
  drawn <- draw_mosaic_plot(wide, gap = 0.2)
  area <- with(drawn$cells, (x2 - x1) * (y2 - y1))

  # The cap leaves three fifths of each axis to the tiles however many levels
  # ask for a gap.
  expect_equal(sum(area), 0.6 * 0.6)
  expect_true(all(area > 0))
})

test_that("the level names on the y axis sit on the reference strip", {
  res <- sa_categorical_fixture()
  drawn <- sa_mosaic_layout(res$cells, 0.015, res$design$null)
  first <- drawn$cells[drawn$cells$row_level == drawn$row_lv[1], ]

  expect_identical(names(drawn$y_at), drawn$col_lv)
  expect_equal(unname(drawn$y_at), (first$y1 + first$y2) / 2)
  # And not on the strip that happens to be widest, which is what would leave
  # the names beside a tile they do not belong to.
  second <- drawn$cells[drawn$cells$row_level == drawn$row_lv[2], ]
  expect_false(isTRUE(all.equal(unname(drawn$y_at),
                                (second$y1 + second$y2) / 2)))
})

test_that("independence expects the same cut in every strip", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_categorical_fixture()
  drawn <- draw_mosaic_plot(res)
  marginal <- colSums(as.table(res)) / res$design$n_used

  expect_identical(dim(drawn$expected_prop), c(2L, 3L))
  for (i in seq_len(nrow(drawn$expected_prop))) {
    expect_equal(drawn$expected_prop[i, ], marginal)
  }
})

test_that("symmetry expects a different cut in each strip", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_categorical_matched_fixture()
  drawn <- draw_mosaic_plot(res)

  expect_identical(drawn$null, "symmetry")
  expect_false(isTRUE(all.equal(drawn$expected_prop[1, ],
                                drawn$expected_prop[2, ])))
  # The expected table under symmetry is its own transpose, so the shares are
  # read off that rather than off the observed margins.
  counts <- unclass(as.table(res))
  expected <- (counts + t(counts)) / 2
  expect_equal(unname(drawn$expected_prop),
               unname(expected / rowSums(expected)))
})

test_that("the standardized residual is refused where it has no value", {
  skip_if_not_installed("withr")
  local_null_device()
  expect_error(
    draw_mosaic_plot(sa_categorical_matched_fixture(),
                     residual = "standardized"),
    "no value under symmetry"
  )
  expect_silent(draw_mosaic_plot(sa_categorical_fixture(),
                                 residual = "standardized"))
})

test_that("a level holding nothing is reported rather than labelled", {
  # Laid out directly: `compare_categorical_groups()` takes its levels from the
  # data, so a level holding no observation is not a table it can be asked for.
  hollow <- as.table(matrix(c(10, 0, 5, 0), nrow = 2,
                            dimnames = list(g = c("a", "b"),
                                            h = c("x", "y"))))
  drawn <- sa_mosaic_layout(sa_categorical_cells(hollow, "independence"),
                            0.015, "independence")

  expect_identical(drawn$empty_levels$row, "b")
  expect_identical(drawn$empty_levels$col, character(0))
  expect_equal(unname(drawn$widths), c(1, 0))
  # A zero-width strip draws no tile, and the gap where it should have been is
  # left standing.
  hollow_tiles <- drawn$cells[drawn$cells$row_level == "b", ]
  expect_true(all(hollow_tiles$x2 - hollow_tiles$x1 == 0))
})

test_that("shading off leaves one fill and no residual to read", {
  skip_if_not_installed("withr")
  local_null_device()
  drawn <- draw_mosaic_plot(sa_categorical_fixture(), shade = FALSE)
  expect_identical(unique(drawn$cells$fill), "white")
  expect_identical(unique(draw_mosaic_plot(sa_categorical_fixture(),
                                           shade = FALSE, dark = TRUE)$cells$fill),
                   "#36454F")
})

test_that("anything but a categorical result is pointed at the comparison", {
  # The null the shading is read under lives on the result, so an input that
  # carries no null is refused rather than held against a guessed one.
  inputs <- list(
    data.frame  = sa_categorical_frame(),
    two_way     = as.table(sa_categorical_fixture()),
    one_way     = table(c("a", "b", "a")),
    plain_matrix = matrix(1:4, 2)
  )
  for (nm in names(inputs)) {
    expect_error(draw_mosaic_plot(inputs[[nm]]),
                 "must be a categorical comparison result", info = nm)
  }
})

test_that("a numeric comparison is pointed at the plots that read it", {
  expect_error(draw_mosaic_plot(sa_two_group_fixture()),
               "numeric comparison")
})

test_that("the reference the comparison settled is the one drawn", {
  skip_if_not_installed("withr")
  local_null_device()
  flipped <- draw_mosaic_plot(
    sa_categorical_fixture(control_label = c(smoker = "y"))
  )
  expect_identical(names(flipped$widths), c("y", "n"))
})


# --- draw_interaction_plot() -------------------------------------------------
# The lines are not compared against a reference image either. What is tested is
# which means the function decided each point was, which panel it put them in,
# and that they are the same numbers the comparison already reported.

sa_inter_fixture <- function() {
  # A third factor so that the pairwise view has something to average away, and
  # a crossover on the last level so the interaction term is real.
  set.seed(20260815)
  lv <- expand.grid(pressure = c("low", "med", "high"), paint = c("f1", "f2"),
                    line = c("l1", "l2"), rep = 1:4, stringsAsFactors = FALSE)
  y <- c(low = 0, med = 9, high = 25)[lv$pressure] +
    ifelse(lv$paint == "f2", 7, 0) +
    ifelse(lv$paint == "f2" & lv$pressure == "high", -9, 0) +
    ifelse(lv$line == "l2", 1, 0) +
    stats::rnorm(nrow(lv), sd = 2)
  compare_factorial_groups(
    data.frame(flaws = y, gloss = -y),
    c("flaws", "gloss"),
    factors   = list(paint = lv$paint, pressure = lv$pressure,
                     line = lv$line),
    factor_lv = list(paint = c("f1", "f2"),
                     pressure = c("low", "med", "high"),
                     line = c("l1", "l2")),
    diagnose = FALSE
  )
}

test_that("the default view is one pair of factors and a panel per feature", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_factorial_fixture()
  drawn <- draw_interaction_plot(res)

  expect_identical(attr(drawn, "view"), "pairwise")
  # Unnamed, the factors go in declaration order, the first one tracing.
  expect_identical(unique(drawn$trace_factor), "wool")
  expect_identical(unique(drawn$x_factor), "tension")
  # One row per point, x varying fastest so the rows of one trace are together.
  expect_identical(nrow(drawn), 6L)
  expect_identical(drawn$x_level, rep(res$design$factor_lv$tension, 2L))
  expect_identical(drawn$trace_level, rep(res$design$factor_lv$wool, each = 3L))
  # Two factors leave nothing to average, so every point is one cell.
  expect_identical(unique(drawn$n_cells), 1L)
  # And no bars asked for is no bars reported, rather than bars of no width.
  expect_true(all(is.na(drawn$lower_conf)))
})

test_that("every point is a marginal mean of the cell table", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_inter_fixture()
  drawn <- draw_interaction_plot(res, x = "pressure", trace = "paint")

  # Three factors and two drawn, so each point averages the two levels of the
  # third: the plot says so on `n_cells` and the means have to follow.
  expect_identical(unique(drawn$n_cells), 2L)
  cells <- res$cells
  for (i in seq_len(nrow(drawn))) {
    at <- cells$features == drawn$features[i] &
      cells$paint == drawn$trace_level[i] &
      cells$pressure == drawn$x_level[i]
    expect_equal(drawn$mean[i], mean(cells$mean[at]))
    expect_equal(drawn$se[i], sqrt(sum(cells$se[at]^2) / sum(at)^2))
  }
})

test_that("a gap along a line is the contrast the comparison reports", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_inter_fixture()
  drawn <- draw_interaction_plot(res, x = "pressure", trace = "paint")

  ph <- res$posthoc$anova_test
  ph <- ph[is.na(ph$stratum) & ph$factor %in% c("paint", "pressure"), ]
  expect_gt(nrow(ph), 0L)
  for (i in seq_len(nrow(ph))) {
    pts <- drawn[drawn$features == ph$features[i], ]
    side <- if (ph$factor[i] == "paint") "trace_level" else "x_level"
    expect_equal(mean(pts$mean[pts[[side]] == ph$group1[i]]) -
                   mean(pts$mean[pts[[side]] == ph$group2[i]]),
                 ph$estimate[i], info = ph$contrast[i])
  }
})

test_that("the matrix view draws the upper triangle of factor pairs", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_inter_fixture()
  # Three factors and nothing named, so "auto" shows every pair rather than
  # picking one silently, and says which feature it settled on.
  expect_message(drawn <- draw_interaction_plot(res), "one feature at a time")

  expect_identical(attr(drawn, "view"), "matrix")
  expect_identical(unique(drawn$features), res$features[1])
  # A pair is drawn once, the earlier factor tracing the later one.
  expect_identical(unique(drawn$panel),
                   c("paint x pressure", "paint x line", "pressure x line"))
  expect_identical(unique(paste(drawn$trace_factor, drawn$x_factor)),
                   c("paint pressure", "paint line", "pressure line"))
})

test_that("the facet view keeps a factor apart instead of averaging it", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_inter_fixture()
  drawn <- draw_interaction_plot(res, x = "pressure", trace = "paint",
                                 facet = "line", feats = "gloss")

  expect_identical(attr(drawn, "view"), "facet")
  expect_identical(unique(drawn$panel), c("line: l1", "line: l2"))
  # All three factors are spoken for, so nothing is left to average away.
  expect_identical(unique(drawn$n_cells), 1L)
  expect_identical(unique(drawn$features), "gloss")
  # Each panel holds one level of the facet, so its points are that level's
  # cells and no others.
  cells <- res$cells[res$cells$features == "gloss", ]
  first <- drawn[drawn$panel == "line: l1", ]
  expect_equal(sort(first$mean), sort(cells$mean[cells$line == "l1"]))
})

test_that("an error bar is one standard error or an interval, as asked", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_inter_fixture()
  se <- draw_interaction_plot(res, x = "pressure", trace = "paint",
                              errorbar = "se")
  ci <- draw_interaction_plot(res, x = "pressure", trace = "paint",
                              errorbar = "ci")

  expect_equal(se$upper_conf - se$mean, se$se)
  # The interval is the same standard error times the multiplier the error
  # degrees of freedom of that feature's own fit give it.
  tbl <- res$tests$anova_test
  df <- tbl$df2[match(ci$features, tbl$features)]
  mult <- stats::qt(1 - (1 - res$parameters$conf_level) / 2, df)
  expect_equal(ci$upper_conf - ci$mean, mult * ci$se)
  expect_true(all(ci$upper_conf > se$upper_conf))
})

test_that("the interaction plot leaves par and the layout as it found them", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_inter_fixture()
  before <- graphics::par(c("mar", "mfrow", "oma", "bg", "fg"))
  invisible(draw_interaction_plot(res, x = "pressure", trace = "paint"))
  expect_identical(graphics::par(c("mar", "mfrow", "oma", "bg", "fg")), before)
  invisible(draw_interaction_plot(res, type = "matrix", dark = TRUE))
  expect_identical(graphics::par(c("mar", "mfrow", "oma", "bg", "fg")), before)
})

test_that("a result with no second factor is pointed at the right function", {
  expect_error(draw_interaction_plot(sa_multi_group_fixture()),
               "factorial comparison result")
  # A result from before the cell means were recorded has to say so rather than
  # fail on a missing column.
  res <- sa_factorial_fixture()
  res$cells <- NULL
  expect_error(draw_interaction_plot(res), "carries no .\\$cells. table")
})

test_that("the factors have to be named once each and to exist", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_inter_fixture()
  expect_error(draw_interaction_plot(res, x = "nope"), "Got nope")
  expect_error(draw_interaction_plot(res, x = c("paint", "line")),
               "single factor name")
  expect_error(draw_interaction_plot(res, x = "paint", trace = "paint"),
               "named twice")
  expect_error(draw_interaction_plot(res, type = "matrix", x = "paint"),
               "nothing for x to choose")
  expect_error(draw_interaction_plot(res, type = "facet", x = "paint"),
               "needs `facet`")
  expect_error(draw_interaction_plot(sa_factorial_fixture(), facet = "wool"),
               "needs at least three factors")
})

test_that("a view whose panels are spent on factors draws one feature", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_inter_fixture()
  expect_error(draw_interaction_plot(res, type = "matrix",
                                     feats = res$features),
               "one feature at a time")
  expect_error(draw_interaction_plot(res, feats = "nope"), "Not found: nope")
  # The pairwise view has a panel per feature, so it takes them all, and
  # `panel_nrow` is the only thing that decides how they are stacked.
  drawn <- draw_interaction_plot(res, x = "pressure", trace = "paint",
                                 panel_nrow = 2)
  expect_identical(unique(drawn$panel), res$features)
  expect_message(draw_interaction_plot(res, type = "matrix", panel_nrow = 2),
                 "not used by the matrix view")
})

# --- the two evaluation plots -------------------------------------------------
# Both draw from the tables the evaluation already holds rather than recomputing
# anything, so what is tested is which rows they decided to draw and that the
# device is left as it was found.

test_that("a regression evaluation panels once there is more than one model", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_perf_reg_fixture()
  drawn <- draw_prediction_plot(res)
  expect_identical(drawn$model, res$models)
  # Two clouds of points on shared axes make a third that belongs to neither,
  # which is why "auto" does not overlay past one model.
  expect_identical(attr(drawn, "view"), "panel")
})

test_that("a single model is overlaid rather than given a grid of one", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_perf_reg_fixture()
  drawn <- draw_prediction_plot(res, models = "wt_only")
  expect_identical(drawn$model, "wt_only")
  expect_identical(attr(drawn, "view"), "overlay")
  # And the view the caller asked for outright is the view they get.
  expect_identical(attr(draw_prediction_plot(res, type = "overlay"), "view"),
                   "overlay")
  expect_identical(attr(draw_prediction_plot(res, models = "wt_only",
                                             type = "panel"), "view"),
                   "panel")
})

test_that("the drawn rows are the metric rows, selected and ordered", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_perf_reg_fixture()
  drawn <- draw_prediction_plot(res, models = c("qsec_only", "baseline"))
  expect_identical(drawn$model, c("qsec_only", "baseline"))
  expect_equal(drawn$rmse,
               res$metrics$rmse[match(c("qsec_only", "baseline"),
                                      res$metrics$model)])
})

test_that("a ROC plot draws every model and reports their AUCs back", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_perf_cls_fixture()
  drawn <- draw_roc_curve(res, anno_auc = TRUE)
  expect_identical(drawn$model, res$models)
  expect_equal(drawn$auc, res$metrics$auc)
  expect_identical(draw_roc_curve(res, models = "width_only")$model,
                   "width_only")
  expect_null(attr(drawn, "view"))
})

test_that("the AUC legend can be sized apart from the model-name legend", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_perf_cls_fixture()
  expect_silent(draw_roc_curve(res, anno_auc = TRUE, cex.anno = 0.9))
  expect_silent(draw_roc_curve(res, cex.legend = 1.2))
  expect_error(draw_roc_curve(res, cex.anno = 0),
               "`cex.anno` must be in \\(0, Inf\\]")
})

test_that("plot() sends each evaluation to its own picture", {
  skip_if_not_installed("withr")
  local_null_device()
  reg <- sa_perf_reg_fixture()
  cls <- sa_perf_cls_fixture()
  # The method carries no logic of its own, so the two entry points cannot draw
  # different rows.
  expect_equal(plot(reg, models = "baseline"),
               draw_prediction_plot(reg, models = "baseline"))
  expect_equal(plot(cls, anno_auc = TRUE), draw_roc_curve(cls, anno_auc = TRUE))
})

test_that("each plot refuses the other one's result by name", {
  skip_if_not_installed("withr")
  local_null_device()
  expect_error(draw_roc_curve(sa_perf_reg_fixture()),
               "Use draw_prediction_plot\\(\\)")
  expect_error(draw_prediction_plot(sa_perf_cls_fixture()),
               "Use draw_roc_curve\\(\\)")
  # A model is one step short of a result, and the message says which step.
  expect_error(draw_roc_curve(sa_perf_cls_parts()$models$both),
               "Score it on held-out rows")
  expect_error(draw_prediction_plot(mtcars), "must be an evaluation result")
})

test_that("the models and the styling are checked before anything is drawn", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_perf_cls_fixture()
  expect_error(draw_roc_curve(res, models = "nope"), "width_only")
  expect_error(draw_roc_curve(res, models = c("baseline", "baseline")),
               "duplicated names")
  expect_error(draw_roc_curve(res, col = c("red", "blue")), "one per drawn")
  expect_error(draw_roc_curve(res, lty = c(1, 2)), "one per drawn")
  expect_error(draw_prediction_plot(sa_perf_reg_fixture(), panel_nrow = 0),
               "`panel_nrow` must be in")
})

test_that("both evaluation plots leave the device as they found it", {
  skip_if_not_installed("withr")
  local_null_device()
  before <- graphics::par(no.readonly = TRUE)
  draw_prediction_plot(sa_perf_reg_fixture(), panel_nrow = 2, dark = TRUE,
                       main = "restored", anno_lm = TRUE)
  after <- graphics::par(no.readonly = TRUE)
  expect_identical(after$bg, before$bg)
  expect_identical(after$fg, before$fg)
  expect_identical(after$mar, before$mar)
  expect_identical(after$mfrow, before$mfrow)
  # A figure title lives in the outer margin, which is the one this plot writes
  # that the forest plot does not.
  expect_identical(after$oma, before$oma)
  # Nothing that carries an absolute size may be written, or the next plot is
  # pinned to the size this one was drawn at.
  expect_identical(after$fin, before$fin)
  expect_identical(after$mai, before$mai)

  draw_roc_curve(sa_perf_cls_fixture(), dark = TRUE)
  after <- graphics::par(no.readonly = TRUE)
  expect_identical(after$bg, before$bg)
  expect_identical(after$mar, before$mar)
  expect_identical(after$fin, before$fin)
})

test_that("a panel grid the caller set up survives both calls", {
  skip_if_not_installed("withr")
  local_null_device()
  withr::local_par(list(mfrow = c(2, 2)))
  draw_prediction_plot(sa_perf_reg_fixture(), type = "overlay")
  expect_identical(graphics::par("mfrow"), c(2L, 2L))
  draw_roc_curve(sa_perf_cls_fixture())
  expect_identical(graphics::par("mfrow"), c(2L, 2L))
})

test_that("a model with no calibration line is drawn without one", {
  skip_if_not_installed("withr")
  local_null_device()
  # An outcome that does not vary leaves the slope NA, and a plot that read it
  # as a number would draw a line through nothing rather than leave it out.
  parts <- sa_perf_reg_parts()
  flat <- parts$test
  flat$mpg <- 20
  suppressWarnings(
    res <- evaluate_regression_models(parts$models$full, newdata = flat)
  )
  expect_true(is.na(res$metrics$calib_slope))
  expect_silent(drawn <- draw_prediction_plot(res, anno_lm = TRUE))
  expect_identical(drawn$model, "baseline")
})

test_that("held-out R-squared can be annotated beside correlation", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_perf_reg_fixture()
  expect_silent(draw_prediction_plot(res, anno_corr = TRUE, anno_rsq = TRUE))
  expect_silent(draw_prediction_plot(res, anno_lm = TRUE, cex.anno = 0.9))
  expect_error(draw_prediction_plot(res, anno_rsq = "yes"),
               "`anno_rsq` must be TRUE or FALSE")
  expect_error(draw_prediction_plot(res, cex.anno = 0),
               "`cex.anno` must be in \\(0, Inf\\]")
})


# draw_dim_reduction_plot() ----------------------------------------------------
#
# The reduction and the clustering are of the same frame on purpose: the whole
# claim of this picture is that a label from one can be painted onto a coordinate
# from the other without either being asked where its rows came from.

sa_scatter_parts <- function() {
  m <- sa_cluster_matrix()
  list(
    m       = m,
    pca     = perform_pca(m),
    kmeans  = cluster_kmeans(m, n_clust = 3, seed = 1),
    group   = factor(paste0("g", sa_cluster_truth()))
  )
}

test_that("the drawn table is the reduction's points and coordinates", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()
  drawn <- draw_dim_reduction_plot(parts$pca)

  expect_identical(drawn$points, parts$pca$points)
  expect_identical(drawn$x, parts$pca$scores$PC1)
  expect_identical(drawn$y, parts$pca$scores$PC2)
  expect_identical(attr(drawn, "view"), "plain")
  # With neither channel used the two columns are absent rather than filled
  # with a placeholder.
  expect_false("cluster" %in% names(drawn))
  expect_false("group" %in% names(drawn))
  expect_identical(length(unique(drawn$col)), 1L)
  expect_identical(unique(drawn$pch), 16L)
})

test_that("colour is the clustering and shape is the group", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()
  drawn <- draw_dim_reduction_plot(parts$pca, group = parts$group,
                                   cluster_result = parts$kmeans)

  expect_identical(attr(drawn, "view"), "both")
  expect_identical(drawn$cluster, parts$kmeans$assignments$cluster)
  expect_identical(levels(drawn$group), levels(parts$group))
  # Each channel is a function of its own source and of nothing else, which is
  # what lets the two be read at once.
  expect_identical(length(unique(drawn$col)), 3L)
  expect_identical(length(unique(drawn$pch)), 3L)
  expect_true(all(tapply(drawn$col, drawn$cluster,
                         function(v) length(unique(v))) == 1L))
  expect_true(all(tapply(drawn$pch, drawn$group,
                         function(v) length(unique(v))) == 1L))
})

test_that("one channel on its own takes the colours", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()

  by_clust <- draw_dim_reduction_plot(parts$pca, cluster_result = parts$kmeans)
  expect_identical(attr(by_clust, "view"), "cluster")
  expect_identical(length(unique(by_clust$col)), 3L)
  expect_identical(unique(by_clust$pch), 16L)

  by_group <- draw_dim_reduction_plot(parts$pca, group = parts$group)
  expect_identical(attr(by_group, "view"), "group")
  expect_identical(length(unique(by_group$col)), 1L)
  expect_identical(length(unique(by_group$pch)), 3L)

  by_group_col <- draw_dim_reduction_plot(parts$pca, group = parts$group,
                                          col = c("red", "green", "blue"))
  expect_identical(length(unique(by_group_col$col)), 3L)
  expect_true(all(tapply(by_group_col$col, by_group_col$group,
                         function(v) length(unique(v))) == 1L))
})

test_that("`col` and `pch` name the group channel when it stands alone", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()
  cols <- c("red", "green", "blue")
  shapes <- c(16L, 17L, 15L)
  drawn <- draw_dim_reduction_plot(parts$pca, group = parts$group,
                                   col = cols, pch = shapes)
  resolved <- sa_scatter_groups(parts$group, NULL, parts$pca$points,
                                parts$pca$design, cols, shapes)
  expect_identical(drawn$pch, resolved$pch)
  expect_identical(drawn$col, resolved$col)
  expect_identical(levels(drawn$group), resolved$levels)
  expect_identical(resolved$pch_lv, shapes)
  expect_identical(resolved$palette, cols)
  expect_error(
    draw_dim_reduction_plot(parts$pca, group = parts$group, col = cols[1:2]),
    "one per group level"
  )
  expect_error(
    draw_dim_reduction_plot(parts$pca, group = parts$group, pch = shapes[1:2]),
    "one per group level"
  )
  expect_error(
    draw_dim_reduction_plot(parts$pca, pch = c(16, 17)),
    "one per point"
  )
})

test_that("`col` stays with the clustering when both channels are given", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()
  cols <- c("red", "green", "blue")
  drawn <- draw_dim_reduction_plot(parts$pca, group = parts$group,
                                   cluster_result = parts$kmeans,
                                   col = cols)
  expect_identical(drawn$col, cols[drawn$cluster])
})

test_that("noise is grey rather than a colour of its own", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()
  # A radius small enough to strand points, so there is noise to look at.
  db <- suppressMessages(cluster_dbscan(parts$m, eps = 0.2, min_pts = 3))
  expect_gt(db$design$n_noise, 0L)

  drawn <- draw_dim_reduction_plot(parts$pca, cluster_result = db)
  noise <- unique(drawn$col[drawn$cluster == 0L])
  expect_length(noise, 1L)
  expect_false(any(noise %in% unique(drawn$col[drawn$cluster > 0L])))
  # `col` names the clusters, so it has nothing to say about the points that
  # joined none of them.
  forced <- draw_dim_reduction_plot(parts$pca, cluster_result = db,
                                    col = "red")
  expect_identical(unique(forced$col[forced$cluster > 0L]), "red")
  expect_false(identical(unique(forced$col[forced$cluster == 0L]), "red"))
})

test_that("only a rotation puts a share of the variance on its axes", {
  parts <- sa_scatter_parts()
  space <- sa_reduction_scatter(parts$pca, c(1L, 2L))
  expect_identical(
    space$xlab,
    paste0("PC1 (", sa_fmt_num(parts$pca$variance$prop_var[1], 3), "%)")
  )

  # An embedding carries no `$variance`, so its axes are the names alone.
  embedding <- parts$pca
  embedding$variance <- NULL
  expect_identical(sa_reduction_scatter(embedding, c(1L, 2L))$xlab, "PC1")
})

test_that("`dims` chooses which two coordinates are drawn", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()
  drawn <- draw_dim_reduction_plot(parts$pca, dims = c(2, 3))

  expect_identical(drawn$x, parts$pca$scores$PC2)
  expect_identical(drawn$y, parts$pca$scores$PC3)
  expect_true(startsWith(sa_reduction_scatter(parts$pca, c(2L, 3L))$xlab,
                         "PC2"))

  expect_error(draw_dim_reduction_plot(parts$pca, dims = c(2, 2)),
               "names 2 twice")
  expect_error(draw_dim_reduction_plot(parts$pca, dims = c(1, 9)),
               "asks for coordinate 9")
  expect_error(draw_dim_reduction_plot(parts$pca, dims = 1),
               "`dims` must be a numeric vector of length 2")
})

test_that("a clustering of the reduction's scores shares its point labels", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()
  on_scores <- cluster_kmeans(parts$pca$scores, feats = c("PC1", "PC2"),
                              n_clust = 2, seed = 2026)
  expect_identical(on_scores$assignments$points, parts$pca$points)
  drawn <- draw_dim_reduction_plot(parts$pca, cluster_result = on_scores,
                                   cluster_lv = c("A", "B"))
  expect_identical(drawn$points, parts$pca$points)
  expect_identical(drawn$cluster, on_scores$assignments$cluster)
})

test_that("a clustering of other points is refused rather than lined up", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()

  other <- cluster_kmeans(sa_reduce_matrix(), n_clust = 2, seed = 1)
  expect_error(
    draw_dim_reduction_plot(parts$pca, cluster_result = other),
    "different points"
  )
  # Same frame, other margin: the two contracts disagree about what a row is.
  by_feat <- perform_pca(parts$m, embedding_scale = "features")
  expect_error(
    draw_dim_reduction_plot(by_feat, cluster_result = parts$kmeans),
    "different points"
  )
  feat_clust <- cluster_kmeans(parts$m, cluster_scale = "features",
                               n_clust = 2, seed = 1)
  expect_s3_class(
    draw_dim_reduction_plot(by_feat, cluster_result = feat_clust),
    "data.frame"
  )
  expect_error(draw_dim_reduction_plot(parts$pca, cluster_result = parts$m),
               "must be a clustering")
})

test_that("a clustering as the first argument names the reduction it needs", {
  parts <- sa_scatter_parts()
  expect_error(draw_dim_reduction_plot(parts$kmeans),
               "gives every point a label and no coordinate")
  expect_error(draw_dim_reduction_plot(parts$kmeans), "perform_pca")
  expect_error(draw_dim_reduction_plot(parts$m), "must be a reduction")
})

test_that("`cluster_lv` names the clusters in the legend", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()
  lv <- c("alpha", "beta", "gamma")
  drawn <- draw_dim_reduction_plot(parts$pca, cluster_result = parts$kmeans,
                                   cluster_lv = lv)
  expect_identical(drawn$cluster, parts$kmeans$assignments$cluster)
  expect_identical(
    sa_scatter_clusters(parts$kmeans, parts$pca$points, NULL, lv,
                        sa_plot_theme(FALSE))$levels,
    lv
  )
  expect_error(
    draw_dim_reduction_plot(parts$pca, cluster_lv = lv),
    "`cluster_lv` names the levels of the clustering"
  )
  expect_error(
    draw_dim_reduction_plot(parts$pca, cluster_result = parts$kmeans,
                            cluster_lv = lv[1:2]),
    "one label per cluster"
  )
  expect_error(
    draw_dim_reduction_plot(parts$pca, cluster_result = parts$kmeans,
                            cluster_lv = c("alpha", "beta", "alpha")),
    "must not repeat a level"
  )
})

test_that("`group_lv` orders the levels and drops no point", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()
  lv <- levels(parts$group)
  drawn <- draw_dim_reduction_plot(parts$pca, group = parts$group,
                                   group_lv = rev(lv))

  # Unlike draw_grouped_boxplot()'s, this argument selects no rows: the
  # reduction has already placed every point and none of them may vanish.
  expect_identical(nrow(drawn), length(parts$pca$points))
  expect_identical(levels(drawn$group), rev(lv))
  expect_error(
    draw_dim_reduction_plot(parts$pca, group = parts$group,
                            group_lv = lv[1:2]),
    "leaves out"
  )
  expect_error(draw_dim_reduction_plot(parts$pca, group = parts$group[-1]),
               "one label per point")
  expect_error(draw_dim_reduction_plot(parts$pca, group_lv = lv),
               "`group_lv` names the levels of `group`")
  expect_error(draw_dim_reduction_plot(parts$pca, group = parts$group,
                                       group_lv = c(lv, lv[1])),
               "must not repeat a level")
})

test_that("a group of one label per original row says which rows went", {
  skip_if_not_installed("withr")
  local_null_device()
  m <- sa_cluster_matrix()
  m[1, 1] <- NA_real_
  res <- suppressMessages(perform_pca(m))
  expect_identical(res$design$n_dropped, 1L)
  # The likeliest mistake is passing the grouping column of the frame that was
  # reduced, which is one longer than the points that survived it.
  expect_error(
    draw_dim_reduction_plot(res, group = factor(sa_cluster_truth())),
    "dropped 1 of the 30 row"
  )
})

test_that("the dim reduction plot leaves the device as it found it", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()
  before <- graphics::par(no.readonly = TRUE)
  draw_dim_reduction_plot(parts$pca, group = parts$group,
                          cluster_result = parts$kmeans, dark = TRUE,
                          anno_points = TRUE, asp = 1)
  after <- graphics::par(no.readonly = TRUE)
  expect_identical(after$bg, before$bg)
  expect_identical(after$fg, before$fg)
  expect_identical(after$col.axis, before$col.axis)
  expect_identical(after$mar, before$mar)
  # The legend sits in a panel of its own, so the layout has to come back too.
  expect_identical(after$mfrow, before$mfrow)
  # Nothing that carries an absolute size may be written, or the next plot is
  # pinned to the size this one was drawn at.
  expect_identical(after$fin, before$fin)
  expect_identical(after$pin, before$pin)
  expect_identical(after$mai, before$mai)

  withr::local_par(list(mfrow = c(2, 2)))
  draw_dim_reduction_plot(parts$pca)
  expect_identical(graphics::par("mfrow"), c(2L, 2L))
})

test_that("plot() on a reduction is this function", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()
  expect_identical(
    plot(parts$pca, group = parts$group, cluster_result = parts$kmeans),
    draw_dim_reduction_plot(parts$pca, group = parts$group,
                            cluster_result = parts$kmeans)
  )
})

test_that("the scatter refuses the arguments it cannot honour", {
  skip_if_not_installed("withr")
  local_null_device()
  parts <- sa_scatter_parts()
  expect_error(draw_dim_reduction_plot(parts$pca, dark = "yes"),
               "`dark` must be TRUE or FALSE")
  expect_error(draw_dim_reduction_plot(parts$pca, anno_points = NA),
               "`anno_points` must be TRUE or FALSE")
  expect_error(draw_dim_reduction_plot(parts$pca, cex = 0),
               "`cex` must be in \\(0, Inf\\]")
  expect_error(draw_dim_reduction_plot(parts$pca, asp = 0),
               "`asp` must be in \\(0, Inf\\]")
  expect_error(draw_dim_reduction_plot(parts$pca, xlim = 1),
               "`xlim` must be NULL or a finite numeric vector of length 2")
  expect_error(draw_dim_reduction_plot(parts$pca,
                                       cluster_result = parts$kmeans,
                                       col = c("red", "blue")),
               "one per cluster")
  # Eleven levels for ten shapes, which is where the shapes stop telling the
  # groups apart and start telling each other apart.
  many <- factor(rep(letters[1:11], length.out = length(parts$pca$points)))
  expect_error(draw_dim_reduction_plot(parts$pca, group = many),
               "shapes to tell them apart")
})


# draw_corrplot() adds three decisions to draw_heatmap() and no drawing of its
# own, so what it has to get right is those three: nothing is standardised, the
# two axes share one order, and a cell is blanked after the tree was built rather
# than before. The rest of the file already covers the cells themselves.

sa_corr_feats <- function() c("mpg", "cyl", "disp", "hp", "drat", "wt", "qsec")

sa_corr_result <- function(...) {
  summarize_association_stats(mtcars, sa_corr_feats(), methods = "pearson", ...)
}

test_that("both axes are drawn in one order, so the diagonal stays diagonal", {
  skip_if_not_installed("withr")
  local_null_device()
  out <- draw_corrplot(sa_corr_result())

  expect_identical(rownames(out$corr), colnames(out$corr))
  expect_identical(rownames(out$matrix), colnames(out$matrix))
  expect_identical(out$corr[!is.na(out$corr)],
                   t(out$corr)[!is.na(t(out$corr))])
  expect_equal(diag(out$corr), stats::setNames(rep(1, 7), rownames(out$corr)))
})

test_that("the order is the tree draw_heatmap() would build on the same terms", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_corr_result()
  out <- draw_corrplot(res, hclust_method = "average")

  # `1 - cor()` is what sa_cluster_dist() means by "correlation", so a corrplot
  # and a heatmap of the same features group them the same way.
  same <- stats::hclust(stats::as.dist(1 - res$pearson$corr), method = "average")
  expect_identical(out$order, same$order)
  expect_identical(out$hclust$order, same$order)
  expect_identical(rownames(out$corr), sa_corr_feats()[out$order])
})

test_that("`cluster = FALSE` keeps the features in the order they arrived", {
  skip_if_not_installed("withr")
  local_null_device()
  out <- draw_corrplot(sa_corr_result(), cluster = FALSE)

  expect_identical(out$order, seq_len(7L))
  expect_identical(rownames(out$corr), sa_corr_feats())
  expect_null(out$hclust)
})

test_that("blanking a cell does not move the tree it was clustered by", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_corr_result()
  # The whole reason the masking happens after the clustering: a cell removed
  # first would change the distance and the picture would no longer be the
  # matrix the reader is being shown.
  bare <- draw_corrplot(res, sig_level = 1)
  masked <- draw_corrplot(res, sig_level = 0.05)

  expect_identical(bare$order, masked$order)
  expect_identical(bare$hclust$merge, masked$hclust$merge)
  expect_identical(bare$n_masked, 0L)
  expect_gt(masked$n_masked, 0L)
})

test_that("the blanked cells are the ones above sig_level, diagonal excepted", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_corr_result()
  out <- draw_corrplot(res, sig_level = 0.05)

  above <- !is.na(out$pvalue) & out$pvalue > 0.05
  diag(above) <- FALSE
  expect_identical(out$n_masked, sum(above))
  expect_true(all(is.na(out$corr[above])))
  # A feature is not tested against itself, so the diagonal has no p-value to
  # fail and survives whatever the level is.
  expect_false(anyNA(diag(out$corr)))
  # Everything else keeps the coefficient it arrived with.
  kept <- !above
  ord <- out$order
  expect_equal(out$corr[kept], res$pearson$corr[ord, ord][kept])
})

test_that("`use_adjusted` chooses which p-value the blanking reads", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_corr_result()
  ord <- draw_corrplot(res)$order

  adjusted <- draw_corrplot(res, use_adjusted = TRUE)
  raw <- draw_corrplot(res, use_adjusted = FALSE)

  expect_equal(adjusted$pvalue, res$pearson$adj_pvalue[ord, ord])
  expect_equal(raw$pvalue, res$pearson$pvalue[ord, ord])
  # An adjusted p-value is never the smaller of the two, so it can only blank
  # more cells, never fewer.
  expect_gte(adjusted$n_masked, raw$n_masked)
})

test_that("a bare matrix is drawn, and a p-value matrix can come beside it", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_corr_result()

  bare <- draw_corrplot(res$pearson$corr)
  expect_null(bare$pvalue)
  expect_identical(bare$n_masked, 0L)
  expect_identical(bare$order, draw_corrplot(res)$order)

  # The two ways in have to reach the same picture, or the convenience of
  # handing over the whole result would be a second behaviour.
  by_hand <- draw_corrplot(res$pearson$corr, pvalue = res$pearson$adj_pvalue)
  by_result <- draw_corrplot(res)
  expect_equal(by_hand$corr, by_result$corr)
  expect_identical(by_hand$n_masked, by_result$n_masked)
})

test_that("`method` names the slot to draw out of a result", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- summarize_association_stats(mtcars, sa_corr_feats(),
                                     methods = c("pearson", "spearman"))
  # NULL takes the first method.
  expect_equal(draw_corrplot(res, cluster = FALSE)$corr,
               draw_corrplot(res, method = "pearson", cluster = FALSE)$corr)
  spearman <- draw_corrplot(res, method = "spearman", cluster = FALSE,
                            sig_level = 1)
  expect_equal(spearman$corr, res$spearman$corr)
})

test_that("nothing is standardised and the colour range is a correlation's", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_corr_result()
  out <- draw_corrplot(res, sig_level = 1)

  # draw_heatmap() z-scores by default, which would replace the coefficient the
  # reader came for with a number about the other coefficients in its row.
  expect_equal(out$matrix, out$corr)
  expect_identical(out$zlim, c(-1, 1))
  # Fixed rather than derived, so the same colour means the same strength from
  # one plot to the next.
  expect_identical(draw_corrplot(res$pearson$corr[1:3, 1:3])$zlim, c(-1, 1))
})

test_that("a feature with no correlation leaves the order alone rather than failing", {
  skip_if_not_installed("withr")
  local_null_device()
  d <- mtcars[sa_corr_feats()]
  d$flat <- 1
  res <- suppressMessages(summarize_association_stats(d, methods = "pearson"))

  expect_message(out <- draw_corrplot(res), "order they arrived")
  expect_identical(out$order, seq_len(8L))
  expect_null(out$hclust)
})

test_that("the corrplot refuses the arguments it cannot honour", {
  skip_if_not_installed("withr")
  local_null_device()
  res <- sa_corr_result()
  corr <- res$pearson$corr

  expect_error(draw_corrplot(corr[, 1:3]), "must be square")
  expect_error(draw_corrplot(corr[1:1, 1:1, drop = FALSE]),
               "at least 2 features")
  expect_error(draw_corrplot(matrix(c(1, 0.5, 0.2, 1), 2, 2)),
               "must be symmetric")
  expect_error(draw_corrplot(matrix(c(1, 2, 2, 1), 2, 2)),
               "outside \\[-1, 1\\]")
  expect_error(draw_corrplot(res, method = "spearman"),
               "must name one of the methods")
  expect_error(draw_corrplot(corr, method = "pearson"), "Leave it NULL")
  expect_error(draw_corrplot(res, pvalue = res$pearson$pvalue),
               "`pvalue` cannot be given")
  expect_error(draw_corrplot(corr, pvalue = corr[1:3, 1:3]),
               "laid out like `cor_matrix`")
  expect_error(draw_corrplot(res, sig_level = 0), "`sig_level` must be in")
  expect_error(draw_corrplot(res, cluster = NA), "`cluster` must be TRUE")
  expect_error(draw_corrplot(res, use_adjusted = "yes"),
               "`use_adjusted` must be TRUE")
  expect_error(draw_corrplot(res, cex.anno = 0), "`cex.anno` must be in")
  expect_error(draw_corrplot(res, cex.axis = 0), "`cex.axis` must be in")
})

test_that("cex.anno and cex.axis are accepted on the signature", {
  skip_if_not_installed("withr")
  local_null_device()
  expect_silent(draw_corrplot(sa_corr_result(), cex.anno = 0.8, cex.axis = 0.7))
})
