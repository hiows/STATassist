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

test_that("the x axis names what a multi-group log2fc compares", {
  multi <- estimate_significance(sa_multi_group_fixture(),
                                 log2fc_cutoff = 0.1)$significance
  lab <- sa_volcano_xlab(multi)
  expect_true(grepl("most extreme level vs setosa",
                    deparse(lab), fixed = TRUE))

  # A two-group verdict and a single contrast both compare two named levels, so
  # neither earns the qualifier.
  two <- estimate_significance(sa_two_group_fixture())$significance
  expect_identical(sa_volcano_xlab(two), expression(log[2] ~ FC))
  ct <- estimate_significance(sa_multi_group_fixture(),
                              by = "contrast")$significance
  expect_identical(sa_volcano_xlab(ct[[1]]), expression(log[2] ~ FC))
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

