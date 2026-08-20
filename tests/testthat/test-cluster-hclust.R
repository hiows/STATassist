# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# A tree is judged on two things. First, that the labels really are the tree cut at
# `n_clust` and not something reconstructed beside it, which is pinned by
# reproducing the `hclust()` and `cutree()` calls directly. Second, that the
# distance is the one `draw_heatmap()` draws its dendrogram from, since the two
# claiming the same clustering while computing two of them is the failure a reader
# has no way to see.
#
# This is also where the shared input reading is tested, since all four clusterings
# read a frame through `sa_cluster_input()` and this is the one with no random
# number in it.

test_that("the labels are the tree cut at `n_clust`", {
  m <- sa_cluster_matrix()
  res <- cluster_hclust(m, n_clust = 3)

  expect_identical(res$analysis, "hclust")
  expect_s3_class(res$fit, "hclust")
  expect_identical(res$fit$method, "average")

  ref <- stats::hclust(stats::dist(base::scale(m)), method = "average")
  expect_equal(res$fit$merge, ref$merge)
  expect_equal(res$fit$height, ref$height)
  expect_identical(res$assignments$cluster,
                   unname(stats::cutree(ref, k = 3)))

  expect_identical(res$parameters$n_clust, 3L)
  expect_identical(res$parameters$dist_method, "euclidean")
  expect_identical(res$parameters$hclust_method, "average")
  expect_identical(res$engine$package, "stats")
  expect_identical(res$engine$method, "hclust")
})


test_that("a tree places every point, so there is never any noise", {
  m <- sa_cluster_matrix()
  for (k in 2:5) {
    res <- cluster_hclust(m, n_clust = k)
    expect_identical(res$design$n_noise, 0L)
    expect_identical(res$design$n_clusters, k)
    expect_false(any(res$assignments$cluster == 0L))
    expect_false(anyNA(res$assignments$silhouette))
  }
})


test_that("the tree does not depend on `n_clust`, so `$fit` can be cut again", {
  m <- sa_cluster_matrix()
  two <- cluster_hclust(m, n_clust = 2)
  three <- cluster_hclust(m, n_clust = 3)

  expect_equal(two$fit$merge, three$fit$merge)
  expect_equal(two$fit$height, three$fit$height)
  expect_identical(two$assignments$cluster,
                   unname(stats::cutree(three$fit, k = 2)))
})


test_that("`cluster_scale` turns the points without transposing the input", {
  m <- sa_cluster_matrix()
  res <- cluster_hclust(m, cluster_scale = "features", n_clust = 2)

  expect_identical(res$design$point_type, "feature")
  expect_identical(res$parameters$cluster_scale, "features")
  expect_identical(res$points, colnames(m))
  expect_identical(nrow(res$assignments), ncol(m))

  # The features are standardised and the transpose is then clustered, which is
  # not the same as clustering the transpose: scale() standardises columns, so on
  # t(m) it would standardise samples and answer a third question. The trees
  # differ, which is the claim; whether a given cut of them differs depends on how
  # coarse the cut is, and on four features a two-way cut is very coarse.
  ref <- stats::hclust(stats::dist(t(base::scale(m))), method = "average")
  expect_identical(res$assignments$cluster, unname(stats::cutree(ref, k = 2)))

  hand <- cluster_hclust(t(m), n_clust = 2)
  expect_equal(hand$fit$height,
               stats::hclust(stats::dist(base::scale(t(m))),
                             method = "average")$height)
  expect_false(isTRUE(all.equal(res$fit$height, hand$fit$height)))

  # `design` describes the input, so its counts do not turn with the argument.
  expect_identical(res$design$n_samples, nrow(m))
  expect_identical(res$design$n_feats, ncol(m))
  expect_identical(res$design$n_used, nrow(m))
})


test_that("the distance is the one draw_heatmap() draws its dendrogram from", {
  m <- sa_cluster_matrix()

  for (dm in c("euclidean", "manhattan", "correlation")) {
    res <- cluster_hclust(m, cluster_scale = "features", n_clust = 2,
                          dist_method = dm, hclust_method = "complete")
    # The heatmap transposes so that features run down the rows, and scales them
    # first, which is the matrix this clusters on the feature scale.
    ref <- sa_heatmap_hclust(t(base::scale(m)), dm, "complete", "feature")
    expect_equal(res$fit$merge, ref$merge)
    expect_equal(res$fit$height, ref$height)
  }
})


test_that("the linkage and the distance are both the caller's", {
  m <- sa_cluster_matrix()

  for (hm in c("average", "complete", "ward.D2")) {
    res <- cluster_hclust(m, n_clust = 3, hclust_method = hm)
    expect_identical(res$fit$method, hm)
    expect_identical(res$parameters$hclust_method, hm)
  }

  manhattan <- cluster_hclust(m, n_clust = 3, dist_method = "manhattan")
  expect_equal(manhattan$fit$height,
               stats::hclust(stats::dist(base::scale(m), method = "manhattan"),
                             method = "average")$height)
  expect_identical(manhattan$parameters$dist_method, "manhattan")
})


test_that("a distance that is not defined everywhere is an error, not a fallback", {
  # A sample with no variance across the features has no correlation with
  # anything. draw_heatmap() keeps the input order and says so, because a picture
  # without a tree is still a picture; here the tree is the whole answer.
  m <- sa_cluster_matrix()
  m[1, ] <- 5
  expect_error(
    cluster_hclust(m, dist_method = "correlation", n_clust = 2,
                   center = FALSE, scale = FALSE),
    "cannot measure the correlation distance"
  )
  # Standardising the features gives that row a profile again, so the same matrix
  # goes through on the defaults. The flatness was in the units, not in the data.
  expect_s3_class(cluster_hclust(m, dist_method = "correlation", n_clust = 2),
                  "sa_cluster")
  # The other two are defined for a flat row, so they are not governed by this.
  expect_s3_class(
    cluster_hclust(m, n_clust = 2, center = FALSE, scale = FALSE),
    "sa_cluster"
  )
})


test_that("`n_clust` is checked against the points, not left to the engine", {
  m <- sa_cluster_matrix()

  expect_error(cluster_hclust(m, n_clust = 1), "must be in \\[2, Inf\\]")
  expect_error(cluster_hclust(m, n_clust = 2.5), "whole number")
  expect_error(cluster_hclust(m, n_clust = nrow(m) + 1L),
               "must not exceed the 30 usable sample")
  # On the feature scale the limit is quoted in features, since that is what a
  # point is there.
  expect_error(cluster_hclust(m, cluster_scale = "features", n_clust = 5),
               "usable feature")
})


test_that("the input is read the way the reductions read it", {
  m <- sa_cluster_matrix()

  # A non-numeric column is left out rather than being an error, so a frame
  # carrying its grouping column can be passed as it arrived.
  frame <- as.data.frame(m)
  frame$group <- rep(c("a", "b", "c"), each = 10)
  expect_message(cluster_hclust(frame, n_clust = 3), "non-numeric column")
  quiet <- suppressMessages(cluster_hclust(frame, n_clust = 3))
  expect_identical(quiet$design$feats, colnames(m))

  # Rows that cannot be measured go before the distance is taken.
  holed <- m
  holed[2, 1] <- NA
  holed[5, 3] <- Inf
  expect_message(cluster_hclust(holed, n_clust = 3), "not complete and finite")
  dropped <- suppressMessages(cluster_hclust(holed, n_clust = 3))
  expect_identical(dropped$design$n_dropped, 2L)
  expect_identical(dropped$design$n_used, nrow(m) - 2L)
  expect_identical(dropped$points, rownames(m)[-c(2, 5)])

  # A feature that never moves cannot be scaled, so it is named and left out.
  flat <- cbind(m, f5 = 1)
  expect_message(cluster_hclust(flat, n_clust = 3), "no variance")
  left_out <- suppressMessages(cluster_hclust(flat, n_clust = 3))
  expect_identical(left_out$design$dropped_feats, "f5")
  expect_identical(left_out$design$n_feats, ncol(m))

  expect_error(cluster_hclust(iris["Species"]), "no numeric column")
  expect_error(cluster_hclust(m[, 1, drop = FALSE]),
               "at least 2 samples and 2 features")
  expect_error(cluster_hclust(list(a = 1)), "data.frame or a matrix")
})


test_that("`center` and `scale` are flags and they change the answer", {
  m <- sa_cluster_matrix()
  m[, 1] <- m[, 1] * 100

  scaled <- cluster_hclust(m, n_clust = 3)
  raw <- cluster_hclust(m, n_clust = 3, scale = FALSE)

  expect_true(scaled$parameters$scale)
  expect_false(raw$parameters$scale)
  # Without scaling the widest feature decides which points are close, which is
  # the whole reason scaling is on by default.
  expect_false(identical(scaled$assignments$cluster, raw$assignments$cluster))

  expect_error(cluster_hclust(m, center = "yes"), "`center` must be")
  expect_error(cluster_hclust(m, scale = NA), "`scale` must be")
})
