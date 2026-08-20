# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# The second density method, and the one whose arguments are easiest to confuse
# with another function's. `eps` here counts shared neighbours and `eps` in
# `cluster_dbscan()` is a radius, so the test that matters most is the one holding
# them apart: the same number passed to both has to mean two different things and
# the wrong kind of value has to be refused rather than quietly accepted.
#
# The silhouette kernel is tested here too, against the definition written out by
# hand on a case small enough to check on paper.

test_that("the labels are sNNclust() run on the scaled matrix", {
  m <- sa_cluster_matrix()
  res <- suppressMessages(cluster_snn(m, k = 6, eps = 3, min_pts = 3))

  expect_identical(res$analysis, "snn")

  ref <- dbscan::sNNclust(base::scale(m), k = 6, eps = 3, minPts = 3)
  expect_identical(res$assignments$cluster, as.integer(ref$cluster))

  expect_identical(res$engine$package, "dbscan")
  expect_identical(res$engine$method, "sNNclust")
  expect_identical(res$parameters$k, 6L)
  expect_identical(res$parameters$eps, 3L)
  expect_identical(res$parameters$min_pts, 3L)
})


test_that("`analysis` matches the function name", {
  res <- suppressMessages(cluster_snn(sa_cluster_matrix()))

  expect_identical(res$analysis, "snn")
})


test_that("the count of clusters is a result rather than an argument", {
  res <- suppressMessages(cluster_snn(sa_cluster_matrix()))

  expect_false("n_clust" %in% names(formals(cluster_snn)))
  expect_identical(res$design$n_clusters, 3L)
  expect_identical(res$assignments$cluster, as.integer(sa_cluster_truth()))
})


test_that("`eps` counts shared neighbours here and is a radius in dbscan", {
  m <- sa_cluster_matrix()

  # 3 is a perfectly good count of shared neighbours and a perfectly good radius,
  # and the two give different answers because they are different quantities.
  shared <- suppressMessages(cluster_snn(m, k = 6, eps = 3, min_pts = 3))
  radius <- suppressMessages(cluster_dbscan(m, eps = 3, min_pts = 3))
  expect_false(identical(shared$assignments$cluster,
                         radius$assignments$cluster))

  # A count cannot be fractional, and it cannot exceed the neighbours each point
  # keeps. A radius is under neither restriction.
  expect_error(cluster_snn(m, k = 6, eps = 2.5), "whole number")
  expect_error(cluster_snn(m, k = 6, eps = 7),
               "cannot share more neighbours than they each keep")
  expect_s3_class(suppressMessages(cluster_dbscan(m, eps = 2.5, min_pts = 3)),
                  "sa_cluster")
})


test_that("the three arguments are derived together and reported", {
  m <- sa_cluster_matrix()
  expect_message(cluster_snn(m), "Pass `k` to set it")

  res <- suppressMessages(cluster_snn(m))
  # 30 points, so k is the square root rounded up, and the other two are half of
  # it.
  expect_identical(res$parameters$k, 6L)
  expect_identical(res$parameters$eps, 3L)
  expect_identical(res$parameters$min_pts, 3L)

  # Naming them takes the same path the derivation would have taken.
  named <- suppressMessages(cluster_snn(m, k = 6, eps = 3, min_pts = 3))
  expect_identical(named$assignments$cluster, res$assignments$cluster)

  expect_error(cluster_snn(m, k = nrow(m)),
               "must not exceed one less than the 30 usable sample")
  expect_error(cluster_snn(m, k = 1), "must be in \\[2, Inf\\]")

  # Checked at the boundary, so a call that is going to fail does not first
  # announce a `k` it never used.
  expect_silent(try(cluster_snn(m, eps = 2.5), silent = TRUE))
})


test_that("noise is cluster 0 and every point being noise is an answer", {
  m <- sa_cluster_matrix()
  res <- suppressMessages(cluster_snn(m, k = 6, eps = 6, min_pts = 6))

  expect_message(cluster_snn(m, k = 6, eps = 6, min_pts = 6),
                 "every one of them is noise")
  expect_identical(res$design$n_clusters, 0L)
  expect_identical(res$design$n_noise, nrow(m))
  expect_identical(nrow(res$clusters), 0L)
  expect_true(all(is.na(res$assignments$silhouette)))
  expect_s3_class(res, "sa_cluster")
})


test_that("`cluster_scale` turns the points", {
  m <- sa_cluster_matrix()
  res <- suppressMessages(cluster_snn(m, cluster_scale = "features", k = 2,
                                      eps = 1, min_pts = 1))

  expect_identical(res$design$point_type, "feature")
  expect_identical(res$points, colnames(m))
  expect_identical(nrow(res$assignments), ncol(m))
  expect_identical(res$design$n_samples, nrow(m))
  expect_identical(res$design$n_feats, ncol(m))
})


test_that("the silhouette is Rousseeuw's definition, including its conventions", {
  # Four points on a line: two at 0 and 1, two at 10 and 11. Small enough that
  # a and b can be worked out on paper.
  d <- stats::dist(matrix(c(0, 1, 10, 11), ncol = 1))
  s <- sa_silhouette(d, c(1L, 1L, 2L, 2L))
  # Point 1: a = 1, b = mean(10, 11) = 10.5, so (10.5 - 1) / 10.5.
  expect_equal(s, c((10.5 - 1) / 10.5, (9.5 - 1) / 9.5,
                    (9.5 - 1) / 9.5, (10.5 - 1) / 10.5))

  # A point alone in its cluster scores 0 rather than 1: it has no `a`, and the
  # alternative would make a singleton the best-placed point in the data.
  expect_identical(sa_silhouette(d, c(1L, 2L, 2L, 2L))[1], 0)

  # Noise scores NA and takes no part in anyone else's a or b.
  with_noise <- sa_silhouette(d, c(1L, 1L, 2L, 0L))
  expect_true(is.na(with_noise[4]))
  expect_equal(with_noise[1], (10 - 1) / 10)

  # A single cluster has nothing to be compared with.
  expect_true(all(is.na(sa_silhouette(d, c(1L, 1L, 1L, 1L)))))
  expect_true(all(is.na(sa_silhouette(d, c(0L, 0L, 0L, 0L)))))

  # Coincident points are a tie rather than a division by zero.
  flat <- stats::dist(matrix(c(0, 0, 0, 0), ncol = 1))
  expect_identical(sa_silhouette(flat, c(1L, 1L, 2L, 2L)), rep(0, 4))

  # A separated fixture scores high, a labelling that ignores it scores low.
  m <- base::scale(sa_cluster_matrix())
  dm <- stats::dist(m)
  expect_gt(mean(sa_silhouette(dm, as.integer(sa_cluster_truth()))), 0.8)
  expect_lt(mean(sa_silhouette(dm, rep(1:3, times = 10))), 0.1)
})
