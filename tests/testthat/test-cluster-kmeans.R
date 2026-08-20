# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# The one clustering here whose answer depends on the random stream, so the two
# things worth pinning are that the arithmetic is `stats::kmeans()`'s own and that
# `seed` makes a run repeatable without leaving the caller's stream disturbed. The
# `nstart` override is checked because it is a default this package changed, and a
# changed default that is not declared is the kind of thing that goes unnoticed
# until someone compares against a textbook.

test_that("the labels are kmeans() run on the scaled matrix", {
  m <- sa_cluster_matrix()
  res <- cluster_kmeans(m, n_clust = 3, seed = 1)

  expect_identical(res$analysis, "kmeans")
  expect_s3_class(res$fit, "kmeans")

  set.seed(1)
  ref <- stats::kmeans(base::scale(m), centers = 3, nstart = 25,
                       iter.max = 100)
  expect_equal(res$fit$centers, ref$centers)
  expect_equal(res$fit$tot.withinss, ref$tot.withinss)
  expect_identical(res$parameters$tot_withinss, ref$tot.withinss)

  expect_identical(res$engine$package, "stats")
  expect_identical(res$engine$method, "kmeans")
})


test_that("the twenty-five starts are an override and are declared as one", {
  res <- cluster_kmeans(sa_cluster_matrix(), n_clust = 3, seed = 1)

  # stats::kmeans() defaults to a single start, which makes the answer a draw
  # rather than a search. Changing it is fine; changing it silently is not.
  expect_identical(res$parameters$n_start, 25L)
  expect_identical(res$engine$overridden, "nstart = 25")
  expect_equal(formals(stats::kmeans)$nstart, 1)
})


test_that("`seed` repeats a run and leaves the caller's stream alone", {
  m <- sa_cluster_matrix()

  expect_identical(cluster_kmeans(m, n_clust = 3, seed = 7)$assignments,
                   cluster_kmeans(m, n_clust = 3, seed = 7)$assignments)

  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  cluster_kmeans(m, n_clust = 3, seed = 7)
  expect_identical(stats::runif(1), before)

  expect_identical(cluster_kmeans(m, n_clust = 3, seed = 7)$parameters$seed, 7)
  expect_null(cluster_kmeans(m, n_clust = 3)$parameters$seed)
})


test_that("every point is given to a centre, so there is never any noise", {
  m <- sa_cluster_matrix()
  for (k in 2:5) {
    res <- cluster_kmeans(m, n_clust = k, seed = 1)
    expect_identical(res$design$n_noise, 0L)
    # k-means always returns the number it was asked for, which is the thing
    # `clusters$silhouette` is there to let a reader doubt.
    expect_identical(res$design$n_clusters, k)
    expect_false(anyNA(res$assignments$silhouette))
  }
})


test_that("`cluster_scale` turns the points", {
  m <- sa_cluster_matrix()
  res <- cluster_kmeans(m, cluster_scale = "features", n_clust = 2, seed = 1)

  expect_identical(res$design$point_type, "feature")
  expect_identical(res$points, colnames(m))
  expect_identical(nrow(res$assignments), ncol(m))

  set.seed(1)
  ref <- stats::kmeans(t(base::scale(m)), centers = 2, nstart = 25,
                       iter.max = 100)
  expect_equal(res$fit$centers, ref$centers)

  expect_identical(res$design$n_samples, nrow(m))
  expect_identical(res$design$n_feats, ncol(m))
})


test_that("the arguments are checked before the engine sees them", {
  m <- sa_cluster_matrix()

  expect_error(cluster_kmeans(m, n_clust = 1), "must be in \\[2, Inf\\]")
  expect_error(cluster_kmeans(m, n_clust = nrow(m) + 1L),
               "must not exceed the 30 usable sample")
  expect_error(cluster_kmeans(m, n_start = 0), "must be in \\[1, Inf\\]")
  expect_error(cluster_kmeans(m, iter_max = 0), "must be in \\[1, Inf\\]")
  expect_error(cluster_kmeans(m, seed = "one"), "single non-missing number")

  # More centres than distinct points is refused here rather than by the engine,
  # so the message names the argument the caller passed.
  dup <- rbind(c(1, 1), c(1, 1), c(2, 2), c(2, 2))
  dimnames(dup) <- list(paste0("s", 1:4), c("f1", "f2"))
  expect_error(cluster_kmeans(dup, n_clust = 3, scale = FALSE),
               "only 2 of the 4 sample\\(s\\).*are distinct")
})


test_that("the silhouette is measured on the distance k-means minimises", {
  res <- cluster_kmeans(sa_cluster_matrix(), n_clust = 3, seed = 1)

  # Not an argument here: a mean is only the centre of its group under squared
  # Euclidean distance, so offering a choice would offer an algorithm that does
  # not exist.
  expect_identical(res$parameters$dist_method, "euclidean")
  expect_false("dist_method" %in% names(formals(cluster_kmeans)))
})
