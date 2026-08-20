# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# The first of the two methods that are not told how many groups to find, so the
# count and the noise are both results and both have to be checked as results.
# What is pinned here beyond the contract is that the derived arguments are
# derived by the rule the documentation states, since a derived argument that
# quietly changed would move every default answer in the package without moving a
# single test that names a number.

test_that("the labels are dbscan() run on the scaled matrix", {
  m <- sa_cluster_matrix()
  res <- suppressMessages(cluster_dbscan(m, eps = 1, min_pts = 4))

  expect_identical(res$analysis, "dbscan")
  expect_s3_class(res$fit, "dbscan_fast")

  ref <- dbscan::dbscan(base::scale(m), eps = 1, minPts = 4)
  expect_identical(res$assignments$cluster, as.integer(ref$cluster))

  expect_identical(res$engine$package, "dbscan")
  expect_identical(res$engine$method, "dbscan")
  expect_identical(res$parameters$eps, 1)
  expect_identical(res$parameters$min_pts, 4L)
  expect_identical(res$parameters$eps_source, "supplied")
})


test_that("the count of clusters is a result rather than an argument", {
  m <- sa_cluster_matrix()
  res <- suppressMessages(cluster_dbscan(m))

  expect_false("n_clust" %in% names(formals(cluster_dbscan)))
  expect_identical(res$design$n_clusters, 3L)
  expect_identical(res$assignments$cluster, as.integer(sa_cluster_truth()))
})


test_that("noise is cluster 0, has no silhouette and no row in `clusters`", {
  m <- sa_cluster_matrix()
  # A radius small enough to strand the edges of every blob but not the middles.
  res <- suppressMessages(cluster_dbscan(m, eps = 0.2, min_pts = 4))

  expect_gt(res$design$n_noise, 0L)
  expect_identical(res$design$n_noise, sum(res$assignments$cluster == 0L))
  expect_false(any(res$clusters$cluster == 0L))
  expect_true(all(is.na(res$assignments$silhouette[
    res$assignments$cluster == 0L])))
  expect_false(any(is.na(res$assignments$silhouette[
    res$assignments$cluster > 0L])))
  expect_identical(sum(res$clusters$size) + res$design$n_noise, nrow(m))
})


test_that("every point being noise is an answer and it says so", {
  m <- sa_cluster_matrix()
  expect_message(cluster_dbscan(m, eps = 0.001, min_pts = 4),
                 "every one of them is noise")

  res <- suppressMessages(cluster_dbscan(m, eps = 0.001, min_pts = 4))
  expect_identical(res$design$n_clusters, 0L)
  expect_identical(res$design$n_noise, nrow(m))
  expect_identical(nrow(res$clusters), 0L)
  expect_true(all(is.na(res$assignments$silhouette)))
  expect_s3_class(res, "sa_cluster")
})


test_that("`eps` is derived as the 95th percentile of the k-distance curve", {
  m <- sa_cluster_matrix()
  # Both derived arguments are announced, so the messages are captured together
  # rather than one at a time, which would let the other one through.
  expect_match(testthat::capture_messages(cluster_dbscan(m)),
               "Pass `eps` to set it", all = FALSE, fixed = TRUE)

  res <- suppressMessages(cluster_dbscan(m))
  ref <- stats::quantile(
    dbscan::kNNdist(base::scale(m), k = res$parameters$min_pts - 1L),
    0.95, names = FALSE
  )
  expect_equal(res$parameters$eps, ref)
  expect_identical(res$parameters$eps_source, "derived")

  # Supplying it takes the same path the derivation would have taken.
  supplied <- suppressMessages(cluster_dbscan(m, eps = ref))
  expect_identical(supplied$parameters$eps_source, "supplied")
  expect_identical(supplied$assignments$cluster, res$assignments$cluster)
})


test_that("`min_pts` is derived from the variables and capped at half the points", {
  m <- sa_cluster_matrix()
  expect_match(testthat::capture_messages(cluster_dbscan(m)),
               "Pass `min_pts` to set it", all = FALSE, fixed = TRUE)

  # 4 variables, so the textbook d + 1 is 5 and neither the floor of 4 nor the
  # cap at half of 30 binds.
  expect_identical(suppressMessages(cluster_dbscan(m))$parameters$min_pts, 5L)

  # Wider than it is tall, which is the shape that makes d + 1 unusable: the cap
  # is what keeps the threshold below the whole data set.
  wide <- cbind(m, matrix(stats::rnorm(nrow(m) * 40), nrow = nrow(m),
                          dimnames = list(NULL, paste0("w", 1:40))))
  expect_identical(suppressMessages(cluster_dbscan(wide))$parameters$min_pts,
                   nrow(m) %/% 2L)

  # Two variables would give 3, which fragments, so the floor of 4 binds.
  narrow <- m[, 1:2]
  expect_identical(suppressMessages(cluster_dbscan(narrow))$parameters$min_pts,
                   4L)
})


test_that("`cluster_scale` turns the points", {
  m <- sa_cluster_matrix()
  res <- suppressMessages(cluster_dbscan(m, cluster_scale = "features",
                                         eps = 1, min_pts = 2))

  expect_identical(res$design$point_type, "feature")
  expect_identical(res$points, colnames(m))
  expect_identical(nrow(res$assignments), ncol(m))

  ref <- dbscan::dbscan(t(base::scale(m)), eps = 1, minPts = 2)
  expect_identical(res$assignments$cluster, as.integer(ref$cluster))

  expect_identical(res$design$n_samples, nrow(m))
  expect_identical(res$design$n_feats, ncol(m))
})


test_that("the two arguments are checked before the engine sees them", {
  m <- sa_cluster_matrix()

  expect_error(cluster_dbscan(m, eps = 0), "must be in \\(0, Inf\\]")
  expect_error(cluster_dbscan(m, eps = -1), "must be in \\(0, Inf\\]")
  expect_error(cluster_dbscan(m, eps = c(1, 2)), "single non-missing number")
  expect_error(cluster_dbscan(m, min_pts = 1), "must be in \\[2, Inf\\]")
  expect_error(cluster_dbscan(m, min_pts = nrow(m) + 1L),
               "must not exceed the 30 usable sample")

  # Checked at the boundary, so a call that is going to fail does not first
  # announce a `min_pts` it never used.
  expect_silent(try(cluster_dbscan(m, eps = 0), silent = TRUE))
})
