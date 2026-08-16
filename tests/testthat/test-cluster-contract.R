# The `sa_cluster` contract, held to on two sides. The guards in `sa_new_cluster()`
# are called directly, since nothing a user can pass reaches them: they fire on a
# mistake inside the package and their job is to fail loudly when a future
# `cluster_*()` assembles the object wrong. The four functions are then checked to
# produce the same shape as each other, which is the promise that lets anything
# downstream read a clustering without being told which method made it.
#
# The whole object is never compared, because `metadata` carries a timestamp.

# Replacement is wholesale rather than by utils::modifyList(), which recurses into
# a data.frame and merges it column by column, so a `clusters` passed here would be
# blended with the default one instead of standing in for it.
sa_cluster_ok <- function(...) {
  args <- list(
    analysis = "kmeans",
    points   = c("a", "b", "c"),
    design   = list(point_type = "sample", n_samples = 3L, n_used = 3L,
                    n_dropped = 0L, n_feats = 2L, feats = c("f1", "f2"),
                    dropped_feats = character(0), n_clusters = 2L,
                    n_noise = 0L),
    parameters  = list(center = TRUE, scale = TRUE, dist_method = "euclidean"),
    assignments = data.frame(points = c("a", "b", "c"),
                             cluster = c(1L, 1L, 2L),
                             silhouette = c(0.5, 0.5, 0),
                             stringsAsFactors = FALSE),
    clusters    = data.frame(cluster = 1:2, size = c(2L, 1L),
                             silhouette = c(0.5, 0), stringsAsFactors = FALSE),
    engine      = list(package = "stats", method = "kmeans", label = "k-means",
                       overridden = character(0)),
    fit         = NULL
  )
  replacements <- list(...)
  for (nm in names(replacements)) args[nm] <- replacements[nm]
  do.call(sa_new_cluster, args)
}


test_that("the constructor accepts a well-formed clustering", {
  res <- sa_cluster_ok()

  expect_s3_class(res, "sa_cluster")
  expect_s3_class(res, "sa_result")
  expect_named(res, c("analysis", "points", "design", "parameters",
                      "assignments", "clusters", "engine", "fit", "metadata"))
})


test_that("the constructor refuses an analysis outside the four", {
  expect_error(sa_cluster_ok(analysis = "gmm"), "internal error")
  expect_error(sa_cluster_ok(analysis = "pca"), "must be one of")
  # The four the contract does cover, so a new one has to be declared to be used.
  expect_identical(sa_cluster_analyses(),
                   c("hclust", "kmeans", "dbscan", "snn"))
})


test_that("the constructor refuses tables that are not aligned with `points`", {
  expect_error(sa_cluster_ok(points = character(0)), "non-empty character")
  expect_error(sa_cluster_ok(points = c("a", "b")), "aligned with")

  shuffled <- data.frame(points = c("c", "b", "a"), cluster = c(1L, 1L, 2L),
                         silhouette = c(0.5, 0.5, 0), stringsAsFactors = FALSE)
  expect_error(sa_cluster_ok(assignments = shuffled), "aligned with")
})


test_that("the constructor refuses a cluster label that is not an integer", {
  as_double <- data.frame(points = c("a", "b", "c"), cluster = c(1, 1, 2),
                          silhouette = c(0.5, 0.5, 0), stringsAsFactors = FALSE)
  expect_error(sa_cluster_ok(assignments = as_double), "must be an integer")

  missing <- data.frame(points = c("a", "b", "c"), cluster = c(1L, 1L, NA_integer_),
                        silhouette = c(0.5, 0.5, 0), stringsAsFactors = FALSE)
  expect_error(sa_cluster_ok(assignments = missing), "no missing value")
})


test_that("the two tables have to be one fact counted twice", {
  # A row for noise is the mistake this is guarding against: noise is not a
  # cluster, and giving it a row would make `nrow(clusters)` the wrong count.
  with_noise <- data.frame(cluster = 0:2, size = c(1L, 2L, 1L),
                           silhouette = c(NA, 0.5, 0), stringsAsFactors = FALSE)
  expect_error(sa_cluster_ok(clusters = with_noise), "never includes noise")

  short <- data.frame(cluster = 1L, size = 2L, silhouette = 0.5,
                      stringsAsFactors = FALSE)
  expect_error(sa_cluster_ok(clusters = short), "but the assignments hold")

  expect_error(
    sa_cluster_ok(design = utils::modifyList(sa_cluster_ok()$design,
                                             list(n_clusters = 5L))),
    "`clusters` has 2 row"
  )
})


test_that("`design$n_noise` has to be the noise that is actually there", {
  noisy <- data.frame(points = c("a", "b", "c"), cluster = c(1L, 1L, 0L),
                      silhouette = c(0.5, 0.5, NA), stringsAsFactors = FALSE)
  one <- data.frame(cluster = 1L, size = 2L, silhouette = 0.5,
                    stringsAsFactors = FALSE)

  expect_error(
    sa_cluster_ok(assignments = noisy, clusters = one,
                  design = utils::modifyList(sa_cluster_ok()$design,
                                             list(n_clusters = 1L))),
    "point\\(s\\) were left unassigned"
  )
  expect_s3_class(
    sa_cluster_ok(assignments = noisy, clusters = one,
                  design = utils::modifyList(sa_cluster_ok()$design,
                                             list(n_clusters = 1L,
                                                  n_noise = 1L))),
    "sa_cluster"
  )
})


test_that("the constructor refuses an engine that does not name itself", {
  for (nm in c("package", "method", "label", "overridden")) {
    bad <- sa_cluster_ok()$engine
    bad[[nm]] <- NULL
    expect_error(sa_cluster_ok(engine = bad), paste0("missing `", nm, "`"))
  }
})


test_that("the constructor refuses a point type that is not one of the two", {
  expect_error(
    sa_cluster_ok(design = utils::modifyList(sa_cluster_ok()$design,
                                             list(point_type = "row"))),
    "must be \"sample\" or \"feature\""
  )
})


test_that("all four functions return the same shape", {
  m <- sa_cluster_matrix()
  fits <- list(
    hclust = cluster_hclust(m, n_clust = 3),
    kmeans = cluster_kmeans(m, n_clust = 3, seed = 1),
    dbscan = suppressMessages(cluster_dbscan(m)),
    snn    = suppressMessages(cluster_snn(m))
  )

  for (nm in names(fits)) {
    res <- fits[[nm]]
    expect_s3_class(res, "sa_cluster")
    expect_named(res, c("analysis", "points", "design", "parameters",
                        "assignments", "clusters", "engine", "fit",
                        "metadata"))
    expect_identical(res$points, rownames(m))
    expect_identical(res$assignments$points, res$points)
    expect_named(res$assignments, c("points", "cluster", "silhouette"))
    expect_named(res$clusters, c("cluster", "size", "silhouette"))
    expect_identical(res$design$n_clusters, nrow(res$clusters))
    expect_identical(sum(res$clusters$size) + res$design$n_noise,
                     length(res$points))
    expect_identical(res$design$point_type, "sample")
    expect_identical(res$design$feats, colnames(m))
    expect_identical(res$design$n_samples, nrow(m))
    expect_identical(res$design$n_feats, ncol(m))
    expect_identical(res$design$n_dropped, 0L)
    # `analysis` names the method, so it is never the name of the family.
    expect_true(res$analysis %in% sa_cluster_analyses())
    expect_false(identical(res$analysis, "cluster"))
  }

  # The engines disagree about everything except this: on a fixture all four can
  # solve, all four solve it the same way.
  truth <- as.integer(sa_cluster_truth())
  for (nm in names(fits)) {
    expect_identical(fits[[nm]]$assignments$cluster, truth)
  }
})


test_that("the print method summarises rather than listing every point", {
  res <- cluster_hclust(sa_cluster_matrix(), n_clust = 3)

  expect_output(print(res), "<sa_cluster> hclust")
  expect_output(print(res), "clusters : 3")
  expect_output(print(res), "linkage")
  expect_output(print(res), "silhouette")
  # The points themselves are in `$assignments`, not in the summary.
  expect_false(any(grepl("s17", utils::capture.output(print(res)))))
  expect_output(expect_invisible(print(res)), "<sa_cluster>", fixed = TRUE)

  noisy <- suppressMessages(cluster_dbscan(sa_cluster_matrix(), eps = 0.1))
  expect_output(print(noisy), "left as noise")
})
